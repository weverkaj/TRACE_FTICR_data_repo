#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
})

get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)))
  }
  getwd()
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

extract_first_int <- function(x) {
  x <- as.character(x)
  has_digit <- grepl("\\d", x)
  out <- rep(NA_integer_, length(x))
  out[has_digit] <- suppressWarnings(as.integer(sub(".*?(\\d+).*", "\\1", x[has_digit])))
  out
}

normalize_treatment <- function(x) {
  y <- tolower(trimws(as.character(x)))
  ifelse(grepl("warm", y), "warmed", ifelse(grepl("control", y), "control", NA_character_))
}

classify_vk <- function(hc, oc, ai_mod) {
  fifelse(
    ai_mod >= 0.67,
    "condensed_aromatic",
    fifelse(
      ai_mod >= 0.5,
      "aromatic",
      fifelse(
        hc >= 1.5 & oc <= 0.3,
        "lipid",
        fifelse(
          hc >= 1.5 & oc > 0.3,
          "carbohydrate",
          fifelse(
            hc < 1.5 & hc >= 0.7 & oc <= 0.67,
            "lignin_like",
            fifelse(hc < 0.7, "condensed_aromatic", "other")
          )
        )
      )
    )
  )
}

transform_metric <- function(metric, value_plot) {
  if (metric == "richness") {
    return(log1p(value_plot))
  }
  if (metric %in% c("frac_lignin", "frac_n_containing")) {
    p <- value_plot / 100
    p <- pmin(pmax(p, 1e-4), 1 - 1e-4)
    return(qlogis(p))
  }
  value_plot
}

standardize <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

jaccard_distance <- function(a, b) {
  if (length(a) == 0L && length(b) == 0L) return(NA_real_)
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0L) return(NA_real_)
  1 - inter / uni
}

bray_curtis <- function(a, b) {
  denom <- sum(a + b, na.rm = TRUE)
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  sum(abs(a - b), na.rm = TRUE) / denom
}

fit_pair_lm <- function(dt, distance_col) {
  mod_dt <- dt[is.finite(get(distance_col)) & is.finite(log2_lag)]
  if (nrow(mod_dt) < 10L || length(unique(mod_dt$log2_lag)) < 2L) {
    return(data.table(n_pairs = nrow(mod_dt), slope_per_lag_doubling = NA_real_, intercept = NA_real_, r2 = NA_real_))
  }
  fit <- lm(mod_dt[[distance_col]] ~ mod_dt$log2_lag)
  data.table(
    n_pairs = nrow(mod_dt),
    slope_per_lag_doubling = unname(coef(fit)[[2]]),
    intercept = unname(coef(fit)[[1]]),
    r2 = summary(fit)$r.squared
  )
}

permute_pair_slope <- function(dt, distance_col, n_perm = 2000L, seed = 20260707L) {
  set.seed(seed)
  mod_dt <- dt[is.finite(get(distance_col)) & lag_days > 0]
  observed <- fit_pair_lm(mod_dt, distance_col)
  if (!is.finite(observed$slope_per_lag_doubling)) {
    observed[, `:=`(
      permutation_p_positive_slope = NA_real_,
      null_slope_mean = NA_real_,
      null_slope_sd = NA_real_
    )]
    return(observed)
  }

  sample_dates <- unique(rbind(
    mod_dt[, .(series_id, sample_num = sample_num_1, date = date_1)],
    mod_dt[, .(series_id, sample_num = sample_num_2, date = date_2)]
  ))

  null_slopes <- rep(NA_real_, n_perm)
  for (b in seq_len(n_perm)) {
    perm_dates <- sample_dates[
      ,
      .(sample_num, perm_date = sample(date, length(date), replace = FALSE)),
      by = series_id
    ]
    p1 <- perm_dates[, .(series_id, sample_num_1 = sample_num, perm_date_1 = perm_date)]
    p2 <- perm_dates[, .(series_id, sample_num_2 = sample_num, perm_date_2 = perm_date)]
    perm_dt <- merge(
      mod_dt[, .(series_id, sample_num_1, sample_num_2, value = get(distance_col))],
      p1,
      by = c("series_id", "sample_num_1"),
      all.x = TRUE
    )
    perm_dt <- merge(
      perm_dt,
      p2,
      by = c("series_id", "sample_num_2"),
      all.x = TRUE
    )
    perm_dt[, lag_days := abs(as.integer(perm_date_2 - perm_date_1))]
    perm_dt <- perm_dt[lag_days > 0]
    perm_dt[, log2_lag := log2(lag_days)]
    setnames(perm_dt, "value", distance_col)
    null_slopes[[b]] <- fit_pair_lm(perm_dt, distance_col)$slope_per_lag_doubling
  }

  observed[, `:=`(
    permutation_p_positive_slope = (1 + sum(null_slopes >= slope_per_lag_doubling, na.rm = TRUE)) /
      (1 + sum(is.finite(null_slopes))),
    null_slope_mean = mean(null_slopes, na.rm = TRUE),
    null_slope_sd = sd(null_slopes, na.rm = TRUE)
  )]
  observed
}

run_pair_tests <- function(pair_dt, distance_cols, group_cols, n_perm = 2000L) {
  groups <- unique(pair_dt[, ..group_cols])
  out <- list()
  idx <- 1L
  for (g_i in seq_len(nrow(groups))) {
    g <- groups[g_i]
    sub <- pair_dt[g, on = group_cols, nomatch = 0]
    for (dist_col in distance_cols) {
      res <- permute_pair_slope(
        sub,
        dist_col,
        n_perm = n_perm,
        seed = 20260707L + g_i + match(dist_col, distance_cols) * 1000L
      )
      res[, distance_metric := dist_col]
      out[[idx]] <- cbind(g, res)
      idx <- idx + 1L
    }
  }
  rbindlist(out, fill = TRUE)
}

script_dir <- get_script_dir()
trace_mcmc_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
docs_dir <- normalizePath(file.path(trace_mcmc_dir, ".."), winslash = "/", mustWork = FALSE)
out_dir <- file.path(trace_mcmc_dir, "fticr_integration", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_data_path <- file.path(out_dir, "fticr_formula_count_temporal_metrics_plot_data.csv")
mol_path <- file.path(
  docs_dir,
  "trace_data_cleanup",
  "data_raw",
  "processed_data",
  "Merged_Output",
  "Test_Processed-Unprocessed_Mol.csv"
)
intensity_path <- file.path(
  docs_dir,
  "trace_data_cleanup",
  "data_raw",
  "processed_data",
  "Merged_Output",
  "Test_Processed-Unprocssed_Data.csv"
)
sample_key_path <- file.path(
  docs_dir,
  "trace_data_cleanup",
  "data_raw",
  "60398_Sihi_Porewater_July82024_DS.csv"
)

required_paths <- c(plot_data_path, mol_path, intensity_path, sample_key_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required input files:\n", paste(missing_paths, collapse = "\n"))
}

message("Loading formula-count metric table...")
long_dt <- fread(plot_data_path)
long_dt[, date := as.IDate(date)]
long_dt[, series_id := interaction(plot, depth_label, drop = TRUE)]
long_dt[, time30 := as.numeric(date - min(date, na.rm = TRUE)) / 30]

metric_info <- unique(long_dt[, .(metric, metric_label, metric_order)])
setorder(metric_info, metric_order)

message("Fitting continuous-time scalar autocorrelation models...")
ctmm_rows <- list()
ctrl <- lmeControl(opt = "optim", msMaxIter = 250, maxIter = 150, returnObject = TRUE)
for (m in metric_info$metric) {
  d <- copy(long_dt[metric == m & is.finite(value_plot)])
  d[, y_transformed := transform_metric(m, value_plot)]
  d[, z := standardize(y_transformed)]
  d <- d[is.finite(z)]
  d[, treatment_label := factor(treatment_label)]
  d[, depth_label := factor(depth_label)]
  d[, series_id := factor(series_id)]

  fit_null <- try(lme(
    z ~ treatment_label * depth_label,
    random = ~1 | series_id,
    data = d,
    method = "ML",
    na.action = na.omit,
    control = ctrl
  ), silent = TRUE)

  fit_car <- try(lme(
    z ~ treatment_label * depth_label,
    random = ~1 | series_id,
    correlation = corCAR1(value = 0.6, form = ~time30 | series_id),
    data = d,
    method = "ML",
    na.action = na.omit,
    control = ctrl
  ), silent = TRUE)

  if (inherits(fit_null, "try-error") || inherits(fit_car, "try-error")) {
    ctmm_rows[[length(ctmm_rows) + 1L]] <- data.table(
      metric = m,
      metric_label = metric_info[metric == m, metric_label],
      metric_order = metric_info[metric == m, metric_order],
      n_samples = nrow(d),
      n_series = uniqueN(d$series_id),
      null_aic = if (!inherits(fit_null, "try-error")) AIC(fit_null) else NA_real_,
      car1_aic = if (!inherits(fit_car, "try-error")) AIC(fit_car) else NA_real_,
      delta_aic_car1_minus_null = NA_real_,
      lrt = NA_real_,
      lrt_p = NA_real_,
      phi_per_30d = NA_real_,
      corr_90d = NA_real_,
      corr_180d = NA_real_,
      corr_365d = NA_real_,
      half_life_days = NA_real_,
      fit_note = paste(
        if (inherits(fit_null, "try-error")) "null_failed" else "null_ok",
        if (inherits(fit_car, "try-error")) "car1_failed" else "car1_ok",
        sep = ";"
      )
    )
    next
  }

  cmp <- anova(fit_null, fit_car)
  phi <- as.numeric(coef(fit_car$modelStruct$corStruct, unconstrained = FALSE))
  half_life <- if (is.finite(phi) && phi > 0 && phi < 1) 30 * log(0.5) / log(phi) else NA_real_
  ctmm_rows[[length(ctmm_rows) + 1L]] <- data.table(
    metric = m,
    metric_label = metric_info[metric == m, metric_label],
    metric_order = metric_info[metric == m, metric_order],
    n_samples = nrow(d),
    n_series = uniqueN(d$series_id),
    null_aic = AIC(fit_null),
    car1_aic = AIC(fit_car),
    delta_aic_car1_minus_null = AIC(fit_car) - AIC(fit_null),
    lrt = cmp$L.Ratio[2],
    lrt_p = cmp$`p-value`[2],
    phi_per_30d = phi,
    corr_90d = phi^(90 / 30),
    corr_180d = phi^(180 / 30),
    corr_365d = phi^(365 / 30),
    half_life_days = half_life,
    fit_note = "ok"
  )
}
ctmm_results <- rbindlist(ctmm_rows, fill = TRUE)
setorder(ctmm_results, metric_order)

message("Loading detected formula sets for pairwise sample similarity...")
mol <- fread(mol_path)
formula_col <- "Molecular Formula"
if (!(formula_col %in% names(mol))) stop("Missing Molecular Formula column in molecular formula table.")
if ("Is Isotopologue" %in% names(mol)) {
  mol <- mol[safe_numeric(`Is Isotopologue`) == 0 | is.na(`Is Isotopologue`)]
}
element_cols <- c("C", "H", "O", "N", "S", "P")
missing_elements <- setdiff(element_cols, names(mol))
if (length(missing_elements) > 0L) {
  stop("Missing element columns in molecular formula table: ", paste(missing_elements, collapse = ", "))
}
for (nm in element_cols) {
  mol[, (nm) := safe_numeric(get(nm))]
}
mol <- mol[!is.na(C) & !is.na(H) & C > 0]
for (nm in c("O", "N", "S", "P")) mol[is.na(get(nm)), (nm) := 0]
mol[, HC := H / C]
mol[, OC := O / C]
mol[, NC := N / C]
mol[, NOSC := 4 - ((4 * C + H - 3 * N - 2 * O - 2 * S + 5 * P) / C)]
mol[, AI_mod_denominator := C - O - S - N - P]
mol[, AI_mod := (1 + C - O - S - 0.5 * (H + N + P)) / AI_mod_denominator]
mol[AI_mod_denominator <= 0 | !is.finite(AI_mod) | AI_mod < 0, AI_mod := 0]
mol[, vk_class := classify_vk(HC, OC, AI_mod)]
mol <- unique(mol, by = formula_col)
setkeyv(mol, formula_col)

intensity <- fread(intensity_path)
if (!(formula_col %in% names(intensity))) stop("Missing Molecular Formula column in intensity matrix.")
setnames(intensity, formula_col, "formula")

blank_cols <- grep("Blank", names(intensity), value = TRUE, ignore.case = TRUE)
field_cols <- grep("^Sihi_60398_[0-9]+_r[0-9]+_", names(intensity), value = TRUE)
if (length(field_cols) == 0L) stop("No field sample columns matched expected Sihi_60398 sample pattern.")

blank_formulas <- character()
if (length(blank_cols) > 0L) {
  blank_detected <- rowSums(!is.na(intensity[, ..blank_cols])) > 0
  blank_formulas <- intensity[blank_detected, formula]
}
valid_formula_set <- setdiff(mol[[formula_col]], blank_formulas)
formula_id <- data.table(formula = valid_formula_set, formula_index = seq_along(valid_formula_set))
formula_class <- merge(
  formula_id,
  mol[, .(formula = get(formula_col), vk_class)],
  by = "formula",
  all.x = TRUE
)

sample_map <- data.table(sample_col = field_cols)
sample_map[, sample_num := suppressWarnings(as.integer(sub("^Sihi_60398_([0-9]+)_r[0-9]+_.*$", "\\1", sample_col)))]
sample_map <- sample_map[!is.na(sample_num)]

key_raw <- fread(sample_key_path, encoding = "Latin-1")
key_raw[, sample_num := extract_first_int(sample_name)]
sample_key <- key_raw[!is.na(sample_num), .(
  sample_num,
  date = as.IDate(collection_date),
  plot = extract_first_int(experimental_factor_other),
  treatment_label = fifelse(
    normalize_treatment(fifelse(
      nzchar(trimws(as.character(experimental_factor))),
      trimws(as.character(experimental_factor)),
      trimws(as.character(climate_environment))
    )) == "control",
    "Control",
    "Warmed"
  ),
  depth = trimws(as.character(depth))
)]
sample_key[, depth_label := fifelse(grepl("^0-10", depth), "0-10 cm", fifelse(grepl("^10-30", depth), "10-30 cm", NA_character_))]
sample_key <- unique(sample_key[!is.na(date) & !is.na(plot) & !is.na(depth_label)], by = "sample_num")

target_samples <- unique(long_dt[, .(sample_num, plot, treatment_label, depth_label, date, series_id)])
sample_formula_rows <- vector("list", nrow(target_samples))
sample_class_rows <- vector("list", nrow(target_samples))
class_levels <- sort(unique(na.omit(formula_class$vk_class)))

for (i in seq_len(nrow(target_samples))) {
  snum <- target_samples$sample_num[[i]]
  cols <- sample_map[sample_num == snum, sample_col]
  if (length(cols) == 0L) next
  detected <- rowSums(!is.na(intensity[, ..cols])) > 0
  detected_formulas <- intersect(intensity[detected, formula], valid_formula_set)
  if (length(detected_formulas) < 5L) next
  ids <- formula_id[J(detected_formulas), on = "formula", nomatch = 0, formula_index]
  classes <- formula_class[J(detected_formulas), on = "formula", nomatch = 0, vk_class]
  class_counts <- table(factor(classes, levels = class_levels))
  class_props <- as.numeric(class_counts) / sum(class_counts)
  sample_formula_rows[[i]] <- data.table(
    sample_num = snum,
    formula_indices = list(sort(ids))
  )
  sample_class_rows[[i]] <- data.table(
    sample_num = snum,
    class_props = list(class_props)
  )
}
sample_formulas <- rbindlist(sample_formula_rows, fill = TRUE)
sample_classes <- rbindlist(sample_class_rows, fill = TRUE)

sample_meta <- merge(target_samples, sample_formulas, by = "sample_num", all.x = TRUE)
sample_meta <- merge(sample_meta, sample_classes, by = "sample_num", all.x = TRUE)
sample_meta <- sample_meta[!sapply(formula_indices, is.null)]

wide_metrics <- dcast(
  long_dt,
  sample_num + plot + treatment_label + depth_label + date + series_id ~ metric,
  value.var = "value_plot"
)
for (m in metric_info$metric) {
  wide_metrics[, paste0(m, "_z") := standardize(transform_metric(m, get(m)))]
}
z_cols <- paste0(metric_info$metric, "_z")
wide_metrics[, metric_vector := split(as.matrix(.SD), seq_len(.N)), .SDcols = z_cols]
sample_meta <- merge(sample_meta, wide_metrics[, .(sample_num, metric_vector)], by = "sample_num", all.x = TRUE)

message("Building pairwise within-plot/depth distances...")
pair_rows <- list()
idx <- 1L
for (grp_name in unique(sample_meta$series_id)) {
  grp <- sample_meta[series_id == grp_name][order(date)]
  if (nrow(grp) < 2L) next
  for (i in seq_len(nrow(grp) - 1L)) {
    for (j in seq.int(i + 1L, nrow(grp))) {
      lag_days <- as.integer(grp$date[[j]] - grp$date[[i]])
      if (!is.finite(lag_days) || lag_days <= 0L) next
      v1 <- as.numeric(grp$metric_vector[[i]])
      v2 <- as.numeric(grp$metric_vector[[j]])
      pair_rows[[idx]] <- data.table(
        series_id = grp$series_id[[1]],
        plot = grp$plot[[1]],
        treatment_label = grp$treatment_label[[1]],
        depth_label = grp$depth_label[[1]],
        sample_num_1 = grp$sample_num[[i]],
        sample_num_2 = grp$sample_num[[j]],
        date_1 = grp$date[[i]],
        date_2 = grp$date[[j]],
        lag_days = lag_days,
        log2_lag = log2(lag_days),
        formula_jaccard_distance = jaccard_distance(grp$formula_indices[[i]], grp$formula_indices[[j]]),
        vk_class_braycurtis = bray_curtis(grp$class_props[[i]], grp$class_props[[j]]),
        metric_vector_distance = sqrt(mean((v2 - v1)^2, na.rm = TRUE))
      )
      idx <- idx + 1L
    }
  }
}
pair_dt <- rbindlist(pair_rows, fill = TRUE)

distance_cols <- c("formula_jaccard_distance", "vk_class_braycurtis", "metric_vector_distance")
message("Running restricted-permutation pairwise distance tests...")
pair_tests <- rbindlist(list(
  run_pair_tests(copy(pair_dt)[, scope_dummy := "overall"], distance_cols, c("scope_dummy"), n_perm = 2000L)[, scope := "overall"],
  run_pair_tests(copy(pair_dt)[, scope_dummy := depth_label], distance_cols, c("scope_dummy"), n_perm = 2000L)[
    ,
    `:=`(scope = "by_depth", depth_label = scope_dummy)
  ],
  run_pair_tests(copy(pair_dt)[, scope_dummy := paste(treatment_label, depth_label, sep = " | ")], distance_cols, c("scope_dummy"), n_perm = 2000L)[
    ,
    c("treatment_label", "depth_label") := tstrsplit(scope_dummy, " \\| ")
  ][, scope := "by_treatment_depth"]
), fill = TRUE)
pair_tests[, scope_dummy := NULL]
setcolorder(pair_tests, c("scope", "distance_metric", setdiff(names(pair_tests), c("scope", "distance_metric"))))
setorder(pair_tests, scope, distance_metric, depth_label, treatment_label)

ctmm_out <- file.path(out_dir, "fticr_q2_scalar_continuous_time_autocorrelation.csv")
pairs_out <- file.path(out_dir, "fticr_q2_pairwise_temporal_distance_pairs.csv")
pair_tests_out <- file.path(out_dir, "fticr_q2_pairwise_temporal_distance_tests.csv")
summary_out <- file.path(out_dir, "fticr_q2_temporal_similarity_summary.txt")

fwrite(ctmm_results, ctmm_out)
fwrite(pair_dt, pairs_out)
fwrite(pair_tests, pair_tests_out)

ctmm_print <- ctmm_results[
  ,
  .(
    metric_label,
    n_samples,
    delta_aic_car1_minus_null,
    lrt_p,
    phi_per_30d,
    corr_90d,
    corr_180d,
    half_life_days
  )
]

pair_overall_print <- pair_tests[
  scope == "overall",
  .(
    distance_metric,
    n_pairs,
    slope_per_lag_doubling,
    r2,
    permutation_p_positive_slope
  )
]

pair_depth_print <- pair_tests[
  scope == "by_depth",
  .(
    distance_metric,
    depth_label,
    n_pairs,
    slope_per_lag_doubling,
    r2,
    permutation_p_positive_slope
  )
][order(distance_metric, depth_label)]

summary_lines <- c(
  "TRACE FTICR Q2 temporal similarity tests",
  "=======================================",
  "",
  "Question:",
  "How temporally variable is FT-ICR molecular composition, and does sample similarity depend on the resampling interval?",
  "",
  "Scalar continuous-time models:",
  "- Each feature was transformed where useful, z-scored, and modeled as treatment * depth with random plot-depth intercepts.",
  "- The temporal model adds continuous AR(1) residual correlation within plot-depth using time in 30-day units.",
  "- Delta AIC < 0 means the temporal correlation model improved fit over the no-correlation model.",
  "- Phi is the estimated 30-day residual correlation; corr_90d/corr_180d are implied correlations at those intervals.",
  capture.output(print(ctmm_print)),
  "",
  "Pairwise dissimilarity tests:",
  "- Within each plot-depth, pairwise dissimilarity was regressed on log2(days between samples).",
  "- Slopes are the expected increase in dissimilarity per doubling of resampling interval.",
  "- P-values are one-sided restricted permutations that shuffle dates within plot-depth series.",
  "",
  "Overall pairwise tests:",
  capture.output(print(pair_overall_print)),
  "",
  "By-depth pairwise tests:",
  capture.output(print(pair_depth_print)),
  "",
  "Outputs:",
  paste("-", ctmm_out),
  paste("-", pairs_out),
  paste("-", pair_tests_out),
  "",
  "Caveats:",
  "- Continuous-time mixed-model likelihood-ratio p-values are approximate because zero temporal correlation is a boundary case.",
  "- Pairwise dissimilarity tests handle non-independence by restricted date permutations, but pair counts still arise from repeated spectra.",
  "- The pairwise tests ask whether dissimilarity increases with time interval; they do not identify environmental drivers of that change."
)
writeLines(summary_lines, summary_out)

message("Wrote: ", ctmm_out)
message("Wrote: ", pairs_out)
message("Wrote: ", pair_tests_out)
message("Wrote: ", summary_out)
