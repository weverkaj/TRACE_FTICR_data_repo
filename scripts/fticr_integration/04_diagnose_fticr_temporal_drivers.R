#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
})

get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)))
  }
  getwd()
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  sd(x)
}

safe_last <- function(x) {
  ok <- which(is.finite(x))
  if (length(ok) == 0L) return(NA_real_)
  x[tail(ok, 1L)]
}

safe_days_since <- function(x, dates, target_date) {
  ok <- which(is.finite(x))
  if (length(ok) == 0L) return(NA_real_)
  as.integer(target_date - dates[tail(ok, 1L)])
}

metric_from_depth <- function(depth) {
  fifelse(depth == 10L, "surface", fifelse(depth == 30L, "subsurface", "unknown"))
}

clean_feature_name <- function(x) {
  y <- make.names(x, unique = TRUE)
  gsub("\\.+", "_", y)
}

response_linear_baseline <- function(sample_dt, response, method = c("linear90", "nearest45")) {
  method <- match.arg(method)
  out <- rep(NA_real_, nrow(sample_dt))
  for (i in seq_len(nrow(sample_dt))) {
    row <- sample_dt[i]
    train <- sample_dt[
      plot == row$plot &
        depth_cm == row$depth_cm &
        sample_id != row$sample_id &
        !is.na(get(response))
    ][order(date)]
    if (nrow(train) == 0L) next
    target <- as.integer(row$date)
    train_days <- as.integer(train$date)
    if (method == "nearest45") {
      nearest_idx <- which.min(abs(train_days - target))
      if (length(nearest_idx) == 1L && abs(train_days[nearest_idx] - target) <= 45L) {
        out[i] <- train[[response]][nearest_idx]
      }
    } else if (method == "linear90") {
      prev_idx <- suppressWarnings(max(which(train_days < target), na.rm = TRUE))
      next_idx <- suppressWarnings(min(which(train_days > target), na.rm = TRUE))
      if (!is.finite(prev_idx) || !is.finite(next_idx)) next
      bracket_gap <- train_days[next_idx] - train_days[prev_idx]
      if (bracket_gap <= 90L) {
        frac <- (target - train_days[prev_idx]) / bracket_gap
        out[i] <- train[[response]][prev_idx] + frac * (train[[response]][next_idx] - train[[response]][prev_idx])
      }
    }
  }
  out
}

group_key <- function(dt, groups) {
  if (length(groups) == 0L) return(rep("all", nrow(dt)))
  x <- as.data.frame(dt[, ..groups])
  for (col in names(x)) {
    x[[col]] <- as.character(x[[col]])
    x[[col]][is.na(x[[col]]) | !nzchar(x[[col]])] <- "missing"
  }
  do.call(paste, c(x, sep = "\r"))
}

response_mean_baseline <- function(sample_dt, response, groups = character()) {
  out <- rep(NA_real_, nrow(sample_dt))
  for (fold in sort(unique(sample_dt$cv_fold))) {
    test_idx <- which(sample_dt$cv_fold == fold)
    train <- sample_dt[cv_fold != fold & is.finite(get(response))]
    if (nrow(train) == 0L) next
    global_mean <- mean(train[[response]], na.rm = TRUE)
    if (length(groups) == 0L) {
      out[test_idx] <- global_mean
      next
    }

    train_key <- group_key(train, groups)
    test_key <- group_key(sample_dt[test_idx], groups)
    group_means <- tapply(train[[response]], train_key, mean, na.rm = TRUE)
    pred <- as.numeric(group_means[test_key])
    pred[!is.finite(pred)] <- global_mean
    out[test_idx] <- pred
  }
  out
}

prepare_predictors <- function(dt, train_idx, test_idx, feature_cols) {
  train <- dt[train_idx, ..feature_cols]
  test <- dt[test_idx, ..feature_cols]

  x_train <- data.frame(row.names = seq_len(nrow(train)))
  x_test <- data.frame(row.names = seq_len(nrow(test)))

  for (col in feature_cols) {
    tr <- train[[col]]
    te <- test[[col]]

    if (is.character(tr) || is.factor(tr)) {
      tr_chr <- as.character(tr)
      te_chr <- as.character(te)
      tr_chr[is.na(tr_chr) | !nzchar(tr_chr)] <- "missing"
      te_chr[is.na(te_chr) | !nzchar(te_chr)] <- "missing"
      levels_all <- sort(unique(c(tr_chr, te_chr, "missing")))
      x_train[[col]] <- factor(tr_chr, levels = levels_all)
      x_test[[col]] <- factor(te_chr, levels = levels_all)
    } else {
      tr_num <- suppressWarnings(as.numeric(tr))
      te_num <- suppressWarnings(as.numeric(te))
      med <- suppressWarnings(stats::median(tr_num, na.rm = TRUE))
      if (!is.finite(med)) med <- 0
      any_missing <- any(!is.finite(tr_num)) || any(!is.finite(te_num))
      if (any_missing) {
        x_train[[paste0(col, "__missing")]] <- as.integer(!is.finite(tr_num))
        x_test[[paste0(col, "__missing")]] <- as.integer(!is.finite(te_num))
      }
      tr_num[!is.finite(tr_num)] <- med
      te_num[!is.finite(te_num)] <- med
      x_train[[col]] <- tr_num
      x_test[[col]] <- te_num
    }
  }

  keep <- vapply(x_train, function(x) length(unique(x)) > 1L, logical(1))
  x_train <- x_train[, keep, drop = FALSE]
  x_test <- x_test[, keep, drop = FALSE]
  list(train = x_train, test = x_test)
}

score_predictions <- function(obs, pred, scale_sd) {
  ok <- is.finite(obs) & is.finite(pred)
  if (!any(ok)) {
    return(data.table(n = 0L, rmse = NA_real_, mae = NA_real_, scaled_mae = NA_real_, r2 = NA_real_))
  }
  err <- pred[ok] - obs[ok]
  sst <- sum((obs[ok] - mean(obs[ok]))^2)
  data.table(
    n = sum(ok),
    rmse = sqrt(mean(err^2)),
    mae = mean(abs(err)),
    scaled_mae = mean(abs(err)) / scale_sd,
    r2 = if (sst > 0) 1 - sum(err^2) / sst else NA_real_
  )
}

script_dir <- get_script_dir()
trace_mcmc_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
out_dir <- file.path(trace_mcmc_dir, "fticr_integration", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fticr_sample_path <- file.path(out_dir, "fticr_sample_features.csv")
chemodiv_path <- file.path(trace_mcmc_dir, "machine_learning", "fticr_chemodiversity_metrics.csv")
flux_path <- file.path(trace_mcmc_dir, "data", "flux_daily_bgc.csv")

required_paths <- c(fticr_sample_path, chemodiv_path, flux_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required input files:\n", paste(missing_paths, collapse = "\n"))
}

message("Loading FTICR response metrics...")
fticr <- fread(fticr_sample_path)
fticr[, date := as.IDate(sample_date)]
fticr[, depth_cm := fifelse(grepl("^0-10", depth), 10L, fifelse(grepl("^10-30", depth), 30L, NA_integer_))]
fticr <- fticr[
  !is.na(plot) & !is.na(date) & !is.na(depth_cm),
  .(
    sample_num,
    plot,
    treatment,
    date,
    depth_cm,
    sample_run_n,
    total_peak_area_mean,
    formula_detected_n_mean,
    cue_wmean,
    nosc_wmean,
    lambda_wmean,
    h_to_c_wmean,
    o_to_c_wmean
  )
]

chemodiv <- fread(chemodiv_path)
chemodiv[, date := as.IDate(date)]
setnames(chemodiv, "mean_nosc", "mean_nosc_pa")
chemodiv <- chemodiv[
  !is.na(plot) & !is.na(date) & !is.na(depth_cm),
  .(
    sample_num,
    plot,
    treatment,
    date,
    depth_cm,
    richness,
    shannon_classes,
    evenness,
    mean_nosc_pa,
    mean_hc,
    mean_oc,
    frac_labile,
    frac_aromatic,
    frac_lignin
  )
]

sample_dt <- merge(
  fticr,
  chemodiv,
  by = c("sample_num", "plot", "treatment", "date", "depth_cm"),
  all = TRUE
)
sample_dt[, sample_id := .I]
sample_dt[, year := as.integer(format(date, "%Y"))]
sample_dt[, doy := as.integer(format(date, "%j"))]
sample_dt[, date_num := as.integer(date)]
sample_dt[, doy_sin := sin(2 * pi * doy / 365.25)]
sample_dt[, doy_cos := cos(2 * pi * doy / 365.25)]
sample_dt[, depth_group := metric_from_depth(depth_cm)]
sample_dt[, treatment := fifelse(is.na(treatment) | !nzchar(treatment), "unknown", treatment)]

response_cols <- c(
  "cue_wmean",
  "nosc_wmean",
  "lambda_wmean",
  "richness",
  "shannon_classes",
  "evenness",
  "frac_labile",
  "frac_aromatic",
  "frac_lignin"
)
response_cols <- response_cols[response_cols %in% names(sample_dt)]

message("Loading daily driver data...")
daily <- fread(flux_path, showProgress = FALSE)
daily[, date := as.IDate(date)]
daily <- daily[!is.na(plot) & !is.na(date)]

climate_cols <- intersect(
  c("flux", "vwc", "temperature", "temp_20_30", "temp_40_50", "vwc_20_30", "vwc_40_50"),
  names(daily)
)

bgc_requested <- c(
  "npoc_ug.mL_lysimiter.surface.chemistry",
  "tn_ug.mL_lysimiter.surface.chemistry",
  "po4_ugP.mL_lysimiter.surface.chemistry",
  "no3_ugN.mL_lysimiter.surface.chemistry",
  "nh4_ugN.mL_lysimiter.surface.chemistry",
  "npoc_mg.l_10cm_pwchemistry",
  "tdn_mg.l_10cm_pwchemistry",
  "npoc_mg.l_30cm_pwchemistry",
  "tdn_mg.l_30cm_pwchemistry",
  "tfvol_ml_plotmean_throughfall",
  "npoc_mgC.L_plotmean_throughfall",
  "tdn_mgN.L_plotmean_throughfall",
  "litterfall_mean_g.m2",
  "rootstock_g.m2_0-10_minirhizotron",
  "rootgrowth_g.m2.day_0-10_minirhizotron",
  "rootmortality_g.m2.day_0-10_minirhizotron",
  "rootstock_g.m2_10-20_minirhizotron",
  "rootgrowth_g.m2.day_10-20_minirhizotron",
  "rootmortality_g.m2.day_10-20_minirhizotron",
  "rootstock_g.m2_20-30_minirhizotron",
  "rootgrowth_g.m2.day_20-30_minirhizotron",
  "rootmortality_g.m2.day_20-30_minirhizotron"
)
bgc_cols <- intersect(bgc_requested, names(daily))
bgc_cols <- bgc_cols[!grepl("fticr", bgc_cols, ignore.case = TRUE)]

message("Building prior-window driver features...")
driver_rows <- vector("list", nrow(sample_dt))
for (i in seq_len(nrow(sample_dt))) {
  row <- sample_dt[i]
  p <- row$plot
  target_date <- row$date
  features <- data.table(sample_id = row$sample_id)

  for (w in c(7L, 30L, 90L)) {
    win <- daily[plot == p & date <= target_date & date > (target_date - w)]
    for (col in climate_cols) {
      nm <- clean_feature_name(paste0(col, "_mean", w, "d"))
      features[, (nm) := safe_mean(win[[col]])]
      nm_sd <- clean_feature_name(paste0(col, "_sd", w, "d"))
      features[, (nm_sd) := safe_sd(win[[col]])]
    }
  }

  for (col in bgc_cols) {
    win180 <- daily[plot == p & date <= target_date & date > (target_date - 180L)]
    win365 <- daily[plot == p & date <= target_date & date > (target_date - 365L)]
    nm_mean <- clean_feature_name(paste0(col, "_mean180d"))
    nm_last <- clean_feature_name(paste0(col, "_last365d"))
    nm_days <- clean_feature_name(paste0(col, "_days_since365d"))
    features[, (nm_mean) := safe_mean(win180[[col]])]
    features[, (nm_last) := safe_last(win365[[col]])]
    features[, (nm_days) := safe_days_since(win365[[col]], win365$date, target_date)]
  }

  driver_rows[[i]] <- features
}

driver_features <- rbindlist(driver_rows, fill = TRUE)
analysis_dt <- merge(sample_dt, driver_features, by = "sample_id", all.x = TRUE)

structure_features <- c("treatment", "depth_group")
calendar_features <- c("year", "date_num", "doy_sin", "doy_cos")
climate_features <- setdiff(names(analysis_dt)[grepl("_(mean|sd)(7|30|90)d$", names(analysis_dt))], character())
bgc_features <- setdiff(
  names(analysis_dt)[grepl("_(mean180d|last365d|days_since365d)$", names(analysis_dt))],
  character()
)

feature_sets <- list(
  calendar = calendar_features,
  structure_calendar = c(structure_features, calendar_features),
  structure_calendar_climate = c(structure_features, calendar_features, climate_features),
  structure_calendar_climate_bgc = c(structure_features, calendar_features, climate_features, bgc_features)
)
feature_sets <- lapply(feature_sets, function(x) intersect(unique(x), names(analysis_dt)))

mean_baselines <- list(
  global_mean = character(),
  depth_mean = "depth_group",
  trt_depth_mean = c("treatment", "depth_group"),
  plot_depth_mean = c("plot", "depth_group")
)

feature_block <- data.table(feature = unique(unlist(feature_sets)))
feature_block[, block := fcase(
  feature %in% structure_features, "experimental_structure",
  feature %in% calendar_features, "calendar_time",
  feature %in% climate_features, "climate_flux",
  default = "bgc_sparse"
)]

set.seed(24601)
folds <- sample(rep(seq_len(5L), length.out = nrow(analysis_dt)))
analysis_dt[, cv_fold := folds]

message("Running cross-validated driver models...")
cv_rows <- list()
pred_rows <- list()
importance_rows <- list()
row_i <- 1L
pred_i <- 1L
imp_i <- 1L

for (response in response_cols) {
  y_all <- analysis_dt[[response]]
  response_scale <- sd(y_all, na.rm = TRUE)
  if (!is.finite(response_scale) || response_scale == 0) next

  for (model_name in names(feature_sets)) {
    feature_cols <- feature_sets[[model_name]]
    cv_pred <- rep(NA_real_, nrow(analysis_dt))

    for (fold in sort(unique(analysis_dt$cv_fold))) {
      test_idx <- which(analysis_dt$cv_fold == fold & is.finite(y_all))
      train_idx <- which(analysis_dt$cv_fold != fold & is.finite(y_all))
      if (length(test_idx) == 0L || length(train_idx) < 20L) next

      prep <- prepare_predictors(analysis_dt, train_idx, test_idx, feature_cols)
      if (ncol(prep$train) == 0L) next

      rf <- ranger(
        x = prep$train,
        y = y_all[train_idx],
        num.trees = 300,
        min.node.size = 5,
        mtry = max(1L, floor(sqrt(ncol(prep$train)))),
        seed = 1000 + fold + nchar(response) + nchar(model_name)
      )
      cv_pred[test_idx] <- predict(rf, data = prep$test)$predictions
    }

    score <- score_predictions(y_all, cv_pred, response_scale)
    score[, response := response]
    score[, model := model_name]
    setcolorder(score, c("response", "model", "n", "rmse", "mae", "scaled_mae", "r2"))
    cv_rows[[row_i]] <- score
    row_i <- row_i + 1L

    pred_rows[[pred_i]] <- data.table(
      sample_id = analysis_dt$sample_id,
      response = response,
      model = model_name,
      observed = y_all,
      prediction = cv_pred
    )
    pred_i <- pred_i + 1L
  }

  # Baselines on the same response scale.
  for (baseline in c("linear90", "nearest45")) {
    pred <- response_linear_baseline(analysis_dt, response, baseline)
    score <- score_predictions(y_all, pred, response_scale)
    score[, response := response]
    score[, model := baseline]
    setcolorder(score, c("response", "model", "n", "rmse", "mae", "scaled_mae", "r2"))
    cv_rows[[row_i]] <- score
    row_i <- row_i + 1L
    pred_rows[[pred_i]] <- data.table(
      sample_id = analysis_dt$sample_id,
      response = response,
      model = baseline,
      observed = y_all,
      prediction = pred
    )
    pred_i <- pred_i + 1L
  }

  # Simple held-out means show whether driver models beat static site structure.
  for (baseline in names(mean_baselines)) {
    pred <- response_mean_baseline(analysis_dt, response, mean_baselines[[baseline]])
    score <- score_predictions(y_all, pred, response_scale)
    score[, response := response]
    score[, model := baseline]
    setcolorder(score, c("response", "model", "n", "rmse", "mae", "scaled_mae", "r2"))
    cv_rows[[row_i]] <- score
    row_i <- row_i + 1L
    pred_rows[[pred_i]] <- data.table(
      sample_id = analysis_dt$sample_id,
      response = response,
      model = baseline,
      observed = y_all,
      prediction = pred
    )
    pred_i <- pred_i + 1L
  }

  # Full-data model for importance.
  feature_cols <- feature_sets[["structure_calendar_climate_bgc"]]
  train_idx <- which(is.finite(y_all))
  prep_full <- prepare_predictors(analysis_dt, train_idx, train_idx, feature_cols)
  if (ncol(prep_full$train) > 0L && length(train_idx) >= 30L) {
    rf_full <- ranger(
      x = prep_full$train,
      y = y_all[train_idx],
      num.trees = 500,
      min.node.size = 5,
      importance = "permutation",
      mtry = max(1L, floor(sqrt(ncol(prep_full$train)))),
      seed = 4242 + nchar(response)
    )
    imp <- data.table(
      response = response,
      feature = names(rf_full$variable.importance),
      importance = as.numeric(rf_full$variable.importance)
    )
    imp[, feature_base := sub("__missing$", "", feature)]
    imp[feature_block, block := i.block, on = .(feature_base = feature)]
    imp[is.na(block), block := fifelse(grepl("__missing$", feature), "missingness_indicator", "derived")]
    importance_rows[[imp_i]] <- imp[order(-importance)]
    imp_i <- imp_i + 1L
  }
}

cv_metrics <- rbindlist(cv_rows, fill = TRUE)
predictions <- rbindlist(pred_rows, fill = TRUE)
importance <- rbindlist(importance_rows, fill = TRUE)

importance_positive <- importance[importance > 0]
block_importance <- importance_positive[
  ,
  .(positive_importance = sum(importance, na.rm = TRUE), top_feature = feature[which.max(importance)]),
  by = .(response, block)
]
block_importance[, block_share := positive_importance / sum(positive_importance), by = response]
setorder(block_importance, response, -block_share)

top_importance <- importance[order(response, -importance)][
  ,
  head(.SD, 15L),
  by = response
]

full_driver_model <- "structure_calendar_climate_bgc"
prediction_wide <- dcast(predictions, sample_id + response + observed ~ model, value.var = "prediction")
comparison_models <- intersect(
  c("linear90", "nearest45", names(mean_baselines), names(feature_sets)),
  names(prediction_wide)
)
overlap_rows <- list()
overlap_i <- 1L
for (anchor in c("linear90", "nearest45")) {
  if (!anchor %in% names(prediction_wide)) next
  for (resp in response_cols) {
    subset_dt <- prediction_wide[response == resp & is.finite(get(anchor))]
    subset_scale <- sd(subset_dt$observed, na.rm = TRUE)
    if (!is.finite(subset_scale) || subset_scale == 0 || nrow(subset_dt) == 0L) next
    for (model in comparison_models) {
      score <- score_predictions(subset_dt$observed, subset_dt[[model]], subset_scale)
      score[, response := resp]
      score[, anchor := anchor]
      score[, model := model]
      setcolorder(score, c("response", "anchor", "model", "n", "rmse", "mae", "scaled_mae", "r2"))
      overlap_rows[[overlap_i]] <- score
      overlap_i <- overlap_i + 1L
    }
  }
}
overlap_metrics <- rbindlist(overlap_rows, fill = TRUE)

message("Writing outputs...")
feature_out <- file.path(out_dir, "fticr_temporal_driver_sample_features.csv")
metrics_out <- file.path(out_dir, "fticr_temporal_driver_cv_metrics.csv")
pred_out <- file.path(out_dir, "fticr_temporal_driver_cv_predictions.csv")
overlap_out <- file.path(out_dir, "fticr_temporal_driver_overlap_metrics.csv")
importance_out <- file.path(out_dir, "fticr_temporal_driver_importance.csv")
block_importance_out <- file.path(out_dir, "fticr_temporal_driver_block_importance.csv")
summary_out <- file.path(out_dir, "fticr_temporal_driver_summary.txt")

fwrite(analysis_dt, feature_out)
fwrite(cv_metrics, metrics_out)
fwrite(predictions, pred_out)
fwrite(overlap_metrics, overlap_out)
fwrite(importance, importance_out)
fwrite(block_importance, block_importance_out)

best_models <- cv_metrics[model %in% names(feature_sets)][
  ,
  .SD[which.min(scaled_mae)],
  by = response
][order(scaled_mae)]

baseline_compare <- merge(
  cv_metrics[model == full_driver_model, .(response, driver_scaled_mae = scaled_mae, driver_r2 = r2, driver_n = n)],
  cv_metrics[model == "linear90", .(response, linear90_scaled_mae = scaled_mae, linear90_n = n)],
  by = "response",
  all = TRUE
)
baseline_compare <- merge(
  baseline_compare,
  cv_metrics[model == "nearest45", .(response, nearest45_scaled_mae = scaled_mae, nearest45_n = n)],
  by = "response",
  all = TRUE
)
baseline_compare[, driver_minus_linear90 := driver_scaled_mae - linear90_scaled_mae]
baseline_compare[, driver_minus_nearest45 := driver_scaled_mae - nearest45_scaled_mae]

summary_lines <- c(
  "FTICR temporal driver diagnostic",
  "================================",
  "",
  paste("Samples:", nrow(analysis_dt)),
  paste("Responses:", paste(response_cols, collapse = ", ")),
  paste("Climate features:", length(climate_features)),
  paste("Sparse BGC features:", length(bgc_features)),
  "",
  "Cross-validated performance (lower scaled_mae is better):",
  capture.output(print(cv_metrics[order(response, scaled_mae)])),
  "",
  "Best driver model by response:",
  capture.output(print(best_models)),
  "",
  "Driver model versus interpolation baselines:",
  capture.output(print(baseline_compare[order(response)])),
  "",
  "Apples-to-apples overlap scores against interpolation baselines:",
  capture.output(print(overlap_metrics[order(response, anchor, scaled_mae)])),
  "",
  "Positive permutation importance by block:",
  capture.output(print(block_importance)),
  "",
  "Top 15 feature importances by response:",
  capture.output(print(top_importance[, .(response, feature, block, importance)][order(response, -importance)])),
  "",
  "Caveats:",
  "- This is an exploratory driver screen, not causal attribution.",
  "- Cross-validation is sample-level 5-fold; out-of-year and out-of-plot validation should be added before treating drivers as stable.",
  "- Sparse BGC predictors are summarized by prior-window means/last values and median-imputed for modeling.",
  "- Existing FTICR-derived daily columns are excluded from predictors to avoid leakage."
)
writeLines(summary_lines, summary_out)

message("Wrote: ", feature_out)
message("Wrote: ", metrics_out)
message("Wrote: ", pred_out)
message("Wrote: ", overlap_out)
message("Wrote: ", importance_out)
message("Wrote: ", block_importance_out)
message("Wrote: ", summary_out)
