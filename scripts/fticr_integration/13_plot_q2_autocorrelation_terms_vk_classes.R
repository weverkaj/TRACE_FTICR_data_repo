#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
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
  if (metric == "richness") return(log1p(value_plot))
  if (metric == "frac_n_containing" || grepl("^vk_frac_", metric)) {
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

clean_metric_label <- function(x) {
  gsub("\\n", " ", x)
}

theme_trace <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(color = "#202124"),
      axis.title = element_text(color = "#202124"),
      plot.margin = margin(5.5, 12, 5.5, 12)
    )
}

save_plot <- function(plot, png_path, pdf_path, width, height) {
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 320, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, units = "in", bg = "white")
  message("Wrote: ", png_path)
  message("Wrote: ", pdf_path)
}

fit_scalar_car1_by_group <- function(long_dt, metric_order_dt) {
  ctrl <- lmeControl(opt = "optim", msMaxIter = 200, maxIter = 100, returnObject = TRUE)
  rows <- list()
  idx <- 1L

  for (metric_name in metric_order_dt$metric) {
    metric_dt <- long_dt[metric == metric_name]
    for (dep in c("0-10 cm", "10-30 cm")) {
      for (trt in c("Control", "Warmed")) {
        d <- copy(metric_dt[depth_label == dep & treatment_label == trt])
        if (nrow(d) < 10L) next
        d[, date := as.IDate(date)]
        d[, series_id := factor(interaction(plot, depth_label, drop = TRUE))]
        d[, time30 := as.numeric(date - min(date, na.rm = TRUE)) / 30]
        d[, z := standardize(transform_metric(metric_name, value_plot))]
        d <- d[is.finite(z)]
        if (nrow(d) < 10L || uniqueN(d$series_id) < 2L) next

        fit0 <- try(lme(
          z ~ 1,
          random = ~1 | series_id,
          data = d,
          method = "ML",
          control = ctrl
        ), silent = TRUE)
        fit1 <- try(lme(
          z ~ 1,
          random = ~1 | series_id,
          correlation = corCAR1(value = 0.5, form = ~time30 | series_id),
          data = d,
          method = "ML",
          control = ctrl
        ), silent = TRUE)

        phi <- if (!inherits(fit1, "try-error")) {
          as.numeric(coef(fit1$modelStruct$corStruct, unconstrained = FALSE))
        } else {
          NA_real_
        }
        rows[[idx]] <- data.table(
          metric = metric_name,
          metric_label = unique(d$metric_label)[[1]],
          metric_order = unique(d$metric_order)[[1]],
          metric_group = unique(d$metric_group)[[1]],
          depth_label = dep,
          treatment_label = trt,
          n_samples = nrow(d),
          n_series = uniqueN(d$series_id),
          phi_per_30d = phi,
          half_life_days = if (is.finite(phi) && phi > 0 && phi < 1) 30 * log(0.5) / log(phi) else NA_real_,
          days_to_corr05 = if (is.finite(phi) && phi > 0 && phi < 1) 30 * log(0.05) / log(phi) else NA_real_,
          delta_aic_car1_minus_null = if (!inherits(fit0, "try-error") && !inherits(fit1, "try-error")) AIC(fit1) - AIC(fit0) else NA_real_,
          fit_ok = !inherits(fit1, "try-error")
        )
        idx <- idx + 1L
      }
    }
  }

  rbindlist(rows, fill = TRUE)
}

script_dir <- get_script_dir()
trace_mcmc_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
docs_dir <- normalizePath(file.path(trace_mcmc_dir, ".."), winslash = "/", mustWork = FALSE)
out_dir <- file.path(trace_mcmc_dir, "fticr_integration", "output")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

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

required_paths <- c(plot_data_path, mol_path, intensity_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required input files:\n", paste(missing_paths, collapse = "\n"))
}

message("Loading existing formula-count metric data...")
long_dt <- fread(plot_data_path)
long_dt[, date := as.IDate(date)]

base_metric_map <- data.table(
  metric = c("richness", "mean_nosc", "mean_hc", "mean_oc", "mean_ai", "frac_n_containing"),
  metric_label = c(
    "Formula richness",
    "Mean NOSC",
    "Mean H/C",
    "Mean O/C",
    "Aromaticity",
    "N-containing formulae (%)"
  ),
  metric_group = c(rep("Bulk molecular metrics", 5), "Elemental formula fractions"),
  metric_order = seq_len(6L)
)
base_dt <- merge(
  long_dt[metric %in% base_metric_map$metric][
    ,
    !c("metric_label", "metric_order"),
    with = FALSE
  ],
  base_metric_map,
  by = "metric",
  all.x = TRUE
)

message("Classifying formulas into van Krevelen classes...")
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
for (nm in element_cols) mol[, (nm) := safe_numeric(get(nm))]
mol <- mol[!is.na(C) & !is.na(H) & C > 0]
for (nm in c("O", "N", "S", "P")) mol[is.na(get(nm)), (nm) := 0]
mol[, HC := H / C]
mol[, OC := O / C]
mol[, AI_mod_denominator := C - O - S - N - P]
mol[, AI_mod := (1 + C - O - S - 0.5 * (H + N + P)) / AI_mod_denominator]
mol[AI_mod_denominator <= 0 | !is.finite(AI_mod) | AI_mod < 0, AI_mod := 0]
mol[, vk_class := classify_vk(HC, OC, AI_mod)]
mol <- unique(mol, by = formula_col)
setkeyv(mol, formula_col)

message("Loading FT-ICR intensity matrix for sample formula detection...")
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
formula_class <- mol[, .(formula = get(formula_col), vk_class)]
formula_class <- formula_class[formula %in% valid_formula_set]
setkey(formula_class, formula)

sample_map <- data.table(sample_col = field_cols)
sample_map[, sample_num := suppressWarnings(as.integer(sub("^Sihi_60398_([0-9]+)_r[0-9]+_.*$", "\\1", sample_col)))]
sample_map <- sample_map[!is.na(sample_num)]

vk_levels <- c("aromatic", "carbohydrate", "condensed_aromatic", "lignin_like", "lipid", "other")
vk_metric_map <- data.table(
  vk_class = vk_levels,
  metric = paste0("vk_frac_", vk_levels),
  metric_label = c(
    "Aromatic formulae (%)",
    "Carbohydrate formulae (%)",
    "Condensed aromatic formulae (%)",
    "Lignin-like formulae (%)",
    "Lipid formulae (%)",
    "Other formulae (%)"
  ),
  metric_group = "Van Krevelen class fractions",
  metric_order = seq.int(max(base_metric_map$metric_order) + 1L, max(base_metric_map$metric_order) + length(vk_levels))
)

message("Computing formula-count van Krevelen class fractions by sample...")
target_samples <- unique(long_dt[, .(sample_num, plot, date, treatment, treatment_label, depth_cm, depth_label)])
vk_rows <- vector("list", nrow(target_samples))
for (i in seq_len(nrow(target_samples))) {
  snum <- target_samples$sample_num[[i]]
  cols <- sample_map[sample_num == snum, sample_col]
  if (length(cols) == 0L) next
  detected <- rowSums(!is.na(intensity[, ..cols])) > 0
  detected_formulas <- intersect(intensity[detected, formula], valid_formula_set)
  if (length(detected_formulas) < 5L) next
  classes <- formula_class[J(detected_formulas), on = "formula", nomatch = 0, vk_class]
  class_counts <- table(factor(classes, levels = vk_levels))
  total <- sum(class_counts)
  if (!is.finite(total) || total == 0) next
  vk_rows[[i]] <- cbind(
    target_samples[i],
    data.table(
      vk_class = vk_levels,
      value = as.numeric(class_counts) / total,
      multiplier = 100
    )
  )
}
vk_dt <- rbindlist(vk_rows, fill = TRUE)
vk_dt <- merge(vk_dt, vk_metric_map, by = "vk_class", all.x = TRUE)
vk_dt[, value_plot := value * multiplier]
vk_dt[, vk_class := NULL]

metric_dt <- rbindlist(list(base_dt, vk_dt), fill = TRUE)
metric_dt <- metric_dt[is.finite(value_plot)]
metric_dt[, date := as.IDate(date)]
metric_order_dt <- unique(metric_dt[, .(metric, metric_label, metric_group, metric_order)])
setorder(metric_order_dt, metric_order)

message("Fitting treatment-depth CAR(1) temporal autocorrelation models...")
autocorr_dt <- fit_scalar_car1_by_group(metric_dt, metric_order_dt)
autocorr_dt[, metric_clean := factor(clean_metric_label(metric_label), levels = rev(metric_order_dt$metric_label))]
autocorr_dt[, depth_label := factor(depth_label, levels = c("0-10 cm", "10-30 cm"))]
autocorr_dt[, treatment_label := factor(treatment_label, levels = c("Control", "Warmed"))]

metric_data_out <- file.path(out_dir, "fticr_q2_autocorrelation_metric_plot_data_vk_classes.csv")
autocorr_out <- file.path(out_dir, "fticr_q2_autocorrelation_terms_by_treatment_depth_vk_classes.csv")
fwrite(metric_dt, metric_data_out)
fwrite(autocorr_dt, autocorr_out)
message("Wrote: ", metric_data_out)
message("Wrote: ", autocorr_out)

treatment_cols <- c("Control" = "#2B6CB0", "Warmed" = "#C45A2A")
depth_shapes <- c("0-10 cm" = 16, "10-30 cm" = 17)
depth_display_labels <- c("0-10 cm" = "10cm", "10-30 cm" = "30cm")
depth_offsets <- c("0-10 cm" = 0.14, "10-30 cm" = -0.14)

axis_dt <- unique(autocorr_dt[, .(metric_clean, metric_order)])
axis_dt[, base_y := max(metric_order, na.rm = TRUE) - metric_order + 1]
plot_dt <- merge(autocorr_dt, axis_dt, by = "metric_clean", all.x = TRUE)
plot_dt[, y_position := base_y + unname(depth_offsets[as.character(depth_label)])]
segment_dt <- dcast(
  plot_dt,
  depth_label + metric_clean + y_position ~ treatment_label,
  value.var = "phi_per_30d"
)

p_autocorr <- ggplot(plot_dt, aes(x = phi_per_30d, y = y_position, color = treatment_label, shape = depth_label)) +
  geom_hline(
    data = axis_dt,
    aes(yintercept = base_y),
    color = "#EEF0F3",
    linewidth = 0.4
  ) +
  geom_vline(xintercept = 0, color = "#6B7280", linewidth = 0.35) +
  geom_segment(
    data = segment_dt[is.finite(Control) & is.finite(Warmed)],
    aes(x = Control, xend = Warmed, y = y_position, yend = y_position),
    inherit.aes = FALSE,
    color = "#B8BDC7",
    linewidth = 0.45
  ) +
  geom_point(size = 2.8, stroke = 1.0, na.rm = TRUE) +
  scale_color_manual(values = treatment_cols) +
  scale_shape_manual(values = depth_shapes, breaks = names(depth_shapes), labels = depth_display_labels) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  scale_y_continuous(
    breaks = axis_dt[order(base_y), base_y],
    labels = as.character(axis_dt[order(base_y), metric_clean]),
    expand = expansion(add = 0.58)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(shape = 16)),
    shape = guide_legend(order = 2, override.aes = list(color = "#4B5563"))
  ) +
  labs(
    x = "30-day temporal autocorrelation",
    y = NULL
  ) +
  theme_trace(base_size = 10.5)

save_plot(
  p_autocorr,
  file.path(fig_dir, "fticr_q2_autocorrelation_terms_vk_classes_by_treatment_depth.png"),
  file.path(fig_dir, "fticr_q2_autocorrelation_terms_vk_classes_by_treatment_depth.pdf"),
  width = 11.0,
  height = 7.8
)
