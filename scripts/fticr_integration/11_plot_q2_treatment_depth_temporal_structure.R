#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
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

transform_metric <- function(metric, value_plot) {
  if (metric == "richness") return(log1p(value_plot))
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

fit_scalar_car1_by_group <- function(long_dt) {
  ctrl <- lmeControl(opt = "optim", msMaxIter = 200, maxIter = 100, returnObject = TRUE)
  rows <- list()
  idx <- 1L

  for (metric_name in sort(unique(long_dt$metric))) {
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

clean_metric_label <- function(x) {
  gsub("\\n\\([^)]*\\)", "", gsub("\\n", " ", x))
}

theme_trace <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = "#202124", size = rel(1.14)),
      plot.subtitle = element_text(color = "#4B5563", size = rel(0.92), margin = margin(b = 8)),
      plot.caption = element_text(color = "#5F6368", size = rel(0.78), hjust = 0, margin = margin(t = 8)),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.text = element_text(face = "bold", color = "#202124"),
      axis.text.y = element_text(color = "#202124"),
      axis.title = element_text(color = "#202124")
    )
}

save_plot <- function(plot, png_path, pdf_path, width, height) {
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 320, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, units = "in", bg = "white")
  message("Wrote: ", png_path)
  message("Wrote: ", pdf_path)
}

script_dir <- get_script_dir()
trace_mcmc_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
out_dir <- file.path(trace_mcmc_dir, "fticr_integration", "output")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pair_test_path <- file.path(out_dir, "fticr_q2_pairwise_temporal_distance_tests.csv")
scalar_slope_path <- file.path(out_dir, "fticr_q2_warming_scalar_semivariance_slope_tests.csv")
random_time_path <- file.path(out_dir, "fticr_q2_effective_randomization_time.csv")
long_metric_path <- file.path(out_dir, "fticr_formula_count_temporal_metrics_plot_data.csv")

required_paths <- c(pair_test_path, scalar_slope_path, random_time_path, long_metric_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required input files:\n", paste(missing_paths, collapse = "\n"))
}

message("Loading Q2 results...")
pair_tests <- fread(pair_test_path)
scalar_slopes <- fread(scalar_slope_path)
random_times <- fread(random_time_path)
long_dt <- fread(long_metric_path)

distance_names <- c(
  formula_jaccard_distance = "Formula identity\nJaccard",
  vk_class_braycurtis = "Chemical class\nBray-Curtis",
  metric_vector_distance = "Scalar metric\nvector"
)
distance_levels <- c("Formula identity\nJaccard", "Chemical class\nBray-Curtis", "Scalar metric\nvector")
treatment_cols <- c("Control" = "#2B6CB0", "Warmed" = "#C45A2A")
treatment_shapes <- c("Control" = 16, "Warmed" = 1)
depth_levels <- c("0-10 cm", "10-30 cm")
depth_display_labels <- c("0-10 cm" = "10cm", "10-30 cm" = "30cm")

message("Preparing distance-slope data...")
distance_slope <- pair_tests[
  scope == "by_treatment_depth",
  .(
    distance_metric,
    treatment_label,
    depth_label,
    n_pairs,
    slope_per_lag_doubling,
    permutation_p_positive_slope
  )
]
distance_slope[, distance_label := factor(distance_names[distance_metric], levels = rev(distance_levels))]
distance_slope[, depth_label := factor(depth_label, levels = depth_levels)]
distance_slope[, treatment_label := factor(treatment_label, levels = c("Control", "Warmed"))]

message("Preparing scalar-slope data...")
scalar_slope_long <- melt(
  scalar_slopes,
  id.vars = c("metric", "metric_label", "metric_order", "depth_label", "n_control_pairs", "n_warmed_pairs"),
  measure.vars = c("control_slope", "warmed_slope"),
  variable.name = "treatment_raw",
  value.name = "slope_per_lag_doubling"
)
scalar_slope_long[, treatment_label := fifelse(treatment_raw == "control_slope", "Control", "Warmed")]
scalar_slope_long[, n_pairs := fifelse(treatment_label == "Control", n_control_pairs, n_warmed_pairs)]
scalar_slope_long[, metric_clean := factor(clean_metric_label(metric_label), levels = rev(clean_metric_label(unique(scalar_slopes[order(metric_order), metric_label]))))]
scalar_slope_long[, depth_label := factor(depth_label, levels = depth_levels)]
scalar_slope_long[, treatment_label := factor(treatment_label, levels = c("Control", "Warmed"))]

message("Preparing whole-sample randomization-time data...")
distance_random <- random_times[scope == "by_treatment_depth"]
distance_random[, distance_label := factor(distance_names[distance_metric], levels = rev(distance_levels))]
distance_random[, depth_label := factor(depth_label, levels = depth_levels)]
distance_random[, treatment_label := factor(treatment_label, levels = c("Control", "Warmed"))]
distance_random[, days_to_median_random := 2^((random_median_distance - intercept) / slope_per_lag_doubling)]
distance_random[!is.finite(days_to_median_random) | days_to_median_random <= 0, days_to_median_random := NA_real_]

message("Preparing common-baseline randomization-time data...")
common_depth_baseline <- random_times[
  scope == "by_depth",
  .(
    distance_metric,
    depth_label,
    common_random_median_distance = random_median_distance,
    common_random_mean_distance = random_mean_distance
  )
]
distance_common_random <- merge(
  distance_random,
  common_depth_baseline,
  by = c("distance_metric", "depth_label"),
  all.x = TRUE
)
distance_common_random[, days_to_common_median_random := 2^((common_random_median_distance - intercept) / slope_per_lag_doubling)]
distance_common_random[
  !is.finite(days_to_common_median_random) | days_to_common_median_random <= 0,
  days_to_common_median_random := NA_real_
]

message("Estimating scalar CAR(1) decorrelation times by treatment and depth...")
scalar_random <- fit_scalar_car1_by_group(long_dt)
scalar_random[, metric_clean := factor(clean_metric_label(metric_label), levels = rev(clean_metric_label(unique(scalar_slopes[order(metric_order), metric_label]))))]
scalar_random[, depth_label := factor(depth_label, levels = depth_levels)]
scalar_random[, treatment_label := factor(treatment_label, levels = c("Control", "Warmed"))]

distance_slope_out <- file.path(out_dir, "fticr_q2_plot_distance_slope_data.csv")
scalar_slope_out <- file.path(out_dir, "fticr_q2_plot_scalar_slope_data.csv")
distance_random_out <- file.path(out_dir, "fticr_q2_plot_distance_randomization_time_data.csv")
distance_common_random_out <- file.path(out_dir, "fticr_q2_plot_distance_common_baseline_randomization_time_data.csv")
scalar_random_out <- file.path(out_dir, "fticr_q2_scalar_car1_treatment_depth_decorrelation.csv")
fwrite(distance_slope, distance_slope_out)
fwrite(scalar_slope_long, scalar_slope_out)
fwrite(distance_random, distance_random_out)
fwrite(distance_common_random, distance_common_random_out)
fwrite(scalar_random, scalar_random_out)
message("Wrote: ", distance_slope_out)
message("Wrote: ", scalar_slope_out)
message("Wrote: ", distance_random_out)
message("Wrote: ", distance_common_random_out)
message("Wrote: ", scalar_random_out)

distance_slope_segments <- dcast(
  distance_slope,
  depth_label + distance_label ~ treatment_label,
  value.var = "slope_per_lag_doubling"
)
scalar_slope_segments <- dcast(
  scalar_slope_long,
  depth_label + metric_clean ~ treatment_label,
  value.var = "slope_per_lag_doubling"
)
distance_random_segments <- dcast(
  distance_random[is.finite(days_to_median_random)],
  depth_label + distance_label ~ treatment_label,
  value.var = "days_to_median_random"
)
distance_common_random_segments <- dcast(
  distance_common_random[is.finite(days_to_common_median_random)],
  depth_label + distance_label ~ treatment_label,
  value.var = "days_to_common_median_random"
)
scalar_random_segments <- dcast(
  scalar_random[is.finite(days_to_corr05)],
  depth_label + metric_clean ~ treatment_label,
  value.var = "days_to_corr05"
)

depth_offsets <- c("0-10 cm" = 0.14, "10-30 cm" = -0.14)
depth_shapes <- c("0-10 cm" = 16, "10-30 cm" = 17)

distance_axis <- data.table(
  distance_label = factor(distance_levels, levels = rev(distance_levels)),
  base_y = rev(seq_along(distance_levels))
)
distance_slope_plot <- merge(distance_slope, distance_axis, by = "distance_label", all.x = TRUE)
distance_slope_plot[, y_position := base_y + unname(depth_offsets[as.character(depth_label)])]
distance_slope_offset_segments <- dcast(
  distance_slope_plot,
  depth_label + distance_label + y_position ~ treatment_label,
  value.var = "slope_per_lag_doubling"
)

metric_axis <- unique(scalar_slope_long[, .(metric_clean, metric_order)])
metric_axis[, base_y := max(metric_order, na.rm = TRUE) - metric_order + 1]
scalar_slope_plot <- merge(scalar_slope_long, metric_axis, by = "metric_clean", all.x = TRUE)
scalar_slope_plot[, y_position := base_y + unname(depth_offsets[as.character(depth_label)])]
scalar_slope_offset_segments <- dcast(
  scalar_slope_plot,
  depth_label + metric_clean + y_position ~ treatment_label,
  value.var = "slope_per_lag_doubling"
)

p_distance_slope <- ggplot(distance_slope_plot, aes(x = slope_per_lag_doubling, y = y_position, color = treatment_label, shape = depth_label)) +
  geom_hline(
    data = distance_axis,
    aes(yintercept = base_y),
    color = "#EEF0F3",
    linewidth = 0.4
  ) +
  geom_vline(xintercept = 0, color = "#6B7280", linewidth = 0.35) +
  geom_segment(
    data = distance_slope_offset_segments[is.finite(Control) & is.finite(Warmed)],
    aes(x = Control, xend = Warmed, y = y_position, yend = y_position),
    inherit.aes = FALSE,
    color = "#B8BDC7",
    linewidth = 0.45
  ) +
  geom_point(size = 3.0, stroke = 1.0) +
  scale_color_manual(values = treatment_cols) +
  scale_shape_manual(values = depth_shapes, breaks = depth_levels, labels = depth_display_labels) +
  scale_x_continuous(
    labels = label_number(accuracy = 0.01),
    limits = c(-0.01, 0.18),
    breaks = seq(0, 0.18, by = 0.03),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_continuous(
    breaks = distance_axis[order(base_y), base_y],
    labels = as.character(distance_axis[order(base_y), distance_label]),
    expand = expansion(add = 0.55)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(shape = 16)),
    shape = guide_legend(order = 2, override.aes = list(color = "#4B5563"))
  ) +
  labs(
    x = "Distance slope per log2(day interval)",
    y = NULL
  ) +
  theme_trace()

p_scalar_slope <- ggplot(scalar_slope_plot, aes(x = slope_per_lag_doubling, y = y_position, color = treatment_label, shape = depth_label)) +
  geom_hline(
    data = metric_axis,
    aes(yintercept = base_y),
    color = "#EEF0F3",
    linewidth = 0.4
  ) +
  geom_vline(xintercept = 0, color = "#6B7280", linewidth = 0.35) +
  geom_segment(
    data = scalar_slope_offset_segments[is.finite(Control) & is.finite(Warmed)],
    aes(x = Control, xend = Warmed, y = y_position, yend = y_position),
    inherit.aes = FALSE,
    color = "#B8BDC7",
    linewidth = 0.45
  ) +
  geom_point(size = 2.7, stroke = 1.0) +
  scale_color_manual(values = treatment_cols) +
  scale_shape_manual(values = depth_shapes, breaks = depth_levels, labels = depth_display_labels) +
  scale_x_continuous(
    labels = label_number(accuracy = 0.01),
    limits = c(-0.05, 0.40),
    breaks = seq(0, 0.40, by = 0.10),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_continuous(
    breaks = metric_axis[order(base_y), base_y],
    labels = as.character(metric_axis[order(base_y), metric_clean]),
    expand = expansion(add = 0.55)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(shape = 16)),
    shape = guide_legend(order = 2, override.aes = list(color = "#4B5563"))
  ) +
  labs(
    x = "Scaled semivariance per log2(day interval)",
    y = NULL
  ) +
  theme_trace()

slope_legend <- get_legend(p_distance_slope + theme(legend.position = "top"))
distance_slope_panel_margin <- margin(5.5, 12, 5.5, 28)
scalar_slope_panel_margin <- margin(5.5, 5.5, 5.5, 56)
p_slope_two_panel <- plot_grid(
  slope_legend,
  plot_grid(
    p_distance_slope + theme(legend.position = "none", plot.margin = distance_slope_panel_margin),
    p_scalar_slope + theme(legend.position = "none", plot.margin = scalar_slope_panel_margin),
    ncol = 2,
    align = "h",
    axis = "tb",
    labels = c("A", "B"),
    label_size = 12,
    label_fontface = "bold",
    rel_widths = c(0.92, 1.22)
  ),
  ncol = 1,
  rel_heights = c(0.08, 1)
)

p_distance_random <- ggplot(distance_random, aes(x = days_to_median_random, y = distance_label, color = treatment_label, shape = treatment_label)) +
  geom_segment(
    data = distance_random_segments[is.finite(Control) & is.finite(Warmed)],
    aes(x = Control, xend = Warmed, y = distance_label, yend = distance_label),
    inherit.aes = FALSE,
    color = "#B8BDC7",
    linewidth = 0.45
  ) +
  geom_point(size = 3.0, stroke = 1.0, na.rm = TRUE) +
  facet_wrap(~depth_label, nrow = 1, labeller = labeller(depth_label = as_labeller(depth_display_labels))) +
  scale_color_manual(values = treatment_cols) +
  scale_shape_manual(values = treatment_shapes) +
  scale_x_log10(
    breaks = c(14, 30, 90, 180, 365, 730, 1095),
    labels = c("14 d", "30 d", "90 d", "180 d", "1 yr", "2 yr", "3 yr"),
    expand = expansion(mult = c(0.08, 0.14))
  ) +
  labs(
    x = "Days to median random-pair distance",
    y = NULL
  ) +
  theme_trace()

p_distance_common_random <- ggplot(distance_common_random, aes(x = days_to_common_median_random, y = distance_label, color = treatment_label, shape = treatment_label)) +
  geom_segment(
    data = distance_common_random_segments[is.finite(Control) & is.finite(Warmed)],
    aes(x = Control, xend = Warmed, y = distance_label, yend = distance_label),
    inherit.aes = FALSE,
    color = "#B8BDC7",
    linewidth = 0.45
  ) +
  geom_point(size = 3.0, stroke = 1.0, na.rm = TRUE) +
  facet_wrap(~depth_label, nrow = 1, labeller = labeller(depth_label = as_labeller(depth_display_labels))) +
  scale_color_manual(values = treatment_cols) +
  scale_shape_manual(values = treatment_shapes) +
  scale_x_log10(
    breaks = c(14, 30, 90, 365, 730, 1825, 3650, 7300, 14600),
    labels = c("14 d", "30 d", "90 d", "1 yr", "2 yr", "5 yr", "10 yr", "20 yr", "40 yr"),
    limits = c(10, 20000),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = "Days to common median random-pair distance",
    y = NULL
  ) +
  theme_trace()

p_scalar_random <- ggplot(scalar_random, aes(x = days_to_corr05, y = metric_clean, color = treatment_label, shape = treatment_label)) +
  geom_segment(
    data = scalar_random_segments[is.finite(Control) & is.finite(Warmed)],
    aes(x = Control, xend = Warmed, y = metric_clean, yend = metric_clean),
    inherit.aes = FALSE,
    color = "#B8BDC7",
    linewidth = 0.45
  ) +
  geom_point(size = 2.7, stroke = 1.0, na.rm = TRUE) +
  facet_wrap(~depth_label, nrow = 1, labeller = labeller(depth_label = as_labeller(depth_display_labels))) +
  scale_color_manual(values = treatment_cols) +
  scale_shape_manual(values = treatment_shapes) +
  scale_x_log10(
    breaks = c(7, 14, 30, 90, 180, 365, 548),
    labels = c("7 d", "14 d", "30 d", "90 d", "180 d", "1 yr", "18 mo"),
    limits = c(5, 650),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = "Days to residual correlation < 0.05",
    y = NULL
  ) +
  theme_trace()

save_plot(
  p_distance_slope,
  file.path(fig_dir, "fticr_q2_distance_slopes_by_treatment_depth.png"),
  file.path(fig_dir, "fticr_q2_distance_slopes_by_treatment_depth.pdf"),
  width = 9.8,
  height = 4.8
)
save_plot(
  p_scalar_slope,
  file.path(fig_dir, "fticr_q2_scalar_semivariance_slopes_by_treatment_depth.png"),
  file.path(fig_dir, "fticr_q2_scalar_semivariance_slopes_by_treatment_depth.pdf"),
  width = 10.2,
  height = 6.0
)
save_plot(
  p_slope_two_panel,
  file.path(fig_dir, "fticr_q2_temporal_slopes_two_panel_by_treatment_depth.png"),
  file.path(fig_dir, "fticr_q2_temporal_slopes_two_panel_by_treatment_depth.pdf"),
  width = 14.0,
  height = 6.6
)
save_plot(
  p_distance_random,
  file.path(fig_dir, "fticr_q2_distance_randomization_time_by_treatment_depth.png"),
  file.path(fig_dir, "fticr_q2_distance_randomization_time_by_treatment_depth.pdf"),
  width = 9.8,
  height = 4.8
)
save_plot(
  p_distance_common_random,
  file.path(fig_dir, "fticr_q2_distance_common_baseline_randomization_time_by_treatment_depth.png"),
  file.path(fig_dir, "fticr_q2_distance_common_baseline_randomization_time_by_treatment_depth.pdf"),
  width = 10.6,
  height = 4.8
)
save_plot(
  p_scalar_random,
  file.path(fig_dir, "fticr_q2_scalar_decorrelation_time_by_treatment_depth.png"),
  file.path(fig_dir, "fticr_q2_scalar_decorrelation_time_by_treatment_depth.pdf"),
  width = 10.2,
  height = 6.0
)
