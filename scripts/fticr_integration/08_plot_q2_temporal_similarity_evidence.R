#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)))
  }
  getwd()
}

script_dir <- get_script_dir()
trace_mcmc_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
out_dir <- file.path(trace_mcmc_dir, "fticr_integration", "output")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ctmm_path <- file.path(out_dir, "fticr_q2_scalar_continuous_time_autocorrelation.csv")
pair_path <- file.path(out_dir, "fticr_q2_pairwise_temporal_distance_tests.csv")
missing_paths <- c(ctmm_path, pair_path)[!file.exists(c(ctmm_path, pair_path))]
if (length(missing_paths) > 0L) {
  stop("Missing required Q2 output files. Run 07_q2_temporal_similarity_tests.R first:\n", paste(missing_paths, collapse = "\n"))
}

ctmm <- fread(ctmm_path)
pair <- fread(pair_path)

metric_names <- c(
  richness = "Formula\nrichness",
  mean_nosc = "Mean\nNOSC",
  mean_hc = "Mean\nH/C",
  mean_oc = "Mean\nO/C",
  mean_ai = "Aromaticity",
  frac_lignin = "Lignin-like\nformulae",
  frac_n_containing = "N-containing\nformulae"
)
ctmm[, metric_short := metric_names[metric]]
ctmm[, metric_short := factor(metric_short, levels = rev(metric_names[ctmm[order(metric_order), metric]]))]
ctmm[, temporal_support := fifelse(delta_aic_car1_minus_null < -2 & lrt_p < 0.05, "Supported", "Weak / absent")]
ctmm[, label := fifelse(
  is.finite(half_life_days) & delta_aic_car1_minus_null < -2,
  paste0("half-life ", round(half_life_days), " d"),
  ""
)]

distance_names <- c(
  formula_jaccard_distance = "Formula identity\nJaccard distance",
  vk_class_braycurtis = "Chemical class\nBray-Curtis",
  metric_vector_distance = "Scalar metric\nvector distance"
)
pair_plot <- pair[scope == "by_depth"]
pair_plot[, distance_label := factor(distance_names[distance_metric], levels = distance_names)]
pair_plot[, depth_label := factor(depth_label, levels = c("0-10 cm", "10-30 cm"))]
pair_plot[, support := fifelse(permutation_p_positive_slope < 0.05, "Permutation p < 0.05", "Weak")]

support_cols <- c("Supported" = "#2B6CB0", "Weak / absent" = "#B8BDC7")
depth_cols <- c("0-10 cm" = "#2B6CB0", "10-30 cm" = "#C45A2A")

p_ctmm <- ggplot(ctmm, aes(x = delta_aic_car1_minus_null, y = metric_short, fill = temporal_support)) +
  geom_vline(xintercept = 0, color = "#6B7280", linewidth = 0.35) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  geom_text(
    aes(label = label),
    hjust = 1.03,
    color = "white",
    size = 3.0,
    fontface = "bold"
  ) +
  scale_fill_manual(values = support_cols, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.12))) +
  labs(
    title = "A. Scalar features: continuous-time autocorrelation",
    subtitle = "Negative Delta AIC favors CAR(1) residual correlation within plot-depth series",
    x = "Delta AIC: CAR(1) model - no temporal correlation",
    y = NULL
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", color = "#202124"),
    plot.subtitle = element_text(color = "#4B5563", size = 9.5),
    axis.text.y = element_text(color = "#202124")
  )

p_pair <- ggplot(pair_plot, aes(x = slope_per_lag_doubling, y = distance_label, color = depth_label, shape = depth_label)) +
  geom_vline(xintercept = 0, color = "#6B7280", linewidth = 0.35) +
  geom_point(size = 3.0, stroke = 0.9) +
  scale_color_manual(values = depth_cols, name = NULL) +
  scale_shape_manual(values = c("0-10 cm" = 16, "10-30 cm" = 1), name = NULL) +
  scale_x_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.08, 0.12))) +
  labs(
    title = "B. Pairwise dissimilarity vs. interval",
    subtitle = "Slope per doubling of days; all shown p < 0.005 by restricted permutation",
    x = "Slope per doubling of resampling interval",
    y = NULL
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", color = "#202124"),
    plot.subtitle = element_text(color = "#4B5563", size = 9.5),
    axis.text.y = element_text(color = "#202124")
  )

combined <- plot_grid(p_ctmm, p_pair, nrow = 1, rel_widths = c(1.08, 1), align = "h", axis = "tb")
title <- ggdraw() +
  draw_label(
    "Evidence that TRACE FT-ICR molecular composition changes with resampling interval",
    x = 0,
    hjust = 0,
    fontface = "bold",
    size = 15,
    color = "#202124"
  ) +
  draw_label(
    "Formula-count metrics and detected-formula composition; n = 195 spectra from 6 plots and 2 depths",
    x = 0,
    y = 0.18,
    hjust = 0,
    size = 10.2,
    color = "#4B5563"
  )
final_plot <- plot_grid(title, combined, ncol = 1, rel_heights = c(0.13, 1))

png_out <- file.path(fig_dir, "fticr_q2_temporal_similarity_evidence.png")
pdf_out <- file.path(fig_dir, "fticr_q2_temporal_similarity_evidence.pdf")
ggsave(png_out, final_plot, width = 13.2, height = 5.8, units = "in", dpi = 320, bg = "white")
ggsave(pdf_out, final_plot, width = 13.2, height = 5.8, units = "in", bg = "white")

message("Wrote: ", png_out)
message("Wrote: ", pdf_out)
