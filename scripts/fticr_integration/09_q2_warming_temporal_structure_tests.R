#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)))
  }
  getwd()
}

plot_pair_id <- function(plot) {
  fifelse(plot %in% c(1L, 2L), "pair_1_2",
    fifelse(plot %in% c(3L, 4L), "pair_3_4",
      fifelse(plot %in% c(5L, 6L), "pair_5_6", NA_character_)
    )
  )
}

calc_slope <- function(dt, response_col) {
  mod_dt <- dt[is.finite(get(response_col)) & is.finite(log2_lag)]
  if (nrow(mod_dt) < 10L || length(unique(mod_dt$log2_lag)) < 2L) return(NA_real_)
  unname(coef(lm(mod_dt[[response_col]] ~ mod_dt$log2_lag))[[2]])
}

calc_treatment_slopes <- function(dt, response_col, group_cols) {
  groups <- unique(dt[, ..group_cols])
  out <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    sub <- dt[g, on = group_cols, nomatch = 0]
    ctrl <- sub[treatment_label == "Control"]
    warm <- sub[treatment_label == "Warmed"]
    out[[i]] <- cbind(
      g,
      data.table(
        response = response_col,
        n_control_pairs = nrow(ctrl),
        n_warmed_pairs = nrow(warm),
        control_slope = calc_slope(ctrl, response_col),
        warmed_slope = calc_slope(warm, response_col)
      )
    )
  }
  ans <- rbindlist(out, fill = TRUE)
  ans[, warmed_minus_control_slope := warmed_slope - control_slope]
  ans
}

exact_plot_pair_permutation <- function(dt, response_col, group_cols) {
  plot_pairs <- sort(na.omit(unique(dt$plot_pair)))
  swap_grid <- as.data.table(expand.grid(rep(list(c(FALSE, TRUE)), length(plot_pairs))))
  setnames(swap_grid, plot_pairs)

  observed <- calc_treatment_slopes(dt, response_col, group_cols)
  observed[, permutation_p_two_sided := NA_real_]
  observed[, permutation_p_warmed_greater := NA_real_]

  groups <- unique(dt[, ..group_cols])
  for (g_i in seq_len(nrow(groups))) {
    g <- groups[g_i]
    sub <- dt[g, on = group_cols, nomatch = 0]
    obs_delta <- observed[g, on = group_cols, warmed_minus_control_slope]
    if (!is.finite(obs_delta)) next

    null_delta <- rep(NA_real_, nrow(swap_grid))
    for (b in seq_len(nrow(swap_grid))) {
      perm_sub <- copy(sub)
      for (pp in plot_pairs) {
        if (isTRUE(swap_grid[[pp]][[b]])) {
          perm_sub[plot_pair == pp & treatment_label == "Control", treatment_perm := "Warmed"]
          perm_sub[plot_pair == pp & treatment_label == "Warmed", treatment_perm := "Control"]
        } else {
          perm_sub[plot_pair == pp, treatment_perm := treatment_label]
        }
      }
      perm_sub[, treatment_original := treatment_label]
      perm_sub[, treatment_label := treatment_perm]
      null_delta[[b]] <- calc_treatment_slopes(perm_sub, response_col, group_cols)$warmed_minus_control_slope
    }

    observed[g, on = group_cols, permutation_p_two_sided := mean(abs(null_delta) >= abs(obs_delta), na.rm = TRUE)]
    observed[g, on = group_cols, permutation_p_warmed_greater := mean(null_delta >= obs_delta, na.rm = TRUE)]
  }

  observed
}

script_dir <- get_script_dir()
trace_mcmc_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
out_dir <- file.path(trace_mcmc_dir, "fticr_integration", "output")

whole_pair_path <- file.path(out_dir, "fticr_q2_pairwise_temporal_distance_pairs.csv")
scalar_pair_path <- file.path(out_dir, "fticr_formula_count_autocorrelation_pairs.csv")
required_paths <- c(whole_pair_path, scalar_pair_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required input files. Run 06 and 07 first:\n", paste(missing_paths, collapse = "\n"))
}

message("Loading pairwise temporal-distance data...")
whole_pairs <- fread(whole_pair_path)
whole_pairs[, plot_pair := plot_pair_id(as.integer(plot))]
whole_pairs <- whole_pairs[!is.na(plot_pair) & lag_days > 0]

scalar_pairs <- fread(scalar_pair_path)
scalar_pairs[, plot_pair := plot_pair_id(as.integer(plot))]
scalar_pairs <- scalar_pairs[!is.na(plot_pair) & lag_days > 0]
scalar_pairs[, log2_lag := log2(lag_days)]

distance_cols <- c("formula_jaccard_distance", "vk_class_braycurtis", "metric_vector_distance")

message("Testing warming differences in whole-sample temporal slopes...")
whole_results <- rbindlist(lapply(distance_cols, function(resp) {
  by_depth <- exact_plot_pair_permutation(whole_pairs, resp, c("depth_label"))
  by_depth[, distance_metric := resp]
  by_depth
}), fill = TRUE)
setcolorder(whole_results, c("distance_metric", setdiff(names(whole_results), "distance_metric")))
setorder(whole_results, distance_metric, depth_label)

message("Testing warming differences in scalar-feature semivariance slopes...")
scalar_results <- rbindlist(lapply(sort(unique(scalar_pairs$metric)), function(metric_name) {
  sub <- scalar_pairs[metric == metric_name]
  by_depth <- exact_plot_pair_permutation(sub, "scaled_semivar", c("metric", "metric_label", "metric_order", "depth_label"))
  by_depth
}), fill = TRUE)
setorder(scalar_results, metric_order, depth_label)

whole_out <- file.path(out_dir, "fticr_q2_warming_pairwise_temporal_slope_tests.csv")
scalar_out <- file.path(out_dir, "fticr_q2_warming_scalar_semivariance_slope_tests.csv")
summary_out <- file.path(out_dir, "fticr_q2_warming_temporal_structure_summary.txt")

fwrite(whole_results, whole_out)
fwrite(scalar_results, scalar_out)

whole_print <- whole_results[
  ,
  .(
    distance_metric,
    depth_label,
    n_control_pairs,
    n_warmed_pairs,
    control_slope,
    warmed_slope,
    warmed_minus_control_slope,
    permutation_p_two_sided,
    permutation_p_warmed_greater
  )
]

scalar_print <- scalar_results[
  ,
  .(
    metric_label,
    depth_label,
    n_control_pairs,
    n_warmed_pairs,
    control_slope,
    warmed_slope,
    warmed_minus_control_slope,
    permutation_p_two_sided,
    permutation_p_warmed_greater
  )
]

summary_lines <- c(
  "TRACE FTICR Q2 warming-dependence of temporal structure",
  "=======================================================",
  "",
  "Question:",
  "Does temporal structure differ under warming?",
  "",
  "Method:",
  "- Response is the slope of dissimilarity or scaled semivariance versus log2(days between samples).",
  "- A positive slope means samples become less similar as resampling interval increases.",
  "- warmed_minus_control_slope > 0 means temporal dissimilarity accumulates faster in warmed plots.",
  "- Treatment comparison uses exact restricted permutations that swap treatment labels only within the three plot pairs: 1/2, 3/4, and 5/6.",
  "- With three plot pairs, the smallest possible two-sided p-value is 0.25 and the smallest one-sided p-value is 0.125; effect sizes matter more than nominal significance here.",
  "",
  "Whole-sample pairwise dissimilarity slopes:",
  capture.output(print(whole_print)),
  "",
  "Scalar-feature scaled semivariance slopes:",
  capture.output(print(scalar_print)),
  "",
  "Outputs:",
  paste("-", whole_out),
  paste("-", scalar_out),
  "",
  "Caveats:",
  "- Pairwise distances are not independent; the exact plot-pair permutation protects the treatment comparison but has low power with n = 3 treatment pairs.",
  "- This tests whether the interval-similarity relationship differs by treatment, not whether warmed and control plots differ in average composition."
)
writeLines(summary_lines, summary_out)

message("Wrote: ", whole_out)
message("Wrote: ", scalar_out)
message("Wrote: ", summary_out)
