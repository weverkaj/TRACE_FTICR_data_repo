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

estimate_crossing <- function(dt, distance_col) {
  mod_dt <- dt[is.finite(get(distance_col)) & is.finite(log2_lag) & lag_days > 0]
  if (nrow(mod_dt) < 10L || length(unique(mod_dt$log2_lag)) < 2L) {
    return(data.table(
      n_pairs = nrow(mod_dt),
      n_series = uniqueN(mod_dt$series_id),
      intercept = NA_real_,
      slope_per_lag_doubling = NA_real_,
      random_mean_distance = NA_real_,
      random_median_distance = NA_real_,
      random_p25_distance = NA_real_,
      random_p75_distance = NA_real_,
      effective_randomization_days = NA_real_,
      fit_note = "too_few_pairs"
    ))
  }

  fit <- lm(mod_dt[[distance_col]] ~ mod_dt$log2_lag)
  intercept <- unname(coef(fit)[[1]])
  slope <- unname(coef(fit)[[2]])
  random_mean <- mean(mod_dt[[distance_col]], na.rm = TRUE)
  random_median <- median(mod_dt[[distance_col]], na.rm = TRUE)
  random_p25 <- as.numeric(quantile(mod_dt[[distance_col]], 0.25, na.rm = TRUE))
  random_p75 <- as.numeric(quantile(mod_dt[[distance_col]], 0.75, na.rm = TRUE))
  crossing <- if (is.finite(slope) && slope > 0) 2^((random_mean - intercept) / slope) else NA_real_

  data.table(
    n_pairs = nrow(mod_dt),
    n_series = uniqueN(mod_dt$series_id),
    intercept = intercept,
    slope_per_lag_doubling = slope,
    random_mean_distance = random_mean,
    random_median_distance = random_median,
    random_p25_distance = random_p25,
    random_p75_distance = random_p75,
    effective_randomization_days = crossing,
    fit_note = "ok"
  )
}

bootstrap_crossing <- function(dt, distance_col, n_boot = 2000L, seed = 20260707L) {
  set.seed(seed)
  series <- sort(unique(dt$series_id))
  if (length(series) < 2L) {
    return(data.table(
      crossing_lwr_days = NA_real_,
      crossing_med_days = NA_real_,
      crossing_upr_days = NA_real_,
      random_mean_lwr = NA_real_,
      random_mean_upr = NA_real_
    ))
  }

  boot_cross <- rep(NA_real_, n_boot)
  boot_random <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    picked <- sample(series, length(series), replace = TRUE)
    pieces <- vector("list", length(picked))
    for (i in seq_along(picked)) {
      x <- copy(dt[series_id == picked[[i]]])
      x[, boot_series := paste0(series_id, "_boot_", i)]
      pieces[[i]] <- x
    }
    boot_dt <- rbindlist(pieces, fill = TRUE)
    boot_dt[, series_id := boot_series]
    est <- estimate_crossing(boot_dt, distance_col)
    boot_cross[[b]] <- est$effective_randomization_days
    boot_random[[b]] <- est$random_mean_distance
  }

  data.table(
    crossing_lwr_days = as.numeric(quantile(boot_cross, 0.025, na.rm = TRUE)),
    crossing_med_days = as.numeric(quantile(boot_cross, 0.5, na.rm = TRUE)),
    crossing_upr_days = as.numeric(quantile(boot_cross, 0.975, na.rm = TRUE)),
    random_mean_lwr = as.numeric(quantile(boot_random, 0.025, na.rm = TRUE)),
    random_mean_upr = as.numeric(quantile(boot_random, 0.975, na.rm = TRUE))
  )
}

run_crossings <- function(pair_dt, distance_cols, group_cols, scope_label, n_boot = 2000L) {
  groups <- if (length(group_cols) == 0L) data.table(.overall = "overall") else unique(pair_dt[, ..group_cols])
  out <- list()
  idx <- 1L

  for (g_i in seq_len(nrow(groups))) {
    g <- groups[g_i]
    sub <- if (length(group_cols) == 0L) pair_dt else pair_dt[g, on = group_cols, nomatch = 0]
    for (dist_col in distance_cols) {
      est <- estimate_crossing(sub, dist_col)
      boot <- bootstrap_crossing(
        sub,
        dist_col,
        n_boot = n_boot,
        seed = 20260707L + g_i + match(dist_col, distance_cols) * 1000L
      )
      row <- cbind(
        if (length(group_cols) == 0L) data.table() else g,
        data.table(scope = scope_label, distance_metric = dist_col),
        est,
        boot
      )
      out[[idx]] <- row
      idx <- idx + 1L
    }
  }

  rbindlist(out, fill = TRUE)
}

script_dir <- get_script_dir()
trace_mcmc_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
out_dir <- file.path(trace_mcmc_dir, "fticr_integration", "output")

pair_path <- file.path(out_dir, "fticr_q2_pairwise_temporal_distance_pairs.csv")
if (!file.exists(pair_path)) {
  stop("Missing pairwise distance table. Run 07_q2_temporal_similarity_tests.R first: ", pair_path)
}

message("Loading pairwise temporal-distance table...")
pairs <- fread(pair_path)
pairs <- pairs[lag_days > 0 & is.finite(log2_lag)]

distance_cols <- c("formula_jaccard_distance", "vk_class_braycurtis", "metric_vector_distance")

message("Estimating effective randomization times...")
results <- rbindlist(list(
  run_crossings(pairs, distance_cols, character(), "overall", n_boot = 2000L),
  run_crossings(pairs, distance_cols, c("depth_label"), "by_depth", n_boot = 2000L),
  run_crossings(pairs, distance_cols, c("treatment_label", "depth_label"), "by_treatment_depth", n_boot = 2000L)
), fill = TRUE)

setcolorder(results, c("scope", "distance_metric", "treatment_label", "depth_label", setdiff(names(results), c("scope", "distance_metric", "treatment_label", "depth_label"))))
setorder(results, scope, distance_metric, depth_label, treatment_label)

out_path <- file.path(out_dir, "fticr_q2_effective_randomization_time.csv")
summary_path <- file.path(out_dir, "fticr_q2_effective_randomization_time_summary.txt")
fwrite(results, out_path)

overall_print <- results[
  scope == "overall",
  .(
    distance_metric,
    n_pairs,
    n_series,
    slope_per_lag_doubling,
    random_mean_distance,
    effective_randomization_days,
    crossing_lwr_days,
    crossing_upr_days
  )
]

depth_print <- results[
  scope == "by_depth",
  .(
    distance_metric,
    depth_label,
    n_pairs,
    n_series,
    slope_per_lag_doubling,
    random_mean_distance,
    effective_randomization_days,
    crossing_lwr_days,
    crossing_upr_days
  )
][order(distance_metric, depth_label)]

treatment_depth_print <- results[
  scope == "by_treatment_depth",
  .(
    distance_metric,
    treatment_label,
    depth_label,
    n_pairs,
    n_series,
    slope_per_lag_doubling,
    random_mean_distance,
    effective_randomization_days,
    crossing_lwr_days,
    crossing_upr_days
  )
][order(distance_metric, depth_label, treatment_label)]

summary_lines <- c(
  "TRACE FTICR effective randomization time",
  "======================================",
  "",
  "Question:",
  "How long does it take for FTICR sample distance to become indistinguishable from random?",
  "",
  "Definition used here:",
  "- Random means a random pair of samples from the same plot-depth series, ignoring resampling interval.",
  "- The fitted distance model is distance ~ log2(days between samples).",
  "- Effective randomization time is the lag where the fitted distance reaches the same-series random-pair mean distance.",
  "- Intervals are bootstrap quantiles from resampling plot-depth series; treatment-depth intervals are especially uncertain because each group has only 3 series.",
  "",
  "Overall estimates:",
  capture.output(print(overall_print)),
  "",
  "By-depth estimates:",
  capture.output(print(depth_print)),
  "",
  "By-treatment-depth estimates:",
  capture.output(print(treatment_depth_print)),
  "",
  "Output:",
  paste("-", out_path),
  "",
  "Caveats:",
  "- This is an effective timescale, not a hard point where samples suddenly lose all information.",
  "- The random baseline depends on the observed sampling window and on restricting random pairs to the same plot-depth series.",
  "- Because the model is log-linear in lag, estimates should be used as approximate order-of-magnitude thresholds."
)
writeLines(summary_lines, summary_path)

message("Wrote: ", out_path)
message("Wrote: ", summary_path)
