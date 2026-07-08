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

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

lag_bin <- function(days) {
  fcase(
    days <= 45, "00-45 d",
    days <= 90, "46-90 d",
    days <= 180, "91-180 d",
    days <= 365, "181-365 d",
    default = ">365 d"
  )
}

make_pairs <- function(dt) {
  pair_rows <- list()
  idx <- 1L

  for (grp in split(dt, by = c("metric", "plot", "depth_label"), keep.by = TRUE)) {
    grp <- grp[order(date)]
    if (nrow(grp) < 2L) next
    for (i in seq_len(nrow(grp) - 1L)) {
      for (j in seq.int(i + 1L, nrow(grp))) {
        pair_rows[[idx]] <- data.table(
          metric = grp$metric[[1]],
          metric_label = grp$metric_label[[1]],
          metric_order = grp$metric_order[[1]],
          plot = grp$plot[[1]],
          treatment_label = grp$treatment_label[[1]],
          depth_label = grp$depth_label[[1]],
          date_1 = grp$date[[i]],
          date_2 = grp$date[[j]],
          lag_days = as.integer(grp$date[[j]] - grp$date[[i]]),
          z_1 = grp$z[[i]],
          z_2 = grp$z[[j]],
          raw_1 = grp$value_plot[[i]],
          raw_2 = grp$value_plot[[j]]
        )
        idx <- idx + 1L
      }
    }
  }

  pairs <- rbindlist(pair_rows, fill = TRUE)
  pairs[, lag_bin := lag_bin(lag_days)]
  pairs[, scaled_semivar := 0.5 * (z_2 - z_1)^2]
  pairs[, abs_raw_change := abs(raw_2 - raw_1)]
  pairs
}

summarize_pairs <- function(pairs, by_cols) {
  out <- pairs[
    ,
    .(
      n_pairs = .N,
      median_lag_days = as.numeric(stats::median(lag_days)),
      mean_scaled_semivar = mean(scaled_semivar, na.rm = TRUE),
      rho_hat = 1 - mean(scaled_semivar, na.rm = TRUE),
      median_abs_raw_change = stats::median(abs_raw_change, na.rm = TRUE)
    ),
    by = by_cols
  ]
  sort_cols <- intersect(c("metric_order", "depth_label", "treatment_label", "lag_bin"), names(out))
  if (length(sort_cols) > 0L) setorderv(out, sort_cols)
  out
}

short_long_test <- function(pairs, group_cols, n_perm = 1000L, seed = 20260707L) {
  set.seed(seed)
  groups <- unique(pairs[, ..group_cols])
  out <- vector("list", nrow(groups))

  for (g_i in seq_len(nrow(groups))) {
    g <- groups[g_i]
    sub <- pairs[g, on = group_cols, nomatch = 0]
    short <- sub[lag_days <= 90]
    long <- sub[lag_days > 180]
    if (nrow(short) == 0L || nrow(long) == 0L) next

    observed_delta <- mean(long$scaled_semivar, na.rm = TRUE) - mean(short$scaled_semivar, na.rm = TRUE)
    observed_short_rho <- 1 - mean(short$scaled_semivar, na.rm = TRUE)
    observed_long_rho <- 1 - mean(long$scaled_semivar, na.rm = TRUE)

    null_delta <- rep(NA_real_, n_perm)
    value_dt <- unique(sub[, .(metric, plot, depth_label, date_1, z_1)])
    setnames(value_dt, c("date_1", "z_1"), c("date", "z"))
    value_dt <- unique(rbind(
      value_dt,
      unique(sub[, .(metric, plot, depth_label, date_2, z_2)])[
        ,
        .(metric, plot, depth_label, date = date_2, z = z_2)
      ]
    ))

    for (b in seq_len(n_perm)) {
      perm_values <- value_dt[
        ,
        .(date, z_perm = sample(z, length(z), replace = FALSE)),
        by = .(metric, plot, depth_label)
      ]
      p1 <- perm_values[, .(metric, plot, depth_label, date_1 = date, z_perm_1 = z_perm)]
      p2 <- perm_values[, .(metric, plot, depth_label, date_2 = date, z_perm_2 = z_perm)]
      perm_sub <- merge(
        sub[, .(metric, plot, depth_label, date_1, date_2, lag_days)],
        p1,
        by = c("metric", "plot", "depth_label", "date_1"),
        all.x = TRUE
      )
      perm_sub <- merge(
        perm_sub,
        p2,
        by = c("metric", "plot", "depth_label", "date_2"),
        all.x = TRUE
      )
      perm_sub[, scaled_semivar := 0.5 * (z_perm_2 - z_perm_1)^2]
      null_delta[[b]] <- mean(perm_sub[lag_days > 180, scaled_semivar], na.rm = TRUE) -
        mean(perm_sub[lag_days <= 90, scaled_semivar], na.rm = TRUE)
    }

    p_one_sided <- (1 + sum(null_delta >= observed_delta, na.rm = TRUE)) / (1 + sum(is.finite(null_delta)))
    out[[g_i]] <- cbind(
      g,
      data.table(
        n_short_pairs = nrow(short),
        n_long_pairs = nrow(long),
        short_scaled_semivar = mean(short$scaled_semivar, na.rm = TRUE),
        long_scaled_semivar = mean(long$scaled_semivar, na.rm = TRUE),
        short_rho_hat = observed_short_rho,
        long_rho_hat = observed_long_rho,
        long_minus_short_semivar = observed_delta,
        permutation_p_short_more_similar = p_one_sided
      )
    )
  }

  rbindlist(out, fill = TRUE)
}

script_dir <- get_script_dir()
trace_mcmc_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
out_dir <- file.path(trace_mcmc_dir, "fticr_integration", "output")
plot_data_path <- file.path(out_dir, "fticr_formula_count_temporal_metrics_plot_data.csv")
if (!file.exists(plot_data_path)) {
  stop("Missing plot-data file. Run 05_plot_fticr_formula_count_timeseries.R first: ", plot_data_path)
}

message("Loading formula-count FTICR plot data...")
dt <- fread(plot_data_path)
dt[, date := as.IDate(date)]
dt <- dt[is.finite(value_plot)]

# Standardize within each plot-depth metric series so pairwise semivariance
# reflects temporal persistence rather than stable plot/depth differences.
dt[
  ,
  `:=`(
    series_mean = mean(value_plot, na.rm = TRUE),
    series_sd = safe_sd(value_plot),
    series_n = sum(is.finite(value_plot))
  ),
  by = .(metric, plot, depth_label)
]
dt <- dt[is.finite(series_sd) & series_sd > 0 & series_n >= 3]
dt[, z := (value_plot - series_mean) / series_sd]

message("Building within-plot/depth sample pairs...")
pairs <- make_pairs(dt)

lag_levels <- c("00-45 d", "46-90 d", "91-180 d", "181-365 d", ">365 d")
pairs[, lag_bin := factor(lag_bin, levels = lag_levels)]

lag_summary <- rbindlist(list(
  summarize_pairs(pairs, c("metric", "metric_label", "metric_order", "lag_bin"))[, scope := "overall"],
  summarize_pairs(pairs, c("metric", "metric_label", "metric_order", "depth_label", "lag_bin"))[, scope := "by_depth"],
  summarize_pairs(pairs, c("metric", "metric_label", "metric_order", "treatment_label", "depth_label", "lag_bin"))[, scope := "by_treatment_depth"]
), fill = TRUE)
setcolorder(lag_summary, c("scope", setdiff(names(lag_summary), "scope")))

message("Running permutation tests for short-lag similarity...")
short_long <- rbindlist(list(
  short_long_test(pairs, c("metric", "metric_label", "metric_order"))[, scope := "overall"],
  short_long_test(pairs, c("metric", "metric_label", "metric_order", "depth_label"))[, scope := "by_depth"],
  short_long_test(pairs, c("metric", "metric_label", "metric_order", "treatment_label", "depth_label"))[, scope := "by_treatment_depth"]
), fill = TRUE)
setcolorder(short_long, c("scope", setdiff(names(short_long), "scope")))
setorder(short_long, scope, metric_order, depth_label, treatment_label)

pairs_out <- file.path(out_dir, "fticr_formula_count_autocorrelation_pairs.csv")
lag_summary_out <- file.path(out_dir, "fticr_formula_count_autocorrelation_lag_bins.csv")
short_long_out <- file.path(out_dir, "fticr_formula_count_autocorrelation_short_long.csv")
summary_out <- file.path(out_dir, "fticr_formula_count_autocorrelation_summary.txt")

fwrite(pairs, pairs_out)
fwrite(lag_summary, lag_summary_out)
fwrite(short_long, short_long_out)

overall <- short_long[scope == "overall"][
  order(metric_order),
  .(
    metric_label,
    n_short_pairs,
    n_long_pairs,
    short_rho_hat,
    long_rho_hat,
    long_minus_short_semivar,
    permutation_p_short_more_similar
  )
]

depth <- short_long[scope == "by_depth"][
  order(metric_order, depth_label),
  .(
    metric_label,
    depth_label,
    n_short_pairs,
    n_long_pairs,
    short_rho_hat,
    long_rho_hat,
    long_minus_short_semivar,
    permutation_p_short_more_similar
  )
]

summary_lines <- c(
  "FTICR formula-count temporal autocorrelation diagnostic",
  "=======================================================",
  "",
  paste("Input:", plot_data_path),
  paste("Samples used after requiring >=3 observations per metric/plot/depth series:", uniqueN(dt$sample_num)),
  paste("Within-plot/depth pairs:", nrow(pairs)),
  "",
  "Method:",
  "- Values are standardized within each metric x plot x depth series.",
  "- Scaled semivariance = 0.5 * (z_t2 - z_t1)^2.",
  "- rho_hat = 1 - scaled semivariance; positive values indicate temporal similarity.",
  "- Short-lag evidence compares <=90 d pairs against >180 d pairs.",
  "- Permutation p-values shuffle values within each metric x plot x depth series.",
  "",
  "Overall short-lag versus long-lag autocorrelation:",
  capture.output(print(overall)),
  "",
  "By-depth short-lag versus long-lag autocorrelation:",
  capture.output(print(depth)),
  "",
  "Caveats:",
  "- This is an irregular-sampling semivariogram diagnostic, not a regular-interval ACF.",
  "- rho_hat is approximate because it pools standardized plot-depth series.",
  "- Positive short-lag rho_hat with positive long-minus-short semivariance is evidence that nearby samples are more similar."
)
writeLines(summary_lines, summary_out)

message("Wrote: ", pairs_out)
message("Wrote: ", lag_summary_out)
message("Wrote: ", short_long_out)
message("Wrote: ", summary_out)
