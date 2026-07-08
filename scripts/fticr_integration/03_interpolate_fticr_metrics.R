#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)))
  }
  getwd()
}

get_arg_value <- function(args, name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  sub(prefix, "", hit[[1]])
}

safe_mean <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

interp_group <- function(obs, grid_dates, metric_cols, max_gap_days) {
  obs <- unique(obs[order(date)], by = "date")
  sample_dates <- obs$date
  sample_days <- as.integer(sample_dates)
  grid_days <- as.integer(grid_dates)

  out <- data.table(
    date = grid_dates,
    fticr_exact_match = grid_dates %in% sample_dates,
    fticr_prev_date = as.IDate(NA),
    fticr_next_date = as.IDate(NA),
    fticr_bracket_gap_days = as.integer(NA),
    fticr_interp_fraction = NA_real_
  )

  for (metric in metric_cols) {
    out[, (metric) := NA_real_]
  }

  if (nrow(obs) == 0L || length(grid_dates) == 0L) {
    return(out)
  }

  exact_idx <- match(grid_dates, sample_dates)
  exact_ok <- !is.na(exact_idx)
  if (any(exact_ok)) {
    out[exact_ok, fticr_prev_date := sample_dates[exact_idx[exact_ok]]]
    out[exact_ok, fticr_next_date := sample_dates[exact_idx[exact_ok]]]
    out[exact_ok, fticr_bracket_gap_days := 0L]
    out[exact_ok, fticr_interp_fraction := 0]
    for (metric in metric_cols) {
      set(out, which(exact_ok), metric, obs[[metric]][exact_idx[exact_ok]])
    }
  }

  interp_rows <- which(!exact_ok)
  if (length(interp_rows) == 0L || length(sample_days) < 2L) {
    return(out)
  }

  prev_idx <- findInterval(grid_days[interp_rows], sample_days)
  next_idx <- prev_idx + 1L
  valid <- prev_idx >= 1L & next_idx <= length(sample_days)
  if (!any(valid)) {
    return(out)
  }

  rows_valid <- interp_rows[valid]
  prev_idx <- prev_idx[valid]
  next_idx <- next_idx[valid]
  bracket_gap <- sample_days[next_idx] - sample_days[prev_idx]
  keep <- bracket_gap <= max_gap_days
  if (!any(keep)) {
    return(out)
  }

  rows_keep <- rows_valid[keep]
  prev_idx <- prev_idx[keep]
  next_idx <- next_idx[keep]
  bracket_gap <- bracket_gap[keep]
  fraction <- (grid_days[rows_keep] - sample_days[prev_idx]) / bracket_gap

  out[rows_keep, fticr_prev_date := sample_dates[prev_idx]]
  out[rows_keep, fticr_next_date := sample_dates[next_idx]]
  out[rows_keep, fticr_bracket_gap_days := as.integer(bracket_gap)]
  out[rows_keep, fticr_interp_fraction := fraction]

  for (metric in metric_cols) {
    y0 <- obs[[metric]][prev_idx]
    y1 <- obs[[metric]][next_idx]
    pred <- y0 + fraction * (y1 - y0)
    pred[is.na(y0) | is.na(y1)] <- NA_real_
    set(out, rows_keep, metric, pred)
  }

  out
}

add_historical_mean_fill <- function(interp_daily, sample_daily, metric_cols) {
  out <- copy(interp_daily)

  for (metric in metric_cols) {
    plot_depth_stats <- sample_daily[
      is.finite(get(metric)),
      .(
        fill_mean = mean(get(metric), na.rm = TRUE),
        fill_sd = stats::sd(get(metric), na.rm = TRUE),
        fill_n = .N
      ),
      by = .(plot, depth_cm)
    ]
    plot_depth_stats[!is.finite(fill_sd), fill_sd := 0]
    setnames(
      plot_depth_stats,
      c("fill_mean", "fill_sd", "fill_n"),
      paste0(metric, c("_plot_depth_mean", "_plot_depth_sd", "_plot_depth_n"))
    )

    depth_stats <- sample_daily[
      is.finite(get(metric)),
      .(
        fill_mean = mean(get(metric), na.rm = TRUE),
        fill_sd = stats::sd(get(metric), na.rm = TRUE),
        fill_n = .N
      ),
      by = depth_cm
    ]
    depth_stats[!is.finite(fill_sd), fill_sd := 0]
    setnames(
      depth_stats,
      c("fill_mean", "fill_sd", "fill_n"),
      paste0(metric, c("_depth_mean", "_depth_sd", "_depth_n"))
    )

    global_values <- sample_daily[is.finite(get(metric)), get(metric)]
    global_mean <- if (length(global_values) > 0L) mean(global_values, na.rm = TRUE) else NA_real_
    global_sd <- if (length(global_values) > 1L) stats::sd(global_values, na.rm = TRUE) else 0
    if (!is.finite(global_sd)) global_sd <- 0
    global_n <- length(global_values)

    out <- merge(out, plot_depth_stats, by = c("plot", "depth_cm"), all.x = TRUE, sort = FALSE)
    out <- merge(out, depth_stats, by = "depth_cm", all.x = TRUE, sort = FALSE)

    x <- as.numeric(out[[metric]])
    plot_depth_mean_col <- paste0(metric, "_plot_depth_mean")
    plot_depth_sd_col <- paste0(metric, "_plot_depth_sd")
    plot_depth_n_col <- paste0(metric, "_plot_depth_n")
    depth_mean_col <- paste0(metric, "_depth_mean")
    depth_sd_col <- paste0(metric, "_depth_sd")
    depth_n_col <- paste0(metric, "_depth_n")

    filled <- x
    method <- rep(NA_character_, nrow(out))
    fill_sd <- rep(NA_real_, nrow(out))
    fill_n <- rep(NA_integer_, nrow(out))

    observed <- is.finite(x)
    exact_row <- !is.na(out$fticr_exact_match) & out$fticr_exact_match
    method[observed & exact_row] <- "exact"
    method[observed & !exact_row] <- "linear_interp"

    need <- !is.finite(filled)
    use_plot_depth <- need & is.finite(out[[plot_depth_mean_col]])
    filled[use_plot_depth] <- out[[plot_depth_mean_col]][use_plot_depth]
    method[use_plot_depth] <- "plot_depth_mean"
    fill_sd[use_plot_depth] <- out[[plot_depth_sd_col]][use_plot_depth]
    fill_n[use_plot_depth] <- as.integer(out[[plot_depth_n_col]][use_plot_depth])

    need <- !is.finite(filled)
    use_depth <- need & is.finite(out[[depth_mean_col]])
    filled[use_depth] <- out[[depth_mean_col]][use_depth]
    method[use_depth] <- "depth_mean"
    fill_sd[use_depth] <- out[[depth_sd_col]][use_depth]
    fill_n[use_depth] <- as.integer(out[[depth_n_col]][use_depth])

    need <- !is.finite(filled)
    use_global <- need & is.finite(global_mean)
    filled[use_global] <- global_mean
    method[use_global] <- "global_mean"
    fill_sd[use_global] <- global_sd
    fill_n[use_global] <- as.integer(global_n)

    method[!is.finite(filled)] <- "missing"

    out[, (paste0(metric, "_filled")) := filled]
    out[, (paste0(metric, "_fill_method")) := method]
    out[, (paste0(metric, "_fill_sd")) := fill_sd]
    out[, (paste0(metric, "_fill_n")) := fill_n]

    out[, c(
      plot_depth_mean_col, plot_depth_sd_col, plot_depth_n_col,
      depth_mean_col, depth_sd_col, depth_n_col
    ) := NULL]
  }

  method_cols <- paste0(metric_cols, "_fill_method")
  out[, fticr_any_historical_mean_fill := Reduce(`|`, lapply(.SD, `%in%`, c("plot_depth_mean", "depth_mean", "global_mean"))), .SDcols = method_cols]
  out
}

args <- commandArgs(trailingOnly = TRUE)
max_gap_days <- as.integer(get_arg_value(args, "max-gap-days", "90"))
if (is.na(max_gap_days) || max_gap_days <= 0L) {
  stop("--max-gap-days must be a positive integer.")
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

message("Loading FTICR sample-level CUE and chemodiversity metrics...")
fticr <- fread(fticr_sample_path)
fticr[, date := as.IDate(sample_date)]
fticr[, depth_cm := fifelse(grepl("^0-10", depth), 10L, fifelse(grepl("^10-30", depth), 30L, NA_integer_))]
fticr <- fticr[
  !is.na(plot) & !is.na(date) & !is.na(depth_cm),
  .(sample_num, plot, treatment, date, depth_cm, cue_wmean)
]

chemodiv <- fread(chemodiv_path)
chemodiv[, date := as.IDate(date)]
chemodiv <- chemodiv[!is.na(plot) & !is.na(date) & !is.na(depth_cm)]

metric_cols <- c(
  "cue_wmean",
  "richness",
  "shannon_classes",
  "evenness",
  "mean_nosc",
  "mean_hc",
  "mean_oc",
  "mean_nc",
  "mean_ai",
  "mean_dbe",
  "frac_labile",
  "frac_aromatic",
  "frac_lignin",
  "frac_other"
)

sample_metrics <- merge(
  fticr,
  chemodiv[, c("sample_num", "plot", "treatment", "date", "depth_cm", setdiff(metric_cols, "cue_wmean")), with = FALSE],
  by = c("sample_num", "plot", "treatment", "date", "depth_cm"),
  all = TRUE
)

for (metric in metric_cols) {
  if (!(metric %in% names(sample_metrics))) {
    sample_metrics[, (metric) := NA_real_]
  }
  sample_metrics[, (metric) := as.numeric(get(metric))]
}

sample_daily <- sample_metrics[
  ,
  lapply(.SD, safe_mean),
  by = .(plot, treatment, date, depth_cm),
  .SDcols = metric_cols
]
setorderv(sample_daily, c("plot", "depth_cm", "date"))

message("Loading daily flux grid...")
flux <- fread(flux_path)
flux[, date := as.IDate(date)]
grid <- unique(flux[!is.na(plot) & !is.na(date), .(plot, date)])

depths <- sort(unique(sample_daily$depth_cm))
interp_parts <- vector("list", length(depths) * length(unique(grid$plot)))
part_i <- 1L
for (depth in depths) {
  for (plot_id in sort(unique(grid$plot))) {
    obs <- sample_daily[plot == plot_id & depth_cm == depth]
    grid_dates <- sort(grid[plot == plot_id, date])
    interp <- interp_group(obs, grid_dates, metric_cols, max_gap_days)
    interp[, plot := plot_id]
    interp[, depth_cm := depth]
    interp_parts[[part_i]] <- interp
    part_i <- part_i + 1L
  }
}

interp_daily <- rbindlist(interp_parts, fill = TRUE)
setcolorder(interp_daily, c(
  "plot",
  "date",
  "depth_cm",
  "fticr_exact_match",
  "fticr_prev_date",
  "fticr_next_date",
  "fticr_bracket_gap_days",
  "fticr_interp_fraction",
  metric_cols
))
setorderv(interp_daily, c("plot", "depth_cm", "date"))

interp_daily <- add_historical_mean_fill(interp_daily, sample_daily, metric_cols)

daily_out_path <- file.path(out_dir, sprintf("fticr_interpolated_daily_metrics_gap%dd.csv", max_gap_days))
fwrite(interp_daily, daily_out_path)

surface_interp <- copy(interp_daily[depth_cm == 10L])
surface_cols <- setdiff(names(surface_interp), c("plot", "date", "depth_cm"))
setnames(surface_interp, surface_cols, paste0("interp_", surface_cols))

flux_join <- merge(flux, surface_interp[, c("plot", "date", paste0("interp_", surface_cols)), with = FALSE],
  by = c("plot", "date"),
  all.x = TRUE,
  sort = FALSE
)
setorderv(flux_join, c("plot", "date"))

flux_out_path <- file.path(out_dir, sprintf("flux_daily_with_fticr_interp10cm_gap%dd.csv", max_gap_days))
fwrite(flux_join, flux_out_path)

coverage <- interp_daily[
  ,
  .(
    rows = .N,
    rows_with_any_metric = sum(rowSums(!is.na(.SD)) > 0L),
    exact_rows = sum(fticr_exact_match, na.rm = TRUE),
    interpolated_rows = sum(!fticr_exact_match & !is.na(fticr_bracket_gap_days), na.rm = TRUE),
    median_bracket_gap_days = as.numeric(median(fticr_bracket_gap_days[!fticr_exact_match], na.rm = TRUE))
  ),
  by = depth_cm,
  .SDcols = metric_cols
]

metric_coverage <- melt(
  interp_daily[, c("depth_cm", metric_cols), with = FALSE],
  id.vars = "depth_cm",
  variable.name = "metric",
  value.name = "value"
)[
  ,
  .(non_missing_rows = sum(!is.na(value)), pct_non_missing = round(100 * mean(!is.na(value)), 1)),
  by = .(depth_cm, metric)
]

filled_metric_cols <- paste0(metric_cols, "_filled")
metric_filled_coverage <- melt(
  interp_daily[, c("depth_cm", filled_metric_cols), with = FALSE],
  id.vars = "depth_cm",
  variable.name = "metric",
  value.name = "value"
)[
  ,
  `:=`(metric = sub("_filled$", "", metric))
][
  ,
  .(non_missing_rows = sum(!is.na(value)), pct_non_missing = round(100 * mean(!is.na(value)), 1)),
  by = .(depth_cm, metric)
]

fill_method_coverage <- rbindlist(lapply(metric_cols, function(metric) {
  method_col <- paste0(metric, "_fill_method")
  ans <- interp_daily[
    ,
    .N,
    by = .(depth_cm, fill_method = get(method_col))
  ]
  ans[, metric := metric]
  setcolorder(ans, c("depth_cm", "metric", "fill_method", "N"))
  ans
}), fill = TRUE)

summary_out_path <- file.path(out_dir, sprintf("fticr_interpolation_summary_gap%dd.txt", max_gap_days))
summary_lines <- c(
  paste("Max interpolation bracket gap days:", max_gap_days),
  paste("Daily interpolated rows:", nrow(interp_daily)),
  paste("Flux daily rows:", nrow(flux)),
  "",
  "Coverage by depth:",
  capture.output(print(coverage)),
  "",
  "Metric coverage by depth:",
  capture.output(print(metric_coverage[order(depth_cm, metric)])),
  "",
  "Filled metric coverage by depth:",
  capture.output(print(metric_filled_coverage[order(depth_cm, metric)])),
  "",
  "Fill methods by depth and metric:",
  capture.output(print(fill_method_coverage[order(depth_cm, metric, fill_method)])),
  "",
  "Interpolation rule:",
  "Linear interpolation is used only between bracketing FTICR samples from the same plot and depth.",
  "Exact sample dates retain observed values.",
  "Dates outside the first/last sample or inside brackets wider than the max gap are left missing.",
  "Unfilled metric values are imputed from plot-depth historical means; depth means and then global means are used only if a plot-depth mean is unavailable.",
  "For each metric, *_filled contains the value passed to models that request imputed covariates, with *_fill_method, *_fill_sd, and *_fill_n documenting the fallback.",
  "The flux join uses 0-10 cm interpolated metrics."
)
writeLines(summary_lines, summary_out_path)

message("FTICR interpolation complete.")
message("Wrote: ", daily_out_path)
message("Wrote: ", flux_out_path)
message("Wrote: ", summary_out_path)
