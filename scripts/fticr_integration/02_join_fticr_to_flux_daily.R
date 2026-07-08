suppressPackageStartupMessages(library(data.table))

get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)))
  }
  getwd()
}

nearest_by_plot <- function(flux_dt, fticr_dt, max_days = 45L) {
  out <- copy(flux_dt)

  # Precreate output columns.
  near_cols <- setdiff(names(fticr_dt), c("plot", "sample_date", "depth", "treatment"))
  for (col in near_cols) {
    out[, (paste0("near_", col)) := NA_real_]
  }
  out[, near_sample_date := as.IDate(NA)]
  out[, near_days_from_sample := as.integer(NA)]

  for (p in sort(unique(out$plot))) {
    i_flux <- which(out$plot == p)
    if (length(i_flux) == 0L) {
      next
    }

    s <- fticr_dt[plot == p][order(sample_date)]
    if (nrow(s) == 0L) {
      next
    }

    d_flux <- as.integer(out$date_id[i_flux])
    d_sample <- as.integer(s$sample_date)

    idx_prev <- findInterval(d_flux, d_sample)
    idx_prev[idx_prev == 0L] <- NA_integer_

    idx_next <- idx_prev + 1L
    idx_next[idx_next > length(d_sample)] <- NA_integer_

    d_prev <- ifelse(is.na(idx_prev), Inf, abs(d_flux - d_sample[idx_prev]))
    d_next <- ifelse(is.na(idx_next), Inf, abs(d_flux - d_sample[idx_next]))

    use_prev <- d_prev <= d_next
    best_idx <- ifelse(use_prev, idx_prev, idx_next)
    best_dist <- pmin(d_prev, d_next)

    # Invalidate matches outside the user window.
    best_idx[best_dist > max_days] <- NA_integer_
    best_dist[best_dist > max_days] <- NA_real_

    out$near_sample_date[i_flux] <- as.IDate(ifelse(is.na(best_idx), NA, d_sample[best_idx]), origin = "1970-01-01")
    out$near_days_from_sample[i_flux] <- as.integer(best_dist)

    valid <- !is.na(best_idx)
    if (any(valid)) {
      for (col in near_cols) {
        target_col <- paste0("near_", col)
        out[[target_col]][i_flux[valid]] <- as.numeric(s[[col]][best_idx[valid]])
      }
    }
  }

  out
}

script_dir <- get_script_dir()
trace_mcmc_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

out_dir <- file.path(trace_mcmc_dir, "fticr_integration", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

flux_path <- file.path(trace_mcmc_dir, "data", "flux_daily_bgc.csv")
fticr_sample_path <- file.path(out_dir, "fticr_sample_features.csv")

required_paths <- c(flux_path, fticr_sample_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop(
    "Missing required input files for FTICR join:\n",
    paste(missing_paths, collapse = "\n"),
    "\nRun 01_build_fticr_features.R first if FTICR sample features are missing."
  )
}

message("Loading daily flux dataset...")
flux_dt <- fread(flux_path)
flux_dt[, date_id := as.IDate(date)]

message("Loading FTICR sample feature dataset...")
fticr_dt <- fread(fticr_sample_path)
fticr_dt[, sample_date := as.IDate(sample_date)]

# Prefer 0-10 cm FTICR depth for surface respiration models.
fticr_surface <- fticr_dt[grepl("^0-10", depth) & !is.na(plot) & !is.na(sample_date)]

feature_cols <- setdiff(
  names(fticr_surface),
  c(
    "sample_num",
    "sample_name_raw",
    "sample_is_extra",
    "sample_date",
    "plot",
    "treatment",
    "treatment_raw",
    "depth"
  )
)

if (length(feature_cols) == 0L) {
  stop("No FTICR feature columns available to join.")
}

fticr_surface_daily <- fticr_surface[
  ,
  c(lapply(.SD, mean, na.rm = TRUE), .(n_fticr_samples = .N)),
  by = .(plot, sample_date, depth, treatment),
  .SDcols = feature_cols
]

exact_join <- merge(
  flux_dt,
  fticr_surface_daily,
  by.x = c("plot", "date_id"),
  by.y = c("plot", "sample_date"),
  all.x = TRUE,
  suffixes = c("", "_fticr")
)
exact_join[, fticr_exact_match := !is.na(n_fticr_samples)]

# Build nearest-date features (by plot) for broader model coverage.
nearest_join <- nearest_by_plot(
  flux_dt = flux_dt,
  fticr_dt = fticr_surface_daily,
  max_days = 45L
)

exact_out_path <- file.path(out_dir, "flux_daily_with_fticr_exact.csv")
nearest_out_path <- file.path(out_dir, "flux_daily_with_fticr_nearest45d.csv")
summary_out_path <- file.path(out_dir, "fticr_flux_join_summary.txt")

fwrite(exact_join, exact_out_path)
fwrite(nearest_join, nearest_out_path)

summary_lines <- c(
  paste("Flux daily rows:", nrow(flux_dt)),
  paste("FTICR 0-10 cm sample rows used:", nrow(fticr_surface)),
  paste("Exact joined rows with FTICR:", sum(exact_join$fticr_exact_match, na.rm = TRUE)),
  paste("Nearest45 rows with FTICR:", sum(!is.na(nearest_join$near_days_from_sample))),
  paste(
    "Nearest45 median days to FTICR sample:",
    as.character(stats::median(nearest_join$near_days_from_sample, na.rm = TRUE))
  ),
  paste(
    "Flux date span:",
    as.character(min(flux_dt$date_id, na.rm = TRUE)),
    "to",
    as.character(max(flux_dt$date_id, na.rm = TRUE))
  )
)
writeLines(summary_lines, summary_out_path)

message("FTICR-to-flux joins complete.")
message("Wrote: ", exact_out_path)
message("Wrote: ", nearest_out_path)
message("Wrote: ", summary_out_path)
