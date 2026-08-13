#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)))
  }
  getwd()
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

metadata_row <- function(output_column, role, analysis_group, source_column = NA_character_,
                         source_file = NA_character_, source_grain = NA_character_,
                         aggregation = NA_character_, description = NA_character_,
                         datasource = NA_character_, measurement = NA_character_,
                         unit = NA_character_, depth = NA_character_,
                         source_summary_stat = NA_character_, source_summary_class = NA_character_,
                         temporal_match = NA_character_, depth_match = NA_character_,
                         notes = NA_character_) {
  data.table(
    output_column = output_column,
    role = role,
    analysis_group = analysis_group,
    source_column = source_column,
    source_file = source_file,
    source_grain = source_grain,
    aggregation = aggregation,
    description = description,
    datasource = datasource,
    measurement = measurement,
    unit = unit,
    depth = depth,
    source_summary_stat = source_summary_stat,
    source_summary_class = source_summary_class,
    temporal_match = temporal_match,
    depth_match = depth_match,
    notes = notes
  )
}

script_dir <- get_script_dir()
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
trace_mcmc_dir <- Sys.getenv("TRACE_MCMC_DIR", unset = "/Users/jrweverk/Documents/TRACE_MCMC")

hourly_path <- file.path(trace_mcmc_dir, "data", "joined_data_flux_final_2_hourly.csv")
daily_path <- file.path(trace_mcmc_dir, "data", "joined_data_flux_final_2_daily.csv")
porewater_path <- file.path(trace_mcmc_dir, "data_cleanup", "data_cleaned", "pwchemistry.csv")
fticr_path <- file.path(repo_dir, "data", "derived", "fticr_chemodiversity_metrics.csv")
monthly_metadata_path <- file.path(
  repo_dir,
  "data",
  "derived",
  "monthly_covariates",
  "fticr_monthly_covariates_metadata.csv"
)
out_dir <- file.path(repo_dir, "data", "derived", "same_day_covariates")
out_path <- file.path(out_dir, "fticr_same_day_covariates.csv")
metadata_out_path <- file.path(out_dir, "fticr_same_day_covariates_metadata.csv")

required_paths <- c(hourly_path, daily_path, porewater_path, fticr_path, monthly_metadata_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required files:\n", paste(missing_paths, collapse = "\n"))
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading FT-ICR sample keys and reduced covariate specification...")
fticr <- fread(fticr_path, na.strings = c("", "NA", "NaN", "-9999"))
required_fticr_cols <- c("sample_num", "plot", "date", "depth_cm", "treatment")
missing_fticr_cols <- setdiff(required_fticr_cols, names(fticr))
if (length(missing_fticr_cols) > 0L) {
  stop("FT-ICR table is missing columns: ", paste(missing_fticr_cols, collapse = ", "))
}

samples <- fticr[, ..required_fticr_cols]
samples[, `:=`(
  sample_num = as.integer(sample_num),
  plot = as.integer(plot),
  date = as.IDate(date),
  depth_cm = as.integer(depth_cm),
  treatment = as.character(treatment)
)]
setorder(samples, sample_num)

if (anyNA(samples[, .(sample_num, plot, date, depth_cm, treatment)])) {
  stop("FT-ICR sample keys contain missing values.")
}
if (anyDuplicated(samples$sample_num)) stop("FT-ICR sample_num is not unique.")
if (anyDuplicated(samples[, .(plot, date, depth_cm)])) {
  stop("FT-ICR plot-date-depth keys are not unique.")
}
if (!all(samples$depth_cm %in% c(10L, 30L))) {
  stop("Exact-depth porewater matching is defined only for 10 and 30 cm FT-ICR samples.")
}

monthly_meta <- fread(monthly_metadata_path, na.strings = c("", "NA"))
required_meta_cols <- c(
  "output_column", "role", "analysis_group", "source_column", "source_file",
  "source_grain", "datasource", "measurement", "unit", "depth",
  "source_summary_stat", "source_summary_class", "notes"
)
missing_meta_cols <- setdiff(required_meta_cols, names(monthly_meta))
if (length(missing_meta_cols) > 0L) {
  stop("Monthly metadata is missing columns: ", paste(missing_meta_cols, collapse = ", "))
}

selected_groups <- c("temperature", "vwc", "porewater")
selected_meta <- monthly_meta[
  role == "covariate" &
    !is.na(source_column) &
    analysis_group %chin% selected_groups
]
hourly_meta <- selected_meta[source_grain == "hourly"]
daily_meta <- selected_meta[source_grain == "daily" & datasource != "pwchemistry"]

if (anyDuplicated(hourly_meta$source_column) || anyDuplicated(daily_meta$source_column)) {
  stop("The reduced covariate specification contains duplicate source columns.")
}

hourly_source_cols <- hourly_meta$source_column
daily_source_cols <- daily_meta$source_column
sample_dates <- unique(samples[, .(plot, date)])

message("Calculating same-day means for hourly temperature and VWC...")
hourly <- fread(
  hourly_path,
  select = c("plot", "date", hourly_source_cols),
  na.strings = c("", "NA", "NaN", "-9999"),
  showProgress = TRUE
)
hourly[, `:=`(plot = as.integer(plot), date = as.IDate(date))]
hourly <- hourly[sample_dates, on = .(plot, date), nomatch = 0L]

hourly_daily <- hourly[
  ,
  lapply(.SD, safe_mean),
  by = .(plot, date),
  .SDcols = hourly_source_cols
]
hourly_output_cols <- paste0(hourly_source_cols, "__same_day_mean")
setnames(hourly_daily, hourly_source_cols, hourly_output_cols)

same_day <- merge(samples, hourly_daily, by = c("plot", "date"), all.x = TRUE, sort = FALSE)
setorder(same_day, sample_num)
rm(hourly, hourly_daily)
invisible(gc())

message("Joining sparse covariates observed on the exact FT-ICR date...")
daily <- fread(
  daily_path,
  select = c("plot", "date", daily_source_cols),
  na.strings = c("", "NA", "NaN", "-9999"),
  showProgress = FALSE
)
daily[, `:=`(plot = as.integer(plot), date = as.IDate(date))]
daily <- daily[sample_dates, on = .(plot, date), nomatch = 0L]
if (anyDuplicated(daily[, .(plot, date)])) {
  stop("The daily joined table contains duplicate plot-date rows for FT-ICR dates.")
}

has_exact_date_data <- vapply(
  daily[, ..daily_source_cols],
  function(x) any(!is.na(x)),
  logical(1)
)
empty_exact_date_cols <- daily_source_cols[!has_exact_date_data]
if (length(empty_exact_date_cols) > 0L) {
  message(
    "Omitting ",
    length(empty_exact_date_cols),
    " source columns with no observations on any FT-ICR date."
  )
}
daily_source_cols <- daily_source_cols[has_exact_date_data]
daily_meta <- daily_meta[match(daily_source_cols, source_column)]
daily <- daily[, c("plot", "date", daily_source_cols), with = FALSE]

daily_output_cols <- paste0(daily_source_cols, "__same_day")
setnames(daily, daily_source_cols, daily_output_cols)
same_day <- merge(same_day, daily, by = c("plot", "date"), all.x = TRUE, sort = FALSE)
setorder(same_day, sample_num)
rm(daily)
invisible(gc())

message("Matching porewater NPOC and TDN to each FT-ICR plot, date, and depth...")
porewater_source_cols <- c(
  "npoc_mg.l_10cm_pwchemistry",
  "npoc_mg.l_30cm_pwchemistry",
  "tdn_mg.l_10cm_pwchemistry",
  "tdn_mg.l_30cm_pwchemistry"
)
porewater <- fread(
  porewater_path,
  select = c("plot", "date", porewater_source_cols),
  na.strings = c("", "NA", "NaN", "-9999"),
  showProgress = FALSE
)
porewater[, `:=`(plot = as.integer(plot), date = as.IDate(date))]
if (anyDuplicated(porewater[, .(plot, date)])) {
  stop("The cleaned porewater table contains duplicate plot-date rows.")
}

porewater_match <- merge(
  samples[, .(sample_num, plot, date, depth_cm)],
  porewater,
  by = c("plot", "date"),
  all.x = TRUE,
  sort = FALSE
)
porewater_match[, porewater_npoc_mg_C_L__exact_sample := fcase(
  depth_cm == 10L, npoc_mg.l_10cm_pwchemistry,
  depth_cm == 30L, npoc_mg.l_30cm_pwchemistry,
  default = NA_real_
)]
porewater_match[, porewater_tdn_mg_N_L__exact_sample := fcase(
  depth_cm == 10L, tdn_mg.l_10cm_pwchemistry,
  depth_cm == 30L, tdn_mg.l_30cm_pwchemistry,
  default = NA_real_
)]
porewater_output_cols <- c(
  "porewater_npoc_mg_C_L__exact_sample",
  "porewater_tdn_mg_N_L__exact_sample"
)
porewater_match <- porewater_match[, c("sample_num", porewater_output_cols), with = FALSE]
same_day <- merge(same_day, porewater_match, by = "sample_num", all.x = TRUE, sort = FALSE)
setorder(same_day, sample_num)

if (anyNA(same_day[[porewater_output_cols[[1L]]]]) || anyNA(same_day[[porewater_output_cols[[2L]]]])) {
  stop("At least one FT-ICR sample lacks an exact plot-date-depth NPOC or TDN match.")
}

key_cols <- c("sample_num", "plot", "date", "depth_cm", "treatment")
setcolorder(same_day, c(key_cols, porewater_output_cols, hourly_output_cols, daily_output_cols))

message("Building exact-date metadata...")
metadata_rows <- list(
  metadata_row(
    "sample_num", "key", "key", description = "Unique FT-ICR spectrum identifier.",
    unit = "identifier", temporal_match = "FT-ICR sample date"
  ),
  metadata_row(
    "plot", "key", "key", description = "TRACE plot identifier (1-6).",
    unit = "identifier", temporal_match = "FT-ICR sample date"
  ),
  metadata_row(
    "date", "key", "key", description = "Porewater collection date for the FT-ICR spectrum.",
    unit = "YYYY-MM-DD", temporal_match = "FT-ICR sample date"
  ),
  metadata_row(
    "depth_cm", "key", "key", description = "Exact FT-ICR porewater sampling depth.",
    unit = "cm", temporal_match = "FT-ICR sample date"
  ),
  metadata_row(
    "treatment", "key", "key", description = "Experimental treatment for the TRACE plot.",
    unit = "control/warmed", temporal_match = "FT-ICR sample date"
  ),
  metadata_row(
    porewater_output_cols[[1L]], "covariate", "matched_porewater_cn",
    source_column = "npoc_mg.l_{depth_cm}cm_pwchemistry",
    source_file = basename(porewater_path), source_grain = "individual porewater sample",
    aggregation = "none", description = "NPOC concentration measured in the porewater sample matching the FT-ICR plot, collection date, and depth.",
    datasource = "pwchemistry", measurement = "NPOC", unit = "mg C L^-1",
    depth = "matched to depth_cm", source_summary_stat = "none",
    temporal_match = "same plot and exact FT-ICR collection date",
    depth_match = "exact FT-ICR depth (10 or 30 cm)",
    notes = "Observed concentration from the matching sample; no interpolation, nearest-date matching, or temporal averaging."
  ),
  metadata_row(
    porewater_output_cols[[2L]], "covariate", "matched_porewater_cn",
    source_column = "tdn_mg.l_{depth_cm}cm_pwchemistry",
    source_file = basename(porewater_path), source_grain = "individual porewater sample",
    aggregation = "none", description = "TDN concentration measured in the porewater sample matching the FT-ICR plot, collection date, and depth.",
    datasource = "pwchemistry", measurement = "TDN", unit = "mg N L^-1",
    depth = "matched to depth_cm", source_summary_stat = "none",
    temporal_match = "same plot and exact FT-ICR collection date",
    depth_match = "exact FT-ICR depth (10 or 30 cm)",
    notes = "Observed concentration from the matching sample; no interpolation, nearest-date matching, or temporal averaging. TDN includes organic and inorganic dissolved nitrogen."
  )
)

for (i in seq_len(nrow(hourly_meta))) {
  src <- hourly_meta[i]
  metadata_rows[[length(metadata_rows) + 1L]] <- metadata_row(
    output_column = hourly_output_cols[[i]],
    role = "covariate",
    analysis_group = src$analysis_group,
    source_column = src$source_column,
    source_file = basename(hourly_path),
    source_grain = "hourly",
    aggregation = "arithmetic mean within exact FT-ICR date",
    description = paste0("Mean of observed hourly ", src$source_column, " values on the FT-ICR collection date within plot."),
    datasource = src$datasource,
    measurement = src$measurement,
    unit = src$unit,
    depth = src$depth,
    source_summary_stat = src$source_summary_stat,
    source_summary_class = src$source_summary_class,
    temporal_match = "same plot and exact FT-ICR collection date",
    depth_match = "source-specific sensor depth; not selected by FT-ICR depth",
    notes = "No interpolation, nearest-date matching, or carry-forward."
  )
}

for (i in seq_len(nrow(daily_meta))) {
  src <- daily_meta[i]
  output_analysis_group <- if (identical(src$datasource, "lysimiter.surface.chemistry")) {
    "surface_lysimeter"
  } else {
    src$analysis_group
  }
  depth_match <- if (identical(src$datasource, "lysimiter.surface.chemistry")) {
    "separate surface lysimeter measurement; not the depth-matched FT-ICR sample"
  } else if (!is.na(src$depth) && nzchar(src$depth)) {
    paste0("source-specific depth ", src$depth, "; not selected by FT-ICR depth")
  } else {
    "plot-level measurement"
  }
  metadata_rows[[length(metadata_rows) + 1L]] <- metadata_row(
    output_column = daily_output_cols[[i]],
    role = "covariate",
    analysis_group = output_analysis_group,
    source_column = src$source_column,
    source_file = basename(daily_path),
    source_grain = "daily or event-level joined observation",
    aggregation = "none across dates; exact-date source value retained",
    description = paste0(src$source_column, " observed for the plot on the FT-ICR collection date."),
    datasource = src$datasource,
    measurement = src$measurement,
    unit = src$unit,
    depth = src$depth,
    source_summary_stat = src$source_summary_stat,
    source_summary_class = src$source_summary_class,
    temporal_match = "same plot and exact FT-ICR collection date",
    depth_match = depth_match,
    notes = paste(
      "No interpolation, nearest-date matching, carry-forward, or averaging across dates.",
      "Any mean or SD encoded in the source variable is the original within-event summary."
    )
  )
}

metadata <- rbindlist(metadata_rows, fill = TRUE)
metadata[, output_index := match(output_column, names(same_day))]
metadata[, output_type := vapply(
  output_column,
  function(column) paste(class(same_day[[column]]), collapse = "/"),
  character(1)
)]
metadata[, output_nonmissing_rows := vapply(
  output_column,
  function(column) sum(!is.na(same_day[[column]])),
  integer(1)
)]
metadata[, output_nonmissing_fraction := output_nonmissing_rows / nrow(same_day)]
metadata[, `:=`(
  fticr_date_start = as.character(min(same_day$date)),
  fticr_date_end = as.character(max(same_day$date)),
  output_grain = "one row per FT-ICR spectrum",
  join_key = "sample_num (also unique by plot + date + depth_cm)",
  missing_value = "blank field means no observation on the exact FT-ICR date"
)]
setorder(metadata, output_index)

if (!identical(metadata$output_column, names(same_day))) {
  stop("Output columns and metadata rows do not match exactly.")
}
if (nrow(same_day) != nrow(samples)) stop("The output row count changed during joins.")
if (anyDuplicated(same_day$sample_num)) stop("The output contains duplicate FT-ICR samples.")
if (any(grepl("__monthly", names(same_day), fixed = TRUE))) {
  stop("Monthly summaries leaked into the same-day dataset.")
}
allowed_analysis_groups <- c("key", "temperature", "vwc", "matched_porewater_cn", "surface_lysimeter")
unexpected_analysis_groups <- setdiff(unique(metadata$analysis_group), allowed_analysis_groups)
if (length(unexpected_analysis_groups) > 0L) {
  stop("Unexpected covariate groups in output: ", paste(unexpected_analysis_groups, collapse = ", "))
}
if (any(vapply(same_day[, setdiff(names(same_day), key_cols), with = FALSE], function(x) all(is.na(x)), logical(1)))) {
  stop("The output contains an entirely empty covariate column.")
}
if (any(names(fticr)[!names(fticr) %chin% required_fticr_cols] %chin% names(same_day))) {
  stop("FT-ICR response metrics leaked into the covariate output.")
}

message("Writing same-day covariate dataset and metadata...")
fwrite(same_day, out_path, na = "")
fwrite(metadata, metadata_out_path, na = "")

message("Wrote: ", out_path)
message("Wrote: ", metadata_out_path)
message("FT-ICR spectra: ", nrow(same_day))
message("Output columns: ", ncol(same_day))
message(
  "Exact-sample porewater NPOC/TDN matches: ",
  sum(complete.cases(same_day[, ..porewater_output_cols])),
  " / ",
  nrow(same_day)
)
