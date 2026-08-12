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
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

mode_nonmissing <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "-9999")] <- NA_character_
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_character_)
  counts <- sort(table(x), decreasing = TRUE)
  names(counts)[[1L]]
}

month_floor <- function(x) {
  as.IDate(format(as.IDate(x), "%Y-%m-01"))
}

normalize_missing_strings <- function(x) {
  if (!is.character(x)) return(x)
  y <- trimws(x)
  y[y %in% c("", "NA", "NaN", "-9999")] <- NA_character_
  y
}

coerce_numeric_like <- function(x, threshold = 0.95) {
  if (!is.character(x)) return(x)
  y <- normalize_missing_strings(x)
  observed <- !is.na(y)
  if (!any(observed)) return(y)
  numeric_y <- suppressWarnings(as.numeric(y))
  parse_fraction <- sum(is.finite(numeric_y[observed])) / sum(observed)
  if (is.finite(parse_fraction) && parse_fraction >= threshold) numeric_y else y
}

add_output_metadata <- function(rows, output_name, role, source_column = NA_character_,
                                source_file = NA_character_, source_grain = NA_character_,
                                aggregation = NA_character_, description = NA_character_,
                                unit_override = NA_character_, notes = NA_character_) {
  rows[[length(rows) + 1L]] <- data.table(
    column_name = output_name,
    role = role,
    source_column = source_column,
    source_file = source_file,
    source_grain = source_grain,
    aggregation = aggregation,
    description = description,
    unit_override = unit_override,
    notes = notes
  )
  rows
}

script_dir <- get_script_dir()
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
trace_mcmc_dir <- Sys.getenv("TRACE_MCMC_DIR", unset = "/Users/jrweverk/Documents/TRACE_MCMC")

hourly_path <- file.path(trace_mcmc_dir, "data", "joined_data_flux_final_2_hourly.csv")
daily_path <- file.path(trace_mcmc_dir, "data", "joined_data_flux_final_2_daily.csv")
source_metadata_path <- file.path(trace_mcmc_dir, "data", "joined_data_flux_final_2_daily_metadata.csv")
fticr_sample_path <- file.path(repo_dir, "data", "derived", "fticr_chemodiversity_metrics.csv")
out_dir <- file.path(repo_dir, "data", "derived", "monthly_covariates")
out_path <- file.path(out_dir, "fticr_monthly_covariates.csv")
metadata_out_path <- file.path(out_dir, "fticr_monthly_covariates_metadata.csv")

required_paths <- c(hourly_path, daily_path, source_metadata_path, fticr_sample_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required files:\n", paste(missing_paths, collapse = "\n"))
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start <- as.IDate(Sys.getenv("FTICR_COVARIATE_START", unset = "2017-05-01"))
analysis_end <- as.IDate(Sys.getenv("FTICR_COVARIATE_END", unset = "2024-04-30"))
if (is.na(analysis_start) || is.na(analysis_end) || analysis_start > analysis_end) {
  stop("Invalid FTICR covariate analysis dates.")
}

message("Loading joined-data metadata...")
source_meta <- fread(source_metadata_path, na.strings = c("", "NA", "-9999"))
required_meta_cols <- c(
  "column_name", "column_index", "observation_count", "datasource",
  "measurement_token", "unit_token", "depth", "summary_stat_token",
  "summary_stat_class"
)
missing_meta_cols <- setdiff(required_meta_cols, names(source_meta))
if (length(missing_meta_cols) > 0L) {
  stop("Source metadata is missing columns: ", paste(missing_meta_cols, collapse = ", "))
}

source_meta[, datasource := fifelse(is.na(datasource), "", datasource)]
source_meta[, column_name := as.character(column_name)]
setkey(source_meta, column_name)

id_cols <- c("timestamp", "plot", "day_hour", "year", "month", "date", "treatment", "season")
fticr_source_cols <- source_meta[datasource == "fticr", column_name]
date_like_source_cols <- source_meta[
  grepl("(^date([._]|$)|installed_date|date[._]time[._]installed)", column_name, ignore.case = TRUE) |
    measurement_token == "date" | unit_token %in% c("date", "installed"),
  column_name
]

high_frequency_cols <- c(
  "temp_cs655_cen_0_10cm",
  "vwc_cs655_cen_0_10cm",
  "temp_cs655_cen_20_30cm",
  "vwc_cs655_cen_20_30cm",
  "temp_cs655_cen_40_50cm",
  "vwc_cs655_cen_40_50cm"
)
monthly_flux_cols <- "flux"
hourly_summary_cols <- c(high_frequency_cols, monthly_flux_cols)
missing_hourly_summary_cols <- setdiff(hourly_summary_cols, source_meta$column_name)
if (length(missing_hourly_summary_cols) > 0L) {
  stop("Missing selected hourly summary columns: ", paste(missing_hourly_summary_cols, collapse = ", "))
}

porewater_datasources <- c("lysimiter.surface.chemistry", "pwchemistry")
excluded_ingrowth_core_datasources <- c("traceingrowthcores", "rootigcores")
core_datasources <- c("nutrient.cores", "rootsoilcores")
root_respiration_cols <- c(
  "q10_rootsoilcores",
  "r25_nmolCO2.g.s_rootsoilcores",
  "plot.vegtemp_c_rootsoilcores",
  "soil.mean.temp_c_surf_rootsoilcores",
  "soil.mean.vwc_cm3_cm3_surf_rootsoilcores",
  "soil.temp_c_deepcm20_rootsoilcores",
  "soil.temp_c_deepcm40_rootsoilcores",
  "soil.moist_c_deepcm20_rootsoilcores",
  "soil.moist_c_deepcm40_rootsoilcores",
  "resp.root.mass_g_rootsoilcores",
  "resp.dead.mass_g_rootsoilcores",
  "resp.rate.at.temp_g.m2_rootsoilcores",
  "err25_mmol.g.d_rootsoilcores",
  "err25_g.C.d_rootsoilcores"
)
missing_root_respiration_cols <- setdiff(root_respiration_cols, source_meta$column_name)
if (length(missing_root_respiration_cols) > 0L) {
  stop("Missing expected root-respiration columns: ", paste(missing_root_respiration_cols, collapse = ", "))
}
volume_cols <- source_meta[
  measurement_token == "volume" | grepl("(^|[._])volume([._]|$)", column_name, ignore.case = TRUE),
  column_name
]
remaining_cols <- source_meta[
  grepl("remaining", column_name, ignore.case = TRUE),
  column_name
]
root_morphology_exclusion_cols <- c(
  "crossings_count_rootsoilcores",
  "forks_count_rootsoilcores",
  "length_cm_rootsoilcores",
  "mean.diam_mm_rootsoilcores",
  "projected.area_cm2_rootsoilcores"
)
missing_root_morphology_cols <- setdiff(root_morphology_exclusion_cols, source_meta$column_name)
if (length(missing_root_morphology_cols) > 0L) {
  stop("Missing expected root-morphology columns: ", paste(missing_root_morphology_cols, collapse = ", "))
}
redundant_rootsoil_microbial_cols <- source_meta[
  datasource == "rootsoilcores" & grepl("ubial", column_name, ignore.case = TRUE),
  column_name
]
if (length(redundant_rootsoil_microbial_cols) == 0L) {
  stop("No redundant root-soil-core microbial biomass columns were identified.")
}
redundant_rootsoil_biogeochemistry_cols <- c(
  "soil.total_carbon_percent_rootsoilcores",
  "soil.total.nitrogen_percent_rootsoilcores",
  "soil.tc.tn_percent_rootsoilcores",
  "soil.NH4_ug.g_rootsoilcores",
  "soil.NO3_ug.g_rootsoilcores",
  "soil.inorganic.nitrogen_ug.g_rootsoilcores",
  "soil.NH4.min.rate_ug.gsoil.d.1_rootsoilcores",
  "soil.NO3.min.rate_ug.gsoil.d.1_rootsoilcores",
  "soil.N.min.rate_ug.gsoil.d.1_rootsoilcores",
  "soil.PO4_ug.g_rootsoilcores",
  "soil.ext_c_ug.g_rootsoilcores",
  "soil.ext.n_ug.g_rootsoilcores"
)
missing_rootsoil_bgc_cols <- setdiff(redundant_rootsoil_biogeochemistry_cols, source_meta$column_name)
if (length(missing_rootsoil_bgc_cols) > 0L) {
  stop("Missing expected redundant root-soil biogeochemistry columns: ", paste(missing_rootsoil_bgc_cols, collapse = ", "))
}
additional_exclusion_cols <- unique(c(
  volume_cols,
  remaining_cols,
  root_morphology_exclusion_cols,
  redundant_rootsoil_microbial_cols,
  redundant_rootsoil_biogeochemistry_cols
))
microbial_biomass_cols <- source_meta[
  grepl("microbial|ubial", column_name, ignore.case = TRUE) &
    !datasource %chin% excluded_ingrowth_core_datasources &
    !column_name %chin% redundant_rootsoil_microbial_cols,
  column_name
]

minirhizotron_meta <- source_meta[datasource == "minirhizotron"]
minirhizotron_meta[, depth_upper_cm := suppressWarnings(as.integer(
  sub("^.*-([0-9]+)cm$", "\\1", depth)
))]
selected_minirhizotron_cols <- c(
  "rootstock_g.m2_0-10_minirhizotron",
  "rootstock_g.m2_10-20_minirhizotron",
  "rootstock_g.m2_20-30_minirhizotron",
  "rootgrowth_g.m2.day_0-10_minirhizotron",
  "rootgrowth_g.m2.day_10-20_minirhizotron",
  "rootgrowth_g.m2.day_20-30_minirhizotron",
  "rootmortality_g.m2.day_0-10_minirhizotron",
  "rootmortality_g.m2.day_10-20_minirhizotron",
  "rootmortality_g.m2.day_20-30_minirhizotron"
)
missing_minirhizotron_cols <- setdiff(selected_minirhizotron_cols, minirhizotron_meta$column_name)
if (length(missing_minirhizotron_cols) > 0L) {
  stop("Missing selected minirhizotron columns: ", paste(missing_minirhizotron_cols, collapse = ", "))
}

month_sequence <- seq(month_floor(analysis_start), month_floor(analysis_end), by = "month")
monthly <- CJ(plot = 1:6, month_start = as.IDate(month_sequence), unique = TRUE)
monthly[, `:=`(
  year = as.integer(format(month_start, "%Y")),
  month = as.integer(format(month_start, "%m")),
  year_month = format(month_start, "%Y-%m"),
  treatment = fifelse(plot %in% c(1L, 3L, 5L), "control", "warmed")
)]
setcolorder(monthly, c("plot", "month_start", "year_month", "year", "month", "treatment"))

metadata_rows <- list()
metadata_rows <- add_output_metadata(
  metadata_rows, "plot", "key", description = "TRACE plot identifier (1-6).", unit_override = "identifier"
)
metadata_rows <- add_output_metadata(
  metadata_rows, "month_start", "key", description = "First day of the calendar month.", unit_override = "YYYY-MM-DD"
)
metadata_rows <- add_output_metadata(
  metadata_rows, "year_month", "key", description = "Calendar month join key for FT-ICR sample dates.", unit_override = "YYYY-MM"
)
metadata_rows <- add_output_metadata(
  metadata_rows, "year", "key", description = "Calendar year.", unit_override = "year"
)
metadata_rows <- add_output_metadata(
  metadata_rows, "month", "key", description = "Calendar month number (1-12).", unit_override = "month"
)
metadata_rows <- add_output_metadata(
  metadata_rows, "treatment", "key", description = "Experimental treatment derived from plot number.", unit_override = "control/warmed"
)

message("Reading the selected hourly temperature, moisture, and flux variables...")
hourly_select <- c("plot", "date", hourly_summary_cols)
hourly <- fread(
  hourly_path,
  select = hourly_select,
  na.strings = c("", "NA", "NaN", "-9999"),
  showProgress = TRUE
)
hourly[, `:=`(plot = as.integer(plot), date = as.IDate(date))]
hourly <- hourly[
  plot %in% 1:6 & date >= analysis_start & date <= analysis_end
]
hourly[, month_start := month_floor(date)]

hourly_source_file <- basename(hourly_path)
for (source_col in hourly_summary_cols) {
  hourly[, value_for_summary := suppressWarnings(as.numeric(get(source_col)))]
  agg <- hourly[
    ,
    .(monthly_mean = safe_mean(value_for_summary)),
    by = .(plot, month_start)
  ]
  output_col <- paste0(source_col, "__monthly_mean")
  setnames(agg, "monthly_mean", output_col)
  monthly <- merge(monthly, agg, by = c("plot", "month_start"), all.x = TRUE, sort = FALSE)

  metadata_rows <- add_output_metadata(
    metadata_rows,
    output_col,
    "covariate",
    source_column = source_col,
    source_file = hourly_source_file,
    source_grain = "hourly",
    aggregation = "mean",
    description = paste0("Calendar-month mean of hourly ", source_col, " within plot."),
    notes = fifelse(
      source_col == "flux",
      "Plot-level soil CO2 flux summarized from hourly observations; no interpolation or imputation.",
      "Central CS655 profile selected to provide one temperature or VWC series at each standard soil depth."
    )
  )
}
hourly[, value_for_summary := NULL]
rm(hourly)
invisible(gc())

message("Reading selected porewater, microbial biomass, nutrient/root-soil core, and minirhizotron covariates...")
daily_header <- names(fread(daily_path, nrows = 0L))
daily_candidate_cols <- unique(c(
  source_meta[datasource %chin% porewater_datasources, column_name],
  source_meta[datasource %chin% core_datasources, column_name],
  microbial_biomass_cols,
  selected_minirhizotron_cols
))
daily_candidate_cols <- setdiff(
  daily_candidate_cols,
  c(
    id_cols,
    date_like_source_cols,
    fticr_source_cols,
    hourly_summary_cols,
    root_respiration_cols,
    additional_exclusion_cols
  )
)
missing_daily_cols <- setdiff(daily_candidate_cols, daily_header)
if (length(missing_daily_cols) > 0L) {
  stop("Selected daily covariates are absent from the daily table: ", paste(missing_daily_cols, collapse = ", "))
}
daily_select <- unique(c("plot", "date", "season", daily_candidate_cols))
daily <- fread(
  daily_path,
  select = daily_select,
  na.strings = c("", "NA", "NaN", "-9999"),
  showProgress = FALSE
)
daily[, `:=`(plot = as.integer(plot), date = as.IDate(date))]
daily <- daily[
  plot %in% 1:6 & date >= analysis_start & date <= analysis_end
]
daily[, month_start := month_floor(date)]

season_monthly <- daily[, .(season = mode_nonmissing(season)), by = .(plot, month_start)]
monthly <- merge(monthly, season_monthly, by = c("plot", "month_start"), all.x = TRUE, sort = FALSE)
metadata_rows <- add_output_metadata(
  metadata_rows,
  "season",
  "key",
  source_column = "season",
  source_file = basename(daily_path),
  source_grain = "daily",
  aggregation = "mode",
  description = "Season assigned to the calendar month in the joined TRACE dataset.",
  unit_override = "category"
)

for (source_col in daily_candidate_cols) {
  daily[, (source_col) := coerce_numeric_like(get(source_col))]
}

is_date_col <- vapply(daily[, ..daily_candidate_cols], inherits, logical(1), what = "Date")
numeric_cols <- daily_candidate_cols[vapply(daily[, ..daily_candidate_cols], is.numeric, logical(1)) & !is_date_col]
categorical_cols <- daily_candidate_cols[vapply(daily[, ..daily_candidate_cols], is.character, logical(1)) & !is_date_col]

numeric_cols <- numeric_cols[vapply(daily[, ..numeric_cols], function(x) any(is.finite(x)), logical(1))]
categorical_cols <- categorical_cols[vapply(
  daily[, ..categorical_cols],
  function(x) any(!is.na(normalize_missing_strings(x))),
  logical(1)
)]

daily_source_file <- basename(daily_path)
for (source_col in numeric_cols) {
  daily[, value_for_summary := suppressWarnings(as.numeric(get(source_col)))]
  agg <- daily[
    ,
    .(monthly_mean = safe_mean(value_for_summary)),
    by = .(plot, month_start)
  ]
  output_col <- paste0(source_col, "__monthly_mean")
  setnames(agg, "monthly_mean", output_col)
  monthly <- merge(monthly, agg, by = c("plot", "month_start"), all.x = TRUE, sort = FALSE)

  metadata_rows <- add_output_metadata(
    metadata_rows,
    output_col,
    "covariate",
    source_column = source_col,
    source_file = daily_source_file,
    source_grain = "daily",
    aggregation = "mean",
    description = paste0("Calendar-month mean of non-missing daily ", source_col, " values within plot."),
    notes = "Observed source values only; no temporal interpolation or imputation. If a plot had multiple observation dates in a month, their arithmetic mean is reported. Existing source-level means and SDs retain their original meaning."
  )
}

for (source_col in categorical_cols) {
  daily[, value_for_summary := normalize_missing_strings(get(source_col))]
  agg <- daily[
    ,
    .(monthly_mode = mode_nonmissing(value_for_summary)),
    by = .(plot, month_start)
  ]
  output_col <- paste0(source_col, "__monthly_mode")
  setnames(agg, "monthly_mode", output_col)
  monthly <- merge(monthly, agg, by = c("plot", "month_start"), all.x = TRUE, sort = FALSE)

  metadata_rows <- add_output_metadata(
    metadata_rows,
    output_col,
    "covariate",
    source_column = source_col,
    source_file = daily_source_file,
    source_grain = "daily",
    aggregation = "mode",
    description = paste0("Most frequent non-missing daily category for ", source_col, " in the plot-month."),
    unit_override = "category",
    notes = "Observed source values only; no temporal interpolation or imputation."
  )
}
if ("value_for_summary" %in% names(daily)) daily[, value_for_summary := NULL]
rm(daily)
invisible(gc())

setorder(monthly, month_start, plot)
setcolorder(
  monthly,
  c("plot", "month_start", "year_month", "year", "month", "treatment", "season", setdiff(names(monthly), c(
    "plot", "month_start", "year_month", "year", "month", "treatment", "season"
  )))
)

if (anyDuplicated(monthly[, .(plot, month_start)])) stop("Duplicate plot-month rows were generated.")
expected_rows <- 6L * length(month_sequence)
if (nrow(monthly) != expected_rows) {
  stop("Expected ", expected_rows, " plot-month rows but generated ", nrow(monthly), ".")
}
if (any(grepl("_fticr($|__)", names(monthly), ignore.case = TRUE))) {
  stop("FT-ICR-derived source columns leaked into the covariate output.")
}

message("Building output metadata...")
metadata <- rbindlist(metadata_rows, fill = TRUE)
if (anyDuplicated(metadata$column_name)) {
  dup <- unique(metadata$column_name[duplicated(metadata$column_name)])
  stop("Duplicate metadata rows for output columns: ", paste(dup, collapse = ", "))
}

selected_source_cols <- unique(metadata[!is.na(source_column), source_column])
if (!setequal(metadata[source_grain == "hourly", source_column], hourly_summary_cols)) {
  stop("The output does not contain exactly the selected temperature, VWC, and flux columns.")
}
if (any(metadata$role == "coverage")) {
  stop("Per-variable coverage columns should not be present in the reduced output.")
}
if (any(source_meta[datasource == "apparency", column_name] %in% selected_source_cols)) {
  stop("Foliar apparency columns leaked into the reduced output.")
}
if (any(source_meta[datasource %chin% excluded_ingrowth_core_datasources, column_name] %in% selected_source_cols)) {
  stop("Ingrowth-core columns leaked into the reduced output.")
}
if (any(root_respiration_cols %in% selected_source_cols)) {
  stop("Root-respiration parameters leaked into the reduced output.")
}
if (any(additional_exclusion_cols %in% selected_source_cols)) {
  stop("An explicitly excluded volume, remaining-mass, root-morphology, or redundant microbial column leaked into the reduced output.")
}
retained_minirhizotron_cols <- intersect(minirhizotron_meta$column_name, selected_source_cols)
if (!setequal(retained_minirhizotron_cols, selected_minirhizotron_cols)) {
  stop("The output does not contain exactly the nine selected minirhizotron variables.")
}
if (!any(microbial_biomass_cols %in% selected_source_cols)) {
  stop("No microbial biomass columns were retained.")
}
if (!any(source_meta[datasource %chin% porewater_datasources, column_name] %in% selected_source_cols)) {
  stop("No porewater columns were retained.")
}
if (!any(source_meta[datasource %chin% core_datasources, column_name] %in% selected_source_cols)) {
  stop("No core columns were retained.")
}

source_backed_metadata <- metadata[!is.na(source_column)]
source_fields <- source_meta[
  source_backed_metadata,
  on = .(column_name = source_column),
  .(
    output_column = i.column_name,
    role = i.role,
    source_column = i.source_column,
    source_file = i.source_file,
    source_grain = i.source_grain,
    aggregation = i.aggregation,
    description = i.description,
    datasource = x.datasource,
    measurement = x.measurement_token,
    unit = fifelse(!is.na(i.unit_override), i.unit_override, x.unit_token),
    depth = x.depth,
    source_summary_stat = x.summary_stat_token,
    source_summary_class = x.summary_stat_class,
    source_observation_count_full_daily_table = x.observation_count,
    notes = i.notes
  )
]

key_rows <- metadata[is.na(source_column)]
if (nrow(key_rows) > 0L) {
  source_fields <- rbind(
    source_fields,
    key_rows[, .(
      output_column = column_name,
      role,
      source_column,
      source_file,
      source_grain,
      aggregation,
      description,
      datasource = NA_character_,
      measurement = NA_character_,
      unit = unit_override,
      depth = NA_character_,
      source_summary_stat = NA_character_,
      source_summary_class = NA_character_,
      source_observation_count_full_daily_table = NA_integer_,
      notes
    )],
    fill = TRUE
  )
}

source_fields[, analysis_group := fcase(
  role == "key", "key",
  source_column %chin% high_frequency_cols & grepl("^temp_", source_column), "temperature",
  source_column %chin% high_frequency_cols & grepl("^vwc_", source_column), "vwc",
  source_column %chin% monthly_flux_cols, "soil_co2_flux",
  source_column %chin% microbial_biomass_cols, "microbial_biomass",
  datasource %chin% porewater_datasources, "porewater",
  datasource %chin% core_datasources, "core",
  datasource == "minirhizotron", "minirhizotron_0_30cm",
  default = "other"
)]
setcolorder(source_fields, c("output_column", "role", "analysis_group"))

source_fields[, output_index := match(output_column, names(monthly))]
source_fields[, output_type := vapply(
  output_column,
  function(nm) paste(class(monthly[[nm]]), collapse = "/"),
  character(1)
)]
source_fields[, output_nonmissing_rows := vapply(
  output_column,
  function(nm) sum(!is.na(monthly[[nm]])),
  integer(1)
)]
source_fields[, output_nonmissing_fraction := output_nonmissing_rows / nrow(monthly)]
source_fields[, `:=`(
  analysis_start = as.character(analysis_start),
  analysis_end = as.character(analysis_end),
  output_grain = "one row per TRACE plot and calendar month",
  join_key = "plot + year_month",
  missing_value = "blank field in CSV"
)]
setorder(source_fields, output_index)

if (!identical(source_fields$output_column, names(monthly))) {
  missing_metadata <- setdiff(names(monthly), source_fields$output_column)
  extra_metadata <- setdiff(source_fields$output_column, names(monthly))
  stop(
    "Output/metadata column mismatch. Missing metadata: ", paste(missing_metadata, collapse = ", "),
    "; extra metadata: ", paste(extra_metadata, collapse = ", ")
  )
}

message("Checking joinability to FT-ICR samples...")
fticr <- fread(fticr_sample_path)
if (!all(c("plot", "date") %in% names(fticr))) {
  stop("FT-ICR sample table must contain plot and date columns.")
}
fticr[, `:=`(plot = as.integer(plot), date = as.IDate(date))]
fticr[, year_month := format(date, "%Y-%m")]
fticr_in_range <- fticr[date >= analysis_start & date <= analysis_end]
matched <- monthly[fticr_in_range, on = .(plot, year_month)]
if (nrow(matched) != nrow(fticr_in_range) || any(is.na(matched$month_start))) {
  stop("Not every in-range FT-ICR sample matched exactly one plot-month covariate row.")
}

message("Writing monthly covariate dataset and metadata...")
fwrite(monthly, out_path, na = "")
fwrite(source_fields, metadata_out_path, na = "")

message("Wrote: ", out_path)
message("Wrote: ", metadata_out_path)
message("Plot-month rows: ", nrow(monthly))
message("Output columns: ", ncol(monthly))
message("Selected source variables: ", length(selected_source_cols))
message("FT-ICR samples matched: ", nrow(fticr_in_range), " / ", nrow(fticr_in_range))
