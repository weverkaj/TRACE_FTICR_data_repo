suppressPackageStartupMessages(library(data.table))

get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)))
  }
  getwd()
}

safe_date <- function(x) {
  out <- suppressWarnings(as.IDate(x))
  out
}

normalize_treatment <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- ifelse(grepl("warm", y), "warmed", ifelse(grepl("control", y), "control", NA_character_))
  out
}

extract_first_int <- function(x) {
  x <- as.character(x)
  has_digit <- grepl("\\d", x)
  out <- rep(NA_integer_, length(x))
  out[has_digit] <- suppressWarnings(as.integer(sub(".*?(\\d+).*", "\\1", x[has_digit])))
  out
}

safe_sd <- function(x) {
  if (length(x) < 2L) {
    return(NA_real_)
  }
  stats::sd(x, na.rm = TRUE)
}

safe_wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & (w > 0)
  if (!any(ok)) {
    return(NA_real_)
  }
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted_mean_cols <- function(values, w_mat) {
  valid <- is.finite(values)
  if (!any(valid)) {
    return(rep(NA_real_, ncol(w_mat)))
  }

  num <- as.numeric(crossprod(values[valid], w_mat[valid, , drop = FALSE]))
  den <- colSums(w_mat[valid, , drop = FALSE], na.rm = TRUE)
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}

script_dir <- get_script_dir()
trace_mcmc_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
docs_dir <- normalizePath(file.path(trace_mcmc_dir, ".."), winslash = "/", mustWork = FALSE)
out_dir <- file.path(trace_mcmc_dir, "fticr_integration", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

merged_intensity_path <- file.path(
  docs_dir,
  "trace_data_cleanup",
  "data_raw",
  "processed_data",
  "Merged_Output",
  "Test_Processed-Unprocssed_Data.csv"
)

merged_mol_path <- file.path(
  docs_dir,
  "trace_data_cleanup",
  "data_raw",
  "processed_data",
  "Merged_Output",
  "Test_Processed-Unprocessed_Mol.csv"
)

thermo_path <- file.path(
  docs_dir,
  "trace_data_cleanup",
  "code",
  "trace_fticr_all.csv"
)

sample_key_path <- file.path(
  docs_dir,
  "trace_data_cleanup",
  "data_raw",
  "60398_Sihi_Porewater_July82024_DS.csv"
)

required_paths <- c(
  merged_intensity_path,
  merged_mol_path,
  thermo_path,
  sample_key_path
)

missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing required input files:\n", paste(missing_paths, collapse = "\n"))
}

message("Loading merged molecular formula table...")
mol_dt <- fread(
  merged_mol_path,
  select = c("Molecular Formula", "C", "H", "N", "O", "P", "S")
)

message("Loading thermodynamic feature table...")
thermo_dt <- fread(
  thermo_path,
  select = c("delGcat", "lambda", "CUE", "nosc", "ne")
)

if (nrow(mol_dt) != nrow(thermo_dt)) {
  stop("Row mismatch: molecular table (", nrow(mol_dt), ") vs thermo table (", nrow(thermo_dt), ").")
}

for (col in c("C", "H", "N", "O", "P", "S")) {
  set(mol_dt, i = which(is.na(mol_dt[[col]])), j = col, value = 0)
}

message("Loading merged intensity matrix...")
int_dt <- fread(merged_intensity_path)
if (!("Molecular Formula" %in% names(int_dt))) {
  stop("Expected `Molecular Formula` column in merged intensity matrix.")
}
if (nrow(int_dt) != nrow(mol_dt)) {
  stop("Row mismatch: intensity matrix (", nrow(int_dt), ") vs molecular table (", nrow(mol_dt), ").")
}

all_signal_cols <- setdiff(names(int_dt), "Molecular Formula")
field_cols <- all_signal_cols[grepl("^Sihi_60398_[0-9]{3}_r[0-9]+_", all_signal_cols)]
qc_cols <- setdiff(all_signal_cols, field_cols)

if (length(field_cols) == 0L) {
  stop("No field sample columns matched expected naming pattern in merged intensity matrix.")
}

W <- as.matrix(int_dt[, ..field_cols])
storage.mode(W) <- "double"
W[is.na(W)] <- 0
rm(int_dt)
invisible(gc())

message("Computing run-level FTICR feature summaries...")
C_num <- mol_dt$C
H_num <- mol_dt$H
N_num <- mol_dt$N
O_num <- mol_dt$O
P_num <- mol_dt$P
S_num <- mol_dt$S

h_to_c <- ifelse(C_num > 0, H_num / C_num, NA_real_)
o_to_c <- ifelse(C_num > 0, O_num / C_num, NA_real_)

is_cho <- as.numeric((N_num == 0) & (P_num == 0) & (S_num == 0))
is_chon <- as.numeric((N_num > 0) & (P_num == 0) & (S_num == 0))
is_chos <- as.numeric((N_num == 0) & (P_num == 0) & (S_num > 0))
is_chop <- as.numeric((N_num == 0) & (P_num > 0) & (S_num == 0))
contains_n <- as.numeric(N_num > 0)
contains_p <- as.numeric(P_num > 0)
contains_s <- as.numeric(S_num > 0)

run_dt <- data.table(
  sample_col = field_cols,
  total_peak_area = colSums(W, na.rm = TRUE),
  formula_detected_n = colSums(W > 0, na.rm = TRUE),
  delGcat_wmean = weighted_mean_cols(thermo_dt$delGcat, W),
  lambda_wmean = weighted_mean_cols(thermo_dt$lambda, W),
  cue_wmean = weighted_mean_cols(thermo_dt$CUE, W),
  nosc_wmean = weighted_mean_cols(thermo_dt$nosc, W),
  ne_wmean = weighted_mean_cols(thermo_dt$ne, W),
  h_to_c_wmean = weighted_mean_cols(h_to_c, W),
  o_to_c_wmean = weighted_mean_cols(o_to_c, W),
  frac_CHO = weighted_mean_cols(is_cho, W),
  frac_CHON = weighted_mean_cols(is_chon, W),
  frac_CHOS = weighted_mean_cols(is_chos, W),
  frac_CHOP = weighted_mean_cols(is_chop, W),
  frac_N_containing = weighted_mean_cols(contains_n, W),
  frac_P_containing = weighted_mean_cols(contains_p, W),
  frac_S_containing = weighted_mean_cols(contains_s, W)
)

run_dt[, sample_num := as.integer(sub("^Sihi_60398_([0-9]{3})_.*$", "\\1", sample_col))]
run_dt[, technical_rep := as.integer(sub("^.*_r([0-9]+)_.*$", "\\1", sample_col))]
run_dt[, run_date_chr := sub("^.*_Fir_([0-9]{2}[A-Za-z]{3}[0-9]{2})_.*$", "\\1", sample_col)]
run_dt[, run_date := as.IDate(run_date_chr, format = "%d%b%y")]
run_dt[, tray_position := as.integer(sub("^.*_p[0-9]{2}_([0-9]+)_.*$", "\\1", sample_col))]

message("Parsing sample key and attaching field metadata...")
key_raw <- as.data.table(
  read.csv(
    sample_key_path,
    fileEncoding = "latin1",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
)

key_raw[, sample_num := extract_first_int(sample_name)]
key_dt <- key_raw[!is.na(sample_num), .(
  sample_num,
  sample_name_raw = as.character(sample_name),
  sample_is_extra = grepl("EXTRA", sample_name, ignore.case = TRUE),
  sample_date = safe_date(collection_date),
  depth = trimws(as.character(depth)),
  treatment_raw = fifelse(
    nzchar(trimws(as.character(experimental_factor))),
    trimws(as.character(experimental_factor)),
    trimws(as.character(climate_environment))
  ),
  plot = extract_first_int(experimental_factor_other)
)]
key_dt[, treatment := normalize_treatment(treatment_raw)]
key_dt <- unique(key_dt, by = "sample_num")

run_out <- merge(run_dt, key_dt, by = "sample_num", all.x = TRUE)

sample_out <- run_out[, .(
  sample_run_n = .N,
  sample_first_run = if (all(is.na(run_date))) as.IDate(NA) else min(run_date, na.rm = TRUE),
  sample_last_run = if (all(is.na(run_date))) as.IDate(NA) else max(run_date, na.rm = TRUE),
  total_peak_area_mean = mean(total_peak_area, na.rm = TRUE),
  total_peak_area_sd = safe_sd(total_peak_area),
  formula_detected_n_mean = mean(formula_detected_n, na.rm = TRUE),
  formula_detected_n_sd = safe_sd(formula_detected_n),
  delGcat_wmean = mean(delGcat_wmean, na.rm = TRUE),
  delGcat_wmean_sd = safe_sd(delGcat_wmean),
  lambda_wmean = mean(lambda_wmean, na.rm = TRUE),
  lambda_wmean_sd = safe_sd(lambda_wmean),
  cue_wmean = mean(cue_wmean, na.rm = TRUE),
  cue_wmean_sd = safe_sd(cue_wmean),
  nosc_wmean = mean(nosc_wmean, na.rm = TRUE),
  nosc_wmean_sd = safe_sd(nosc_wmean),
  ne_wmean = mean(ne_wmean, na.rm = TRUE),
  h_to_c_wmean = mean(h_to_c_wmean, na.rm = TRUE),
  o_to_c_wmean = mean(o_to_c_wmean, na.rm = TRUE),
  frac_CHO = mean(frac_CHO, na.rm = TRUE),
  frac_CHON = mean(frac_CHON, na.rm = TRUE),
  frac_CHOS = mean(frac_CHOS, na.rm = TRUE),
  frac_CHOP = mean(frac_CHOP, na.rm = TRUE),
  frac_N_containing = mean(frac_N_containing, na.rm = TRUE),
  frac_P_containing = mean(frac_P_containing, na.rm = TRUE),
  frac_S_containing = mean(frac_S_containing, na.rm = TRUE)
), by = .(
  sample_num,
  sample_name_raw,
  sample_is_extra,
  sample_date,
  plot,
  treatment,
  treatment_raw,
  depth
)]

surface_sample_out <- sample_out[grepl("^0-10", depth)]
subsurface_sample_out <- sample_out[grepl("^10-30", depth)]

run_out_path <- file.path(out_dir, "fticr_run_features.csv")
sample_out_path <- file.path(out_dir, "fticr_sample_features.csv")
surface_out_path <- file.path(out_dir, "fticr_sample_features_0_10cm.csv")
subsurface_out_path <- file.path(out_dir, "fticr_sample_features_10_30cm.csv")
summary_out_path <- file.path(out_dir, "fticr_feature_build_summary.txt")

fwrite(run_out, run_out_path)
fwrite(sample_out, sample_out_path)
fwrite(surface_sample_out, surface_out_path)
fwrite(subsurface_sample_out, subsurface_out_path)

summary_lines <- c(
  paste("Field run columns in merged matrix:", length(field_cols)),
  paste("QC/blank columns excluded:", length(qc_cols)),
  paste("Unique field samples in key output:", nrow(sample_out)),
  paste("0-10 cm sample rows:", nrow(surface_sample_out)),
  paste("10-30 cm sample rows:", nrow(subsurface_sample_out)),
  paste("Sample dates span:", as.character(min(sample_out$sample_date, na.rm = TRUE)), "to", as.character(max(sample_out$sample_date, na.rm = TRUE))),
  paste("Plots represented:", paste(sort(unique(sample_out$plot[!is.na(sample_out$plot)])), collapse = ", "))
)
writeLines(summary_lines, summary_out_path)

message("FTICR feature build complete.")
message("Wrote: ", run_out_path)
message("Wrote: ", sample_out_path)
message("Wrote: ", surface_out_path)
message("Wrote: ", subsurface_out_path)
message("Wrote: ", summary_out_path)
