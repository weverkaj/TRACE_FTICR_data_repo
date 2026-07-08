#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

extract_first_int <- function(x) {
  x <- as.character(x)
  has_digit <- grepl("\\d", x)
  out <- rep(NA_integer_, length(x))
  out[has_digit] <- suppressWarnings(as.integer(sub(".*?(\\d+).*", "\\1", x[has_digit])))
  out
}

normalize_treatment <- function(x) {
  y <- tolower(trimws(as.character(x)))
  ifelse(grepl("warm", y), "warmed", ifelse(grepl("control", y), "control", NA_character_))
}

safe_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

script_dir <- get_script_dir()
trace_mcmc_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
docs_dir <- normalizePath(file.path(trace_mcmc_dir, ".."), winslash = "/", mustWork = FALSE)
out_dir <- file.path(trace_mcmc_dir, "fticr_integration", "output")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

chemodiv_path <- file.path(trace_mcmc_dir, "machine_learning", "fticr_chemodiversity_metrics.csv")
mol_path <- file.path(
  docs_dir,
  "trace_data_cleanup",
  "data_raw",
  "processed_data",
  "Merged_Output",
  "Test_Processed-Unprocessed_Mol.csv"
)
intensity_path <- file.path(
  docs_dir,
  "trace_data_cleanup",
  "data_raw",
  "processed_data",
  "Merged_Output",
  "Test_Processed-Unprocssed_Data.csv"
)
sample_key_path <- file.path(
  docs_dir,
  "trace_data_cleanup",
  "data_raw",
  "60398_Sihi_Porewater_July82024_DS.csv"
)

required_paths <- c(chemodiv_path, mol_path, intensity_path, sample_key_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required input files:\n", paste(missing_paths, collapse = "\n"))
}

message("Loading formula-count FTICR metrics...")
metrics <- fread(chemodiv_path)
metrics[, date := as.IDate(date)]
metrics[, treatment := normalize_treatment(treatment)]

message("Computing formula-count N-containing fraction...")
mol <- fread(mol_path)
formula_col <- "Molecular Formula"
if (!(formula_col %in% names(mol))) stop("Missing Molecular Formula column in molecular formula table.")

if ("Is Isotopologue" %in% names(mol)) {
  mol <- mol[safe_numeric(`Is Isotopologue`) == 0 | is.na(`Is Isotopologue`)]
}

element_cols <- c("C", "H", "O", "N", "S", "P")
missing_elements <- setdiff(element_cols, names(mol))
if (length(missing_elements) > 0L) {
  stop("Missing element columns in molecular formula table: ", paste(missing_elements, collapse = ", "))
}
for (nm in element_cols) {
  mol[, (nm) := safe_numeric(get(nm))]
}
mol <- mol[!is.na(C) & !is.na(H) & C > 0]
for (nm in c("O", "N", "S", "P")) {
  mol[is.na(get(nm)), (nm) := 0]
}
mol <- unique(mol, by = formula_col)
setkeyv(mol, formula_col)

intensity <- fread(intensity_path)
if (!(formula_col %in% names(intensity))) stop("Missing Molecular Formula column in intensity matrix.")
setnames(intensity, formula_col, "formula")

blank_cols <- grep("Blank", names(intensity), value = TRUE, ignore.case = TRUE)
field_cols <- grep("^Sihi_60398_[0-9]+_r[0-9]+_", names(intensity), value = TRUE)
if (length(field_cols) == 0L) stop("No field sample columns matched expected Sihi_60398 sample pattern.")

blank_formulas <- character()
if (length(blank_cols) > 0L) {
  blank_detected <- rowSums(!is.na(intensity[, ..blank_cols])) > 0
  blank_formulas <- intensity[blank_detected, formula]
}
valid_formula_set <- setdiff(mol[[formula_col]], blank_formulas)

sample_map <- data.table(sample_col = field_cols)
sample_map[, sample_num := suppressWarnings(as.integer(sub("^Sihi_60398_([0-9]+)_r[0-9]+_.*$", "\\1", sample_col)))]
sample_map <- sample_map[!is.na(sample_num)]

key_raw <- fread(sample_key_path, encoding = "Latin-1")
key_raw[, sample_num := extract_first_int(sample_name)]
sample_key <- key_raw[!is.na(sample_num), .(
  sample_num,
  date = as.IDate(collection_date),
  plot = extract_first_int(experimental_factor_other),
  treatment = normalize_treatment(fifelse(
    nzchar(trimws(as.character(experimental_factor))),
    trimws(as.character(experimental_factor)),
    trimws(as.character(climate_environment))
  )),
  depth = trimws(as.character(depth))
)]
sample_key[, depth_cm := fifelse(grepl("^0-10", depth), 10L, fifelse(grepl("^10-30", depth), 30L, NA_integer_))]
sample_key <- unique(sample_key[!is.na(date) & !is.na(plot) & !is.na(depth_cm)], by = "sample_num")

n_records <- vector("list", length(unique(sample_map$sample_num)))
sample_nums <- sort(unique(sample_map$sample_num))
for (i in seq_along(sample_nums)) {
  snum <- sample_nums[[i]]
  cols <- sample_map[sample_num == snum, sample_col]
  detected <- rowSums(!is.na(intensity[, ..cols])) > 0
  formulas <- intersect(intensity[detected, formula], valid_formula_set)
  if (length(formulas) < 5L) next

  meta_row <- sample_key[sample_num == snum]
  if (nrow(meta_row) == 0L) next

  props <- mol[J(formulas), nomatch = 0]
  if (nrow(props) < 5L) next

  n_records[[i]] <- data.table(
    sample_num = snum,
    plot = meta_row$plot[[1]],
    treatment = meta_row$treatment[[1]],
    date = meta_row$date[[1]],
    depth_cm = meta_row$depth_cm[[1]],
    n_containing_formulae = sum(props$N > 0, na.rm = TRUE),
    frac_n_containing = mean(props$N > 0, na.rm = TRUE)
  )
}
n_metrics <- rbindlist(n_records, fill = TRUE)

plot_dt <- merge(
  metrics,
  n_metrics,
  by = c("sample_num", "plot", "treatment", "date", "depth_cm"),
  all.x = TRUE
)

plot_dt[, depth_label := fifelse(depth_cm == 10L, "0-10 cm", fifelse(depth_cm == 30L, "10-30 cm", paste(depth_cm, "cm")))]
plot_dt[, treatment_label := fifelse(treatment == "control", "Control", fifelse(treatment == "warmed", "Warmed", treatment))]

metric_map <- data.table(
  metric = c(
    "richness",
    "mean_nosc",
    "mean_hc",
    "mean_oc",
    "mean_ai",
    "frac_lignin",
    "frac_n_containing"
  ),
  metric_label = c(
    "Formula richness\n(detected formulae)",
    "Mean NOSC\n(formula-count mean)",
    "Mean H/C\n(formula-count mean)",
    "Mean O/C\n(formula-count mean)",
    "Aromaticity\n(mean AI_mod)",
    "Lignin-like formulae\n(%)",
    "N-containing formulae\n(%)"
  ),
  multiplier = c(1, 1, 1, 1, 1, 100, 100),
  metric_order = seq_len(7L)
)

plot_dt[, (metric_map$metric) := lapply(.SD, as.numeric), .SDcols = metric_map$metric]

long_dt <- melt(
  plot_dt,
  id.vars = c("sample_num", "plot", "date", "treatment", "treatment_label", "depth_cm", "depth_label"),
  measure.vars = metric_map$metric,
  variable.name = "metric",
  value.name = "value"
)
long_dt <- merge(long_dt, metric_map, by = "metric", all.x = TRUE)
long_dt[, value_plot := value * multiplier]
long_dt <- long_dt[is.finite(value_plot)]
long_dt[, metric_label := factor(metric_label, levels = metric_map[order(metric_order), metric_label])]
long_dt[, depth_label := factor(depth_label, levels = c("0-10 cm", "10-30 cm"))]
long_dt[, treatment_label := factor(treatment_label, levels = c("Control", "Warmed"))]

summary_dt <- long_dt[
  ,
  .(
    n_plots = uniqueN(plot),
    mean_value = mean(value_plot, na.rm = TRUE),
    se_value = safe_se(value_plot)
  ),
  by = .(date, treatment_label, depth_label, metric_label)
]
summary_dt[, ymin := mean_value - se_value]
summary_dt[, ymax := mean_value + se_value]
summary_dt[!is.finite(ymin), `:=`(ymin = mean_value, ymax = mean_value)]

plot_data_out <- file.path(out_dir, "fticr_formula_count_temporal_metrics_plot_data.csv")
summary_data_out <- file.path(out_dir, "fticr_formula_count_temporal_metrics_summary.csv")
png_out <- file.path(fig_dir, "fticr_formula_count_temporal_metrics_by_treatment_depth.png")
pdf_out <- file.path(fig_dir, "fticr_formula_count_temporal_metrics_by_treatment_depth.pdf")

fwrite(long_dt[order(metric_order, depth_cm, treatment_label, date, plot)], plot_data_out)
fwrite(summary_dt[order(metric_label, depth_label, treatment_label, date)], summary_data_out)

treatment_colors <- c("Control" = "#2B6CB0", "Warmed" = "#C45A2A")
treatment_linetypes <- c("Control" = "solid", "Warmed" = "22")

p <- ggplot() +
  geom_point(
    data = long_dt,
    aes(x = date, y = value_plot, color = treatment_label),
    alpha = 0.22,
    size = 1.05,
    stroke = 0
  ) +
  geom_linerange(
    data = summary_dt,
    aes(x = date, ymin = ymin, ymax = ymax, color = treatment_label),
    alpha = 0.42,
    linewidth = 0.35
  ) +
  geom_line(
    data = summary_dt,
    aes(x = date, y = mean_value, color = treatment_label, linetype = treatment_label, group = treatment_label),
    linewidth = 0.72
  ) +
  geom_point(
    data = summary_dt,
    aes(x = date, y = mean_value, color = treatment_label),
    size = 1.55,
    stroke = 0
  ) +
  facet_grid(metric_label ~ depth_label, scales = "free_y", switch = "y") +
  scale_color_manual(values = treatment_colors, name = NULL) +
  scale_linetype_manual(values = treatment_linetypes, name = NULL) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.015, 0.035))
  ) +
  labs(
    title = "FT-ICR formula-count metrics through time by treatment and depth",
    subtitle = paste0(
      "Points are plot-level samples; lines connect observed treatment means by sample date; bars are +/- 1 SE. ",
      "Means and fractions use detected formula counts, not peak intensity."
    ),
    x = NULL,
    y = "Metric value",
    caption = paste0(
      "Source: fticr_chemodiversity_metrics.csv plus formula-count N-containing fraction recomputed from detected formulae. ",
      "No interpolation; n = ", uniqueN(long_dt$sample_num), " samples."
    )
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "#202124", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10.5, color = "#4A5568", margin = margin(b = 10)),
    plot.caption = element_text(size = 8.4, color = "#5F6368", hjust = 0, margin = margin(t = 8)),
    legend.position = "top",
    legend.justification = "left",
    legend.box.margin = margin(0, 0, 4, 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#E6E8EB", linewidth = 0.25),
    panel.grid.major.y = element_line(color = "#E6E8EB", linewidth = 0.25),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(color = "#3C4043"),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", color = "#202124", lineheight = 0.95),
    strip.text.x = element_text(face = "bold", color = "#202124"),
    panel.spacing.x = unit(0.7, "lines"),
    panel.spacing.y = unit(0.82, "lines"),
    plot.margin = margin(12, 18, 10, 12)
  )

ggsave(png_out, p, width = 13.5, height = 15.5, units = "in", dpi = 320, bg = "white")
ggsave(pdf_out, p, width = 13.5, height = 15.5, units = "in", bg = "white")

message("Wrote: ", plot_data_out)
message("Wrote: ", summary_data_out)
message("Wrote: ", png_out)
message("Wrote: ", pdf_out)
