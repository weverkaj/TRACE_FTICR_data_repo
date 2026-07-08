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

fit_and_summarize <- function(dt, response_col, scope, treatment_label = NA_character_, depth_label = NA_character_) {
  mod_dt <- dt[is.finite(get(response_col)) & is.finite(log2_lag) & lag_days > 0]
  if (nrow(mod_dt) < 3L || uniqueN(mod_dt$log2_lag) < 2L) {
    return(list(
      summary_text = paste("Not enough observations for", response_col, scope),
      coeff = data.table()
    ))
  }

  fit <- lm(mod_dt[[response_col]] ~ mod_dt$log2_lag)
  coef_tab <- as.data.table(coef(summary(fit)), keep.rownames = "term")
  setnames(coef_tab, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"), c("estimate", "std_error", "t_value", "p_value"))
  coef_tab[, `:=`(
    scope = scope,
    distance_metric = response_col,
    treatment_label = treatment_label,
    depth_label = depth_label,
    n_pairs = nrow(mod_dt),
    r_squared = summary(fit)$r.squared,
    adj_r_squared = summary(fit)$adj.r.squared,
    residual_standard_error = summary(fit)$sigma,
    df_residual = fit$df.residual
  )]
  setcolorder(coef_tab, c(
    "scope",
    "distance_metric",
    "treatment_label",
    "depth_label",
    "n_pairs",
    "term",
    "estimate",
    "std_error",
    "t_value",
    "p_value",
    "r_squared",
    "adj_r_squared",
    "residual_standard_error",
    "df_residual"
  ))

  header <- paste0(
    "\n\n",
    strrep("=", 90), "\n",
    "Scope: ", scope, "\n",
    "Distance metric: ", response_col, "\n",
    "Treatment: ", ifelse(is.na(treatment_label), "all", treatment_label), "\n",
    "Depth: ", ifelse(is.na(depth_label), "all", depth_label), "\n",
    "Model: ", response_col, " ~ log2_lag\n",
    "n pairwise observations: ", nrow(mod_dt), "\n",
    strrep("=", 90), "\n"
  )

  list(
    summary_text = paste0(header, paste(capture.output(summary(fit)), collapse = "\n"), "\n"),
    coeff = coef_tab
  )
}

script_dir <- get_script_dir()
trace_mcmc_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
out_dir <- file.path(trace_mcmc_dir, "fticr_integration", "output")

pairs_path <- file.path(out_dir, "fticr_q2_pairwise_temporal_distance_pairs.csv")
if (!file.exists(pairs_path)) {
  stop("Missing pairwise distance table: ", pairs_path)
}

pair_dt <- fread(pairs_path)
distance_cols <- c("formula_jaccard_distance", "vk_class_braycurtis", "metric_vector_distance")

results <- list()
idx <- 1L

for (response_col in distance_cols) {
  results[[idx]] <- fit_and_summarize(pair_dt, response_col, scope = "overall")
  idx <- idx + 1L

  for (dep in c("0-10 cm", "10-30 cm")) {
    results[[idx]] <- fit_and_summarize(
      pair_dt[depth_label == dep],
      response_col,
      scope = "by_depth",
      depth_label = dep
    )
    idx <- idx + 1L
  }

  for (dep in c("0-10 cm", "10-30 cm")) {
    for (trt in c("Control", "Warmed")) {
      results[[idx]] <- fit_and_summarize(
        pair_dt[depth_label == dep & treatment_label == trt],
        response_col,
        scope = "by_treatment_depth",
        treatment_label = trt,
        depth_label = dep
      )
      idx <- idx + 1L
    }
  }
}

summary_out <- file.path(out_dir, "fticr_q2_lm_model_summaries.txt")
coeff_out <- file.path(out_dir, "fticr_q2_lm_model_coefficients.csv")

summary_header <- paste(
  "TRACE FTICR pairwise temporal-distance linear model summaries",
  "",
  "Models are ordinary least-squares fits of distance ~ log2(lag_days).",
  "These summaries are descriptive; p-values assume independent pairwise observations.",
  "For inference, use the restricted-permutation results in fticr_q2_pairwise_temporal_distance_tests.csv.",
  sep = "\n"
)
writeLines(c(summary_header, vapply(results, `[[`, character(1), "summary_text")), summary_out)
fwrite(rbindlist(lapply(results, `[[`, "coeff"), fill = TRUE), coeff_out)

message("Wrote: ", summary_out)
message("Wrote: ", coeff_out)
