#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(emmeans)
  library(ggplot2)
  library(nlme)
  library(scales)
})

get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)))
  }
  getwd()
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

normalize_treatment <- function(x) {
  y <- tolower(trimws(as.character(x)))
  ifelse(grepl("warm", y), "Warmed", ifelse(grepl("control", y), "Control", NA_character_))
}

depth_display <- function(x) {
  fifelse(as.integer(x) == 10L, "10cm", fifelse(as.integer(x) == 30L, "30cm", NA_character_))
}

block_from_plot <- function(x) {
  fifelse(
    as.integer(x) %in% c(1L, 2L),
    "A",
    fifelse(as.integer(x) %in% c(3L, 4L), "B", fifelse(as.integer(x) %in% c(5L, 6L), "C", NA_character_))
  )
}

read_fticr_csv <- function(path, ...) {
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(fread(path, ...))
  }
  if (requireNamespace("R.utils", quietly = TRUE)) {
    return(fread(path, ...))
  }
  cmd <- paste("gzip -dc", shQuote(path))
  fread(cmd = cmd, ...)
}

first_existing_path <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) {
    stop("None of the expected input files exists:\n", paste(paths, collapse = "\n"))
  }
  hit[[1]]
}

clamp_probability <- function(x, eps = 1e-8) {
  pmin(pmax(x, eps), 1 - eps)
}

metric_to_model_scale <- function(metric, value, richness) {
  if (metric == "richness") {
    return(log(value))
  }
  if (grepl("^vk_frac_", metric)) {
    proportion <- clamp_probability(value / 100)
    count <- proportion * richness
    return(qlogis((count + 0.5) / (richness + 1)))
  }
  value
}

metric_from_model_scale <- function(metric, value) {
  if (metric == "richness") {
    return(exp(value))
  }
  if (grepl("^vk_frac_", metric)) {
    return(100 * plogis(value))
  }
  value
}

align_model_matrix <- function(fit, newdata) {
  x <- model.matrix(delete.response(terms(fit)), newdata)
  beta_names <- names(fixed.effects(fit))
  missing_cols <- setdiff(beta_names, colnames(x))
  if (length(missing_cols) > 0L) {
    x <- cbind(x, matrix(0, nrow = nrow(x), ncol = length(missing_cols), dimnames = list(NULL, missing_cols)))
  }
  x[, beta_names, drop = FALSE]
}

extract_car_phi <- function(fit) {
  if (is.null(fit$modelStruct$corStruct)) {
    return(NA_real_)
  }
  as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
}

fit_lme_with_fallback <- function(d) {
  ctrl <- lmeControl(
    opt = "optim",
    msMaxIter = 300,
    maxIter = 200,
    tolerance = 1e-6,
    returnObject = TRUE
  )

  fit_car <- try(
    lme(
      y ~ time_year * treatment_label * depth_label,
      random = list(
        plot = pdDiag(~1),
        series_id = pdDiag(~time_year)
      ),
      correlation = corCAR1(value = 0.4, form = ~time30 | plot/series_id),
      data = d,
      method = "REML",
      na.action = na.omit,
      control = ctrl
    ),
    silent = TRUE
  )
  if (!inherits(fit_car, "try-error")) {
    return(list(fit = fit_car, fit_method = "lme_random_plot_and_plot_depth_slope_car1", fit_note = "ok"))
  }

  fit_no_car <- try(
    lme(
      y ~ time_year * treatment_label * depth_label,
      random = list(
        plot = pdDiag(~1),
        series_id = pdDiag(~time_year)
      ),
      data = d,
      method = "REML",
      na.action = na.omit,
      control = ctrl
    ),
    silent = TRUE
  )
  if (!inherits(fit_no_car, "try-error")) {
    return(list(
      fit = fit_no_car,
      fit_method = "lme_random_plot_and_plot_depth_slope_no_car1",
      fit_note = paste("CAR1 failed:", as.character(fit_car))
    ))
  }

  fit_intercept_car <- try(
    lme(
      y ~ time_year * treatment_label * depth_label,
      random = ~1 | plot/depth_label,
      correlation = corCAR1(value = 0.4, form = ~time30 | plot/depth_label),
      data = d,
      method = "REML",
      na.action = na.omit,
      control = ctrl
    ),
    silent = TRUE
  )
  if (!inherits(fit_intercept_car, "try-error")) {
    return(list(
      fit = fit_intercept_car,
      fit_method = "lme_random_plot_depth_intercept_car1_fallback",
      fit_note = paste(
        "Random-slope CAR1 failed:", as.character(fit_car),
        "| random-slope no CAR1 failed:", as.character(fit_no_car)
      )
    ))
  }

  list(
    fit = NULL,
    fit_method = "failed",
    fit_note = paste(
      "Random-slope CAR1:", as.character(fit_car),
      "| random-slope no CAR1:", as.character(fit_no_car),
      "| random-intercept CAR1:", as.character(fit_intercept_car)
    )
  )
}

natural_annual_change <- function(fit, metric, treatment_value, depth_value, time_mid) {
  new0 <- data.frame(
    time_year = time_mid,
    treatment_label = factor(treatment_value, levels = c("Control", "Warmed")),
    depth_label = factor(depth_value, levels = c("10cm", "30cm")),
    plot = factor(levels(fit$data$plot)[[1]], levels = levels(fit$data$plot)),
    time30 = 0
  )
  new1 <- new0
  new1$time_year <- new1$time_year + 1

  x0 <- align_model_matrix(fit, new0)[1, ]
  x1 <- align_model_matrix(fit, new1)[1, ]
  xs <- x1 - x0
  beta <- fixed.effects(fit)
  v <- vcov(fit)
  eta <- sum(x0 * beta)
  slope <- sum(xs * beta)

  if (metric == "richness") {
    estimate <- 100 * (exp(slope) - 1)
    gradient <- 100 * exp(slope) * xs
    units <- "percent change per year"
  } else if (grepl("^vk_frac_", metric)) {
    p0 <- plogis(eta)
    p1 <- plogis(eta + slope)
    estimate <- 100 * (p1 - p0)
    d_eta <- 100 * (p1 * (1 - p1) - p0 * (1 - p0))
    d_slope <- 100 * p1 * (1 - p1)
    gradient <- d_eta * x0 + d_slope * xs
    units <- "percentage points per year at group midpoint"
  } else {
    estimate <- slope
    gradient <- xs
    units <- "metric units per year"
  }

  se <- sqrt(as.numeric(t(gradient) %*% v %*% gradient))
  data.table(
    annual_change_natural = estimate,
    annual_change_natural_se = se,
    annual_change_natural_lower95 = estimate - qnorm(0.975) * se,
    annual_change_natural_upper95 = estimate + qnorm(0.975) * se,
    annual_change_units = units,
    group_midpoint_time_year = time_mid
  )
}

exact_lme_treatment_slope_test <- function(d, observed_fit) {
  swap_grid <- CJ(swap_a = c(FALSE, TRUE), swap_b = c(FALSE, TRUE), swap_c = c(FALSE, TRUE))
  swap_grid[, permutation_id := .I]
  null_rows <- list()

  for (perm_i in seq_len(nrow(swap_grid))) {
    row <- swap_grid[perm_i]
    swap_blocks <- c("A", "B", "C")[as.logical(unlist(row[, .(swap_a, swap_b, swap_c)]))]
    if (perm_i == 1L) {
      fit <- observed_fit
      fit_method <- "observed_fit"
    } else {
      perm_d <- copy(d)
      perm_d[
        ,
        treatment_label := factor(
          permute_treatment(as.character(treatment_label), block, swap_blocks),
          levels = c("Control", "Warmed")
        )
      ]
      perm_fit_result <- fit_lme_with_fallback(perm_d)
      fit <- perm_fit_result$fit
      fit_method <- perm_fit_result$fit_method
    }

    if (is.null(fit)) {
      null_rows[[perm_i]] <- data.table(
        permutation_id = row$permutation_id,
        swap_blocks = paste(swap_blocks, collapse = ","),
        depth_label = c("10cm", "30cm"),
        warmed_minus_control_slope = NA_real_,
        fit_method = fit_method
      )
      next
    }

    trend_emm <- emtrends(fit, ~treatment_label * depth_label, var = "time_year")
    contrast_dt <- as.data.table(summary(pairs(trend_emm, by = "depth_label", reverse = TRUE)))
    null_rows[[perm_i]] <- contrast_dt[
      ,
      .(
        permutation_id = row$permutation_id,
        swap_blocks = paste(swap_blocks, collapse = ","),
        depth_label = as.character(depth_label),
        warmed_minus_control_slope = estimate,
        fit_method = fit_method
      )
    ]
  }

  null_dt <- rbindlist(null_rows, fill = TRUE)
  observed <- null_dt[permutation_id == 1L]
  tests <- observed[
    ,
    {
      null <- null_dt[depth_label == .BY$depth_label, warmed_minus_control_slope]
      obs <- warmed_minus_control_slope
      list(
        observed_warmed_minus_control_slope = obs,
        n_exact_permutations = sum(is.finite(null)),
        exact_p_warmed_slope_lower = mean(null <= obs + 1e-12, na.rm = TRUE),
        exact_p_warmed_slope_higher = mean(null >= obs - 1e-12, na.rm = TRUE),
        exact_p_two_sided = mean(abs(null) >= abs(obs) - 1e-12, na.rm = TRUE),
        exact_null_min = min(null, na.rm = TRUE),
        exact_null_median = median(null, na.rm = TRUE),
        exact_null_max = max(null, na.rm = TRUE)
      )
    },
    by = depth_label
  ]

  list(null_distribution = null_dt, tests = tests)
}

fit_metric_models <- function(metric_dt, metric_info, analysis_window, window_start = as.IDate(NA)) {
  source_dt <- copy(metric_dt)
  if (!is.na(window_start)) {
    source_dt <- source_dt[date >= window_start]
  }

  slope_rows <- list()
  contrast_rows <- list()
  anova_rows <- list()
  coefficient_rows <- list()
  diagnostic_rows <- list()
  prediction_rows <- list()
  exact_test_rows <- list()
  exact_null_rows <- list()
  summary_lines <- character()
  out_i <- 1L

  for (metric_name in metric_info$metric) {
    info <- metric_info[metric == metric_name]
    d <- copy(source_dt[metric == metric_name & is.finite(value_plot) & is.finite(richness_count)])
    d[, treatment_label := factor(treatment_label, levels = c("Control", "Warmed"))]
    d[, depth_label := factor(depth_label, levels = c("10cm", "30cm"))]
    d[, plot := factor(plot)]
    d[, series_id := factor(interaction(plot, depth_label, drop = TRUE))]
    d[, block := block_from_plot(as.integer(as.character(plot)))]
    d[, time_year := as.numeric(date - as.IDate("2021-01-01")) / 365.25]
    d[, time30 := as.numeric(date - min(date)) / 30]
    d[, y := metric_to_model_scale(metric_name, value_plot, richness_count)]
    d <- d[is.finite(y)]
    setorder(d, plot, depth_label, date)

    fit_result <- fit_lme_with_fallback(d)
    fit <- fit_result$fit
    diagnostic_rows[[out_i]] <- data.table(
      analysis_window = analysis_window,
      metric = metric_name,
      metric_label = info$metric_label,
      n_samples = nrow(d),
      n_plots = uniqueN(d$plot),
      n_plot_depth_series = uniqueN(interaction(d$plot, d$depth_label)),
      first_date = min(d$date),
      last_date = max(d$date),
      fit_method = fit_result$fit_method,
      fit_note = fit_result$fit_note,
      car1_phi_per_30d = if (!is.null(fit)) extract_car_phi(fit) else NA_real_,
      aic = if (!is.null(fit)) AIC(fit) else NA_real_,
      bic = if (!is.null(fit)) BIC(fit) else NA_real_
    )

    if (is.null(fit)) {
      out_i <- out_i + 1L
      next
    }

    trend_emm <- emtrends(fit, ~treatment_label * depth_label, var = "time_year")
    trend_dt <- as.data.table(summary(trend_emm, infer = c(TRUE, TRUE)))
    trend_col <- grep("\\.trend$", names(trend_dt), value = TRUE)[[1]]
    setnames(
      trend_dt,
      c(trend_col, "lower.CL", "upper.CL"),
      c("model_scale_slope_per_year", "model_scale_lower95", "model_scale_upper95")
    )
    trend_dt[, `:=`(
      analysis_window = analysis_window,
      metric = metric_name,
      metric_label = info$metric_label,
      metric_group = info$metric_group,
      metric_order = info$metric_order,
      fit_method = fit_result$fit_method,
      n_samples = nrow(d),
      car1_phi_per_30d = extract_car_phi(fit)
    )]

    natural_rows <- trend_dt[
      ,
      {
        group_d <- d[treatment_label == .BY$treatment_label & depth_label == .BY$depth_label]
        midpoint <- mean(range(group_d$time_year))
        natural_annual_change(
          fit,
          metric_name,
          as.character(.BY$treatment_label),
          as.character(.BY$depth_label),
          midpoint
        )
      },
      by = .(treatment_label, depth_label)
    ]
    trend_dt <- merge(trend_dt, natural_rows, by = c("treatment_label", "depth_label"), all.x = TRUE)
    slope_rows[[out_i]] <- trend_dt

    treatment_contrasts <- as.data.table(summary(pairs(trend_emm, by = "depth_label", reverse = TRUE), infer = c(TRUE, TRUE)))
    setnames(treatment_contrasts, c("estimate", "lower.CL", "upper.CL"), c(
      "warmed_minus_control_model_slope",
      "lower95",
      "upper95"
    ))
    treatment_contrasts[, `:=`(
      analysis_window = analysis_window,
      metric = metric_name,
      metric_label = info$metric_label,
      contrast_type = "Warmed minus Control time slope",
      fit_method = fit_result$fit_method
    )]

    exact_test <- exact_lme_treatment_slope_test(d, fit)
    exact_test$tests[, `:=`(
      analysis_window = analysis_window,
      metric = metric_name,
      metric_label = info$metric_label
    )]
    exact_test$null_distribution[, `:=`(
      analysis_window = analysis_window,
      metric = metric_name,
      metric_label = info$metric_label
    )]
    exact_test_rows[[out_i]] <- exact_test$tests
    exact_null_rows[[out_i]] <- exact_test$null_distribution

    depth_contrasts <- as.data.table(summary(pairs(trend_emm, by = "treatment_label", reverse = TRUE), infer = c(TRUE, TRUE)))
    setnames(depth_contrasts, c("estimate", "lower.CL", "upper.CL"), c(
      "depth30_minus_depth10_model_slope",
      "lower95",
      "upper95"
    ))
    depth_contrasts[, `:=`(
      analysis_window = analysis_window,
      metric = metric_name,
      metric_label = info$metric_label,
      contrast_type = "30cm minus 10cm time slope",
      fit_method = fit_result$fit_method
    )]
    contrast_rows[[out_i]] <- rbindlist(list(treatment_contrasts, depth_contrasts), fill = TRUE)

    fit_anova <- as.data.table(anova(fit, type = "marginal"), keep.rownames = "term")
    fit_anova[, `:=`(
      analysis_window = analysis_window,
      metric = metric_name,
      metric_label = info$metric_label
    )]
    anova_rows[[out_i]] <- fit_anova

    coef_dt <- as.data.table(summary(fit)$tTable, keep.rownames = "term")
    coef_dt[, `:=`(
      analysis_window = analysis_window,
      metric = metric_name,
      metric_label = info$metric_label
    )]
    coefficient_rows[[out_i]] <- coef_dt

    summary_lines <- c(
      summary_lines,
      paste0("\n\nMETRIC: ", info$metric_label, " | WINDOW: ", analysis_window, "\n"),
      capture.output(summary(fit)),
      "\nMarginal fixed-effect tests:\n",
      capture.output(anova(fit, type = "marginal"))
    )

    if (analysis_window == "All available") {
      group_ranges <- d[, .(date_min = min(date), date_max = max(date)), by = .(treatment_label, depth_label)]
      grid_parts <- lapply(seq_len(nrow(group_ranges)), function(i) {
        g <- group_ranges[i]
        dates <- as.IDate(round(seq(as.numeric(g$date_min), as.numeric(g$date_max), length.out = 90L)))
        unique(data.table(
          date = dates,
          treatment_label = factor(as.character(g$treatment_label), levels = c("Control", "Warmed")),
          depth_label = factor(as.character(g$depth_label), levels = c("10cm", "30cm"))
        ))
      })
      pred_grid <- rbindlist(grid_parts)
      pred_grid[, time_year := as.numeric(date - as.IDate("2021-01-01")) / 365.25]
      pred_grid[, time30 := as.numeric(date - min(d$date)) / 30]
      pred_grid[, plot := factor(levels(d$plot)[[1]], levels = levels(d$plot))]

      x <- align_model_matrix(fit, pred_grid)
      beta <- fixed.effects(fit)
      v <- vcov(fit)
      pred_grid[, fit_model_scale := as.numeric(x %*% beta)]
      pred_grid[, se_model_scale := sqrt(rowSums((x %*% v) * x))]
      pred_grid[, `:=`(
        fit = metric_from_model_scale(metric_name, fit_model_scale),
        lower95 = metric_from_model_scale(metric_name, fit_model_scale - qnorm(0.975) * se_model_scale),
        upper95 = metric_from_model_scale(metric_name, fit_model_scale + qnorm(0.975) * se_model_scale),
        metric = metric_name,
        metric_label = info$metric_label,
        metric_group = info$metric_group,
        metric_order = info$metric_order
      )]
      prediction_rows[[out_i]] <- pred_grid
    }

    out_i <- out_i + 1L
  }

  list(
    slopes = rbindlist(slope_rows, fill = TRUE),
    contrasts = rbindlist(contrast_rows, fill = TRUE),
    anova = rbindlist(anova_rows, fill = TRUE),
    coefficients = rbindlist(coefficient_rows, fill = TRUE),
    diagnostics = rbindlist(diagnostic_rows, fill = TRUE),
    predictions = rbindlist(prediction_rows, fill = TRUE),
    exact_treatment_tests = rbindlist(exact_test_rows, fill = TRUE),
    exact_treatment_null = rbindlist(exact_null_rows, fill = TRUE),
    summary_lines = summary_lines
  )
}

jaccard_distance <- function(a, b) {
  union_n <- length(union(a, b))
  if (union_n == 0L) {
    return(NA_real_)
  }
  1 - length(intersect(a, b)) / union_n
}

bray_curtis <- function(a, b) {
  denominator <- sum(a + b)
  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  sum(abs(a - b)) / denominator
}

build_cross_plot_pair_distances <- function(sample_meta) {
  pair_rows <- list()
  out_i <- 1L
  groups <- unique(sample_meta[, .(date, depth_label)])

  for (g_i in seq_len(nrow(groups))) {
    g <- groups[g_i]
    d <- sample_meta[g, on = c("date", "depth_label"), nomatch = 0]
    if (nrow(d) < 2L) {
      next
    }
    combinations <- combn(seq_len(nrow(d)), 2L)
    for (pair_i in seq_len(ncol(combinations))) {
      i <- combinations[1L, pair_i]
      j <- combinations[2L, pair_i]
      if (d$block[[i]] == d$block[[j]]) {
        next
      }
      pair_rows[[out_i]] <- data.table(
        date = d$date[[i]],
        depth_label = d$depth_label[[i]],
        sample_num_1 = d$sample_num[[i]],
        sample_num_2 = d$sample_num[[j]],
        plot_1 = d$plot[[i]],
        plot_2 = d$plot[[j]],
        block_1 = d$block[[i]],
        block_2 = d$block[[j]],
        block_pair = paste(sort(c(d$block[[i]], d$block[[j]])), collapse = "-"),
        observed_treatment_1 = as.character(d$treatment_label[[i]]),
        observed_treatment_2 = as.character(d$treatment_label[[j]]),
        formula_jaccard_distance = jaccard_distance(
          d$formula_indices[[i]],
          d$formula_indices[[j]]
        ),
        vk_class_braycurtis = bray_curtis(
          d$vk_class_proportions[[i]],
          d$vk_class_proportions[[j]]
        )
      )
      out_i <- out_i + 1L
    }
  }

  rbindlist(pair_rows, fill = TRUE)
}

permute_treatment <- function(observed, block, swap_blocks) {
  swap <- block %in% swap_blocks
  fifelse(
    !swap,
    observed,
    fifelse(observed == "Control", "Warmed", fifelse(observed == "Warmed", "Control", NA_character_))
  )
}

aggregate_pair_distances <- function(pair_dt, swap_blocks = character()) {
  d <- copy(pair_dt)
  d[, treatment_1 := permute_treatment(observed_treatment_1, block_1, swap_blocks)]
  d[, treatment_2 := permute_treatment(observed_treatment_2, block_2, swap_blocks)]
  d <- d[treatment_1 == treatment_2]
  d[, treatment_label := treatment_1]

  long <- melt(
    d,
    id.vars = c(
      "date", "depth_label", "sample_num_1", "sample_num_2", "plot_1", "plot_2",
      "block_1", "block_2", "block_pair", "treatment_label"
    ),
    measure.vars = c("formula_jaccard_distance", "vk_class_braycurtis"),
    variable.name = "distance_metric",
    value.name = "distance"
  )
  long <- long[is.finite(distance)]

  group_means <- long[
    ,
    .(
      mean_distance = mean(distance),
      n_pairs = .N,
      n_plots = uniqueN(c(plot_1, plot_2)),
      distance_sd = if (.N > 1L) sd(distance) else NA_real_
    ),
    by = .(date, depth_label, treatment_label, distance_metric)
  ]
  list(pair_long = long, group_means = group_means)
}

fit_gls_trend <- function(d, response_col) {
  d <- copy(d[is.finite(get(response_col))])
  setorder(d, date)
  if (nrow(d) < 6L || uniqueN(d$date) < 6L) {
    return(list(fit = NULL, result = data.table(
      n_dates = nrow(d),
      slope_per_year = NA_real_,
      se = NA_real_,
      df = NA_real_,
      lower95 = NA_real_,
      upper95 = NA_real_,
      t_value = NA_real_,
      p_value = NA_real_,
      p_one_sided_convergence = NA_real_,
      car1_phi_per_30d = NA_real_,
      fit_method = "insufficient_dates"
    )))
  }

  formula_text <- paste(response_col, "~ time_year")
  formula_obj <- as.formula(formula_text)
  fit_car <- try(
    gls(
      formula_obj,
      correlation = corCAR1(value = 0.3, form = ~time30),
      data = d,
      method = "REML",
      control = glsControl(opt = "optim", msMaxIter = 250, returnObject = TRUE)
    ),
    silent = TRUE
  )
  fit_method <- "gls_car1"
  fit <- fit_car

  if (inherits(fit_car, "try-error")) {
    fit <- lm(formula_obj, data = d)
    fit_method <- "lm_no_car1_fallback"
  }

  tt <- if (inherits(fit, "gls")) summary(fit)$tTable else summary(fit)$coefficients
  slope <- tt["time_year", "Value"]
  se <- tt["time_year", "Std.Error"]
  df <- if (inherits(fit, "gls")) fit$dims$N - length(coef(fit)) else df.residual(fit)
  t_value <- slope / se
  p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
  p_one_sided <- if (slope <= 0) p_value / 2 else 1 - p_value / 2

  list(
    fit = fit,
    result = data.table(
      n_dates = nrow(d),
      slope_per_year = slope,
      se = se,
      df = df,
      lower95 = slope - qt(0.975, df) * se,
      upper95 = slope + qt(0.975, df) * se,
      t_value = t_value,
      p_value = p_value,
      p_one_sided_convergence = p_one_sided,
      car1_phi_per_30d = if (inherits(fit, "gls")) as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE)) else NA_real_,
      fit_method = fit_method
    )
  )
}

fit_gls_mean <- function(d, response_col) {
  d <- copy(d[is.finite(get(response_col))])
  setorder(d, date)
  if (nrow(d) < 5L) {
    return(data.table(
      n_dates = nrow(d),
      mean_difference = NA_real_,
      se = NA_real_,
      df = NA_real_,
      lower95 = NA_real_,
      upper95 = NA_real_,
      t_value = NA_real_,
      p_value = NA_real_,
      car1_phi_per_30d = NA_real_,
      fit_method = "insufficient_dates"
    ))
  }

  formula_obj <- as.formula(paste(response_col, "~ 1"))
  fit_car <- try(
    gls(
      formula_obj,
      correlation = corCAR1(value = 0.3, form = ~time30),
      data = d,
      method = "REML",
      control = glsControl(opt = "optim", msMaxIter = 250, returnObject = TRUE)
    ),
    silent = TRUE
  )
  fit_method <- "gls_car1"
  fit <- fit_car
  if (inherits(fit_car, "try-error")) {
    fit <- lm(formula_obj, data = d)
    fit_method <- "lm_no_car1_fallback"
  }

  tt <- if (inherits(fit, "gls")) summary(fit)$tTable else summary(fit)$coefficients
  estimate <- tt["(Intercept)", "Value"]
  se <- tt["(Intercept)", "Std.Error"]
  df <- if (inherits(fit, "gls")) fit$dims$N - length(coef(fit)) else df.residual(fit)
  t_value <- estimate / se

  data.table(
    n_dates = nrow(d),
    mean_difference = estimate,
    se = se,
    df = df,
    lower95 = estimate - qt(0.975, df) * se,
    upper95 = estimate + qt(0.975, df) * se,
    t_value = t_value,
    p_value = 2 * pt(abs(t_value), df = df, lower.tail = FALSE),
    car1_phi_per_30d = if (inherits(fit, "gls")) as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE)) else NA_real_,
    fit_method = fit_method
  )
}

fit_convergence_models <- function(group_means) {
  d <- copy(group_means)
  d[, time_year := as.numeric(date - as.IDate("2021-01-01")) / 365.25]
  d[, time30 := as.numeric(date - min(date)) / 30]

  group_results <- d[
    ,
    fit_gls_trend(.SD, "mean_distance")$result,
    by = .(distance_metric, treatment_label, depth_label)
  ]

  paired <- dcast(
    d,
    distance_metric + depth_label + date + time_year + time30 ~ treatment_label,
    value.var = "mean_distance"
  )
  paired <- paired[is.finite(Control) & is.finite(Warmed)]
  paired[, warmed_minus_control_distance := Warmed - Control]
  paired_results <- paired[
    ,
    fit_gls_trend(.SD, "warmed_minus_control_distance")$result,
    by = .(distance_metric, depth_label)
  ]
  paired_results[, interpretation := "Negative slope means warmed plots converge faster than controls."]
  paired_level_results <- paired[
    ,
    fit_gls_mean(.SD, "warmed_minus_control_distance"),
    by = .(distance_metric, depth_label)
  ]
  paired_level_results[, interpretation := "Negative mean difference means warmed plots are more similar on matched dates."]

  prediction_rows <- list()
  out_i <- 1L
  for (metric_name in unique(d$distance_metric)) {
    for (depth_value in unique(d$depth_label)) {
      for (treatment_value in c("Control", "Warmed")) {
        sub <- d[
          distance_metric == metric_name &
            depth_label == depth_value &
            treatment_label == treatment_value
        ]
        fit_result <- fit_gls_trend(sub, "mean_distance")
        fit <- fit_result$fit
        if (is.null(fit)) {
          next
        }
        grid <- data.table(
          date = as.IDate(round(seq(as.numeric(min(sub$date)), as.numeric(max(sub$date)), length.out = 80L)))
        )
        grid <- unique(grid)
        grid[, time_year := as.numeric(date - as.IDate("2021-01-01")) / 365.25]
        grid[, time30 := as.numeric(date - min(sub$date)) / 30]
        x <- model.matrix(~time_year, grid)
        beta <- coef(fit)
        v <- vcov(fit)
        grid[, fit := as.numeric(x %*% beta)]
        grid[, se := sqrt(rowSums((x %*% v) * x))]
        grid[, `:=`(
          lower95 = fit - qnorm(0.975) * se,
          upper95 = fit + qnorm(0.975) * se,
          distance_metric = metric_name,
          depth_label = depth_value,
          treatment_label = treatment_value
        )]
        prediction_rows[[out_i]] <- grid
        out_i <- out_i + 1L
      }
    }
  }

  list(
    group_results = group_results,
    paired_dates = paired,
    paired_results = paired_results,
    paired_level_results = paired_level_results,
    predictions = rbindlist(prediction_rows, fill = TRUE)
  )
}

ordinary_delta_slope <- function(group_means, metric_name, depth_value) {
  d <- group_means[distance_metric == metric_name & depth_label == depth_value]
  d[, time_year := as.numeric(date - as.IDate("2021-01-01")) / 365.25]
  paired <- dcast(d, date + time_year ~ treatment_label, value.var = "mean_distance")
  paired <- paired[is.finite(Control) & is.finite(Warmed)]
  if (nrow(paired) < 5L) {
    return(NA_real_)
  }
  paired[, delta := Warmed - Control]
  unname(coef(lm(delta ~ time_year, data = paired))[["time_year"]])
}

ordinary_delta_mean <- function(group_means, metric_name, depth_value) {
  d <- group_means[distance_metric == metric_name & depth_label == depth_value]
  paired <- dcast(d, date ~ treatment_label, value.var = "mean_distance")
  paired <- paired[is.finite(Control) & is.finite(Warmed)]
  if (nrow(paired) < 5L) {
    return(NA_real_)
  }
  mean(paired$Warmed - paired$Control)
}

exact_block_randomization <- function(pair_dt) {
  swap_grid <- CJ(swap_a = c(FALSE, TRUE), swap_b = c(FALSE, TRUE), swap_c = c(FALSE, TRUE))
  swap_grid[, permutation_id := .I]
  metric_names <- c("formula_jaccard_distance", "vk_class_braycurtis")
  depth_values <- c("10cm", "30cm")
  null_rows <- list()
  out_i <- 1L

  for (perm_i in seq_len(nrow(swap_grid))) {
    row <- swap_grid[perm_i]
    swap_blocks <- c("A", "B", "C")[as.logical(unlist(row[, .(swap_a, swap_b, swap_c)]))]
    aggregated <- aggregate_pair_distances(pair_dt, swap_blocks = swap_blocks)$group_means
    for (metric_name in metric_names) {
      for (depth_value in depth_values) {
        null_rows[[out_i]] <- data.table(
          permutation_id = row$permutation_id,
          swap_blocks = paste(swap_blocks, collapse = ","),
          distance_metric = metric_name,
          depth_label = depth_value,
          warmed_minus_control_slope = ordinary_delta_slope(aggregated, metric_name, depth_value),
          warmed_minus_control_mean = ordinary_delta_mean(aggregated, metric_name, depth_value)
        )
        out_i <- out_i + 1L
      }
    }
  }

  null_dt <- rbindlist(null_rows)
  observed <- null_dt[permutation_id == 1L]
  tests <- observed[
    ,
    {
      null <- null_dt[
        distance_metric == .BY$distance_metric &
          depth_label == .BY$depth_label,
        warmed_minus_control_slope
      ]
      null_mean <- null_dt[
        distance_metric == .BY$distance_metric &
          depth_label == .BY$depth_label,
        warmed_minus_control_mean
      ]
      obs <- warmed_minus_control_slope
      obs_mean <- warmed_minus_control_mean
      list(
        observed_warmed_minus_control_slope = obs,
        n_exact_permutations = sum(is.finite(null)),
        p_one_sided_warmed_converges_faster = mean(null <= obs, na.rm = TRUE),
        p_two_sided = mean(abs(null) >= abs(obs), na.rm = TRUE),
        null_min = min(null, na.rm = TRUE),
        null_median = median(null, na.rm = TRUE),
        null_max = max(null, na.rm = TRUE),
        observed_warmed_minus_control_mean = obs_mean,
        p_one_sided_warmed_more_similar = mean(null_mean <= obs_mean, na.rm = TRUE),
        p_mean_two_sided = mean(abs(null_mean) >= abs(obs_mean), na.rm = TRUE),
        null_mean_min = min(null_mean, na.rm = TRUE),
        null_mean_median = median(null_mean, na.rm = TRUE),
        null_mean_max = max(null_mean, na.rm = TRUE)
      )
    },
    by = .(distance_metric, depth_label)
  ]

  list(null_distribution = null_dt, tests = tests)
}

theme_trace <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E7E9ED", linewidth = 0.28),
      strip.text = element_text(face = "bold", color = "#202124"),
      strip.text.y.left = element_text(angle = 0, hjust = 1, lineheight = 0.95),
      strip.placement = "outside",
      axis.text = element_text(color = "#31353B"),
      axis.title = element_text(color = "#202124"),
      plot.margin = margin(7, 10, 7, 10)
    )
}

save_plot_pair <- function(plot, stem, fig_dir, width, height) {
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 320, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, units = "in", bg = "white")
  message("Wrote: ", png_path)
  message("Wrote: ", pdf_path)
}

format_p <- function(x) {
  ifelse(!is.finite(x), "NA", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}

script_dir <- get_script_dir()
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_dir, "data")
out_dir <- file.path(repo_dir, "output", "fticr_integration")
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

metric_data_path <- file.path(out_dir, "fticr_q2_autocorrelation_metric_plot_data_vk_classes.csv")
mol_path <- file.path(data_dir, "raw", "merged_output", "Test_Processed-Unprocessed_Mol.csv")
intensity_path <- first_existing_path(c(
  file.path(data_dir, "raw", "merged_output", "Test_Processed-Unprocssed_Data.csv"),
  file.path(data_dir, "raw", "merged_output", "Test_Processed-Unprocssed_Data.csv.gz")
))

required_paths <- c(metric_data_path, mol_path, intensity_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required input files:\n", paste(missing_paths, collapse = "\n"))
}

message("Loading count-based sample metrics...")
all_metric_dt <- fread(metric_data_path)
all_metric_dt[, date := as.IDate(date)]
all_metric_dt[, treatment_label := normalize_treatment(treatment_label)]
all_metric_dt[, depth_label := depth_display(depth_cm)]

target_metric_info <- data.table(
  metric = c(
    "mean_nosc",
    "richness",
    "mean_ai",
    "vk_frac_aromatic",
    "vk_frac_carbohydrate",
    "vk_frac_condensed_aromatic",
    "vk_frac_lignin_like",
    "vk_frac_lipid",
    "vk_frac_other"
  ),
  metric_label = c(
    "Mean NOSC",
    "Formula richness",
    "Mean AI_mod",
    "Aromatic formulae",
    "Carbohydrate formulae",
    "Condensed aromatic formulae",
    "Lignin-like formulae",
    "Lipid formulae",
    "Other formulae"
  ),
  metric_group = c(rep("Scalar metrics", 3L), rep("Van Krevelen fractions", 6L)),
  metric_order = seq_len(9L)
)

missing_metrics <- setdiff(target_metric_info$metric, unique(all_metric_dt$metric))
if (length(missing_metrics) > 0L) {
  stop("Missing requested metrics from count-based metric table: ", paste(missing_metrics, collapse = ", "))
}

richness_dt <- all_metric_dt[
  metric == "richness",
  .(sample_num, richness_count = value_plot)
]
metric_dt <- merge(
  all_metric_dt[metric %in% target_metric_info$metric][
    ,
    .(sample_num, plot, date, treatment_label, depth_label, metric, value_plot)
  ],
  target_metric_info,
  by = "metric",
  all.x = TRUE
)
metric_dt <- merge(metric_dt, richness_dt, by = "sample_num", all.x = TRUE)
setorder(metric_dt, metric_order, depth_label, treatment_label, plot, date)

sample_grain <- unique(metric_dt[, .(sample_num, plot, date, treatment_label, depth_label)])
duplicate_sample_keys <- sample_grain[, .N, by = .(sample_num)][N > 1L]
if (nrow(duplicate_sample_keys) > 0L) {
  stop("Sample numbers are not unique at the intended sample grain.")
}

class_check <- metric_dt[grepl("^vk_frac_", metric), .(class_sum = sum(value_plot)), by = sample_num]
class_sum_max_abs_error <- max(abs(class_check$class_sum - 100), na.rm = TRUE)

message("Fitting calendar-time trend models...")
trend_all <- fit_metric_models(metric_dt, target_metric_info, analysis_window = "All available")
trend_2018 <- fit_metric_models(
  metric_dt,
  target_metric_info,
  analysis_window = "2018 onward",
  window_start = as.IDate("2018-01-01")
)

trend_slopes <- rbindlist(list(trend_all$slopes, trend_2018$slopes), fill = TRUE)
trend_slopes[, p_fdr_within_window := p.adjust(p.value, method = "BH"), by = analysis_window]
trend_contrasts <- rbindlist(list(trend_all$contrasts, trend_2018$contrasts), fill = TRUE)
trend_contrasts[, p_fdr_within_contrast_window := p.adjust(p.value, method = "BH"), by = .(analysis_window, contrast_type)]
trend_exact_tests <- rbindlist(
  list(trend_all$exact_treatment_tests, trend_2018$exact_treatment_tests),
  fill = TRUE
)
trend_exact_null <- rbindlist(
  list(trend_all$exact_treatment_null, trend_2018$exact_treatment_null),
  fill = TRUE
)
trend_contrasts <- merge(
  trend_contrasts,
  trend_exact_tests,
  by = c("analysis_window", "metric", "metric_label", "depth_label"),
  all.x = TRUE
)
trend_anova <- rbindlist(list(trend_all$anova, trend_2018$anova), fill = TRUE)
trend_coefficients <- rbindlist(list(trend_all$coefficients, trend_2018$coefficients), fill = TRUE)
trend_diagnostics <- rbindlist(list(trend_all$diagnostics, trend_2018$diagnostics), fill = TRUE)

message("Loading formula detections for exact-date cross-plot distances...")
mol_header <- fread(mol_path, nrows = 0L)
mol_select <- intersect(c("Molecular Formula", "Is Isotopologue"), names(mol_header))
mol <- fread(mol_path, select = mol_select)
if (!("Molecular Formula" %in% names(mol))) {
  stop("Molecular formula annotation table lacks `Molecular Formula`.")
}
if ("Is Isotopologue" %in% names(mol)) {
  mol <- mol[safe_numeric(`Is Isotopologue`) == 0 | is.na(`Is Isotopologue`)]
}
valid_annotated_formulas <- unique(mol[["Molecular Formula"]])

intensity_header <- read_fticr_csv(intensity_path, nrows = 0L)
intensity_names <- names(intensity_header)
formula_col <- "Molecular Formula"
if (!(formula_col %in% intensity_names)) {
  stop("Intensity matrix lacks `Molecular Formula`.")
}
field_cols <- grep("^Sihi_60398_[0-9]+_r[0-9]+_", intensity_names, value = TRUE)
blank_cols <- grep("Blank", intensity_names, value = TRUE, ignore.case = TRUE)
sample_map <- data.table(sample_col = field_cols)
sample_map[
  ,
  sample_num := suppressWarnings(as.integer(sub("^Sihi_60398_([0-9]+)_r[0-9]+_.*$", "\\1", sample_col)))
]
sample_map <- sample_map[!is.na(sample_num) & sample_num %in% sample_grain$sample_num]
selected_intensity_cols <- unique(c(formula_col, blank_cols, sample_map$sample_col))

message("Reading ", length(selected_intensity_cols), " formula-matrix columns...")
intensity <- read_fticr_csv(intensity_path, select = selected_intensity_cols)
setnames(intensity, formula_col, "formula")

blank_formulas <- character()
if (length(blank_cols) > 0L) {
  blank_matrix <- as.matrix(intensity[, ..blank_cols])
  storage.mode(blank_matrix) <- "double"
  blank_detected <- rowSums(is.finite(blank_matrix) & blank_matrix > 0) > 0
  blank_formulas <- intensity[blank_detected, formula]
}
valid_formula_set <- setdiff(intersect(intensity$formula, valid_annotated_formulas), blank_formulas)
formula_id <- data.table(formula = valid_formula_set, formula_index = seq_along(valid_formula_set))
setkey(formula_id, formula)

sample_formula_rows <- vector("list", nrow(sample_grain))
for (i in seq_len(nrow(sample_grain))) {
  sample_number <- sample_grain$sample_num[[i]]
  sample_cols <- sample_map[sample_num == sample_number, sample_col]
  if (length(sample_cols) == 0L) {
    next
  }
  sample_matrix <- as.matrix(intensity[, ..sample_cols])
  storage.mode(sample_matrix) <- "double"
  detected <- rowSums(is.finite(sample_matrix) & sample_matrix > 0) > 0
  formulas <- intersect(intensity[detected, formula], valid_formula_set)
  indices <- formula_id[J(formulas), nomatch = 0L, formula_index]
  sample_formula_rows[[i]] <- data.table(
    sample_num = sample_number,
    formula_indices = list(sort(indices))
  )
}
sample_formulas <- rbindlist(sample_formula_rows, fill = TRUE)

vk_metrics <- target_metric_info[metric_group == "Van Krevelen fractions", metric]
vk_wide <- dcast(
  metric_dt[metric %in% vk_metrics],
  sample_num + plot + date + treatment_label + depth_label ~ metric,
  value.var = "value_plot"
)
for (metric_name in vk_metrics) {
  set(vk_wide, j = metric_name, value = vk_wide[[metric_name]] / 100)
}
vk_wide[
  ,
  vk_class_proportions := split(as.matrix(.SD), seq_len(.N)),
  .SDcols = vk_metrics
]

sample_meta <- merge(
  sample_grain,
  sample_formulas,
  by = "sample_num",
  all.x = TRUE
)
sample_meta <- merge(
  sample_meta,
  vk_wide[, .(sample_num, vk_class_proportions)],
  by = "sample_num",
  all.x = TRUE
)
sample_meta[, block := block_from_plot(plot)]
sample_meta <- sample_meta[
  !sapply(formula_indices, is.null) &
    !sapply(vk_class_proportions, is.null) &
    !is.na(block)
]

formula_coverage <- nrow(sample_meta) / nrow(sample_grain)
if (formula_coverage < 0.95) {
  warning("Only ", round(100 * formula_coverage, 1), "% of sample metadata rows matched formula detections.")
}

message("Computing exact-date cross-plot distances...")
candidate_pair_dt <- build_cross_plot_pair_distances(sample_meta)
observed_aggregate <- aggregate_pair_distances(candidate_pair_dt)
observed_pair_long <- observed_aggregate$pair_long
observed_group_means <- observed_aggregate$group_means

convergence_models <- fit_convergence_models(observed_group_means)
randomization <- exact_block_randomization(candidate_pair_dt)

distance_labels <- c(
  formula_jaccard_distance = "Formula Jaccard distance",
  vk_class_braycurtis = "Chemical-class Bray-Curtis distance"
)
observed_pair_long[, distance_label := distance_labels[distance_metric]]
observed_group_means[, distance_label := distance_labels[distance_metric]]
convergence_models$group_results[, distance_label := distance_labels[distance_metric]]
convergence_models$paired_results[, distance_label := distance_labels[distance_metric]]
convergence_models$paired_level_results[, distance_label := distance_labels[distance_metric]]
convergence_models$predictions[, distance_label := distance_labels[distance_metric]]
randomization$tests[, distance_label := distance_labels[distance_metric]]

qa_dt <- rbindlist(list(
  data.table(
    check = c(
      "sample_count",
      "sample_count_2018_onward",
      "first_sample_date",
      "last_sample_date",
      "duplicate_sample_numbers",
      "max_abs_vk_fraction_sum_error_percent",
      "formula_detection_sample_coverage",
      "blank_formula_count_removed"
    ),
    value = c(
      as.character(nrow(sample_grain)),
      as.character(sample_grain[date >= as.IDate("2018-01-01"), .N]),
      as.character(min(sample_grain$date)),
      as.character(max(sample_grain$date)),
      as.character(nrow(duplicate_sample_keys)),
      format(class_sum_max_abs_error, scientific = TRUE),
      sprintf("%.6f", formula_coverage),
      as.character(length(unique(blank_formulas)))
    )
  ),
  observed_group_means[
    ,
    .(
      check = paste0(
        "exact_date_groups_",
        distance_metric,
        "_",
        treatment_label,
        "_",
        depth_label
      ),
      value = as.character(.N)
    ),
    by = .(distance_metric, treatment_label, depth_label)
  ][, c("distance_metric", "treatment_label", "depth_label") := NULL]
), fill = TRUE)

message("Writing analysis tables...")
output_tables <- list(
  fticr_calendar_time_trend_slopes = trend_slopes,
  fticr_calendar_time_trend_contrasts = trend_contrasts,
  fticr_calendar_time_trend_exact_treatment_tests = trend_exact_tests,
  fticr_calendar_time_trend_exact_treatment_null = trend_exact_null,
  fticr_calendar_time_trend_marginal_tests = trend_anova,
  fticr_calendar_time_trend_coefficients = trend_coefficients,
  fticr_calendar_time_trend_diagnostics = trend_diagnostics,
  fticr_calendar_time_trend_predictions = trend_all$predictions,
  fticr_cross_plot_composition_pair_distances_exact_date = observed_pair_long,
  fticr_cross_plot_composition_group_means_exact_date = observed_group_means,
  fticr_cross_plot_convergence_slopes = convergence_models$group_results,
  fticr_cross_plot_convergence_paired_treatment_slopes = convergence_models$paired_results,
  fticr_cross_plot_convergence_paired_treatment_levels = convergence_models$paired_level_results,
  fticr_cross_plot_convergence_paired_dates = convergence_models$paired_dates,
  fticr_cross_plot_convergence_predictions = convergence_models$predictions,
  fticr_cross_plot_convergence_exact_randomization_tests = randomization$tests,
  fticr_cross_plot_convergence_exact_randomization_null = randomization$null_distribution,
  fticr_calendar_time_trends_data_quality = qa_dt
)
for (output_name in names(output_tables)) {
  output_path <- file.path(out_dir, paste0(output_name, ".csv"))
  fwrite(output_tables[[output_name]], output_path)
  message("Wrote: ", output_path)
}

model_summary_path <- file.path(out_dir, "fticr_calendar_time_trend_model_summaries.txt")
writeLines(
  c(
    "TRACE FTICR calendar-time trend mixed models",
    "============================================",
    "",
    "Response transformations:",
    "- Formula richness: natural log.",
    "- Van Krevelen fractions: empirical logit using formula richness as the count denominator.",
    "- Mean NOSC and mean AI_mod: untransformed.",
    "",
    "Fixed effects: time_year * treatment * depth.",
    "Random effects: plot intercepts plus plot-depth intercepts and calendar-time slopes.",
    "Residual correlation: continuous-time CAR(1) within plot-depth when the model converged.",
    "",
    trend_all$summary_lines,
    "\n\n2018-ONWARD SENSITIVITY MODELS\n",
    trend_2018$summary_lines
  ),
  model_summary_path
)
message("Wrote: ", model_summary_path)

treatment_cols <- c("Control" = "#2B6CB0", "Warmed" = "#C45A2A")
treatment_linetypes <- c("Control" = "solid", "Warmed" = "22")

plot_observed <- metric_dt[
  metric %in% target_metric_info$metric &
    is.finite(value_plot)
]
plot_observed[, treatment_label := factor(treatment_label, levels = c("Control", "Warmed"))]
plot_observed[, depth_label := factor(depth_label, levels = c("10cm", "30cm"))]
trend_predictions <- copy(trend_all$predictions)
trend_predictions[, treatment_label := factor(treatment_label, levels = c("Control", "Warmed"))]
trend_predictions[, depth_label := factor(depth_label, levels = c("10cm", "30cm"))]

scalar_metrics <- target_metric_info[metric_group == "Scalar metrics", metric]
scalar_levels <- target_metric_info[metric %in% scalar_metrics, metric_label]
scalar_obs <- plot_observed[metric %in% scalar_metrics]
scalar_pred <- trend_predictions[metric %in% scalar_metrics]
scalar_obs[, metric_label := factor(metric_label, levels = scalar_levels)]
scalar_pred[, metric_label := factor(metric_label, levels = scalar_levels)]

p_scalar <- ggplot() +
  geom_ribbon(
    data = scalar_pred,
    aes(x = date, ymin = lower95, ymax = upper95, fill = treatment_label, group = treatment_label),
    alpha = 0.14,
    color = NA
  ) +
  geom_point(
    data = scalar_obs,
    aes(x = date, y = value_plot, color = treatment_label),
    alpha = 0.25,
    size = 1.1
  ) +
  geom_line(
    data = scalar_pred,
    aes(x = date, y = fit, color = treatment_label, linetype = treatment_label),
    linewidth = 0.82
  ) +
  facet_grid(metric_label ~ depth_label, scales = "free_y", switch = "y") +
  scale_color_manual(values = treatment_cols) +
  scale_fill_manual(values = treatment_cols) +
  scale_linetype_manual(values = treatment_linetypes) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = NULL, y = "Observed and model-estimated value") +
  theme_trace() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot_pair(
  p_scalar,
  "fticr_calendar_time_trends_scalar_metrics",
  fig_dir,
  width = 10.5,
  height = 8.0
)

fraction_metrics <- target_metric_info[metric_group == "Van Krevelen fractions", metric]
fraction_levels <- target_metric_info[metric %in% fraction_metrics, metric_label]
fraction_obs <- plot_observed[metric %in% fraction_metrics]
fraction_pred <- trend_predictions[metric %in% fraction_metrics]
fraction_obs[, metric_label := factor(metric_label, levels = fraction_levels)]
fraction_pred[, metric_label := factor(metric_label, levels = fraction_levels)]

p_fractions <- ggplot() +
  geom_ribbon(
    data = fraction_pred,
    aes(x = date, ymin = lower95, ymax = upper95, fill = treatment_label, group = treatment_label),
    alpha = 0.14,
    color = NA
  ) +
  geom_point(
    data = fraction_obs,
    aes(x = date, y = value_plot, color = treatment_label),
    alpha = 0.23,
    size = 0.95
  ) +
  geom_line(
    data = fraction_pred,
    aes(x = date, y = fit, color = treatment_label, linetype = treatment_label),
    linewidth = 0.78
  ) +
  facet_grid(metric_label ~ depth_label, scales = "free_y", switch = "y") +
  scale_color_manual(values = treatment_cols) +
  scale_fill_manual(values = treatment_cols) +
  scale_linetype_manual(values = treatment_linetypes) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = NULL, y = "Formula-count fraction (%)") +
  theme_trace(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot_pair(
  p_fractions,
  "fticr_calendar_time_trends_vk_fractions",
  fig_dir,
  width = 10.5,
  height = 12.2
)

convergence_obs <- copy(observed_group_means)
convergence_pred <- copy(convergence_models$predictions)
distance_level_order <- unname(distance_labels)
distance_plot_labels <- c(
  "Formula Jaccard distance" = "Formula Jaccard\ndistance",
  "Chemical-class Bray-Curtis distance" = "Chemical-class\nBray-Curtis distance"
)
convergence_obs[, distance_plot_label := factor(
  distance_plot_labels[as.character(distance_label)],
  levels = unname(distance_plot_labels[distance_level_order])
)]
convergence_pred[, distance_plot_label := factor(
  distance_plot_labels[as.character(distance_label)],
  levels = unname(distance_plot_labels[distance_level_order])
)]
convergence_obs[, depth_label := factor(depth_label, levels = c("10cm", "30cm"))]
convergence_pred[, depth_label := factor(depth_label, levels = c("10cm", "30cm"))]
convergence_obs[, treatment_label := factor(treatment_label, levels = c("Control", "Warmed"))]
convergence_pred[, treatment_label := factor(treatment_label, levels = c("Control", "Warmed"))]

p_convergence <- ggplot() +
  geom_ribbon(
    data = convergence_pred,
    aes(x = date, ymin = lower95, ymax = upper95, fill = treatment_label, group = treatment_label),
    alpha = 0.14,
    color = NA
  ) +
  geom_point(
    data = convergence_obs,
    aes(x = date, y = mean_distance, color = treatment_label, size = n_plots),
    alpha = 0.75
  ) +
  geom_line(
    data = convergence_pred,
    aes(x = date, y = fit, color = treatment_label, linetype = treatment_label),
    linewidth = 0.82
  ) +
  facet_grid(distance_plot_label ~ depth_label, scales = "free_y", switch = "y") +
  scale_color_manual(values = treatment_cols) +
  scale_fill_manual(values = treatment_cols) +
  scale_linetype_manual(values = treatment_linetypes) +
  scale_size_continuous(range = c(1.9, 3.0), breaks = c(2, 3), guide = "none") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = c(0.02, 0.03))) +
  labs(x = NULL, y = "Mean same-date distance among plots") +
  theme_trace() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot_pair(
  p_convergence,
  "fticr_cross_plot_composition_convergence",
  fig_dir,
  width = 10.5,
  height = 6.8
)

primary_slopes <- trend_slopes[analysis_window == "All available"]
primary_treatment_contrasts <- trend_contrasts[
  analysis_window == "All available" &
    contrast_type == "Warmed minus Control time slope"
]
sensitivity_direction_check <- merge(
  trend_exact_tests[
    analysis_window == "All available",
    .(metric, depth_label, all_available_slope = observed_warmed_minus_control_slope)
  ],
  trend_exact_tests[
    analysis_window == "2018 onward",
    .(metric, depth_label, sensitivity_slope = observed_warmed_minus_control_slope)
  ],
  by = c("metric", "depth_label")
)
n_sensitivity_same_direction <- sensitivity_direction_check[
  sign(all_available_slope) == sign(sensitivity_slope),
  .N
]
n_sensitivity_contrasts <- nrow(sensitivity_direction_check)

trend_summary_lines <- primary_slopes[
  order(metric_order, depth_label, treatment_label),
  sprintf(
    "- %s, %s, %s: %+.4f %s (95%% CI %+.4f to %+.4f; model-scale p = %s).",
    metric_label,
    depth_label,
    as.character(treatment_label),
    annual_change_natural,
    annual_change_units,
    annual_change_natural_lower95,
    annual_change_natural_upper95,
    format_p(p.value)
  )
]

contrast_summary_lines <- primary_treatment_contrasts[
  order(metric, depth_label),
  sprintf(
    "- %s, %s: warmed-minus-control model-scale slope = %+.4f (95%% CI %+.4f to %+.4f; exact paired-plot p = %.3f two-sided; mixed-model p = %s).",
    metric_label,
    as.character(depth_label),
    warmed_minus_control_model_slope,
    lower95,
    upper95,
    exact_p_two_sided,
    format_p(p.value)
  )
]

convergence_summary_lines <- convergence_models$group_results[
  order(distance_metric, depth_label, treatment_label),
  sprintf(
    "- %s, %s, %s: slope %+.4f distance units/year (95%% CI %+.4f to %+.4f; n = %d exact dates; p = %s).",
    distance_labels[distance_metric],
    depth_label,
    treatment_label,
    slope_per_year,
    lower95,
    upper95,
    n_dates,
    format_p(p_value)
  )
]

randomization_summary_lines <- randomization$tests[
  order(distance_metric, depth_label),
  sprintf(
    "- %s, %s: warmed-minus-control slope %+.4f/year; exact one-sided p = %.3f, two-sided p = %.3f.",
    distance_labels[distance_metric],
    depth_label,
    observed_warmed_minus_control_slope,
    p_one_sided_warmed_converges_faster,
    p_two_sided
  )
]

level_summary_lines <- convergence_models$paired_level_results[
  order(distance_metric, depth_label),
  sprintf(
    "- %s, %s: warmed minus control = %+.4f distance units (95%% CI %+.4f to %+.4f; exact one-sided p for warmed being more similar = %.3f).",
    distance_labels[distance_metric],
    depth_label,
    mean_difference,
    lower95,
    upper95,
    randomization$tests[
      distance_metric == .BY$distance_metric & depth_label == .BY$depth_label,
      p_one_sided_warmed_more_similar
    ]
  ),
  by = .(distance_metric, depth_label)
]$V1

summary_md_path <- file.path(out_dir, "fticr_calendar_time_trends_and_convergence_summary.md")
writeLines(
  c(
    "# TRACE FT-ICR calendar-time trends and cross-plot convergence",
    "",
    "## Questions",
    "",
    "1. Do mean NOSC, formula richness, mean AI_mod, and formula-count Van Krevelen fractions change through calendar time depending on treatment and depth?",
    "2. Do warmed plots become more compositionally similar to one another through time?",
    "",
    "## Data",
    "",
    sprintf(
      "- The packaged table contains %d spectra from %s through %s; %d spectra occur on or after 2018-01-01.",
      nrow(sample_grain),
      min(sample_grain$date),
      max(sample_grain$date),
      sample_grain[date >= as.IDate("2018-01-01"), .N]
    ),
    "- The primary trend models use all available spectra. A complete 2018-onward sensitivity set is saved in the same output tables.",
    "- All molecular averages and fractions are based on detected formula counts, not intensity.",
    "- The six fractions are aromatic, carbohydrate, condensed aromatic, lignin-like, lipid, and other formulae.",
    "",
    "## Calendar-time trend model",
    "",
    "- Each metric is modeled as `time * treatment * depth` with a random plot intercept and random plot-depth intercepts and time slopes.",
    "- Irregular-time residual dependence is represented with a continuous-time CAR(1) term within plot-depth when that model converges.",
    "- Richness is log transformed. Fractions use an empirical logit based on formula counts. NOSC and AI_mod remain in their raw units.",
    "- Reported natural-scale fraction slopes are percentage-point changes per year evaluated at the midpoint of each treatment-depth series.",
    "",
    "### Treatment-depth slopes",
    "",
    trend_summary_lines,
    "",
    "### Warming differences in temporal slope",
    "",
    contrast_summary_lines,
    "",
    "- The exact paired-plot p-value is the primary treatment inference. It swaps treatment labels within plot pairs 1/2, 3/4, and 5/6 and therefore has only eight possible assignments.",
    sprintf(
      "- Omitting the 2017 spectra preserved the direction of %d of %d warmed-minus-control slope contrasts.",
      n_sensitivity_same_direction,
      n_sensitivity_contrasts
    ),
    "",
    "## Cross-plot convergence",
    "",
    "- Spectra are compared only when plots share the exact same collection date and depth; no interpolation is used.",
    "- Formula Jaccard distance tests convergence in formula identity. Six-class Bray-Curtis distance tests convergence in broad chemical composition.",
    "- A negative within-treatment slope means plots become more similar through time.",
    "- The paired treatment comparison uses dates observed in both treatments. The exact test also swaps treatment labels within the three randomized plot pairs, yielding only eight possible assignments.",
    "",
    "### Within-treatment slopes",
    "",
    convergence_summary_lines,
    "",
    "### Average matched-date treatment difference",
    "",
    level_summary_lines,
    "",
    "- These are level differences, not convergence rates. A stable negative difference means warmed plots are more similar on average without becoming progressively more similar.",
    "",
    "### Exact treatment-randomization tests",
    "",
    randomization_summary_lines,
    "",
    "## Interpretation limits",
    "",
    "- Treatment is replicated at only three paired plots, so treatment-by-time inference is intrinsically low power. The smallest attainable one-sided exact p-value is 0.125.",
    "- Mixed-model p-values use repeated spectra and are retained as descriptive model output; treatment conclusions should be based on the exact paired-plot test and effect sizes.",
    "- Exact-date convergence coverage ends on 2022-11-16 at 10cm and 2023-05-30 at 30cm because later dates do not include at least two plots in both treatments.",
    "- Formula Jaccard level differences can reflect both formula overlap and formula richness; chemical-class Bray-Curtis uses class proportions and is less sensitive to total formula count.",
    "- Separate fraction models are compositionally dependent because the six fractions sum to 100%. Interpret the pattern across classes rather than treating the tests as independent discoveries.",
    "- Linear slopes summarize average directional change. They do not establish that trajectories are strictly linear or that warming caused any calendar-time pattern.",
    "",
    "## Reproducibility",
    "",
    "- Script: `scripts/fticr_integration/14_calendar_time_trends_and_cross_plot_convergence.R`",
    "- Full mixed-model summaries: `output/fticr_integration/fticr_calendar_time_trend_model_summaries.txt`",
    "- Model-ready tables and exact-randomization results are in `output/fticr_integration/`.",
    "- Figures are in `output/fticr_integration/figures/` and contain no embedded titles or subtitles."
  ),
  summary_md_path
)
message("Wrote: ", summary_md_path)

message("Calendar-time trend and cross-plot convergence analysis complete.")
