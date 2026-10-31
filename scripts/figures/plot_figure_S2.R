#!/usr/bin/env Rscript
# Figure S2: shared-frailty bias (a) and corrected coverage (b).

argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(dirname(normalizePath(script))), "utils", "bootstrap_regression.R"))
SetRegressionWd()

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(optparse)
})

SourceRegressionUtil("ground_truth_io.R")

# -----------------------------------------------------------------------------
# Plot helpers
# -----------------------------------------------------------------------------

GroundTruthLongFrailty <- function(gt) {
  out_mcf <- gt %>%
    dplyr::select(
      frailty_variance,
      beta_name,
      beta_value,
      tau,
      truth = beta_A_MCF
    ) %>%
    dplyr::mutate(term = "A", type_gt = "MCF")
  out_auc <- gt %>%
    dplyr::select(
      frailty_variance,
      beta_name,
      beta_value,
      tau,
      truth = beta_A_AUMCF
    ) %>%
    dplyr::mutate(term = "A", type_gt = "AUMCF")
  out <- dplyr::bind_rows(out_mcf, out_auc)
  return(out)
}

#' Wilson score (interval) bounds for binomial proportion; vectorized on k
#' successes and n trials.
WilsonCoverageCi <- function(n_cover, n_rep, conf_level = 0.95) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  k <- as.numeric(n_cover)
  n <- as.numeric(n_rep)
  p_hat <- ifelse(n > 0, k / n, NA_real_)
  z2 <- z^2
  denom <- 1 + z2 / n
  centre <- (p_hat + z2 / (2 * n)) / denom
  half <- z * sqrt(p_hat * (1 - p_hat) / n + z2 / (4 * n^2)) / denom
  lower <- centre - half
  upper <- centre + half
  lower <- pmax(0, lower)
  upper <- pmin(1, upper)
  bad <- !is.finite(n) | n <= 0
  lower[bad] <- NA_real_
  upper[bad] <- NA_real_
  out <- list(lower = lower, upper = upper)
  return(out)
}

LoadFrailtyData <- function(dir_in) {
  files <- list.files(
    path = dir_in,
    pattern = "\\.parquet$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(files) == 0L) {
    stop("No parquet files found in ", dir_in)
  }
  out <- dplyr::bind_rows(lapply(files, arrow::read_parquet))
  ValidateStudyInference(out, "frailty")
  req <- c("estimate", "lower", "upper", "term", "beta_name", "beta_value", "n", "tau", "type", "frailty_variance")
  miss <- setdiff(req, names(out))
  if (length(miss) > 0L) {
    stop("Frailty simulation data missing columns: ", paste(miss, collapse = ", "))
  }
  return(out)
}

#' Join simulations to frailty ground truth and build plot tibbles.
#'
#' @return A list with `coverage` (aggregated proportions + Wilson CIs) and
#'   `bias` (replicate-level bias for boxplots).
PrepareFrailtyPlotData <- function(data, gt_long) {
  joined <- data %>%
    dplyr::mutate(
      type = toupper(type),
      type_gt = dplyr::if_else(type == "AUC", "AUMCF", "MCF")
    ) %>%
    dplyr::left_join(
      gt_long,
      by = c("frailty_variance", "beta_name", "beta_value", "tau", "term", "type_gt")
    ) %>%
    dplyr::mutate(
      bias = estimate - truth,
      covers = (truth >= lower) & (truth <= upper)
    )
  if (anyNA(joined$truth)) {
    stop("Some frailty rows are missing matching ground truth.", call. = FALSE)
  }

  tau_levels <- paste0("tau == ", sort(unique(joined$tau)))
  n_levels <- sort(unique(joined$n))
  frailty_levels <- sort(unique(joined$frailty_variance))

  out_cov <- joined %>%
    dplyr::group_by(type, n, tau, frailty_variance) %>%
    dplyr::summarise(
      n_rep = sum(!is.na(covers)),
      n_cover = sum(covers, na.rm = TRUE),
      coverage = dplyr::if_else(n_rep > 0, n_cover / n_rep, NA_real_),
      .groups = "drop"
    )
  ci <- WilsonCoverageCi(out_cov$n_cover, out_cov$n_rep, conf_level = 0.95)
  out_cov <- out_cov %>%
    dplyr::mutate(
      coverage_lower = ci$lower,
      coverage_upper = ci$upper,
      type_label = dplyr::recode(type, AUC = "AUMCF", MCF = "MCF"),
      tau_label = factor(paste0("tau == ", tau), levels = tau_levels),
      N = factor(n, levels = n_levels),
      frailty_f = factor(frailty_variance, levels = frailty_levels)
    )

  out_bias <- joined %>%
    dplyr::filter(is.finite(bias)) %>%
    dplyr::mutate(
      type_label = dplyr::recode(type, AUC = "AUMCF", MCF = "MCF"),
      tau_label = factor(paste0("tau == ", tau), levels = tau_levels),
      N = factor(n, levels = n_levels),
      frailty_f = factor(frailty_variance, levels = frailty_levels)
    )

  out <- list(coverage = out_cov, bias = out_bias)
  return(out)
}

BuildCoveragePlot <- function(data) {
  dodge_width <- 0.9
  ggplot(data, aes(x = frailty_f, y = coverage, fill = N, group = N)) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    geom_hline(
      yintercept = 0.95,
      linetype = "dashed",
      linewidth = 0.35,
      color = "gray35"
    ) +
    geom_col(position = position_dodge(width = dodge_width), width = 0.82) +
    geom_errorbar(
      aes(ymin = coverage_lower, ymax = coverage_upper),
      position = position_dodge(width = dodge_width),
      width = 0.22,
      linewidth = 0.35
    ) +
    facet_grid(
      rows = vars(type_label),
      cols = vars(tau_label),
      labeller = labeller(tau_label = label_parsed)
    ) +
    scale_fill_brewer(name = "N per Arm", palette = "Blues") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Frailty Variance", y = "Coverage")
}

BuildBiasPlot <- function(data) {
  ggplot(data, aes(x = frailty_f, y = bias, fill = N)) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "gray50") +
    geom_boxplot(
      position = position_dodge2(width = 0.85, preserve = "single"),
      width = 0.78,
      linewidth = 0.35,
      outlier.size = 0.45,
      alpha = 0.92
    ) +
    facet_grid(
      rows = vars(type_label),
      cols = vars(tau_label),
      labeller = labeller(tau_label = label_parsed)
    ) +
    scale_fill_brewer(name = "N per Arm", palette = "Blues") +
    labs(x = "Frailty Variance", y = "Bias")
}

opt_list <- list(
  make_option(c("--dir"), type = "character", default = NULL),
  make_option(c("--ground_truth"), type = "character", default = NULL),
  make_option(c("--out_dir"), type = "character", default = NULL)
)

parsed <- parse_args(OptionParser(option_list = opt_list))
reg_dir <- RegressionDir()

source(file.path(reg_dir, "scripts", "utils", "ground_truth_io.R"), local = TRUE)

dir_in <- parsed$dir
if (is.null(dir_in)) {
  dir_in <- file.path(reg_dir, "simulations", "pseudoreg", "frailty")
}
gt_path <- parsed$ground_truth
if (is.null(gt_path)) {
  gt_path <- file.path(reg_dir, "simulations", "frailty_ground_truth")
}
out_dir <- parsed$out_dir
if (is.null(out_dir)) {
  out_dir <- FigureDir("frailty")
}

raw <- LoadFrailtyData(dir_in)
gt_bundle <- LoadGroundTruthRawAndAveraged(gt_path)
gt_long <- GroundTruthLongFrailty(gt_bundle$averaged)
plot_bundle <- PrepareFrailtyPlotData(raw, gt_long)

p_cov <- BuildCoveragePlot(plot_bundle$coverage)
p_bias <- BuildBiasPlot(plot_bundle$bias)

SaveRegressionFigure(
  plot = p_cov,
  path_pdf = file.path(out_dir, "figS2b.pdf"),
  width = 10,
  height = 6
)
SaveRegressionFigure(
  plot = p_bias,
  path_pdf = file.path(out_dir, "figS2a.pdf"),
  width = 10,
  height = 6
)
message("Wrote frailty plots to ", out_dir)
