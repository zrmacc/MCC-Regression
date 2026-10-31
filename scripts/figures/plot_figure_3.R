#!/usr/bin/env Rscript
# Figure 3: uniform QQ plot of corrected normal-Wald p-values.

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
SourceRegressionUtil("t1e_sim.R")

#' Load canonical type I error parquet p-values.
#'
#' @param dir_in Canonical directory containing parquet outputs.
#' @return Data frame of combined simulation results.
LoadPvalues <- function(dir_in) {
  out <- LoadCanonicalT1eData(dir_in)
  req_cols <- c("pvalue_arm", "n", "n_per_arm", "tau", "type", "source")
  miss_cols <- setdiff(req_cols, names(out))
  if (length(miss_cols) > 0L) {
    stop(
      "Combined data are missing required columns: ",
      paste(miss_cols, collapse = ", ")
    )
  }
  if (any(out$n_per_arm == 25L) || any(out$n == 25L)) {
    stop("Type I QQ input must not contain 25-per-arm rows.", call. = FALSE)
  }
  if (!all(out$n == out$n_per_arm)) {
    stop("Type I QQ input n must equal n_per_arm.", call. = FALSE)
  }
  if (!all(out$n_per_arm %in% c(50L, 100L, 200L))) {
    stop("Type I QQ input n_per_arm must be 50, 100, or 200.", call. = FALSE)
  }
  counts <- out %>%
    dplyr::count(.data$n_per_arm, .data$tau, .data$type, name = "n_rep")
  if (nrow(counts) != 18L || any(counts$n_rep != 100000L)) {
    stop("Type I QQ input must contain 18 complete 100,000-replicate cells.", call. = FALSE)
  }
  out$n <- out$n_per_arm
  return(out)
}

#' Prepare QQ plotting data for one group.
#'
#' @param pvals Numeric vector of p-values.
#' @param p_max Maximum p-value to include in the QQ plot.
#' @return Data frame with expected and observed -log10(p).
PrepQQDF <- function(
    pvals,
    p_max = 0.05
) {
  pvals <- pvals[is.finite(pvals)]
  n_orig <- length(pvals)
  pvals <- pvals[pvals > 0 & pvals <= p_max]
  n_keep <- length(pvals)

  if (n_keep == 0L) {
    out <- data.frame(
      exp = numeric(0),
      obs = numeric(0)
    )
    return(out)
  }

  obs <- -log10(sort(pvals))
  exp <- -log10(seq_len(n_keep) / (n_orig + 1))
  out <- data.frame(
    exp = exp,
    obs = obs
  )
  return(out)
}

#' Prepare QQ confidence-band data for one group.
#'
#' @param n_sim Number of simulation replicates in the group.
#' @param p_max Maximum p-value to include in the QQ plot.
#' @return Data frame with expected values and pointwise beta bands.
PrepBandsDF <- function(
    n_sim,
    p_max = 0.05
) {
  if (n_sim <= 0L) {
    out <- data.frame(
      exp = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    )
    return(out)
  }

  coord <- seq_len(n_sim)
  exp <- -log10(coord / (n_sim + 1))
  lower <- -log10(
    stats::qbeta(
      0.975,
      coord,
      n_sim - coord + 1
    )
  )
  upper <- -log10(
    stats::qbeta(
      0.025,
      coord,
      n_sim - coord + 1
    )
  )

  out <- data.frame(
    exp = exp,
    lower = lower,
    upper = upper
  )
  out <- out %>%
    dplyr::filter(exp >= -log10(p_max))
  return(out)
}

#' Add plotting labels and factors.
#'
#' @param data Data frame with n, tau, and type columns.
#' @return Data frame with plotting factors added.
AddPlotLabels <- function(data) {
  out <- data %>%
    dplyr::mutate(
      type = toupper(type),
      type_label = dplyr::recode(type, AUC = "AUMCF", MCF = "MCF"),
      type_label = factor(
        type_label,
        levels = c("AUMCF", "MCF")
      ),
      tau_label = factor(
        paste0("tau == ", tau),
        levels = c("tau == 2", "tau == 3", "tau == 4")
      ),
      N = factor(
        n,
        levels = c(50, 100, 200)
      )
    )
  return(out)
}

#' Prepare QQ point and band data.
#'
#' @param data Combined simulation results.
#' @param p_max Maximum p-value to include in the QQ plot.
#' @return List with QQ point data and confidence-band data.
PreparePlotData <- function(
    data,
    p_max = 0.05
) {
  data <- data %>%
    dplyr::filter(
      n %in% c(50, 100, 200),
      tau %in% c(2, 3, 4),
      toupper(type) %in% c("AUC", "MCF")
    )

  if (nrow(data) == 0L) {
    stop("No rows remain after filtering to n in {50,100,200}, tau in {2,3,4}, and type in {AUC,MCF}.")
  }

  qq_df <- data %>%
    dplyr::group_by(
      n,
      tau,
      type
    ) %>%
    dplyr::group_modify(
      ~ PrepQQDF(
        pvals = .x$pvalue_arm,
        p_max = p_max
      )
    ) %>%
    dplyr::ungroup()

  band_sizes <- data %>%
    dplyr::group_by(
      n,
      tau,
      type
    ) %>%
    dplyr::summarise(
      n_sim = dplyr::n(),
      .groups = "drop"
    )

  band_list <- lapply(
    seq_len(nrow(band_sizes)),
    function(i) {
      band_df <- PrepBandsDF(
        n_sim = band_sizes$n_sim[i],
        p_max = p_max
      )
      out <- data.frame(
        n = band_sizes$n[i],
        tau = band_sizes$tau[i],
        type = band_sizes$type[i],
        exp = band_df$exp,
        lower = band_df$lower,
        upper = band_df$upper
      )
      return(out)
    }
  )

  bands_df <- dplyr::bind_rows(band_list)

  qq_df <- AddPlotLabels(qq_df)
  bands_df <- AddPlotLabels(bands_df)

  out <- list(
    qq_df = qq_df,
    bands_df = bands_df
  )
  return(out)
}

#' Build the faceted QQ plot.
#'
#' @param qq_df QQ point data.
#' @param bands_df QQ confidence-band data.
#' @return ggplot object.
BuildQQPlot <- function(
    qq_df,
    bands_df
) {
  gg_theme <- theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      axis.line = element_line(color = "black")
    )

  out <- ggplot(
    data = qq_df,
    mapping = aes(
      x = exp,
      y = obs,
      color = N
    )
  ) +
    gg_theme +
    facet_grid(
      rows = vars(tau_label),
      cols = vars(type_label),
      labeller = labeller(tau_label = label_parsed)
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed"
    ) +
    geom_ribbon(
      data = bands_df,
      mapping = aes(
        x = exp,
        ymin = lower,
        ymax = upper,
        group = N
      ),
      inherit.aes = FALSE,
      alpha = 0.15,
      fill = "#deebf7"
    ) +
    geom_point(
      size = 0.8,
      alpha = 0.7
    ) +
    scale_color_brewer(
      name = "N per Arm",
      palette = "Blues"
    ) +
    scale_x_continuous(
      name = expression(Expected ~ -log[10](p))
    ) +
    scale_y_continuous(
      name = expression(Observed ~ -log[10](p))
    )

  return(out)
}

opt_list <- list(
  make_option(
    c("--dir"),
    type = "character",
    help = "Directory containing parquet outputs (default: simulations/pseudoreg/t1e/null)",
    default = NULL
  ),
  make_option(
    c("--out"),
    type = "character",
    help = "Output PDF path (default: figures/fig3.pdf)",
    default = NULL
  ),
  make_option(
    c("--p_max"),
    type = "double",
    help = "Maximum p-value to include in QQ plot (default: 0.05)",
    default = 0.05
  ),
  make_option(
    c("--width"),
    type = "double",
    help = "Plot width in inches (default: 10)",
    default = 10
  ),
  make_option(
    c("--height"),
    type = "double",
    help = "Plot height in inches (default: 10)",
    default = 10
  )
)

parsed <- parse_args(
  OptionParser(option_list = opt_list)
)
reg_dir <- RegressionDir()

dir_in <- parsed$dir
if (is.null(dir_in)) {
  dir_in <- file.path(
    reg_dir,
    "simulations",
    "pseudoreg",
    "t1e",
    "null"
  )
}
if (!dir.exists(dir_in)) {
  stop("Input directory does not exist: ", dir_in)
}

out_file <- parsed$out
if (is.null(out_file)) {
  out_file <- file.path(
    FigureDir("t1e"),
    "fig3.pdf"
  )
}

p_max <- parsed$p_max
if (!is.finite(p_max) || p_max <= 0 || p_max >= 1) {
  stop("--p_max must be between 0 and 1.")
}

width <- parsed$width
height <- parsed$height

data <- LoadPvalues(dir_in)
plot_data <- PreparePlotData(
  data = data,
  p_max = p_max
)

counts_df <- data %>%
  dplyr::filter(
    n %in% c(50, 100, 200),
    tau %in% c(2, 3, 4),
    toupper(type) %in% c("AUC", "MCF")
  ) %>%
  dplyr::group_by(
    n,
    tau,
    type
  ) %>%
  dplyr::summarise(
    n_rep = dplyr::n(),
    .groups = "drop"
  )
print(counts_df)

qq_plot <- BuildQQPlot(
  qq_df = plot_data$qq_df,
  bands_df = plot_data$bands_df
)

SaveRegressionFigure(
  plot = qq_plot,
  path_pdf = out_file,
  width = width,
  height = height,
  also_png = TRUE
)
