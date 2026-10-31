#!/usr/bin/env Rscript
# Figure 4: corrected normal-Wald power; no ancillary noncentrality plot.

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

# -----------------------------------------------------------------------------
# Helpers
#' Parsed facet label for beta_A: beta_{A} = value.
#'
#' @param beta_A Numeric.
#' @return Character (expression for label_parsed).
BetaAFacetLabel <- function(beta_A) {
  rr <- exp(beta_A)
  if (abs(rr - 0.6) < 0.01) {
    return("beta[A] == ln(0.6)")
  }
  if (abs(rr - 0.7) < 0.01) {
    return("beta[A] == ln(0.7)")
  }
  if (abs(rr - 0.8) < 0.01) {
    return("beta[A] == ln(0.8)")
  }
  if (abs(rr - 0.9) < 0.01) {
    return("beta[A] == ln(0.9)")
  }
  return(paste0("beta[A] == ln(", sprintf("%.2f", rr), ")"))
}

#' Legend label for beta_A (simpler format for color legend).
#'
#' @param beta_A Numeric.
#' @return Character.
BetaALegendLabel <- function(beta_A) {
  rr <- exp(beta_A)
  if (abs(rr - 0.6) < 0.01) {
    return("ln(0.6)")
  }
  if (abs(rr - 0.7) < 0.01) {
    return("ln(0.7)")
  }
  if (abs(rr - 0.8) < 0.01) {
    return("ln(0.8)")
  }
  if (abs(rr - 0.9) < 0.01) {
    return("ln(0.9)")
  }
  return(paste0("ln(", sprintf("%.2f", rr), ")"))
}

#' Load power simulation data and aggregate.
#'
#' @param dir_in Directory to search for power_*.parquet recursively.
#' @return Data frame with n, tau, beta_A, type, power, ncp_proxy, n_rep, and 95% CIs.
LoadPowerData <- function(dir_in) {
  files <- list.files(
    dir_in,
    pattern = "power_.*\\.parquet$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(files) == 0L) {
    stop("No power_*.parquet files in ", dir_in)
  }

  combined <- dplyr::bind_rows(lapply(files, read_parquet))
  ValidateStudyInference(combined, "power")
  req_cols <- c("pvalue", "chi2", "n", "tau", "beta_A", "type")
  miss <- setdiff(req_cols, names(combined))
  if (length(miss) > 0L) {
    stop("Data missing columns: ", paste(miss, collapse = ", "))
  }

  summary_df <- combined %>%
    dplyr::group_by(n, tau, beta_A, type) %>%
    dplyr::summarise(
      power = mean(pvalue <= 0.05),
      ncp_proxy = mean(chi2),
      ncp_se = sd(chi2) / sqrt(dplyr::n()),
      n_rep = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      power_se = sqrt(power * (1 - power) / n_rep),
      power_lower = pmax(0, power - 1.96 * power_se),
      power_upper = pmin(1, power + 1.96 * power_se),
      ncp_lower = ncp_proxy - 1.96 * ncp_se,
      ncp_upper = ncp_proxy + 1.96 * ncp_se
    )

  beta_A_vals <- sort(unique(summary_df$beta_A))
  beta_A_labels <- vapply(beta_A_vals, BetaAFacetLabel, character(1L))
  beta_A_legend_labels <- vapply(beta_A_vals, BetaALegendLabel, character(1L))

  summary_df %>%
    dplyr::mutate(
      N = n,
      tau_label = factor(
        paste0("tau == ", tau),
        levels = c("tau == 2", "tau == 3", "tau == 4")
      ),
      beta_A_label = factor(
        vapply(beta_A, function(b) {
          i <- which.min(abs(beta_A_vals - b))
          beta_A_labels[i]
        }, character(1L)),
        levels = beta_A_labels
      ),
      beta_A_legend = factor(
        vapply(beta_A, function(b) {
          i <- which.min(abs(beta_A_vals - b))
          beta_A_legend_labels[i]
        }, character(1L)),
        levels = beta_A_legend_labels
      ),
      type_label = dplyr::recode(
        toupper(type),
        AUC = "AUMCF",
        MCF = "MCF"
      )
    )
}

#' Build a power or NCP plot with MCF and AUMCF as row facets.
#'
#' All beta_A curves on a single plot per tau, colored by beta_A.
#' Rows = type (MCF, AUMCF), Columns = tau.
#'
#' @param data Plot data (all types).
#' @param metric "power" or "ncp_proxy".
#' @return ggplot object.
BuildPowerPlot <- function(data, metric) {
  gg_theme <- theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      axis.line = element_line(color = "black")
    )

  y_var <- if (metric == "power") "power" else "ncp_proxy"
  y_lab <- if (metric == "power") "Power" else "NCP Proxy (Mean Chi-Squared)"

  p <- ggplot(
    data,
    aes(
      x = N,
      y = .data[[y_var]],
      color = beta_A_legend,
      group = beta_A_legend
    )
  ) +
    gg_theme +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.5) +
    facet_grid(
      rows = vars(type_label),
      cols = vars(tau_label),
      labeller = labeller(tau_label = label_parsed)
    ) +
    scale_color_viridis_d(
      option = "plasma",
      end = 0.85,
      direction = -1,
      name = expression(alpha[A])
    ) +
    labs(
      x = "N per Arm",
      y = y_lab
    )

  if (metric == "power") {
    p <- p +
      geom_linerange(
        aes(ymin = power_lower, ymax = power_upper),
        linewidth = 0.4
      ) +
      geom_hline(
        yintercept = 0.05,
        linetype = "dashed",
        linewidth = 0.3,
        color = "grey50"
      ) +
      scale_y_continuous(limits = c(0, 1))
  } else {
    p <- p +
      geom_errorbar(
        aes(ymin = ncp_lower, ymax = ncp_upper),
        width = 5,
        linewidth = 0.3
      )
  }

  return(p)
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

opt_list <- list(
  make_option(
    c("--dir"),
    type = "character",
    help = "Directory with power_*.parquet (default: simulations/pseudoreg/power)",
    default = NULL
  ),
  make_option(
    c("--out_dir"),
    type = "character",
    help = "Output directory for PDFs (default: figures)",
    default = NULL
  )
)
parsed <- parse_args(OptionParser(option_list = opt_list))
reg_dir <- RegressionDir()
dir_in <- parsed$dir
out_dir <- parsed$out_dir

if (is.null(dir_in)) {
  dir_in <- file.path(reg_dir, "simulations", "pseudoreg", "power")
}
if (is.null(out_dir)) {
  out_dir <- FigureDir("power")
}

plot_data <- LoadPowerData(dir_in)

metrics <- "power"

for (metric in metrics) {
  p <- BuildPowerPlot(plot_data, metric)

  out_file <- file.path(out_dir, "fig4.pdf")

  SaveRegressionFigure(
    plot = p,
    path_pdf = out_file,
    width = 9,
    height = 5
  )
}

message("Done.")
