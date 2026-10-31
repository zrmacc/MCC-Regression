#!/usr/bin/env Rscript
# Figure S1: exact-jackknife versus empirical-influence-function pseudo-values.

argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(dirname(normalizePath(script))), "utils", "bootstrap_regression.R"))
SetRegressionWd()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(optparse)
})
#' Default RDS written by scripts/sims/sim_pseudovalue.R.
#'
#' @return Character path.
DefaultPseudovalRds <- function() {
  file.path(RegressionDir(), "simulations", "pseudovalue", "pseudovalue_calibration.rds")
}

#' Shared theme aligned with other Regression figures.
#'
#' @return ggplot2 theme.
PseudoCalibTheme <- function() {
  out <- theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      axis.line = element_line(color = "black")
    )
  return(out)
}

#' Format Pearson \eqn{R^2} for \code{geom_text(parse = TRUE)} (two significant digits).
#'
#' Uses \code{paste()} so the numeric part is shown as typed (e.g. \code{1.00}, not plotmath
#' shortening to \code{1}).
#'
#' @param r2 Scalar; squared correlation of jackknife vs influence-function pseudo-values.
#' @return Character plotmath string, or \code{""} if not finite.
.FormatR2ConcordanceLabel <- function(r2) {
  if (!is.finite(r2)) {
    return("")
  }
  r2s <- signif(r2, 2)
  if (r2s == 1) {
    r2_chr <- "1.00"
  } else {
    r2_chr <- as.character(r2s)
  }
  out <- paste0("paste(italic(R)^2, \" = \", \"", r2_chr, "\")")
  return(out)
}

#' Shared axis limits for concordance (identical x and y so \eqn{y = x} is at 45 degrees).
#'
#' @param df Data with \code{pseudo_jk} and \code{pseudo_inf}.
#' @param pad_mult Padding as a fraction of the observed span.
#' @return Named numeric: \code{lo}, \code{hi}.
PseudoCalibSquareAxisLimits <- function(df, pad_mult = 0.02) {
  rng <- range(c(df$pseudo_jk, df$pseudo_inf), finite = TRUE)
  span <- unname(diff(rng))
  if (!is.finite(span) || span < 1e-12) {
    mid <- mean(rng, na.rm = TRUE)
    if (!is.finite(mid)) {
      mid <- 0
    }
    pad <- pmax(abs(mid) * 0.02, 1e-6)
    lo_lim <- mid - pad
    hi_lim <- mid + pad
  } else {
    pad <- pad_mult * span
    lo_lim <- rng[1L] - pad
    hi_lim <- rng[2L] + pad
  }
  out <- c(lo = lo_lim, hi = hi_lim)
  return(out)
}

#' Per-facet \eqn{R^2} annotation tibble for concordance plots.
#'
#' Each facet gets its own \eqn{R^2} string; \code{label_x} and \code{label_y} are
#' shared (relative 0.2 / 0.8 on the square limits used for both axes).
#'
#' @param df Data frame with \code{tau_label}, \code{N_arm}, \code{pseudo_jk},
#'   \code{pseudo_inf}.
#' @param label_x Scalar x position in data units.
#' @param label_y Scalar y position in data units.
#' @return Data frame with facet columns plus \code{label} for \code{parse = TRUE}.
PseudoCalibR2Labels <- function(df, label_x, label_y) {
  lx <- label_x
  ly <- label_y
  ann <- df %>%
    dplyr::group_by(tau_label, N_arm) %>%
    dplyr::summarise(
      r2 = stats::cor(pseudo_jk, pseudo_inf)^2,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = vapply(
        r2,
        .FormatR2ConcordanceLabel,
        character(1L)
      ),
      label_x = lx,
      label_y = ly
    ) %>%
    dplyr::filter(nzchar(label)) %>%
    dplyr::select(
      tau_label,
      N_arm,
      label_x,
      label_y,
      label
    )
  return(ann)
}

#' Per-cell \eqn{R^2} for mean reporting (same cells as facets: \code{type}, \code{tau}, \code{n}).
#'
#' @param plot_df Data with \code{type}, \code{tau}, \code{n}, \code{pseudo_jk}, \code{pseudo_inf}.
#' @return Data frame with one \eqn{R^2} per group.
PseudoCalibCellR2 <- function(plot_df) {
  out <- plot_df %>%
    dplyr::group_by(type, tau, n) %>%
    dplyr::summarise(
      r2 = stats::cor(pseudo_jk, pseudo_inf)^2,
      .groups = "drop"
    )
  return(out)
}

#' Build one calibration plot for a single pseudo-value type.
#'
#' @param df Rows for one \code{type} (AUC or MCF); columns include pseudo_jk,
#'   pseudo_inf, tau, n.
#' @return ggplot object.
BuildPseudoCalibrationPlot <- function(df) {
  tau_levels <- sort(unique(df$tau))
  tau_label_levels <- paste0("tau == ", tau_levels)
  df <- df %>%
    dplyr::mutate(
      tau_label = factor(
        paste0("tau == ", tau),
        levels = tau_label_levels
      ),
      N_arm = factor(
        n,
        levels = sort(unique(n))
      )
    )

  lim <- PseudoCalibSquareAxisLimits(df)
  lo_lim <- lim[["lo"]]
  hi_lim <- lim[["hi"]]
  span_lim <- hi_lim - lo_lim
  label_x_g <- lo_lim + 0.2 * span_lim
  label_y_g <- lo_lim + 0.8 * span_lim
  r2_ann <- PseudoCalibR2Labels(
    df = df,
    label_x = label_x_g,
    label_y = label_y_g
  )

  n_labeller <- function(value) {
    paste0("n = ", value)
  }

  out <- ggplot(data = df) +
    PseudoCalibTheme() +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "gray40"
    ) +
    geom_point(
      mapping = aes(
        x = pseudo_jk,
        y = pseudo_inf,
        color = N_arm
      ),
      alpha = 0.75,
      size = 0.9
    ) +
    geom_text(
      data = r2_ann,
      mapping = aes(
        x = label_x,
        y = label_y,
        label = label
      ),
      inherit.aes = FALSE,
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      color = "gray25",
      show.legend = FALSE
    ) +
    scale_color_brewer(
      name = "N per Arm",
      palette = "Blues"
    ) +
    facet_grid(
      rows = vars(tau_label),
      cols = vars(N_arm),
      labeller = labeller(
        tau_label = label_parsed,
        N_arm = n_labeller
      )
    ) +
    scale_x_continuous(
      limits = c(lo_lim, hi_lim),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(lo_lim, hi_lim),
      expand = c(0, 0)
    ) +
    coord_fixed(ratio = 1) +
    labs(
      x = "Jackknife Pseudovalue",
      y = "Influence Function\nPseudovalue"
    )
  return(out)
}

#' Load and validate RDS from Development simulations.
#'
#' @param path Path to \code{pseudovalue_calibration.rds}.
#' @return Data frame.
ReadPseudovalCalibration <- function(path) {
  if (!file.exists(path)) {
    stop(
      "RDS not found: ", path, "\n",
      "Run scripts/sims/sim_pseudovalue.R to create it, ",
      "or pass --data to an existing file."
    )
  }
  data <- readRDS(path)
  req <- c("pseudo_jk", "pseudo_inf", "type", "tau", "n")
  miss <- setdiff(req, names(data))
  if (length(miss) > 0L) {
    stop(
      "RDS missing columns: ",
      paste(miss, collapse = ", ")
    )
  }
  return(data)
}

opt_list <- list(
  make_option(
    c("--data"),
    type = "character",
    help = paste0(
      "Path to pseudovalue_calibration.rds ",
      "(default: simulations/pseudovalue/pseudovalue_calibration.rds)"
    ),
    default = NULL
  ),
  make_option(
    c("--out_dir"),
    type = "character",
    help = "Directory for PDF outputs (default: figures)",
    default = NULL
  ),
  make_option(
    c("--width"),
    type = "double",
    help = "Figure width in inches (default: 10)",
    default = 10
  ),
  make_option(
    c("--height"),
    type = "double",
    help = "Figure height in inches (default: 11)",
    default = 11
  )
)

parsed <- parse_args(OptionParser(option_list = opt_list))
reg_dir <- RegressionDir()

data_path <- parsed$data
if (is.null(data_path)) {
  data_path <- DefaultPseudovalRds()
}

out_dir <- parsed$out_dir
if (is.null(out_dir)) {
  out_dir <- FigureDir("pseudoval")
}

width <- parsed$width
height <- parsed$height

raw <- ReadPseudovalCalibration(data_path)

plot_df <- raw %>%
  dplyr::filter(
    n %in% c(50, 100, 200),
    tau %in% c(2, 3, 4),
    toupper(type) %in% c("AUC", "MCF")
  )

if (nrow(plot_df) == 0L) {
  stop(
    "No rows after filtering to n in {50, 100, 200}, tau in {2, 3, 4}, ",
    "and type in {AUC, MCF}."
  )
}

r2_by_cell <- PseudoCalibCellR2(plot_df)
mean_r2 <- mean(r2_by_cell$r2)
message(
  "Mean R^2 across conditions (unrounded numeric): ",
  sprintf("%.17g", mean_r2)
)

plotspec <- list(
  list(
    type_key = "AUC",
    base_name = "figS1b"
  ),
  list(
    type_key = "MCF",
    base_name = "figS1a"
  )
)

for (spec in plotspec) {
  type_df <- plot_df %>%
    dplyr::filter(toupper(type) == spec$type_key)
  if (nrow(type_df) == 0L) {
    warning(
      "No rows for type ",
      spec$type_key,
      "; skipping ",
      spec$base_name
    )
    next
  }
  p <- BuildPseudoCalibrationPlot(df = type_df)
  SaveRegressionFigure(
    plot = p,
    path_pdf = file.path(out_dir, paste0(spec$base_name, ".pdf")),
    width = width,
    height = height,
    also_png = TRUE
  )
}

message("Done.")
