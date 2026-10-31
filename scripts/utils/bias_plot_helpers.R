# Shared plotting functions; no command-line execution or output side effects.
#' Convert a beta name to a parsed facet label (MCF: beta, AUMCF: gamma).
#'
#' @param beta_name Character scalar (beta_A, beta_X, or beta_AX from the DGP).
#' @param type Pseudo-value type: \code{AUC} (AUMCF) uses gamma, \code{MCF} uses beta.
#' @return Character scalar.
BetaNameLabel <- function(beta_name, type) {
  use_gamma <- toupper(type) == "AUC"
  prefix <- if (use_gamma) {
    "gamma"
  } else {
    "beta"
  }
  if (beta_name == "beta_A") {
    out <- paste0(prefix, "[A]")
  } else if (beta_name == "beta_X") {
    out <- paste0(prefix, "[Z]")
  } else if (beta_name == "beta_AX") {
    out <- paste0(prefix, "[AZ]")
  } else {
    out <- beta_name
  }
  return(out)
}

#' Convert a beta value to a discrete axis label.
#'
#' @param beta_value Numeric scalar on the log scale.
#' @return Character scalar.
BetaValueLabel <- function(beta_value) {
  rr <- exp(beta_value)
  if (abs(rr - 0.75) < 1e-8) {
    out <- "ln(0.75)"
  } else if (abs(rr - 1.00) < 1e-8) {
    out <- "ln(1.0)"
  } else if (abs(rr - 1.25) < 1e-8) {
    out <- "ln(1.25)"
  } else {
    rr_lab <- format(signif(rr, 3), trim = TRUE)
    out <- paste0("ln(", rr_lab, ")")
  }
  return(out)
}

#' Load combined estimation simulation (or legacy bias) results.
#'
#' Expects columns: rep, term, estimate, se, lower, upper, pvalue, n, tau, type,
#'   beta_name, beta_value. If a "bias" column is present, it is used; otherwise
#'   ground truth must be supplied and bias is computed at plot time.
#'
#' @param dir_in Directory containing parquet outputs.
#' @return Data frame.
LoadBiasData <- function(dir_in) {
  files <- ListMergedSimulationParquets(dir_in)
  if (length(files) == 0L) {
    stop("No parquet files found in ", dir_in)
  }

  out <- bind_rows(
    lapply(
      files,
      read_parquet
    )
  )

  req_cols <- c(
    "estimate",
    "term",
    "beta_name",
    "beta_value",
    "n",
    "tau",
    "type"
  )
  miss_cols <- setdiff(req_cols, names(out))
  if (length(miss_cols) > 0L) {
    stop(
      "Combined data are missing required columns: ",
      paste(miss_cols, collapse = ", ")
    )
  }

  return(out)
}

#' Reshape ground truth table to long (term, type_gt, truth).
#'
#' Ground truth has one row per (beta_name, beta_value, tau) and columns
#' beta_A_MCF, beta_A_AUMCF, beta_X_MCF, beta_X_AUMCF, beta_AX_MCF, beta_AX_AUMCF.
#'
#' @param gt Data frame with ground truth schema.
#' @return Data frame with columns beta_name, beta_value, tau, term, type_gt, truth.
GroundTruthLong <- function(gt) {
  spec <- list(
    beta_A_MCF = c("A", "MCF"),
    beta_A_AUMCF = c("A", "AUMCF"),
    beta_X_MCF = c("X", "MCF"),
    beta_X_AUMCF = c("X", "AUMCF"),
    beta_AX_MCF = c("A:X", "MCF"),
    beta_AX_AUMCF = c("A:X", "AUMCF")
  )
  out_list <- lapply(
    names(spec),
    function(col_nm) {
      gt %>%
        dplyr::select(beta_name, beta_value, tau, truth = dplyr::all_of(col_nm)) %>%
        dplyr::mutate(
          term = spec[[col_nm]][1L],
          type_gt = spec[[col_nm]][2L]
        )
    }
  )
  out <- bind_rows(out_list)
  return(out)
}

#' Join estimation simulation output to ground truth and compute bias.
#'
#' Adds type_gt (AUMCF for type AUC, MCF for type MCF), joins to GT long,
#' and sets bias = estimate - truth.
#'
#' @param data Estimation simulation output with term, type, beta_name, beta_value, tau.
#' @param gt_long Ground truth in long form from GroundTruthLong().
#' @return Data frame with bias column.
JoinGroundTruthAndBias <- function(data, gt_long) {
  data <- data %>%
    dplyr::mutate(
      type_gt = dplyr::if_else(toupper(.data$type) == "AUC", "AUMCF", "MCF"),
      beta_value_key = GroundTruthJoinKey(.data$beta_value)
    )
  gt_long <- gt_long %>%
    dplyr::mutate(beta_value_key = GroundTruthJoinKey(.data$beta_value)) %>%
    dplyr::select(-"beta_value")
  out <- data %>%
    dplyr::left_join(
      gt_long,
      by = c("beta_name", "beta_value_key", "tau", "term", "type_gt")
    ) %>%
    dplyr::mutate(bias = .data$estimate - .data$truth)
  return(out)
}

#' Filter to matching term.
#'
#' Keeps only rows where term matches beta_name (e.g., beta_A -> term A).
#'
#' @param data Data frame with bias results.
#' @return Filtered data frame.
FilterMatchingTerm <- function(data) {
  term_map <- c(
    beta_A = "A",
    beta_X = "X",
    beta_AX = "A:X"
  )
  
  out <- data %>%
    dplyr::filter(term == term_map[beta_name])
  return(out)
}

#' Add plotting labels and factors.
#'
#' @param data Data frame with bias results.
#' @return Data frame.
AddPlotLabels <- function(data) {
  beta_value_levels <- sort(unique(data$beta_value))
  beta_value_labels <- vapply(
    beta_value_levels,
    BetaValueLabel,
    character(1)
  )

  out <- data %>%
    dplyr::mutate(
      type = toupper(type),
      type_label = dplyr::recode(
        type,
        AUC = "AUMCF",
        MCF = "MCF"
      ),
      tau_label = factor(
        paste0("tau == ", tau),
        levels = c("tau == 2", "tau == 3", "tau == 4")
      ),
      beta_name_label = factor(
        mapply(
          BetaNameLabel,
          beta_name,
          type,
          SIMPLIFY = TRUE,
          USE.NAMES = FALSE
        ),
        levels = c(
          "beta[A]", "beta[Z]", "beta[AZ]",
          "gamma[A]", "gamma[Z]", "gamma[AZ]"
        )
      ),
      beta_value_label = factor(
        vapply(beta_value, BetaValueLabel, character(1)),
        levels = beta_value_labels
      ),
      N = factor(
        n,
        levels = c(50, 100, 200)
      )
    )
  return(out)
}

#' Build a faceted boxplot for one pseudo-value type.
#'
#' @param data Plotting data for one type.
#' @return ggplot object.
BuildBiasPlot <- function(data) {
  gg_theme <- theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      axis.line = element_line(color = "black")
    )

  out <- ggplot(data = data) +
    gg_theme +
    geom_boxplot(
      mapping = aes(
        x = beta_value_label,
        y = bias,
        fill = N
      ),
      alpha = 0.6,
      outlier.size = 0.5,
      size = 0.3
    ) +
    facet_grid(
      rows = vars(tau_label),
      cols = vars(beta_name_label),
      labeller = labeller(
        tau_label = label_parsed,
        beta_name_label = label_parsed
      )
    ) +
    scale_fill_brewer(
      name = "N per Arm",
      palette = "Blues"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.2
    ) +
    scale_x_discrete(name = "Rate Parameter") +
    scale_y_continuous(name = "Bias")

  return(out)
}

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------

