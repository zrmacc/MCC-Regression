# Shared plotting functions; no command-line execution or output side effects.
#' Convert a beta name to a parsed facet label for either estimand.
#'
#' @param beta_name Character scalar (beta_A, beta_X, or beta_AX from the DGP).
#' @param type Pseudo-value type, retained for compatibility; both use beta.
#' @return Character scalar.
BetaNameLabel <- function(beta_name, type) {
  prefix <- "beta"
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

#' Load corrected estimation simulation results.
#'
#' Expects columns: term, lower, upper, estimate, n, tau, type, beta_name,
#'   beta_value. If "covers" is present it is used; otherwise ground truth
#'   must be supplied and covers is computed at plot time.
#'
#' @param dir_in Directory containing parquet outputs.
#' @return Data frame.
LoadCoverageData <- function(dir_in) {
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
    "term",
    "lower",
    "upper",
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

  ValidateCorrectedInferenceData(out)
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

#' Join estimation simulation output to ground truth and compute coverage.
#'
#' Adds type_gt (AUMCF for type AUC, MCF for type MCF), joins to GT long,
#' and sets covers = (truth >= lower) & (truth <= upper).
#'
#' @param data Estimation simulation output with term, type, lower, upper, beta_name, beta_value, tau.
#' @param gt_long Ground truth in long form from GroundTruthLong().
#' @return Data frame with covers column.
JoinGroundTruthAndCoverage <- function(data, gt_long) {
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
    dplyr::mutate(
      covers = (.data$truth >= .data$lower) & (.data$truth <= .data$upper)
    )
  return(out)
}

#' Aggregate coverage data.
#'
#' Computes coverage proportion per (type, n, tau, beta_name, beta_value).
#' Filters to the term matching beta_name (e.g., beta_A -> term A).
#'
#' @param data Raw coverage data.
#' @return Aggregated data frame with coverage column.
AggregateCoverage <- function(data) {
  # Map beta_name to the corresponding term
  term_map <- c(
    beta_A = "A",
    beta_X = "X",
    beta_AX = "A:X"
  )
  
  # Filter to matching term only
  data <- data %>%
    dplyr::filter(term == term_map[beta_name])
  
  out <- data %>%
    dplyr::group_by(
      type,
      n,
      tau,
      beta_name,
      beta_value
    ) %>%
    dplyr::summarise(
      coverage = mean(covers),
      n_reps = dplyr::n(),
      se = sqrt(coverage * (1 - coverage) / n_reps),
      .groups = "drop"
    )
  return(out)
}

#' Add plotting labels and factors.
#'
#' @param data Data frame with aggregated coverage results.
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
          "beta[A]", "beta[Z]", "beta[AZ]"
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

#' Build a faceted bar plot for one pseudo-value type.
#'
#' @param data Plotting data for one type.
#' @param base_size Base font size in points.
#' @param n_breaks Retained for compatibility; coverage ticks are fixed at 5 percentage points.
#' @return ggplot object.
BuildCoveragePlot <- function(data, base_size = 11, n_breaks = 5) {
  gg_theme <- theme_bw(base_size = base_size) +
    theme(
      legend.position = "top",
      panel.spacing.y = grid::unit(1, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      axis.line = element_line(color = "black")
    )

  out <- ggplot(data = data) +
    gg_theme +
    geom_col(
      mapping = aes(
        x = beta_value_label,
        y = coverage,
        fill = N
      ),
      position = position_dodge(width = 0.8),
      width = 0.7,
      alpha = 0.8
    ) +
    geom_errorbar(
      mapping = aes(
        x = beta_value_label,
        ymin = coverage - 1.96 * se,
        ymax = coverage + 1.96 * se,
        group = N
      ),
      position = position_dodge(width = 0.8),
      width = 0.2,
      linewidth = 0.3
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
      yintercept = 0.95,
      linetype = "dashed",
      linewidth = 0.3,
      color = "grey50"
    ) +
    scale_x_discrete(name = "Rate Parameter") +
    scale_y_continuous(
      name = "Coverage",
      breaks = seq(0.80, 1.00, by = 0.05),
      labels = function(x) paste0(round(100 * x), "%"),
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(ylim = c(0.80, 1.00))

  return(out)
}

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------

