# Shared data-generating-process helpers for pseudo-value regression simulations.
# Requires MCC to be loaded by the calling script.

#' Build the recurrent-event coefficient vector for one-at-a-time DGP sweeps.
#'
#' @param beta_name One of \code{beta_A}, \code{beta_X}, \code{beta_AX}.
#' @param beta_value Log-rate ratio applied to that coefficient.
#' @return Named numeric vector \code{c(intercept, A, X, AX)}.
MakeBetaVector <- function(beta_name, beta_value) {
  out <- c(intercept = 0, A = 0, X = 0, AX = 0)
  if (beta_name == "beta_A") {
    out["A"] <- beta_value
  } else if (beta_name == "beta_X") {
    out["X"] <- beta_value
  } else if (beta_name == "beta_AX") {
    out["AX"] <- beta_value
  } else {
    stop("Unknown beta_name: ", beta_name, call. = FALSE)
  }
  return(out)
}

#' Balanced two-arm design with standard normal covariate and interaction.
#'
#' @param n_per_arm Sample size per treatment arm.
#' @return Data frame with columns \code{idx}, \code{A}, \code{X}, \code{AX}.
MakeDesignData <- function(n_per_arm) {
  total_n <- 2L * n_per_arm
  A <- c(rep(0L, n_per_arm), rep(1L, n_per_arm))
  X <- stats::rnorm(total_n, mean = 0, sd = 1)
  AX <- A * X
  out <- data.frame(
    idx = seq_len(total_n),
    A = A,
    X = X,
    AX = AX
  )
  return(out)
}

#' Map a varied DGP coefficient name to the matching \code{lm} term label.
#'
#' @param beta_name One of \code{beta_A}, \code{beta_X}, \code{beta_AX}.
#' @return Character scalar: \code{A}, \code{X}, or \code{A:X}.
BetaNameToTerm <- function(beta_name) {
  term_map <- c(beta_A = "A", beta_X = "X", beta_AX = "A:X")
  out <- term_map[beta_name]
  if (is.na(out)) {
    stop("Unknown beta_name: ", beta_name, call. = FALSE)
  }
  return(unname(out))
}

#' Confirm that the loaded MCC exposes the required PseudoReg inference.
#'
#' @return Invisibly, the expected variance-method label.
ValidatePseudoReg <- function() {
  if (!exists("PseudoReg", where = asNamespace("MCC"), inherits = FALSE)) {
    stop(
      "The loaded MCC package does not contain PseudoReg(). ",
      "Install the current MCC source before running this simulation.",
      call. = FALSE
    )
  }
  expected_method <- CorrectedVarianceMethod()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) saved_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv), add = TRUE)
  set.seed(2026103001L)
  arm <- rep(0:1, each = 10L)
  data <- MCC::GenData(covariates = cbind(intercept = 1, arm = arm), beta_event = c(0, .1), beta_death = c(0, 0), tau = 3)
  fit <- MCC::PseudoReg(~ arm, data = data, tau = 1, covariates = data.frame(idx = seq_along(arm), arm = arm))
  if (!identical(fit$variance_method, expected_method)) stop("Install MCC with HC3 + second order inference; the version number alone is insufficient.", call. = FALSE)
  if (any(!is.finite(fit$coefficients$std_error)) || any(fit$coefficients$std_error <= 0)) stop("MCC capability fit returned invalid SEs.", call. = FALSE)
  bread <- solve(fit$Ahat)
  manual <- bread %*% crossprod(fit$h11_hc3 + fit$h12) %*% t(bread) / length(arm)^2
  if (max(abs(manual - fit$vcov)) > 1e-10) stop("MCC capability covariance check failed.", call. = FALSE)
  return(invisible(expected_method))
}

#' Extract one coefficient row from a PseudoReg fit.
#'
#' @param fit Object returned by \code{MCC::PseudoReg()}.
#' @param term Coefficient term to extract.
#' @return One-row data frame with the stable Regression simulation schema.
ExtractPseudoRegCoefficient <- function(fit, term) {
  expected_method <- "HC3 + second order"
  if (!inherits(fit, "PseudoReg")) {
    stop("fit must inherit from PseudoReg.", call. = FALSE)
  }
  if (!identical(fit$variance_method, expected_method)) {
    stop(
      "Unexpected PseudoReg variance method: ",
      fit$variance_method,
      ". Expected: ",
      expected_method,
      ".",
      call. = FALSE
    )
  }
  coefficient_row <- fit$coefficients[
    fit$coefficients$term == term,
    ,
    drop = FALSE
  ]
  if (nrow(coefficient_row) != 1L) {
    stop("Expected exactly one PseudoReg coefficient row for: ", term, call. = FALSE)
  }
  required <- c(
    "estimate",
    "std_error",
    "lower",
    "upper",
    "p_value"
  )
  if (!all(required %in% names(coefficient_row))) {
    stop("PseudoReg coefficient table has an unexpected schema.", call. = FALSE)
  }
  values <- unlist(coefficient_row[1L, required], use.names = FALSE)
  if (any(!is.finite(values))) {
    stop("PseudoReg returned non-finite inference for term: ", term, call. = FALSE)
  }
  out <- data.frame(
    term = term,
    estimate = coefficient_row$estimate,
    se = coefficient_row$std_error,
    lower = coefficient_row$lower,
    upper = coefficient_row$upper,
    pvalue = coefficient_row$p_value,
    variance_method = fit$variance_method,
    check.names = FALSE
  )
  return(out)
}

#' Derive a deterministic seed for one simulation replicate.
#'
#' @param base_seed Non-negative integer base seed.
#' @param config_index One-based configuration index.
#' @param rep_index One-based global replicate index.
#' @return Integer seed accepted by \code{set.seed()}.
RegressionRepSeed <- function(base_seed, config_index, rep_index) {
  if (
    length(base_seed) != 1L ||
      !is.finite(base_seed) ||
      base_seed < 0 ||
      length(config_index) != 1L ||
      !is.finite(config_index) ||
      config_index < 1 ||
      length(rep_index) != 1L ||
      !is.finite(rep_index) ||
      rep_index < 1
  ) {
    stop("Seed inputs must be finite non-negative or positive scalars.", call. = FALSE)
  }
  modulus <- .Machine$integer.max
  raw_seed <- (
    as.double(base_seed) +
      as.double(config_index) * 1000003 +
      as.double(rep_index) * 1009
  ) %% modulus
  out <- as.integer(max(1, raw_seed))
  return(out)
}

#' Create a filesystem-safe representation of a numeric value.
#'
#' @param value Numeric scalar.
#' @return Character scalar.
NumericSlug <- function(value) {
  out <- format(signif(value, 15), scientific = FALSE, trim = TRUE)
  out <- gsub("-", "m", out, fixed = TRUE)
  out <- gsub(".", "p", out, fixed = TRUE)
  return(out)
}

#' Apply a simulation function using deterministic forked workers.
#'
#' Replicate functions must set their own deterministic seed. Forking is used
#' on Unix-like platforms; Windows falls back to sequential evaluation.
#'
#' @param items Vector or list of work items.
#' @param fun Function applied to each item.
#' @param workers Number of worker processes.
#' @return List of results in the same order as \code{items}.
ParallelSimulationLapply <- function(items, fun, workers) {
  if (!is.finite(workers) || workers < 1L) {
    stop("workers must be a positive integer.", call. = FALSE)
  }
  workers <- as.integer(workers)
  if (workers == 1L || .Platform$OS.type == "windows") {
    out <- lapply(items, fun)
  } else {
    out <- parallel::mclapply(
      items,
      fun,
      mc.cores = workers,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  }
  failed <- vapply(out, inherits, logical(1L), what = "try-error")
  if (any(failed)) {
    messages <- unique(as.character(out[failed]))
    stop(
      "Parallel simulation worker failure: ",
      paste(messages, collapse = " | "),
      call. = FALSE
    )
  }
  return(out)
}
