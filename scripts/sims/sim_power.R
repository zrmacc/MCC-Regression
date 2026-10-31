#!/usr/bin/env Rscript
# Supports Figure 4: corrected normal-Wald power.

argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(dirname(normalizePath(script))), "utils", "bootstrap_regression.R"))
SetRegressionWd()

suppressPackageStartupMessages({
  library(MCC)
  library(arrow)
  library(dplyr)
  library(glue)
  library(optparse)
})
SourceRegressionUtil("regression_dgp.R")

# -----------------------------------------------------------------------------
# Command-line helpers
#' Parse the type option.
#'
#' @param type_string Comma-separated type string.
#' @return Character vector with values in {"AUC", "MCF"}.
ParseType <- function(type_string) {
  out <- trimws(strsplit(type_string, ",")[[1L]])
  out <- toupper(out)
  out <- out[out %in% c("AUC", "MCF")]
  if (length(out) == 0L) {
    stop("--type must include at least one of: AUC, MCF")
  }
  return(out)
}

#' Parse a numeric comma-separated option.
#'
#' @param value_string Comma-separated numeric string.
#' @return Numeric vector.
ParseNumericOption <- function(value_string) {
  out <- as.numeric(trimws(strsplit(value_string, ",")[[1L]]))
  if (anyNA(out)) {
    stop("Failed to parse a numeric option value.")
  }
  return(out)
}

# -----------------------------------------------------------------------------
# Simulation helpers
# -----------------------------------------------------------------------------

#' Build the event-model coefficient vector for power study.
#'
#' Fixed: beta_X = ln(1.25), beta_AX = 0. beta_A is varied.
#'
#' @param beta_A Numeric value for treatment effect.
#' @return Named numeric vector of coefficients.
MakePowerBetaVector <- function(beta_A) {
  out <- c(
    intercept = 0,
    A = beta_A,
    X = log(1.25),
    AX = 0
  )
  return(out)
}

#' Create a file-safe beta tag.
#'
#' @param beta_value Numeric beta value on the log scale.
#' @return Character scalar.
BetaTag <- function(beta_value) {
  rr <- exp(beta_value)
  rr_lab <- sprintf("%.2f", rr)
  rr_lab <- gsub("\\.", "p", rr_lab)
  out <- paste0("rr_", rr_lab)
  return(out)
}

#' Generate subject-level design variables.
#'
#' @param n_per_arm Sample size per arm.
#' @return Data frame with one row per subject.
MakeDesignData <- function(n_per_arm) {
  total_n <- 2L * n_per_arm
  A <- c(rep(0L, n_per_arm), rep(1L, n_per_arm))
  X <- stats::rnorm(total_n, mean = 0, sd = 1)
  AX <- A * X
  out <- data.frame(idx = seq_len(total_n), A = A, X = X, AX = AX)
  return(out)
}

#' Run a single simulation replicate.
#'
#' Fits PseudoReg(~ A + X + A:X) and extracts the p-value for coefficient on A.
#' Converts p-value to chi-squared statistic as NCP proxy.
#'
#' @param n Sample size per arm.
#' @param tau Pseudo-value truncation time.
#' @param tau_gen Data-generation truncation time.
#' @param type Pseudo-value type (AUC or MCF).
#' @param beta_A Treatment effect (log rate ratio).
#' @param seed Deterministic replicate seed.
#' @return Data frame with one row: pvalue, chi2.
SimInst <- function(
  n,
  tau,
  tau_gen,
  type,
  beta_A,
  seed
) {
  set.seed(seed)
  design_df <- MakeDesignData(n_per_arm = n)
  covariates <- as.matrix(
    data.frame(
      intercept = 1,
      A = design_df$A,
      X = design_df$X,
      AX = design_df$AX
    )
  )

  beta_event <- MakePowerBetaVector(beta_A)
  beta_death <- rep(0, 4)

  data <- GenData(
    covariates = covariates,
    beta_event = beta_event,
    beta_death = beta_death,
    tau = tau_gen
  )

  fit <- MCC::PseudoReg(
    ~ A + X + A:X,
    data = data,
    tau = tau,
    type = type,
    covariates = design_df
  )
  inference <- ExtractPseudoRegCoefficient(fit = fit, term = "A")
  pval <- inference$pvalue
  if (is.na(pval)) {
    stop("P-value for coefficient on A is NA.")
  }

  chi2 <- stats::qchisq(pval, df = 1, lower.tail = FALSE)

  out <- data.frame(
    pvalue = pval,
    chi2 = chi2,
    seed = seed,
    variance_method = inference$variance_method
  )
  return(out)
}

#' Run many simulation replicates.
#'
#' @param n Sample size per arm.
#' @param tau Pseudo-value truncation time.
#' @param tau_gen Data-generation truncation time.
#' @param type Pseudo-value type.
#' @param beta_A Treatment effect.
#' @param reps Number of replicates.
#' @param config_index One-based configuration index.
#' @param base_seed Base seed.
#' @param workers Number of parallel worker processes.
#' @param testing Logical; cap reps at 100 when TRUE.
#' @return Data frame with reps rows.
Simulation <- function(
  n,
  tau,
  tau_gen,
  type,
  beta_A,
  reps,
  config_index,
  base_seed,
  workers,
  testing = FALSE
) {
  if (testing) {
    reps <- min(reps, 100L)
  }

  message(glue("Sample size per arm n = {n}"))
  message(glue("Tau (pseudo) = {tau}, tau_gen = {tau_gen}, type = {type}"))
  message(glue("beta_A = {signif(beta_A, 6)} (RR = {signif(exp(beta_A), 3)})"))
  message(glue("Replicates = {reps}"))
  if (testing) {
    message("Testing mode (reps capped at 100)")
  }

  pct_marks <- seq(0.1, 1, by = 0.1)
  milestones <- unique(c(round(pct_marks * reps), reps))
  tic <- Sys.time()

  data_list <- ParallelSimulationLapply(
    seq_len(reps),
    function(i) {
      if (i %in% milestones) {
        pct <- round(100 * i / reps, 1)
        toc <- Sys.time()
        delta <- round(as.numeric(difftime(toc, tic, units = "secs")), 2)
        message(glue("Instance {i}/{reps} ({pct}%). Elapsed: {delta}s"))
        tic <<- toc
      }

      out <- SimInst(
        n = n,
        tau = tau,
        tau_gen = tau_gen,
        type = type,
        beta_A = beta_A,
        seed = RegressionRepSeed(
          base_seed = base_seed,
          config_index = config_index,
          rep_index = i
        )
      )
      out$rep <- i
      return(out)
    },
    workers = workers
  )

  out <- dplyr::bind_rows(data_list)
  return(out)
}

# -----------------------------------------------------------------------------
# Command-line options
# -----------------------------------------------------------------------------

default_beta_A <- paste(
  signif(log(c(0.6, 0.7, 0.8, 0.9)), 15),
  collapse = ","
)

opt_list <- list(
  make_option(
    c("--n"),
    type = "character",
    help = "Sample size(s) per arm, comma-separated",
    default = "50,75,100,125,150,175,200"
  ),
  make_option(
    c("--tau"),
    type = "character",
    help = "Pseudo-value truncation time(s), comma-separated",
    default = "2,3,4"
  ),
  make_option(
    c("--tau_gen"),
    type = "numeric",
    help = "Truncation time for data generation (fixed; default 4)",
    default = 4
  ),
  make_option(
    c("--type"),
    type = "character",
    help = "Pseudo-value type(s), comma-separated: AUC, MCF",
    default = "AUC,MCF"
  ),
  make_option(
    c("--beta_A"),
    type = "character",
    help = "Treatment effect value(s), comma-separated (log scale)",
    default = default_beta_A
  ),
  make_option(
    c("--reps"),
    type = "integer",
    help = "Number of simulation replicates",
    default = 500L
  ),
  make_option(
    c("--testing", "-t"),
    action = "store_true",
    help = "Short run (cap reps at 100)",
    default = FALSE
  ),
  make_option(
    c("--batch"),
    type = "character",
    help = "Batch ID for this run (default: timestamp)",
    default = NULL
  ),
  make_option(
    c("--out_dir"),
    type = "character",
    help = "Output directory under Regression (default: simulations/pseudoreg/power)",
    default = NULL
  ),
  make_option(
    c("--seed"),
    type = "integer",
    help = "Base seed",
    default = 20260904L
  ),
  make_option(
    c("--resume"),
    action = "store_true",
    help = "Skip output files that already exist",
    default = FALSE
  ),
  make_option(
    c("--workers"),
    type = "integer",
    help = "Parallel worker processes",
    default = 1L
  )
)

parsed <- parse_args(OptionParser(option_list = opt_list))
n_vec <- as.integer(ParseNumericOption(parsed$n))
tau_vec <- ParseNumericOption(parsed$tau)
type_vec <- ParseType(parsed$type)
beta_A_vec <- ParseNumericOption(parsed$beta_A)
tau_gen <- parsed$tau_gen
reps <- parsed$reps
testing <- parsed$testing
batch <- parsed$batch
out_dir <- parsed$out_dir
base_seed <- parsed$seed
resume <- isTRUE(parsed$resume)
workers <- parsed$workers

if (testing) {
  reps <- min(reps, 100L)
}

configs <- expand.grid(
  n = n_vec,
  tau = tau_vec,
  type = type_vec,
  beta_A = beta_A_vec,
  stringsAsFactors = FALSE
)
configs <- AttachCanonicalConfigIndex(configs, "power")

reg_dir <- RegressionDir()
ValidatePseudoReg()

if (is.null(batch)) {
  batch <- format(Sys.time(), "%Y%m%d_%H%M%S")
}

if (is.null(out_dir)) {
  out_dir <- file.path(reg_dir, "simulations", "pseudoreg", "power")
} else {
  out_dir <- ResolveRegressionPath(out_dir)
}

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------

message(glue("tau_gen (data generation) = {tau_gen} (fixed)"))
message(glue("Configurations: {nrow(configs)}"))
message(glue("Replicates per config: {reps}"))
if (testing) {
  message("Testing mode (reps capped at 100)")
}

for (k in seq_len(nrow(configs))) {
  expected_config <- configs[k, , drop = FALSE]
  expected_config$tau_gen <- tau_gen
  n_k <- configs$n[k]
  tau_k <- configs$tau[k]
  type_k <- configs$type[k]
  beta_A_k <- configs$beta_A[k]
  n_fmt <- format(n_k, scientific = FALSE, trim = TRUE)
  tau_fmt <- format(tau_k, scientific = FALSE, trim = TRUE)
  config_subdir <- glue("n_{n_fmt}_tau_{tau_fmt}")
  config_dir <- file.path(out_dir, config_subdir)
  dir.create(config_dir, showWarnings = FALSE, recursive = TRUE)
  beta_tag <- BetaTag(beta_A_k)
  out_file <- file.path(
    config_dir,
    glue("power_{type_k}_{beta_tag}_{batch}.parquet")
  )
  if (file.exists(out_file)) {
    if (resume) {
      ValidateResumedInference(read_parquet(out_file), expected_config, reps, base_seed)
      message("Skipping completed config: ", out_file)
      next
    }
    stop("Output already exists; use --resume or a new --batch: ", out_file)
  }

  message(
    glue(
      "--- Config {k}/{nrow(configs)}: n = {n_k}, tau = {tau_k}, ",
      "type = {type_k}, beta_A = {signif(beta_A_k, 6)} ---"
    )
  )

  results <- Simulation(
    n = n_k,
    tau = tau_k,
    tau_gen = tau_gen,
    type = type_k,
    beta_A = beta_A_k,
    reps = reps,
    config_index = configs$config_index[k],
    base_seed = base_seed,
    workers = workers,
    testing = testing
  )

  results$n <- n_k
  results$tau <- tau_k
  results$type <- type_k
  results$beta_A <- beta_A_k
  results$tau_gen <- tau_gen
  results$batch <- batch
  results$config_index <- configs$config_index[k]
  ValidateResumedInference(results, expected_config, reps, base_seed)
  write_parquet(results, out_file)
  message(glue("Wrote {out_file}"))
}

message(glue("Done. Outputs under {out_dir}"))
