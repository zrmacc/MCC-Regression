#!/usr/bin/env Rscript
# Supports Figures 1 and 2: shared bias and corrected-coverage simulations.

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
SourceRegressionUtil("regression_options.R")

# -----------------------------------------------------------------------------
# Simulation helpers
# -----------------------------------------------------------------------------

#' Run a single simulation replicate.
#'
#' Fits the full model \code{~ A + X + A:X} with \code{MCC::PseudoReg()} but records only the
#' coefficient that matches beta_name (the one varied in the DGP), so we
#' never compare an estimate with the wrong ground truth.
#'
#' @param n Sample size per arm.
#' @param tau Pseudo-value truncation time.
#' @param tau_gen Data-generation truncation time.
#' @param type Pseudo-value type (AUC or MCF).
#' @param beta_name Name of coefficient varied in the DGP.
#' @param beta_value Numeric value of that coefficient.
#' @param seed Deterministic replicate seed.
#' @return Data frame with one row (the term matching beta_name).
SimInst <- function(
  n,
  tau,
  tau_gen,
  type,
  beta_name,
  beta_value,
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
  beta_event <- MakeBetaVector(beta_name = beta_name, beta_value = beta_value)
  beta_death <- rep(0, 4)

  data <- GenData(
    covariates = covariates,
    beta_event = beta_event,
    beta_death = beta_death,
    tau = tau_gen
  )
  term <- BetaNameToTerm(beta_name)
  fit <- MCC::PseudoReg(
    ~ A + X + A:X,
    data = data,
    tau = tau,
    type = type,
    covariates = design_df
  )
  inference <- ExtractPseudoRegCoefficient(fit = fit, term = term)
  metadata <- data.frame(
    rep = 1L,
    beta_name = beta_name,
    beta_value = beta_value,
    seed = seed
  )
  out <- cbind(metadata, inference)
  return(out)
}

#' Run many simulation replicates for one config.
#'
#' @param n Sample size per arm.
#' @param tau Pseudo-value truncation time.
#' @param tau_gen Data-generation truncation time.
#' @param type Pseudo-value type.
#' @param beta_name Name of coefficient varied in the DGP.
#' @param beta_value Value of the varied coefficient.
#' @param reps Number of simulation replicates.
#' @param config_index One-based configuration index.
#' @param base_seed Base seed for deterministic replicate seeds.
#' @param workers Number of parallel worker processes.
#' @param testing If TRUE, cap replicates at 100.
#' @return Data frame of replicate-level inference.
Simulation <- function(
  n,
  tau,
  tau_gen,
  type,
  beta_name,
  beta_value,
  reps,
  config_index,
  base_seed,
  workers,
  testing = FALSE
) {
  if (testing) {
    reps <- min(reps, 100L)
  }
  message(glue("n = {n}, tau = {tau}, type = {type}, beta_name = {beta_name}, reps = {reps}"))
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
        beta_name = beta_name,
        beta_value = beta_value,
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
  out <- do.call(rbind, data_list)
  return(out)
}

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

default_beta_values <- paste(
  signif(log(c(0.75, 1.0, 1.25)), 15),
  collapse = ","
)

opt_list <- list(
  make_option(
    c("--n"),
    type = "character",
    help = "Sample size(s) per arm, comma-separated",
    default = "50,100,200"
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
    help = "Truncation time for data generation",
    default = 4
  ),
  make_option(
    c("--type"),
    type = "character",
    help = "Pseudo-value type(s): AUC, MCF",
    default = "AUC,MCF"
  ),
  make_option(
    c("--beta_name"),
    type = "character",
    help = "Coefficient(s) to vary: beta_A, beta_X, beta_AX",
    default = "beta_A,beta_X,beta_AX"
  ),
  make_option(
    c("--beta_value"),
    type = "character",
    help = "Log-scale coefficient value(s), comma-separated",
    default = default_beta_values
  ),
  make_option(
    c("--reps"),
    type = "integer",
    help = "Number of simulation replicates",
    default = 2000L
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
    help = "Batch ID (default: timestamp)",
    default = NULL
  ),
  make_option(
    c("--outdir"),
    type = "character",
    help = paste0(
      "Output directory under Regression relative to Regression root ",
      "(default: simulations/pseudoreg/estimation; use this name for consistency with plot scripts)"
    ),
    default = NULL
  ),
  make_option(
    c("--seed"),
    type = "integer",
    help = "Base seed for deterministic per-replicate seeds",
    default = 20260904L
  ),
  make_option(
    c("--resume"),
    action = "store_true",
    help = "Skip design-cell files that already exist",
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
beta_name_vec <- ParseBetaName(parsed$beta_name)
beta_value_vec <- ParseNumericOption(parsed$beta_value)
tau_gen <- parsed$tau_gen
reps <- parsed$reps
testing <- parsed$testing
batch <- parsed$batch
outdir <- parsed$outdir
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
  beta_name = beta_name_vec,
  beta_value = beta_value_vec,
  stringsAsFactors = FALSE
)
configs <- AttachCanonicalConfigIndex(configs, "estimation")

reg_dir <- RegressionDir()
ValidatePseudoReg()

if (is.null(batch)) {
  batch <- format(Sys.time(), "%Y%m%d_%H%M%S")
}

if (is.null(outdir)) {
  out_dir <- file.path(reg_dir, "simulations", "pseudoreg", "estimation")
} else {
  out_dir <- ResolveRegressionPath(outdir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message(glue("tau_gen = {tau_gen} (fixed)"))
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
  beta_name_k <- configs$beta_name[k]
  beta_value_k <- configs$beta_value[k]
  config_slug <- paste(
    sprintf("config_%03d", configs$config_index[k]),
    paste0("n_", n_k),
    paste0("tau_", NumericSlug(tau_k)),
    type_k,
    beta_name_k,
    paste0("value_", NumericSlug(beta_value_k)),
    batch,
    sep = "_"
  )
  out_file <- file.path(out_dir, paste0(config_slug, ".parquet"))
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
      "type = {type_k}, beta_name = {beta_name_k}, ",
      "beta_value = {signif(beta_value_k, 6)} ---"
    )
  )

  results <- Simulation(
    n = n_k,
    tau = tau_k,
    tau_gen = tau_gen,
    type = type_k,
    beta_name = beta_name_k,
    beta_value = beta_value_k,
    reps = reps,
    config_index = configs$config_index[k],
    base_seed = base_seed,
    workers = workers,
    testing = testing
  )
  results$n <- n_k
  results$tau <- tau_k
  results$type <- type_k
  results$tau_gen <- tau_gen
  results$batch <- batch
  results$config_index <- configs$config_index[k]
  ValidateResumedInference(results, expected_config, reps, base_seed)
  write_parquet(results, out_file)
  message(glue("Wrote {out_file}"))
}

message(glue("Done. Output directory: {out_dir}"))
