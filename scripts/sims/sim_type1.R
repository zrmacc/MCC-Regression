#!/usr/bin/env Rscript
# Supports Figure 3 and Table 1: corrected normal-Wald type I error, n per arm.

argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(dirname(normalizePath(script))), "utils", "bootstrap_regression.R"))
SetRegressionWd()

suppressPackageStartupMessages({
  library(MCC)
  library(glue)
  library(arrow)
  library(optparse)
  library(dplyr)
})
SourceRegressionUtil("regression_dgp.R")
SourceRegressionUtil("t1e_sim.R")

# -----------------------------------------------------------------------------
# Command-line and paths
# -----------------------------------------------------------------------------

opt_list <- list(
  make_option(
    c("--n"),
    type = "character",
    help = "Sample size(s) per arm, comma-separated (e.g. 50,100,200)",
    default = "50,100,200"
  ),
  make_option(
    c("--tau"),
    type = "character",
    help = "Pseudo-value truncation time(s), comma-separated (e.g. 2,3,4)",
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
    help = "Pseudo-value type(s), comma-separated: AUC, MCF (e.g. AUC,MCF)",
    default = "AUC,MCF"
  ),
  make_option(
    c("--reps"),
    type = "integer",
    help = "Number of simulation replicates in this invocation",
    default = 10000L
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
    c("--outdir"),
    type = "character",
    help = "Output directory under Regression (default: simulations/pseudoreg/t1e/null)",
    default = NULL
  ),
  make_option(
    c("--seed"),
    type = "integer",
    help = "Base seed (production 200-per-arm default 2026092101)",
    default = 2026092101L
  ),
  make_option(
    c("--rep_start"),
    type = "integer",
    help = "First global replicate ID",
    default = 1L
  ),
  make_option(
    c("--resume"),
    action = "store_true",
    help = "Skip output files that already exist and validate as complete",
    default = FALSE
  ),
  make_option(
    c("--workers"),
    type = "integer",
    help = "Must be 1; R-native multicore execution is not used",
    default = 1L
  )
)
parsed <- parse_args(OptionParser(option_list = opt_list))
n_vec <- as.integer(trimws(strsplit(parsed$n, ",")[[1L]]))
tau_vec <- as.numeric(trimws(strsplit(parsed$tau, ",")[[1L]]))
type_vec <- trimws(strsplit(parsed$type, ",")[[1L]])
type_vec <- type_vec[tolower(type_vec) %in% c("auc", "mcf")]
if (length(type_vec) == 0L) {
  stop("--type must include at least one of: AUC, MCF")
}
if (any(!is.finite(n_vec) | n_vec < 1L)) {
  stop("--n must be a comma-separated list of positive integers (per arm).", call. = FALSE)
}
tau_gen <- parsed$tau_gen
reps <- parsed$reps
testing <- parsed$testing
batch <- parsed$batch
outdir <- parsed$outdir
base_seed <- parsed$seed
rep_start <- parsed$rep_start
resume <- isTRUE(parsed$resume)
workers <- parsed$workers
if (!is.finite(rep_start) || rep_start < 1L) {
  stop("--rep_start must be a positive integer.", call. = FALSE)
}
if (!identical(as.integer(workers), 1L)) {
  stop("--workers must be 1 for this sequential type I driver.", call. = FALSE)
}
if (testing) {
  reps <- min(reps, 100L)
}
canonical <- CanonicalT1eConfigs()
configs <- FilterCanonicalT1eConfigs(
  configs = canonical,
  n_vec = n_vec,
  tau_vec = tau_vec,
  type_vec = type_vec
)
if (nrow(configs) == 0L) {
  stop("No canonical type I cells remain after applying CLI filters.", call. = FALSE)
}
reg_dir <- RegressionDir()
ValidatePseudoReg()

if (is.null(batch)) {
  batch <- format(Sys.time(), "%Y%m%d_%H%M%S")
}
if (is.null(outdir)) {
  outdir <- file.path(reg_dir, "simulations", "pseudoreg", "t1e", "null")
} else if (!grepl("^/", outdir)) {
  outdir <- file.path(reg_dir, outdir)
}

# -----------------------------------------------------------------------------
# Simulation functions
# -----------------------------------------------------------------------------

#' Single replicate: generate data under H0 and fit PseudoReg(~ arm).
#'
#' @param n Sample size per arm.
#' @param tau Truncation time for pseudo-value construction.
#' @param tau_gen Truncation time for data generation (fixed across configs).
#' @param type Pseudo-value type: "AUC" or "MCF".
#' @param seed Deterministic replicate seed.
#' @param max_attempts Maximum deterministic attempts to obtain an estimable
#'   risk set at \code{tau}.
#' @return Data.frame with one row: pvalue_arm, seed, attempts, arm counts.
SimInst <- function(
  n,
  tau,
  tau_gen = 4,
  type = "AUC",
  seed,
  max_attempts = 100L
) {
  n_per_arm <- as.integer(n)
  if (!is.finite(n_per_arm) || n_per_arm < 1L) {
    stop("n (per arm) must be a positive integer.", call. = FALSE)
  }
  n_total <- as.integer(2L * n_per_arm)
  arm <- c(rep(0L, n_per_arm), rep(1L, n_per_arm))
  if (length(arm) != n_total) {
    stop("Realized total sample size does not equal 2 * n_per_arm.", call. = FALSE)
  }
  n_arm0 <- as.integer(sum(arm == 0L))
  n_arm1 <- as.integer(sum(arm == 1L))
  if (!identical(n_arm0, n_per_arm) || !identical(n_arm1, n_per_arm)) {
    stop("Realized arm counts must each equal n_per_arm.", call. = FALSE)
  }
  covariates <- matrix(c(rep(1, n_total), arm), ncol = 2L)
  colnames(covariates) <- c("intercept", "arm")
  beta_event <- c(0, 0)
  beta_death <- c(0, 0)
  subject <- data.frame(
    idx = seq_len(n_total),
    arm = arm
  )
  fit <- NULL
  attempt_seed <- seed
  attempt <- NA_integer_
  for (attempt in seq_len(max_attempts)) {
    attempt_seed <- if (attempt == 1L) {
      seed
    } else {
      RegressionRepSeed(
        base_seed = seed,
        config_index = attempt - 1L,
        rep_index = 1L
      )
    }
    set.seed(attempt_seed)
    data <- GenData(
      covariates = covariates,
      beta_event = beta_event,
      beta_death = beta_death,
      tau = tau_gen
    )
    fit_or_error <- tryCatch(
      MCC::PseudoReg(
        ~ arm,
        data = data,
        tau = tau,
        type = type,
        covariates = subject
      ),
      error = identity
    )
    if (!inherits(fit_or_error, "error")) {
      fit <- fit_or_error
      break
    }
    expected_error <- grepl(
      paste(
        "within the observed time range",
        "No subject remains at risk at tau",
        sep = "|"
      ),
      conditionMessage(fit_or_error)
    )
    if (!expected_error) {
      stop(fit_or_error)
    }
  }
  if (is.null(fit)) {
    stop(
      "Could not generate an estimable risk set at tau after ",
      max_attempts,
      " deterministic attempts.",
      call. = FALSE
    )
  }
  inference <- ExtractPseudoRegCoefficient(fit = fit, term = "arm")
  out <- data.frame(
    rep = 1L,
    pvalue_arm = inference$pvalue,
    seed = attempt_seed,
    attempts = as.integer(attempt),
    variance_method = inference$variance_method,
    n_per_arm = n_per_arm,
    n_total = n_total,
    n_arm0 = n_arm0,
    n_arm1 = n_arm1,
    stringsAsFactors = FALSE
  )
  return(out)
}

#' Run many replicates and return a data frame with one row per replicate.
#'
#' @param n Sample size per arm.
#' @param tau Pseudo-value truncation time.
#' @param tau_gen Truncation time for data generation.
#' @param type Pseudo-value type: "AUC" or "MCF".
#' @param reps Number of replicates.
#' @param rep_start First global replicate ID.
#' @param config_index One-based canonical configuration index.
#' @param base_seed Base seed.
#' @param workers Must be 1.
#' @param testing If TRUE, cap reps at 100.
#' @return Data.frame with one row per replicate.
Simulation <- function(
  n,
  tau,
  tau_gen,
  type,
  reps,
  rep_start,
  config_index,
  base_seed,
  workers,
  testing = FALSE
) {
  if (testing) {
    reps <- min(reps, 100L)
  }
  message(glue("Sample size per arm n = {n} (total {2L * as.integer(n)})"))
  message(glue("Tau (pseudo) = {tau}, tau_gen = {tau_gen}, type = {type}"))
  message(glue("Replicates = {reps} starting at {rep_start}"))
  if (testing) {
    message("Testing mode (reps capped at 100)")
  }
  pct_marks <- seq(0.1, 1, by = 0.1)
  milestones <- unique(c(round(pct_marks * reps), reps))
  tic <- Sys.time()
  rep_ids <- seq.int(from = rep_start, length.out = reps)
  data_list <- ParallelSimulationLapply(
    seq_along(rep_ids),
    function(i) {
      rep_id <- rep_ids[i]
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
        seed = RegressionRepSeed(
          base_seed = base_seed,
          config_index = config_index,
          rep_index = rep_id
        )
      )
      out$rep <- rep_id
      return(out)
    },
    workers = workers
  )
  out <- do.call(rbind, data_list)
  return(out)
}

# -----------------------------------------------------------------------------
# Run: loop over filtered canonical cells, keeping global config_index
# -----------------------------------------------------------------------------

message(glue("tau_gen (data generation) = {tau_gen} (fixed)"))
message(glue(
  "Configurations: {nrow(configs)} (n per arm = {paste(n_vec, collapse = ', ')}, tau = {paste(tau_vec, collapse = ', ')}, type = {paste(type_vec, collapse = ', ')})"
))
message(glue("Canonical indices: {paste(configs$config_index, collapse = ', ')}"))
message(glue("Replicates per config: {reps} from {rep_start}"))
if (testing) {
  message("Testing mode (reps capped at 100)")
}

expected_rep_ids <- seq.int(from = rep_start, length.out = reps)
rep_end <- rep_start + reps - 1L

for (k in seq_len(nrow(configs))) {
  n_k <- as.integer(configs$n_per_arm[k])
  tau_k <- as.numeric(configs$tau[k])
  type_k <- as.character(configs$type[k])
  config_index_k <- as.integer(configs$config_index[k])
  config_dir <- file.path(
    outdir,
    T1eCellDirName(n_per_arm = n_k, tau = tau_k)
  )
  dir.create(config_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(
    config_dir,
    glue("pvalues_{type_k}_{batch}_rep{rep_start}-{rep_end}.parquet")
  )
  cell_config <- list(
    n_per_arm = n_k,
    tau = tau_k,
    type = type_k,
    tau_gen = tau_gen,
    config_index = config_index_k,
    source = "new_per_arm_simulation"
  )
  if (file.exists(out_file)) {
    if (
      resume &&
        IsCompleteT1eParquet(
          path = out_file,
          expected_rep_ids = expected_rep_ids,
          config = cell_config,
          allow_source = "new_per_arm_simulation"
        )
    ) {
      message("Skipping completed config batch: ", out_file)
      next
    }
    if (resume) {
      message("Malformed or incomplete file; deleting and recomputing: ", out_file)
      DeleteT1eParquet(out_file)
    } else {
      stop("Output already exists; use --resume or a new --batch: ", out_file)
    }
  }
  message(glue(
    "--- Config {k}/{nrow(configs)}: n_per_arm = {n_k}, tau = {tau_k}, type = {type_k}, config_index = {config_index_k} ---"
  ))
  results <- Simulation(
    n = n_k,
    tau = tau_k,
    tau_gen = tau_gen,
    type = type_k,
    reps = reps,
    rep_start = rep_start,
    config_index = config_index_k,
    base_seed = base_seed,
    workers = workers,
    testing = testing
  )
  results$n <- n_k
  results$tau <- tau_k
  results$type <- type_k
  results$tau_gen <- tau_gen
  results$batch <- batch
  results$config_index <- config_index_k
  results$source <- "new_per_arm_simulation"
  results$source_n_legacy <- NA_integer_
  results$source_batch <- batch
  results$source_config_index <- config_index_k
  results <- results[, T1eRequiredColumns(), drop = FALSE]
  AtomicWriteT1eParquet(
    data = results,
    path = out_file,
    expected_rep_ids = expected_rep_ids,
    config = cell_config,
    allow_source = "new_per_arm_simulation"
  )
  message(glue("Wrote {out_file}"))
}

message(glue("Done. Outputs under {outdir}"))
