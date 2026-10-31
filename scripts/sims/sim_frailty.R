#!/usr/bin/env Rscript
# Supports Figure S2: shared-frailty bias and corrected coverage.

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
ParseType <- function(type_string) {
  out <- trimws(strsplit(type_string, ",")[[1L]])
  out <- toupper(out)
  out <- out[out %in% c("AUC", "MCF")]
  if (length(out) == 0L) {
    stop("--type must include at least one of: AUC, MCF")
  }
  return(out)
}

ParseNumericOption <- function(s) {
  out <- as.numeric(trimws(strsplit(s, ",")[[1L]]))
  if (anyNA(out)) {
    stop("Failed to parse numeric option.")
  }
  return(out)
}

MakeDesignData <- function(n_per_arm) {
  total_n <- 2L * n_per_arm
  A <- c(rep(0L, n_per_arm), rep(1L, n_per_arm))
  X <- stats::rnorm(total_n, mean = 0, sd = 1)
  AX <- A * X
  out <- data.frame(idx = seq_len(total_n), A = A, X = X, AX = AX)
  return(out)
}

SimInst <- function(
  n,
  tau,
  tau_gen,
  type,
  beta_a,
  frailty_variance,
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
  beta_event <- c(intercept = 0, A = beta_a, X = 0, AX = 0)
  beta_death <- c(intercept = 0, A = 0, X = 0, AX = 0)

  data <- GenData(
    covariates = covariates,
    beta_event = beta_event,
    beta_death = beta_death,
    frailty_variance = frailty_variance,
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
  metadata <- data.frame(
    rep = 1L,
    beta_name = "beta_A",
    beta_value = beta_a,
    frailty_variance = frailty_variance,
    seed = seed
  )
  out <- cbind(metadata, inference)
  return(out)
}

Simulation <- function(
  n,
  tau,
  tau_gen,
  type,
  beta_a,
  frailty_variance,
  reps,
  config_index,
  base_seed,
  workers,
  testing = FALSE
) {
  if (testing) {
    reps <- min(reps, 100L)
  }
  message(glue(
    "n = {n}, tau = {tau}, type = {type}, frailty_variance = {frailty_variance}, reps = {reps}"
  ))
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
        beta_a = beta_a,
        frailty_variance = frailty_variance,
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

opt_list <- list(
  make_option(c("--n"), type = "character", default = "50,100,200"),
  make_option(c("--tau"), type = "character", default = "2,3,4"),
  make_option(c("--tau_gen"), type = "numeric", default = 4),
  make_option(c("--type"), type = "character", default = "AUC,MCF"),
  make_option(c("--frailty_variance"), type = "character", default = "0,1,2,3"),
  make_option(c("--beta_a"), type = "numeric", default = log(0.75)),
  make_option(c("--reps"), type = "integer", default = 2000L),
  make_option(c("--testing", "-t"), action = "store_true", default = FALSE),
  make_option(c("--batch"), type = "character", default = NULL),
  make_option(c("--outdir"), type = "character", default = NULL),
  make_option(
    c("--seed"),
    type = "integer",
    default = 20260904L
  ),
  make_option(
    c("--resume"),
    action = "store_true",
    default = FALSE
  ),
  make_option(
    c("--workers"),
    type = "integer",
    default = 1L
  )
)

parsed <- parse_args(OptionParser(option_list = opt_list))
n_vec <- as.integer(ParseNumericOption(parsed$n))
tau_vec <- ParseNumericOption(parsed$tau)
type_vec <- ParseType(parsed$type)
frailty_variance_vec <- ParseNumericOption(parsed$frailty_variance)
tau_gen <- parsed$tau_gen
beta_a <- parsed$beta_a
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
  frailty_variance = frailty_variance_vec,
  stringsAsFactors = FALSE
)
configs <- AttachCanonicalConfigIndex(configs, "frailty")

reg_dir <- RegressionDir()
ValidatePseudoReg()
if (is.null(batch)) {
  batch <- format(Sys.time(), "%Y%m%d_%H%M%S")
}
if (is.null(outdir)) {
  out_dir <- file.path(reg_dir, "simulations", "pseudoreg", "frailty")
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
  cfg <- configs[k, ]
  expected_config <- configs[k, , drop = FALSE]
  expected_config$tau_gen <- tau_gen
  expected_config$beta_value <- beta_a
  config_slug <- paste(
    sprintf("config_%03d", cfg$config_index),
    paste0("n_", cfg$n),
    paste0("tau_", NumericSlug(cfg$tau)),
    cfg$type,
    paste0("frailty_", NumericSlug(cfg$frailty_variance)),
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
      "--- Config {k}/{nrow(configs)}: n = {cfg$n}, tau = {cfg$tau}, ",
      "type = {cfg$type}, frailty_variance = {cfg$frailty_variance} ---"
    )
  )
  results <- Simulation(
    n = cfg$n,
    tau = cfg$tau,
    tau_gen = tau_gen,
    type = cfg$type,
    beta_a = beta_a,
    frailty_variance = cfg$frailty_variance,
    reps = reps,
    config_index = cfg$config_index,
    base_seed = base_seed,
    workers = workers,
    testing = testing
  )
  results$n <- cfg$n
  results$tau <- cfg$tau
  results$type <- cfg$type
  results$tau_gen <- tau_gen
  results$batch <- batch
  results$config_index <- cfg$config_index
  ValidateResumedInference(results, expected_config, reps, base_seed)
  write_parquet(results, out_file)
  message(glue("Wrote {out_file}"))
}
message(glue("Done. Output directory: {out_dir}"))
