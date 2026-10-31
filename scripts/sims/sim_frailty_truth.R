#!/usr/bin/env Rscript
# Supports Figure S2: large-sample shared-frailty regression truth.

argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(dirname(normalizePath(script))), "utils", "bootstrap_regression.R"))
SetRegressionWd()

suppressPackageStartupMessages({
  library(MCC)
  library(arrow)
  library(dplyr)
  library(optparse)
})
SourceRegressionUtil("regression_dgp.R")
ParseNumericOption <- function(s) {
  out <- as.numeric(trimws(strsplit(s, ",")[[1L]]))
  if (anyNA(out)) {
    stop("Failed to parse a numeric option value.")
  }
  return(out)
}

FormatElapsedHMS <- function(seconds) {
  total_seconds <- as.integer(round(seconds))
  hh <- total_seconds %/% 3600
  mm <- (total_seconds %% 3600) %/% 60
  ss <- total_seconds %% 60
  out <- sprintf("%02d:%02d:%02d", hh, mm, ss)
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

MCFAtTau <- function(mcf, tau) {
  w <- which(mcf$time <= tau)
  if (length(w) == 0L) {
    return(0)
  }
  return(mcf$mcf[max(w)])
}

ComputeFrailtyGroundTruthOne <- function(
    tau_vec,
    tau_gen,
    n_per_arm,
    beta_a,
    frailty_variance,
    censoring_rate,
    seed
) {
  set.seed(seed)
  design_df <- MakeDesignData(n_per_arm = n_per_arm)
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
    censoring_rate = censoring_rate,
    tau = tau_gen
  )

  arm0 <- data %>% dplyr::filter(A == 0)
  arm1 <- data %>% dplyr::filter(A == 1)
  mcf0 <- CalcMCF(idx = arm0$idx, status = arm0$status, time = arm0$time, calc_var = FALSE)
  mcf1 <- CalcMCF(idx = arm1$idx, status = arm1$status, time = arm1$time, calc_var = FALSE)

  out_list <- lapply(
    tau_vec,
    function(tau) {
      m0_tau <- MCFAtTau(mcf0, tau)
      m1_tau <- MCFAtTau(mcf1, tau)
      aumcf0 <- AUC(times = mcf0$time, values = mcf0$mcf, tau = tau)
      aumcf1 <- AUC(times = mcf1$time, values = mcf1$mcf, tau = tau)

      pseudo_mcf <- GenPseudo(data = data, tau = tau, type = "MCF")
      pseudo_auc <- GenPseudo(data = data, tau = tau, type = "AUC")
      reg_mcf <- pseudo_mcf %>%
        dplyr::select(idx, pseudo) %>%
        dplyr::inner_join(design_df, by = "idx")
      reg_auc <- pseudo_auc %>%
        dplyr::select(idx, pseudo) %>%
        dplyr::inner_join(design_df, by = "idx")
      fit_mcf <- stats::lm(pseudo ~ A + X + A:X, data = reg_mcf)
      fit_auc <- stats::lm(pseudo ~ A + X + A:X, data = reg_auc)

      out <- data.frame(
        tau = tau,
        frailty_variance = frailty_variance,
        beta_name = "beta_A",
        beta_value = beta_a,
        delta_MCF = m1_tau - m0_tau,
        delta_AUMCF = aumcf1 - aumcf0,
        beta_A_MCF = unname(stats::coef(fit_mcf)["A"]),
        beta_A_AUMCF = unname(stats::coef(fit_auc)["A"])
      )
      return(out)
    }
  )
  out <- dplyr::bind_rows(out_list)
  return(out)
}

BuildFrailtyGroundTruthTable <- function(
  tau_vec,
  tau_gen,
  frailty_variance_vec,
  n_per_arm,
  beta_a,
  censoring_rate,
  seed,
  rep_ids,
  job_start_time,
  checkpoint_dir,
  resume,
  workers
) {
  total_reps <- length(rep_ids)
  milestone_idx <- if (total_reps <= 10L) {
    seq_len(total_reps)
  } else {
    sort(unique(pmin(
      total_reps,
      pmax(1L, as.integer(round(seq(0.1, 1, length.out = 10L) * total_reps)))
    )))
  }
  out_list <- lapply(
    frailty_variance_vec,
    function(fv) {
      rep_blocks <- ParallelSimulationLapply(
        seq_len(total_reps),
        function(j) {
          r <- rep_ids[j]
          frailty_tag <- gsub(".", "p", format(fv, trim = TRUE), fixed = TRUE)
          part_file <- file.path(
            checkpoint_dir,
            sprintf("frailty_%s_rep_%03d.parquet", frailty_tag, r)
          )
          if (resume && file.exists(part_file)) {
            one <- arrow::read_parquet(part_file)
            if (nrow(one) != length(tau_vec)) {
              stop("Invalid frailty truth checkpoint: ", part_file, call. = FALSE)
            }
            return(one)
          }
          if (j %in% milestone_idx) {
            elapsed <- as.numeric(difftime(
              Sys.time(),
              job_start_time,
              units = "secs"
            ))
            message(sprintf(
              "Start replicate %d/%d (elapsed %s)",
              j,
              total_reps,
              FormatElapsedHMS(elapsed)
            ))
          }
          sim_seed <- seed + as.integer(round(fv * 1e5)) + r
          one <- ComputeFrailtyGroundTruthOne(
            tau_vec = tau_vec,
            tau_gen = tau_gen,
            n_per_arm = n_per_arm,
            beta_a = beta_a,
            frailty_variance = fv,
            censoring_rate = censoring_rate,
            seed = sim_seed
          )
          one$gt_rep <- r
          dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
          arrow::write_parquet(one, part_file)
          message("Checkpointed ", part_file)
          return(one)
        },
        workers = workers
      )
      out <- dplyr::bind_rows(rep_blocks)
      return(out)
    }
  )
  out <- dplyr::bind_rows(out_list) %>%
    dplyr::select(
      gt_rep,
      frailty_variance,
      beta_name,
      beta_value,
      tau,
      delta_MCF,
      delta_AUMCF,
      beta_A_MCF,
      beta_A_AUMCF
    )
  return(out)
}

opt_list <- list(
  make_option(c("--tau"), type = "character", default = "2,3,4"),
  make_option(c("--tau_gen"), type = "numeric", default = 4),
  make_option(c("--frailty_variance"), type = "character", default = "0,1,2,3"),
  make_option(c("--beta_a"), type = "numeric", default = log(0.75)),
  make_option(c("--n_ground"), type = "integer", default = 10000L),
  make_option(c("--reps"), type = "integer", default = 1000L),
  make_option(
    c("--gt_reps"),
    type = "integer",
    help = "Deprecated alias for --reps (kept for compatibility).",
    default = NULL
  ),
  make_option(c("--censoring_rate"), type = "numeric", default = 0),
  make_option(c("--seed"), type = "integer", default = 0L),
  make_option(c("--out"), type = "character", default = NULL),
  make_option(c("--overwrite"), action = "store_true", default = FALSE),
  make_option(c("--resume"), action = "store_true", default = FALSE),
  make_option(c("--workers"), type = "integer", default = 1L)
)

parsed <- parse_args(OptionParser(option_list = opt_list))
reg_dir <- RegressionDir()

tau_vec <- ParseNumericOption(parsed$tau)
tau_gen <- parsed$tau_gen
frailty_variance_vec <- ParseNumericOption(parsed$frailty_variance)
beta_a <- parsed$beta_a
n_per_arm <- parsed$n_ground
reps <- parsed$reps
if (!is.null(parsed$gt_reps)) {
  reps <- parsed$gt_reps
}
censoring_rate <- parsed$censoring_rate
seed <- parsed$seed
overwrite <- isTRUE(parsed$overwrite)
resume <- isTRUE(parsed$resume)
workers <- parsed$workers

if (!is.finite(reps) || reps < 1L) {
  stop("--reps must be >= 1.")
}
rep_ids <- seq_len(reps)

out_file <- parsed$out
if (is.null(out_file)) {
  out_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_file <- file.path(
    reg_dir,
    "simulations",
    "frailty_ground_truth",
    sprintf(
      "frailty_ground_truth_n%d_nrep%d_%s.parquet",
      n_per_arm,
      reps,
      out_stamp
    )
  )
} else {
  out_file <- ResolveRegressionPath(out_file)
  if (file.exists(out_file) && !overwrite && !resume) {
    stop("Output file already exists; use --overwrite.")
  }
}
checkpoint_dir <- paste0(out_file, ".parts")

message(paste0(
  "Computing frailty ground truth: n_per_arm = ",
  n_per_arm,
  ", reps = ",
  reps,
  ", frailty_variance = ",
  paste(frailty_variance_vec, collapse = ","),
  ", tau = ",
  paste(tau_vec, collapse = ","),
  ", tau_gen = ",
  tau_gen
))
tic <- Sys.time()
gt <- BuildFrailtyGroundTruthTable(
  tau_vec = tau_vec,
  tau_gen = tau_gen,
  frailty_variance_vec = frailty_variance_vec,
  n_per_arm = n_per_arm,
  beta_a = beta_a,
  censoring_rate = censoring_rate,
  seed = seed,
  rep_ids = rep_ids,
  job_start_time = tic,
  checkpoint_dir = checkpoint_dir,
  resume = resume,
  workers = workers
)
toc <- Sys.time()
elapsed_sec <- as.numeric(difftime(toc, tic, units = "secs"))
message(sprintf("Frailty ground truth table: %d rows", nrow(gt)))
message(sprintf(
  "Total wall time (full job, all frailty cells and replicates): %s",
  FormatElapsedHMS(elapsed_sec)
))
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write_parquet(gt, out_file)
message(sprintf("Wrote %s", out_file))
