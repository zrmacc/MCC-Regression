#!/usr/bin/env Rscript
# Supports Figures 1 and 2 and Table S1: large-sample recurrent-event-DGP regression truth.

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
# Paths and parsing
# -----------------------------------------------------------------------------

FormatElapsedHMS <- function(seconds) {
  total_seconds <- as.integer(round(seconds))
  hh <- total_seconds %/% 3600
  mm <- (total_seconds %% 3600) %/% 60
  ss <- total_seconds %% 60
  out <- sprintf("%02d:%02d:%02d", hh, mm, ss)
  return(out)
}

# -----------------------------------------------------------------------------
# Ground truth computation (shared logic with bias sim)
# -----------------------------------------------------------------------------

MCFAtTau <- function(mcf, tau) {
  w <- which(mcf$time <= tau)
  if (length(w) == 0L) {
    return(0)
  }
  return(mcf$mcf[max(w)])
}

#' Compute true marginal treatment effect for one DGP at given taus.
ComputeGroundTruthOne <- function(
    beta_name,
    beta_value,
    tau_gen,
    tau_vec,
    n_per_arm,
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
  beta_event <- MakeBetaVector(beta_name = beta_name, beta_value = beta_value)
  beta_death <- rep(0, 4)

  data <- GenData(
    covariates = covariates,
    beta_event = beta_event,
    beta_death = beta_death,
    censoring_rate = censoring_rate,
    tau = tau_gen
  )

  arm0 <- data %>% dplyr::filter(A == 0)
  arm1 <- data %>% dplyr::filter(A == 1)

  mcf0 <- CalcMCF(
    idx = arm0$idx,
    status = arm0$status,
    time = arm0$time,
    calc_var = FALSE
  )
  mcf1 <- CalcMCF(
    idx = arm1$idx,
    status = arm1$status,
    time = arm1$time,
    calc_var = FALSE
  )

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

      data.frame(
        tau = tau,
        delta_MCF = m1_tau - m0_tau,
        delta_AUMCF = aumcf1 - aumcf0,
        beta_A_MCF = unname(stats::coef(fit_mcf)["A"]),
        beta_A_AUMCF = unname(stats::coef(fit_auc)["A"]),
        beta_X_MCF = unname(stats::coef(fit_mcf)["X"]),
        beta_X_AUMCF = unname(stats::coef(fit_auc)["X"]),
        beta_AX_MCF = unname(stats::coef(fit_mcf)["A:X"]),
        beta_AX_AUMCF = unname(stats::coef(fit_auc)["A:X"])
      )
    }
  )
  dplyr::bind_rows(out_list)
}

BuildGroundTruthTable <- function(
    beta_name_vec,
    beta_value_vec,
    tau_vec,
    tau_gen,
    n_per_arm,
    censoring_rate,
    seed,
    rep_ids
) {
  grid <- expand.grid(
    beta_name = beta_name_vec,
    beta_value = beta_value_vec,
    stringsAsFactors = FALSE
  )
  value_cols <- c(
    "delta_MCF",
    "delta_AUMCF",
    "beta_A_MCF",
    "beta_A_AUMCF",
    "beta_X_MCF",
    "beta_X_AUMCF",
    "beta_AX_MCF",
    "beta_AX_AUMCF"
  )
  out_list <- lapply(
    seq_len(nrow(grid)),
    function(i) {
      rep_blocks <- lapply(
        rep_ids,
        function(r) {
          sim_seed <- seed + (i - 1L) * 1000000L + r
          one <- ComputeGroundTruthOne(
            beta_name = grid$beta_name[i],
            beta_value = grid$beta_value[i],
            tau_gen = tau_gen,
            tau_vec = tau_vec,
            n_per_arm = n_per_arm,
            censoring_rate = censoring_rate,
            seed = sim_seed
          )
          one$gt_rep <- r
          one$beta_name <- grid$beta_name[i]
          one$beta_value <- grid$beta_value[i]
          return(one)
        }
      )
      dplyr::bind_rows(rep_blocks)
    }
  )
  dplyr::bind_rows(out_list) %>%
    dplyr::select(
      gt_rep,
      beta_name,
      beta_value,
      tau,
      dplyr::all_of(value_cols)
    )
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

default_beta_values <- paste(
  signif(log(c(0.75, 1.0, 1.25)), 15),
  collapse = ","
)

opt_list <- list(
  make_option(
    c("--tau"),
    type = "character",
    help = "Truncation time(s), comma-separated (e.g. 2,3,4)",
    default = "2,3,4"
  ),
  make_option(
    c("--tau_gen"),
    type = "numeric",
    help = "Truncation time for data generation (default 4)",
    default = 4
  ),
  make_option(
    c("--beta_name"),
    type = "character",
    help = "Coefficient(s): beta_A, beta_X, beta_AX",
    default = "beta_A,beta_X,beta_AX"
  ),
  make_option(
    c("--beta_value"),
    type = "character",
    help = "Log-scale coefficient value(s), comma-separated",
    default = default_beta_values
  ),
  make_option(
    c("--n_ground"),
    type = "integer",
    help = "Sample size per arm in each ground-truth simulation (default 10000)",
    default = 10000L
  ),
  make_option(
    c("--reps"),
    type = "integer",
    help = paste0(
      "Number of independent ground-truth replicates to save (one row block per ",
      "replicate; default 1000)"
    ),
    default = 1000L
  ),
  make_option(
    c("--gt_reps"),
    type = "integer",
    help = "Deprecated alias for --reps (kept for compatibility).",
    default = NULL
  ),
  make_option(
    c("--gt_rep_start"),
    type = "integer",
    help = "Global replicate index to start from (default 1)",
    default = 1L
  ),
  make_option(
    c("--gt_rep_end"),
    type = "integer",
    help = "Global replicate index to end at (default: reps)",
    default = NA
  ),
  make_option(
    c("--censoring_rate"),
    type = "numeric",
    help = "Censoring rate for ground-truth data (default 0: everyone censored at tau_gen)",
    default = 0
  ),
  make_option(
    c("--seed"),
    type = "integer",
    help = "Random seed for ground-truth data generation",
    default = 0L
  ),
  make_option(
    c("--out"),
    type = "character",
    help = paste0(
      "Output parquet path. If omitted, a new file is created under ",
      "simulations/ground_truth/ with a timestamp suffix to avoid overwrites."
    ),
    default = NULL
  ),
  make_option(
    c("--overwrite"),
    action = "store_true",
    default = FALSE,
    help = paste0(
      "Allow replacing an existing file when --out is set (ignored for default ",
      "timestamped paths)"
    )
  )
)

parsed <- parse_args(OptionParser(option_list = opt_list))
reg_dir <- RegressionDir()

tau_vec <- ParseNumericOption(parsed$tau)
tau_gen <- parsed$tau_gen
beta_name_vec <- ParseBetaName(parsed$beta_name)
beta_value_vec <- ParseNumericOption(parsed$beta_value)
n_per_arm <- parsed$n_ground
reps <- parsed$reps
if (!is.null(parsed$gt_reps)) {
  reps <- parsed$gt_reps
}
gt_rep_start <- parsed$gt_rep_start
gt_rep_end <- parsed$gt_rep_end
censoring_rate <- parsed$censoring_rate
seed <- parsed$seed
overwrite <- isTRUE(parsed$overwrite)

if (!is.finite(gt_rep_start) || gt_rep_start < 1L) {
  stop("--gt_rep_start must be >= 1.")
}
if (is.na(gt_rep_end)) {
  gt_rep_end <- reps
}
if (!is.finite(gt_rep_end) || gt_rep_end < gt_rep_start) {
  stop("--gt_rep_end must be >= --gt_rep_start.")
}
if (gt_rep_end > reps) {
  stop("--gt_rep_end cannot exceed --reps.")
}
rep_ids <- seq.int(gt_rep_start, gt_rep_end)

out_file <- parsed$out
user_specified_out <- !is.null(out_file)
if (!user_specified_out) {
  out_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_file <- file.path(
    reg_dir,
    "simulations",
    "ground_truth",
    glue("ground_truth_n{n_per_arm}_nrep{length(rep_ids)}_r{gt_rep_start}-{gt_rep_end}_{out_stamp}.parquet")
  )
} else {
  out_file <- ResolveRegressionPath(out_file)
  if (file.exists(out_file) && !overwrite) {
    stop(
      "Output file already exists: ", out_file, "\n",
      "Choose a different --out path, or rerun with --overwrite to replace it.",
      call. = FALSE
    )
  }
}

message(glue(
  "Computing ground truth: n_per_arm = {n_per_arm}, reps = {reps}, ",
  "rep_range = {gt_rep_start}:{gt_rep_end} ",
  "(stacked replicates, no averaging), tau = {paste(tau_vec, collapse = ',')}, ",
  "tau_gen = {tau_gen}"
))
tic <- Sys.time()
ground_truth_table <- BuildGroundTruthTable(
  beta_name_vec = beta_name_vec,
  beta_value_vec = beta_value_vec,
  tau_vec = tau_vec,
  tau_gen = tau_gen,
  n_per_arm = n_per_arm,
  censoring_rate = censoring_rate,
  seed = seed,
  rep_ids = rep_ids
)
toc <- Sys.time()
elapsed_sec <- as.numeric(difftime(toc, tic, units = "secs"))
message(glue("Ground truth table: {nrow(ground_truth_table)} rows"))
message(
  glue(
    "Total wall time (full job, all DGP cells and replicates): {FormatElapsedHMS(elapsed_sec)}"
  )
)

dir.create(
  dirname(out_file),
  recursive = TRUE,
  showWarnings = FALSE
)
write_parquet(ground_truth_table, out_file)
message(glue("Wrote {out_file}"))
