#!/usr/bin/env Rscript
# Supports Figures S3 and S4: large-sample terminal-event-DGP regression truth.

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

# -----------------------------------------------------------------------------
# Paths and parsing
ParseBetaName <- function(s) {
  out <- trimws(strsplit(s, ",")[[1L]])
  out <- out[out %in% c("beta_A", "beta_X", "beta_AX")]
  if (length(out) == 0L) {
    stop("--beta_name must include at least one of: beta_A, beta_X, beta_AX")
  }
  return(out)
}

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

# -----------------------------------------------------------------------------
# Ground truth computation
# -----------------------------------------------------------------------------

MakeBetaVector <- function(beta_name, beta_value) {
  out <- c(intercept = 0, A = 0, X = 0, AX = 0)
  if (beta_name == "beta_A") {
    out["A"] <- beta_value
  } else if (beta_name == "beta_X") {
    out["X"] <- beta_value
  } else if (beta_name == "beta_AX") {
    out["AX"] <- beta_value
  } else {
    stop("Unknown beta_name: ", beta_name)
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

MCFAtTau <- function(mcf, tau) {
  w <- which(mcf$time <= tau)
  if (length(w) == 0L) {
    return(0)
  }
  return(mcf$mcf[max(w)])
}

#' Compute ground truth for one DGP at given taus.
#'
#' Varies the death-rate model and keeps the recurrent-event model at baseline.
ComputeTerminalEventGroundTruthOne <- function(
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

  beta_event <- c(intercept = 0, A = 0, X = 0, AX = 0)
  beta_death <- MakeBetaVector(beta_name = beta_name, beta_value = beta_value)

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

      out <- data.frame(
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
      return(out)
    }
  )
  out <- dplyr::bind_rows(out_list)
  return(out)
}

BuildTerminalEventGroundTruthTable <- function(
    beta_name_vec,
    beta_value_vec,
    tau_vec,
    tau_gen,
    n_per_arm,
    censoring_rate,
    seed,
    rep_ids,
    job_start_time
) {
  grid <- expand.grid(
    beta_name = beta_name_vec,
    beta_value = beta_value_vec,
    stringsAsFactors = FALSE
  )
  total_reps <- length(rep_ids)
  milestone_idx <- if (total_reps <= 10L) {
    seq_len(total_reps)
  } else {
    sort(unique(pmin(
      total_reps,
      pmax(1L, as.integer(round(seq(0.1, 1, length.out = 10L) * total_reps)))
    )))
  }
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
      message(glue(
        "--- DGP cell {i}/{nrow(grid)}: beta_death = {grid$beta_name[i]} ",
        "at {signif(grid$beta_value[i], 6)}, taus = {paste(tau_vec, collapse = ',')}, ",
        "tau_gen = {tau_gen} ---"
      ))
      rep_blocks <- vector(mode = "list", length = total_reps)
      for (j in seq_len(total_reps)) {
        if (j %in% milestone_idx) {
          elapsed <- as.numeric(difftime(Sys.time(), job_start_time, units = "secs"))
          message(glue(
            "Start replicate {j}/{total_reps} for {grid$beta_name[i]}={signif(grid$beta_value[i], 6)} ",
            "(taus {paste(tau_vec, collapse = ',')}; elapsed {FormatElapsedHMS(elapsed)})"
          ))
        }
        r <- rep_ids[j]
          sim_seed <- seed + (i - 1L) * 1000000L + r
          one <- ComputeTerminalEventGroundTruthOne(
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
          rep_blocks[[j]] <- one
      }
      dplyr::bind_rows(rep_blocks)
    }
  )
  out <- dplyr::bind_rows(out_list) %>%
    dplyr::select(
      gt_rep,
      beta_name,
      beta_value,
      tau,
      dplyr::all_of(value_cols)
    )
  return(out)
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
    help = "Coefficient(s) varied in beta_death: beta_A, beta_X, beta_AX",
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
      "simulations/terminal_event_ground_truth/ with a timestamp suffix to avoid overwrites."
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
censoring_rate <- parsed$censoring_rate
seed <- parsed$seed
overwrite <- isTRUE(parsed$overwrite)

if (!is.finite(reps) || reps < 1L) {
  stop("--reps must be >= 1.")
}
rep_ids <- seq_len(reps)

out_file <- parsed$out
user_specified_out <- !is.null(out_file)
if (!user_specified_out) {
  out_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_file <- file.path(
    reg_dir,
    "simulations",
    "terminal_event_ground_truth",
    glue("terminal_event_ground_truth_n{n_per_arm}_nrep{length(rep_ids)}_{out_stamp}.parquet")
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
  "Computing terminal-event ground truth: n_per_arm = {n_per_arm}, reps = {reps}, ",
  "tau = {paste(tau_vec, collapse = ',')}, tau_gen = {tau_gen}"
))
tic <- Sys.time()
ground_truth_table <- BuildTerminalEventGroundTruthTable(
  beta_name_vec = beta_name_vec,
  beta_value_vec = beta_value_vec,
  tau_vec = tau_vec,
  tau_gen = tau_gen,
  n_per_arm = n_per_arm,
  censoring_rate = censoring_rate,
  seed = seed,
  rep_ids = rep_ids,
  job_start_time = tic
)
toc <- Sys.time()
elapsed_sec <- as.numeric(difftime(toc, tic, units = "secs"))
message(glue("Terminal-event ground truth table: {nrow(ground_truth_table)} rows"))
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

