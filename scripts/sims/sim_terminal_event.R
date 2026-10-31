#!/usr/bin/env Rscript
# Supports Figures S3 and S4: terminal-event covariate effects, corrected bias/coverage.

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
SourceRegressionUtil("terminal_event_sim.R")

# -----------------------------------------------------------------------------
# Simulation helpers
# -----------------------------------------------------------------------------

#' Run a single terminal-event simulation replicate.
#'
#' Holds recurrent-event coefficients at zero and varies one death-rate
#' coefficient. Fits \code{~ A + X + A:X} with \code{MCC::PseudoReg()} and
#' records only the matching term.
#'
#' @param n Sample size per arm.
#' @param tau Pseudo-value truncation time.
#' @param tau_gen Data-generation truncation time.
#' @param type Pseudo-value type (AUC or MCF).
#' @param beta_name Name of the death-rate coefficient varied in the DGP.
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
  beta_event <- c(intercept = 0, A = 0, X = 0, AX = 0)
  beta_death <- MakeBetaVector(beta_name = beta_name, beta_value = beta_value)
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
    beta_name = beta_name,
    beta_value = beta_value,
    stringsAsFactors = FALSE
  )
  out <- cbind(metadata, inference)
  return(out)
}

#' Attempt one replicate up to \code{max_attempts} times with derived seeds.
#'
#' @param n Sample size per arm.
#' @param tau Pseudo-value truncation time.
#' @param tau_gen Data-generation truncation time.
#' @param type Pseudo-value type.
#' @param beta_name Varied death-rate coefficient name.
#' @param beta_value Varied death-rate coefficient value.
#' @param base_seed Base seed for \code{RegressionRepSeed()}.
#' @param config_index Global configuration index.
#' @param rep_index One-based replicate ID to store.
#' @param max_attempts Maximum attempts including the original seed.
#' @return A one-row data frame, or \code{NULL} if every attempt failed.
SimInstWithRetry <- function(
  n,
  tau,
  tau_gen,
  type,
  beta_name,
  beta_value,
  base_seed,
  config_index,
  rep_index,
  max_attempts
) {
  last_error <- NULL
  last_seed <- NA_integer_
  for (attempt in seq_len(max_attempts)) {
    seed <- RegressionRepSeed(
      base_seed = base_seed,
      config_index = config_index,
      rep_index = rep_index + (attempt - 1L) * 1000003L
    )
    last_seed <- seed
    result <- tryCatch(
      {
        row <- SimInst(
          n = n,
          tau = tau,
          tau_gen = tau_gen,
          type = type,
          beta_name = beta_name,
          beta_value = beta_value,
          seed = seed
        )
        row$rep <- as.integer(rep_index)
        row$seed <- as.integer(seed)
        row$attempt <- as.integer(attempt)
        row
      },
      error = function(e) e
    )
    if (!inherits(result, "error")) {
      return(result)
    }
    if (IsSystemicPseudoRegError(result)) {
      stop(conditionMessage(result), call. = FALSE)
    }
    last_error <- result
  }
  attr(last_error, "seed") <- last_seed
  return(last_error)
}

#' Attach configuration metadata shared by every output row.
#'
#' @param data Replicate-level data frame.
#' @param config Named list of cell metadata.
#' @param shard Shard index, or 1 when the cell is unsharded.
#' @return Data frame with metadata columns.
AttachTerminalEventMetadata <- function(data, config, shard) {
  data$config_index <- as.integer(config$config_index)
  data$block_id <- as.integer(config$block_id)
  data$shard <- as.integer(shard)
  data$n <- as.integer(config$n)
  data$tau <- as.numeric(config$tau)
  data$tau_gen <- as.numeric(config$tau_gen)
  data$type <- as.character(config$type)
  data$beta_name <- as.character(config$beta_name)
  data$beta_value <- as.numeric(config$beta_value)
  data$batch <- as.character(config$batch)
  out <- data[, TerminalEventRequiredColumns(), drop = FALSE]
  return(out)
}

#' Run assigned replicates for one shard, retrying failed fits.
#'
#' @param config Named list of cell metadata.
#' @param rep_ids Integer vector of replicate IDs.
#' @param shard Shard index.
#' @param base_seed Base seed.
#' @param max_attempts Maximum attempts per replicate.
#' @param failure_jsonl JSONL path for exhausted replicates.
#' @return List with \code{data} (possibly NULL) and \code{failed} records.
RunShardReplicates <- function(
  config,
  rep_ids,
  shard,
  base_seed,
  max_attempts,
  failure_jsonl
) {
  rows <- vector("list", length(rep_ids))
  failed <- list()
  n_rep <- length(rep_ids)
  pct_marks <- seq(0.1, 1, by = 0.1)
  milestones <- unique(c(round(pct_marks * n_rep), n_rep))
  tic <- Sys.time()
  for (i in seq_along(rep_ids)) {
    if (i %in% milestones) {
      pct <- round(100 * i / n_rep, 1)
      toc <- Sys.time()
      delta <- round(as.numeric(difftime(toc, tic, units = "secs")), 2)
      message(
        glue(
          "Config {config$config_index} shard {shard}: ",
          "instance {i}/{n_rep} ({pct}%). Elapsed: {delta}s"
        )
      )
      tic <- toc
    }
    result <- SimInstWithRetry(
      n = config$n,
      tau = config$tau,
      tau_gen = config$tau_gen,
      type = config$type,
      beta_name = config$beta_name,
      beta_value = config$beta_value,
      base_seed = base_seed,
      config_index = config$config_index,
      rep_index = as.integer(rep_ids[i]),
      max_attempts = max_attempts
    )
    if (inherits(result, "error")) {
      seed_used <- attr(result, "seed")
      if (is.null(seed_used)) {
        seed_used <- NA_integer_
      }
      record <- list(
        timestamp = format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ"),
        block_id = as.integer(config$block_id),
        config_index = as.integer(config$config_index),
        shard = as.integer(shard),
        rep = as.integer(rep_ids[i]),
        attempts = as.integer(max_attempts),
        seed = as.integer(seed_used),
        n = as.integer(config$n),
        tau = as.numeric(config$tau),
        type = as.character(config$type),
        beta_name = as.character(config$beta_name),
        beta_value = as.numeric(config$beta_value),
        error = conditionMessage(result)
      )
      AppendJsonlRecord(path = failure_jsonl, record = record)
      failed[[length(failed) + 1L]] <- record
      message(
        glue(
          "Config {config$config_index} rep {rep_ids[i]} failed after ",
          "{max_attempts} attempts: {conditionMessage(result)}"
        )
      )
      next
    }
    rows[[i]] <- result
  }
  ok <- !vapply(rows, is.null, logical(1L))
  data <- NULL
  if (any(ok)) {
    data <- do.call(rbind, rows[ok])
    data <- AttachTerminalEventMetadata(
      data = data,
      config = config,
      shard = shard
    )
  }
  out <- list(data = data, failed = failed)
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
    c("--reps"),
    type = "integer",
    help = "Number of simulation replicates per cell",
    default = 2000L
  ),
  make_option(
    c("--testing", "-t"),
    action = "store_true",
    default = FALSE,
    help = "Testing mode (cap reps at 100)"
  ),
  make_option(
    c("--batch"),
    type = "character",
    help = "Batch label for output filenames (default: timestamp)",
    default = NULL
  ),
  make_option(
    c("--outdir"),
    type = "character",
    help = "Output directory relative to Regression (default: simulations/pseudoreg/terminal_event)",
    default = NULL
  ),
  make_option(
    c("--seed"),
    type = "integer",
    help = "Base seed for deterministic per-replicate seeds",
    default = 2026092001L
  ),
  make_option(
    c("--resume"),
    action = "store_true",
    default = FALSE,
    help = "Skip valid shards/cells; delete and recompute malformed files"
  ),
  make_option(
    c("--block_id"),
    type = "integer",
    help = "Required block identifier (1 through n_blocks)",
    default = NULL
  ),
  make_option(
    c("--n_blocks"),
    type = "integer",
    help = "Number of independent execution blocks",
    default = 6L
  ),
  make_option(
    c("--shard_size"),
    type = "integer",
    help = "Replicates per shard (cell file is written directly if shard_size >= reps)",
    default = 400L
  ),
  make_option(
    c("--max_attempts"),
    type = "integer",
    help = "Maximum attempts per replicate, including the original seed",
    default = 5L
  )
)

parsed <- parse_args(OptionParser(option_list = opt_list))
n_vec <- as.integer(ParseNumericOption(parsed$n))
tau_vec <- ParseNumericOption(parsed$tau)
type_vec <- ParseType(parsed$type)
beta_name_vec <- ParseBetaName(parsed$beta_name)
beta_value_vec <- ParseNumericOption(parsed$beta_value)
tau_gen <- parsed$tau_gen
reps <- as.integer(parsed$reps)
testing <- isTRUE(parsed$testing)
batch <- parsed$batch
outdir <- parsed$outdir
base_seed <- as.integer(parsed$seed)
resume <- isTRUE(parsed$resume)
block_id <- parsed$block_id
n_blocks <- as.integer(parsed$n_blocks)
shard_size <- as.integer(parsed$shard_size)
max_attempts <- as.integer(parsed$max_attempts)

if (is.null(block_id)) {
  stop("--block_id is required (integer from 1 through --n_blocks).", call. = FALSE)
}
block_id <- as.integer(block_id)
if (!is.finite(n_blocks) || n_blocks < 1L) {
  stop("--n_blocks must be >= 1.", call. = FALSE)
}
if (!is.finite(block_id) || block_id < 1L || block_id > n_blocks) {
  stop("--block_id must satisfy 1 <= block_id <= n_blocks.", call. = FALSE)
}
if (!is.finite(reps) || reps < 1L) {
  stop("--reps must be a positive integer.", call. = FALSE)
}
if (!is.finite(shard_size) || shard_size < 1L) {
  stop("--shard_size must be a positive integer.", call. = FALSE)
}
if (!is.finite(max_attempts) || max_attempts < 1L) {
  stop("--max_attempts must be a positive integer.", call. = FALSE)
}

if (testing) {
  reps <- min(reps, 100L)
}

canonical <- CanonicalTerminalEventConfigs()
if (nrow(canonical) != 162L) {
  stop("Canonical terminal-event grid must contain 162 cells.", call. = FALSE)
}
filtered <- FilterCanonicalTerminalEventConfigs(
  configs = canonical,
  n_vec = n_vec,
  tau_vec = tau_vec,
  type_vec = type_vec,
  beta_name_vec = beta_name_vec,
  beta_value_vec = beta_value_vec
)
filtered$block_id <- TerminalEventBlockId(
  config_index = filtered$config_index,
  n_blocks = n_blocks
)
block_configs <- filtered[filtered$block_id == block_id, , drop = FALSE]

reg_dir <- RegressionDir()
ValidatePseudoReg()

if (is.null(batch)) {
  batch <- format(Sys.time(), "%Y%m%d_%H%M%S")
}
if (is.null(outdir)) {
  out_dir <- file.path(reg_dir, "simulations", "pseudoreg", "terminal_event")
} else {
  out_dir <- ResolveRegressionPath(outdir)
}
shard_dir <- file.path(out_dir, "shards")
fail_dir <- file.path(out_dir, "failures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(shard_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fail_dir, recursive = TRUE, showWarnings = FALSE)
failure_jsonl <- file.path(fail_dir, sprintf("block_%02d.jsonl", block_id))

use_direct_cell <- shard_size >= reps
shard_plan <- TerminalEventShardPlan(reps = reps, shard_size = shard_size)

message(glue("tau_gen = {tau_gen} (fixed)"))
message(glue("Block {block_id}/{n_blocks}"))
message(glue("Canonical cells: {nrow(canonical)}"))
message(glue("Cells after design filters: {nrow(filtered)}"))
message(glue("Cells assigned to this block: {nrow(block_configs)}"))
message(glue("Replicates per cell: {reps}"))
message(glue("Shard size: {shard_size}; direct cell write: {use_direct_cell}"))
message(glue("Max attempts per replicate: {max_attempts}"))
if (testing) {
  message("Testing mode (reps capped at 100)")
}

n_failed_cells <- 0L
n_complete_cells <- 0L
n_skipped_cells <- 0L

for (row_id in seq_len(nrow(block_configs))) {
  cell <- block_configs[row_id, , drop = FALSE]
  config <- list(
    config_index = as.integer(cell$config_index),
    block_id = as.integer(cell$block_id),
    n = as.integer(cell$n),
    tau = as.numeric(cell$tau),
    type = as.character(cell$type),
    beta_name = as.character(cell$beta_name),
    beta_value = as.numeric(cell$beta_value),
    tau_gen = tau_gen,
    batch = batch
  )
  slug <- TerminalEventConfigSlug(
    config_index = config$config_index,
    n = config$n,
    tau = config$tau,
    type = config$type,
    beta_name = config$beta_name,
    beta_value = config$beta_value,
    batch = config$batch
  )
  cell_file <- file.path(out_dir, paste0(slug, ".parquet"))
  failed_marker <- file.path(
    fail_dir,
    sprintf("config_%03d.failed.json", config$config_index)
  )
  expected_cell_ids <- seq_len(reps)

  message(
    glue(
      "--- Config {config$config_index} (block {config$block_id}): ",
      "n = {config$n}, tau = {config$tau}, type = {config$type}, ",
      "beta_name = {config$beta_name}, ",
      "beta_value = {signif(config$beta_value, 6)} ---"
    )
  )

  if (TryValidateTerminalEventParquet(
    path = cell_file,
    expected_rep_ids = expected_cell_ids,
    config = config
  )) {
    if (resume) {
      message("Skipping completed cell: ", cell_file)
      leftover <- list.files(
        shard_dir,
        pattern = paste0("^", slug, "_shard_"),
        full.names = TRUE
      )
      if (length(leftover) > 0L) {
        unlink(leftover)
      }
      n_skipped_cells <- n_skipped_cells + 1L
      n_complete_cells <- n_complete_cells + 1L
      next
    }
    stop("Output already exists; use --resume or a new --batch: ", cell_file)
  }
  if (file.exists(cell_file)) {
    if (!resume) {
      stop("Output already exists; use --resume or a new --batch: ", cell_file)
    }
    message("Deleting malformed cell file: ", cell_file)
    UnlinkTerminalEventPath(cell_file)
  }
  if (CellFailedMarkerExists(failed_marker)) {
    if (resume) {
      message("Skipping finished-failed cell: ", failed_marker)
      n_skipped_cells <- n_skipped_cells + 1L
      n_failed_cells <- n_failed_cells + 1L
      next
    }
    stop("Finished-failed marker exists; use --resume or a new --batch: ", failed_marker)
  }

  if (use_direct_cell) {
    shard_result <- RunShardReplicates(
      config = config,
      rep_ids = expected_cell_ids,
      shard = 1L,
      base_seed = base_seed,
      max_attempts = max_attempts,
      failure_jsonl = failure_jsonl
    )
    if (length(shard_result$failed) > 0L) {
      WriteFailedCellMarker(
        path = failed_marker,
        record = list(
          timestamp = format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ"),
          config_index = config$config_index,
          block_id = config$block_id,
          n_failed_reps = length(shard_result$failed),
          n_success_reps = if (is.null(shard_result$data)) {
            0L
          } else {
            nrow(shard_result$data)
          },
          error = shard_result$failed[[1L]]$error
        )
      )
      n_failed_cells <- n_failed_cells + 1L
      next
    }
    cell_data <- shard_result$data
    cell_data$shard <- 1L
    AtomicWriteTerminalEventParquet(
      data = cell_data,
      path = cell_file,
      expected_rep_ids = expected_cell_ids,
      config = config
    )
    message("Wrote ", cell_file)
    n_complete_cells <- n_complete_cells + 1L
    next
  }

  shard_files <- character(nrow(shard_plan))
  cell_failed <- FALSE
  for (s in seq_len(nrow(shard_plan))) {
    plan_row <- shard_plan[s, , drop = FALSE]
    shard_file <- file.path(
      shard_dir,
      sprintf(
        "%s_shard_%02d_of_%02d.parquet",
        slug,
        plan_row$shard,
        plan_row$n_shards
      )
    )
    shard_files[s] <- shard_file
    shard_ids <- seq.int(plan_row$rep_start, plan_row$rep_end)
    if (TryValidateTerminalEventParquet(
      path = shard_file,
      expected_rep_ids = shard_ids,
      config = config
    )) {
      if (resume) {
        message("Skipping completed shard: ", shard_file)
        next
      }
      stop("Shard already exists; use --resume or a new --batch: ", shard_file)
    }
    if (file.exists(shard_file)) {
      if (!resume) {
        stop("Shard already exists; use --resume or a new --batch: ", shard_file)
      }
      message("Deleting malformed shard: ", shard_file)
      UnlinkTerminalEventPath(shard_file)
    }
    shard_result <- RunShardReplicates(
      config = config,
      rep_ids = shard_ids,
      shard = as.integer(plan_row$shard),
      base_seed = base_seed,
      max_attempts = max_attempts,
      failure_jsonl = failure_jsonl
    )
    if (length(shard_result$failed) > 0L) {
      WriteFailedCellMarker(
        path = failed_marker,
        record = list(
          timestamp = format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ"),
          config_index = config$config_index,
          block_id = config$block_id,
          shard = as.integer(plan_row$shard),
          n_failed_reps = length(shard_result$failed),
          n_success_reps = if (is.null(shard_result$data)) {
            0L
          } else {
            nrow(shard_result$data)
          },
          error = shard_result$failed[[1L]]$error
        )
      )
      cell_failed <- TRUE
      n_failed_cells <- n_failed_cells + 1L
      break
    }
    AtomicWriteTerminalEventParquet(
      data = shard_result$data,
      path = shard_file,
      expected_rep_ids = shard_ids,
      config = config
    )
    message("Wrote ", shard_file)
  }

  if (cell_failed) {
    next
  }

  shard_tables <- lapply(
    seq_along(shard_files),
    function(s) {
      arrow::read_parquet(shard_files[s])
    }
  )
  cell_data <- do.call(rbind, shard_tables)
  AtomicWriteTerminalEventParquet(
    data = cell_data,
    path = cell_file,
    expected_rep_ids = expected_cell_ids,
    config = config
  )
  unlink(shard_files)
  message("Wrote ", cell_file)
  n_complete_cells <- n_complete_cells + 1L
}

message(
  glue(
    "Block {block_id} done. complete={n_complete_cells}, ",
    "failed={n_failed_cells}, skipped={n_skipped_cells}, ",
    "assigned={nrow(block_configs)}"
  )
)
message(glue("Output directory: {out_dir}"))
