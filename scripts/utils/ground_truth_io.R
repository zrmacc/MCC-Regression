# Shared helpers for bias ground-truth parquet I/O and averaging over replicates.
# Load from driver scripts via SourceRegressionUtil("ground_truth_io.R") after the
# bootstrap block (see README.md).

suppressPackageStartupMessages({
  library(dplyr)
})

#' Names of numeric ground-truth columns averaged across replicates.
#'
#' @return Character vector of column names.
GroundTruthValueColNames <- function() {
  out <- c(
    "delta_MCF",
    "delta_AUMCF",
    "beta_A_MCF",
    "beta_A_AUMCF",
    "beta_X_MCF",
    "beta_X_AUMCF",
    "beta_AX_MCF",
    "beta_AX_AUMCF"
  )
  return(out)
}

#' Names of numeric frailty ground-truth columns averaged across replicates.
#'
#' @return Character vector of column names.
FrailtyGroundTruthValueColNames <- function() {
  out <- c(
    "delta_MCF",
    "delta_AUMCF",
    "beta_A_MCF",
    "beta_A_AUMCF"
  )
  return(out)
}

#' Detect the value-column schema present in a ground-truth table.
#'
#' @param gt Ground-truth data frame.
#' @return Character vector of value columns.
DetectGroundTruthValueColNames <- function(gt) {
  full_cols <- GroundTruthValueColNames()
  frailty_cols <- FrailtyGroundTruthValueColNames()
  if (all(full_cols %in% names(gt))) {
    return(full_cols)
  }
  if (all(frailty_cols %in% names(gt))) {
    return(frailty_cols)
  }
  stop(
    "Ground truth table does not match known schemas. Missing required value columns.",
    call. = FALSE
  )
}

#' Resolve a path to one or more ground-truth parquet files.
#'
#' If \code{path} is a directory, returns all \code{*.parquet} in that directory
#' (non-recursive). If it is a file, returns that path as a length-one vector.
#'
#' @param path File or directory path.
#' @return Character vector of parquet file paths.
GroundTruthParquetFiles <- function(path) {
  if (dir.exists(path)) {
    out <- list.files(
      path = path,
      pattern = "\\.parquet$",
      full.names = TRUE,
      recursive = FALSE
    )
    if (length(out) == 0L) {
      stop("No .parquet files in directory: ", path)
    }
  } else {
    if (!file.exists(path)) {
      stop("Ground truth path not found: ", path)
    }
    out <- path
  }
  return(out)
}

#' Prepare one ground-truth table for stacking and averaging.
#'
#' Coerces key and value columns to numeric/character as expected downstream.
#' If \code{gt_rep} is absent (legacy one draw per cell), sets \code{gt_rep = 1}
#' so the table stacks cleanly with multi-replicate outputs. Drops extra
#' columns so all chunks match.
#'
#' @param gt Data frame read from a ground-truth parquet file.
#' @return Data frame with columns \code{gt_rep}, \code{beta_name},
#'   \code{beta_value}, \code{tau}, and the value columns from
#'   \code{GroundTruthValueColNames()}.
PrepareGroundTruthChunk <- function(gt) {
  value_cols <- DetectGroundTruthValueColNames(gt)
  key_cols <- c("beta_name", "beta_value", "tau")
  if ("frailty_variance" %in% names(gt)) {
    key_cols <- c(key_cols, "frailty_variance")
  }
  req <- c(key_cols, value_cols)
  miss <- setdiff(req, names(gt))
  if (length(miss) > 0L) {
    stop(
      "Ground truth table missing columns: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  gt_rep_vec <- if ("gt_rep" %in% names(gt)) {
    as.integer(gt$gt_rep)
  } else {
    rep(1L, nrow(gt))
  }
  out <- gt %>%
    dplyr::select(dplyr::all_of(req))
  out$beta_name <- as.character(out$beta_name)
  out$beta_value <- as.numeric(out$beta_value)
  out$tau <- as.numeric(out$tau)
  if ("frailty_variance" %in% names(out)) {
    out$frailty_variance <- as.numeric(out$frailty_variance)
  }
  for (v in value_cols) {
    out[[v]] <- as.numeric(out[[v]])
  }
  out$gt_rep <- gt_rep_vec
  out <- out %>%
    dplyr::select(gt_rep, dplyr::all_of(req))
  return(out)
}

#' Read and row-bind ground-truth parquet files.
#'
#' Each file is passed through \code{PrepareGroundTruthChunk()} so legacy
#' parquets (without \code{gt_rep}) combine correctly with newer outputs.
#'
#' @param files Character vector of file paths.
#' @return Combined data frame.
ReadGroundTruthParquetStack <- function(files) {
  chunks <- lapply(
    files,
    function(f) {
      raw <- arrow::read_parquet(f)
      PrepareGroundTruthChunk(raw)
    }
  )
  out <- dplyr::bind_rows(chunks)
  return(out)
}

#' Average ground-truth numeric columns over replicates per design cell.
#'
#' Groups by \code{beta_name}, \code{beta_value}, and \code{tau}. Legacy tables
#' with one row per cell are unchanged after averaging.
#'
#' @param gt_raw Stacked ground-truth data frame.
#' @return Data frame with one row per \code{(beta_name, beta_value, tau)}.
AverageGroundTruthOverReplicates <- function(gt_raw) {
  value_cols <- DetectGroundTruthValueColNames(gt_raw)
  miss <- setdiff(value_cols, names(gt_raw))
  if (length(miss) > 0L) {
    stop(
      "Ground truth data missing columns: ",
      paste(miss, collapse = ", ")
    )
  }
  req <- c("beta_name", "beta_value", "tau")
  if ("frailty_variance" %in% names(gt_raw)) {
    req <- c(req, "frailty_variance")
  }
  miss_req <- setdiff(req, names(gt_raw))
  if (length(miss_req) > 0L) {
    stop(
      "Ground truth data missing columns: ",
      paste(miss_req, collapse = ", ")
    )
  }
  out <- gt_raw %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(req))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(value_cols), mean),
      .groups = "drop"
    )
  return(out)
}

#' Load ground truth from a file or directory, then average over replicates.
#'
#' @param path Path to a parquet file or a directory of parquet files.
#' @return List with components \code{raw}, \code{averaged}, and \code{files}.
LoadGroundTruthRawAndAveraged <- function(path) {
  files <- GroundTruthParquetFiles(path)
  if (length(files) > 1L) {
    message(
      "Combining ",
      length(files),
      " ground-truth parquet file(s): ",
      paste(basename(files), collapse = ", ")
    )
  }
  raw <- ReadGroundTruthParquetStack(files)
  averaged <- AverageGroundTruthOverReplicates(raw)
  out <- list(
    raw = raw,
    averaged = averaged,
    files = files
  )
  return(out)
}

#' Round a log-scale coefficient for stable ground-truth joins.
#'
#' @param value Numeric vector.
#' @return Numeric vector rounded to 12 decimal places.
GroundTruthJoinKey <- function(value) {
  out <- round(as.numeric(value), 12)
  return(out)
}

#' List merged simulation parquet files, ignoring shards, logs, and temps.
#'
#' Plotting and validation must not ingest per-cell shards, failure sidecars,
#' logs, or in-progress \code{.tmp} writes. Matches files ending in
#' \code{.parquet} while excluding \code{_shard_} names and known sidecar
#' directories.
#'
#' @param dir_in Directory containing simulation parquet outputs.
#' @return Character vector of parquet paths.
ListMergedSimulationParquets <- function(dir_in) {
  if (!dir.exists(dir_in)) {
    stop("Input directory does not exist: ", dir_in, call. = FALSE)
  }
  files <- list.files(
    path = dir_in,
    pattern = "\\.parquet$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(files) == 0L) {
    return(character())
  }
  files <- files[!grepl("/(shards|logs|failures)(/|$)", files)]
  files <- files[!grepl("_shard_[0-9]+_of_[0-9]+", basename(files))]
  files <- files[!grepl("\\.tmp$", files)]
  return(files)
}
