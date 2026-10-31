# Shared helpers for the type I error per-arm study: canonical grid, schema,
# parquet validation, and atomic writes.

#' Canonical 18-cell type I error grid with stable indices.
#'
#' \code{n_per_arm} varies fastest, then \code{tau}, then \code{type}, matching
#' \code{expand.grid} column order. Indices stay attached when a run filters to
#' \code{n_per_arm = 200} only.
#'
#' @return Data frame with 18 rows and a one-based \code{config_index}.
CanonicalT1eConfigs <- function() {
  configs <- expand.grid(
    n_per_arm = as.integer(c(50, 100, 200)),
    tau = c(2, 3, 4),
    type = c("AUC", "MCF"),
    stringsAsFactors = FALSE
  )
  configs$n_per_arm <- as.integer(configs$n_per_arm)
  configs$n <- configs$n_per_arm
  configs$n_total <- as.integer(2L * configs$n_per_arm)
  configs$config_index <- seq_len(nrow(configs))
  return(configs)
}

#' Subset the canonical type I error grid while keeping global indices.
#'
#' @param configs Canonical grid from \code{CanonicalT1eConfigs()}.
#' @param n_vec Per-arm sample sizes to retain.
#' @param tau_vec Truncation times to retain.
#' @param type_vec Pseudo-value types to retain.
#' @return Filtered configuration data frame.
FilterCanonicalT1eConfigs <- function(configs, n_vec, tau_vec, type_vec) {
  out <- configs[
    configs$n_per_arm %in% as.integer(n_vec) &
      configs$tau %in% as.numeric(tau_vec) &
      configs$type %in% as.character(type_vec),
    ,
    drop = FALSE
  ]
  rownames(out) <- NULL
  return(out)
}

#' Replicate slices for the six durable 200-per-arm production blocks.
#'
#' @return Data frame with \code{block_id}, \code{rep_start}, and \code{reps}.
T1eN200BlockRanges <- function() {
  out <- data.frame(
    block_id = 1:6,
    rep_start = as.integer(c(1, 16668, 33335, 50002, 66669, 83336)),
    reps = as.integer(c(16667, 16667, 16667, 16667, 16667, 16665))
  )
  return(out)
}

#' Required columns in canonical type I error output.
#'
#' @return Character vector.
T1eRequiredColumns <- function() {
  out <- c(
    "rep",
    "pvalue_arm",
    "seed",
    "attempts",
    "variance_method",
    "n",
    "n_per_arm",
    "n_total",
    "n_arm0",
    "n_arm1",
    "tau",
    "type",
    "tau_gen",
    "batch",
    "config_index",
    "source",
    "source_n_legacy",
    "source_batch",
    "source_config_index"
  )
  return(out)
}

#' Allowed provenance labels for canonical rows.
#'
#' @return Character vector.
T1eAllowedSources <- function() {
  out <- c("legacy_total_n_relabel", "new_per_arm_simulation")
  return(out)
}

#' Directory for one canonical cell.
#'
#' @param n_per_arm Subjects per arm.
#' @param tau Truncation time.
#' @return Character scalar such as \code{n_per_arm_50_tau_2}.
T1eCellDirName <- function(n_per_arm, tau) {
  out <- paste0(
    "n_per_arm_",
    format(as.integer(n_per_arm), scientific = FALSE, trim = TRUE),
    "_tau_",
    format(as.numeric(tau), scientific = FALSE, trim = TRUE)
  )
  return(out)
}

#' List canonical type I error parquets without crossing into the archive.
#'
#' @param dir_in Canonical directory (\code{simulations/pseudoreg/t1e/null}).
#' @return Character vector of parquet paths.
ListCanonicalT1eParquets <- function(dir_in) {
  if (!dir.exists(dir_in)) {
    stop("Input directory does not exist: ", dir_in, call. = FALSE)
  }
  norm <- normalizePath(dir_in, mustWork = TRUE)
  if (grepl("archive_legacy_total_n", norm)) {
    stop(
      "Refusing to read the immutable legacy archive as canonical input: ",
      dir_in,
      call. = FALSE
    )
  }
  files <- list.files(
    path = dir_in,
    pattern = "\\.parquet$",
    full.names = TRUE,
    recursive = TRUE
  )
  files <- files[!grepl("archive_legacy_total_n", files)]
  files <- files[!grepl("\\.tmp$", files)]
  return(files)
}

#' Validate a canonical type I error table.
#'
#' @param data Data frame.
#' @param expected_rep_ids Required replicate IDs.
#' @param config Optional named list of expected metadata scalars.
#' @param allow_source Character vector of allowed \code{source} values.
#' @return Invisibly, \code{data}.
ValidateT1eTable <- function(
  data,
  expected_rep_ids,
  config = NULL,
  allow_source = T1eAllowedSources()
) {
  req <- T1eRequiredColumns()
  miss <- setdiff(req, names(data))
  if (length(miss) > 0L) {
    stop("Type I output missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (nrow(data) != length(expected_rep_ids)) {
    stop(
      "Type I output has ",
      nrow(data),
      " rows; expected ",
      length(expected_rep_ids),
      ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(data$rep) > 0L) {
    stop("Duplicate replicate IDs in type I output.", call. = FALSE)
  }
  if (!identical(sort(as.integer(data$rep)), sort(as.integer(expected_rep_ids)))) {
    stop("Type I replicate IDs do not match the expected set.", call. = FALSE)
  }
  if (any(!is.finite(data$pvalue_arm)) || any(data$pvalue_arm < 0) || any(data$pvalue_arm > 1)) {
    stop("Type I p-values must be finite and in [0, 1].", call. = FALSE)
  }
  if (any(!is.finite(data$seed))) {
    stop("Type I seeds must be finite.", call. = FALSE)
  }
  if (!all(data$variance_method == "HC3 + second order")) {
    stop("Type I variance_method must be HC3 + second order.", call. = FALSE)
  }
  if (!all(data$source %in% allow_source)) {
    stop("Unexpected type I source label.", call. = FALSE)
  }
  if (!all(data$n == data$n_per_arm)) {
    stop("Type I n must equal n_per_arm.", call. = FALSE)
  }
  if (!all(data$n_total == 2L * data$n_per_arm)) {
    stop("Type I n_total must equal 2 * n_per_arm.", call. = FALSE)
  }
  if (!all(data$n_arm0 == data$n_per_arm) || !all(data$n_arm1 == data$n_per_arm)) {
    stop("Type I arm counts must equal n_per_arm.", call. = FALSE)
  }
  if (any(data$n_per_arm == 25L)) {
    stop("Canonical type I output must not contain 25-per-arm rows.", call. = FALSE)
  }
  if (!is.null(config)) {
    CheckScalar <- function(col, expected) {
      if (!all(data[[col]] == expected)) {
        stop("Type I column ", col, " does not match the expected configuration.", call. = FALSE)
      }
      return(invisible(NULL))
    }
    if (!is.null(config$n_per_arm)) {
      CheckScalar("n_per_arm", as.integer(config$n_per_arm))
    }
    if (!is.null(config$tau)) {
      CheckScalar("tau", as.numeric(config$tau))
    }
    if (!is.null(config$type)) {
      CheckScalar("type", as.character(config$type))
    }
    if (!is.null(config$tau_gen)) {
      CheckScalar("tau_gen", as.numeric(config$tau_gen))
    }
    if (!is.null(config$config_index)) {
      CheckScalar("config_index", as.integer(config$config_index))
    }
    if (!is.null(config$source)) {
      CheckScalar("source", as.character(config$source))
    }
  }
  return(invisible(data))
}

#' Test whether an existing parquet is complete for a replicate slice.
#'
#' @param path Parquet path.
#' @param expected_rep_ids Required replicate IDs.
#' @param config Optional expected metadata.
#' @param allow_source Allowed source labels.
#' @return Logical scalar.
IsCompleteT1eParquet <- function(
  path,
  expected_rep_ids,
  config = NULL,
  allow_source = T1eAllowedSources()
) {
  ok <- tryCatch(
    {
      data <- arrow::read_parquet(path)
      ValidateT1eTable(
        data = data,
        expected_rep_ids = expected_rep_ids,
        config = config,
        allow_source = allow_source
      )
      TRUE
    },
    error = function(e) FALSE
  )
  return(isTRUE(ok))
}

#' Write a type I parquet atomically after validation.
#'
#' @param data Data frame to write.
#' @param path Final parquet path.
#' @param expected_rep_ids Required replicate IDs.
#' @param config Optional expected metadata.
#' @param allow_source Allowed source labels.
#' @return Invisibly, \code{path}.
AtomicWriteT1eParquet <- function(
  data,
  path,
  expected_rep_ids,
  config = NULL,
  allow_source = T1eAllowedSources()
) {
  ValidateT1eTable(
    data = data,
    expected_rep_ids = expected_rep_ids,
    config = config,
    allow_source = allow_source
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(path, ".tmp")
  if (file.exists(tmp_path)) {
    unlink(tmp_path)
  }
  arrow::write_parquet(data, tmp_path)
  checked <- arrow::read_parquet(tmp_path)
  ValidateT1eTable(
    data = checked,
    expected_rep_ids = expected_rep_ids,
    config = config,
    allow_source = allow_source
  )
  if (file.exists(path)) {
    unlink(path)
  }
  renamed <- file.rename(tmp_path, path)
  if (!isTRUE(renamed)) {
    stop("Failed to rename ", tmp_path, " to ", path, ".", call. = FALSE)
  }
  return(invisible(path))
}

#' Wilson score interval for a binomial proportion.
#'
#' @param k Number of successes.
#' @param n Number of trials.
#' @param conf_level Confidence level (default 0.95).
#' @return Data frame with \code{lower} and \code{upper}.
WilsonProportionCi <- function(k, n, conf_level = 0.95) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  k <- as.numeric(k)
  n <- as.numeric(n)
  p_hat <- ifelse(n > 0, k / n, NA_real_)
  z2 <- z^2
  denom <- 1 + z2 / n
  centre <- (p_hat + z2 / (2 * n)) / denom
  half <- z * sqrt(p_hat * (1 - p_hat) / n + z2 / (4 * n^2)) / denom
  lower <- pmax(0, centre - half)
  upper <- pmin(1, centre + half)
  bad <- !is.finite(n) | n <= 0
  lower[bad] <- NA_real_
  upper[bad] <- NA_real_
  out <- data.frame(lower = lower, upper = upper)
  return(out)
}

#' Load canonical type I error parquets from a directory.
#'
#' @param dir_in Canonical directory.
#' @return Combined data frame.
LoadCanonicalT1eData <- function(dir_in) {
  files <- ListCanonicalT1eParquets(dir_in)
  if (length(files) == 0L) {
    stop("No canonical parquet files found in ", dir_in, call. = FALSE)
  }
  out <- dplyr::bind_rows(lapply(files, arrow::read_parquet))
  if (!all(T1eRequiredColumns() %in% names(out))) stop("Incomplete type I schema.")
  keys <- paste(out$n_per_arm, out$tau, out$type, out$rep, sep="/")
  if (anyDuplicated(keys)) stop("Duplicate type I design/replicate keys.")
  if (!all(out$variance_method == "HC3 + second order")) stop("Incorrect type I covariance method.")
  ValidateStudyInference(out, "type1")
  return(out)
}

#' Delete a parquet and any in-progress \code{.tmp} sibling.
#'
#' @param path Final parquet path.
#' @return Invisibly, \code{NULL}.
DeleteT1eParquet <- function(path) {
  if (file.exists(path)) {
    unlink(path)
  }
  tmp_path <- paste0(path, ".tmp")
  if (file.exists(tmp_path)) {
    unlink(tmp_path)
  }
  return(invisible(NULL))
}
