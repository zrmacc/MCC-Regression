# Shared helpers for the terminal-event PseudoReg simulation driver and
# validator. Requires regression_dgp.R (NumericSlug, RegressionRepSeed).

#' Canonical 162-cell terminal-event production grid with stable indices.
#'
#' Uses the same \code{expand.grid} column order as
#' \code{scripts/sims/estimation/sim_estimation.R}: \code{n} varies fastest.
#'
#' @return Data frame with 162 rows and a one-based \code{config_index}.
CanonicalTerminalEventConfigs <- function() {
  configs <- expand.grid(
    n = as.integer(c(50, 100, 200)),
    tau = c(2, 3, 4),
    type = c("AUC", "MCF"),
    beta_name = c("beta_A", "beta_X", "beta_AX"),
    beta_value = log(c(0.75, 1, 1.25)),
    stringsAsFactors = FALSE
  )
  configs$n <- as.integer(configs$n)
  configs$config_index <- seq_len(nrow(configs))
  return(configs)
}

#' Round a log-scale coefficient for stable membership tests.
#'
#' @param value Numeric vector.
#' @return Numeric vector rounded to 12 decimal places.
RoundBetaValue <- function(value) {
  out <- round(as.numeric(value), 12)
  return(out)
}

#' Subset the canonical grid by CLI filters while keeping global indices.
#'
#' @param configs Canonical grid from \code{CanonicalTerminalEventConfigs()}.
#' @param n_vec Sample sizes per arm to retain.
#' @param tau_vec Truncation times to retain.
#' @param type_vec Pseudo-value types to retain.
#' @param beta_name_vec Coefficient names to retain.
#' @param beta_value_vec Log-scale coefficient values to retain.
#' @return Filtered data frame, possibly with zero rows.
FilterCanonicalTerminalEventConfigs <- function(
  configs,
  n_vec,
  tau_vec,
  type_vec,
  beta_name_vec,
  beta_value_vec
) {
  keep_n <- as.integer(configs$n) %in% as.integer(n_vec)
  keep_tau <- configs$tau %in% as.numeric(tau_vec)
  keep_type <- toupper(configs$type) %in% toupper(type_vec)
  keep_name <- configs$beta_name %in% beta_name_vec
  keep_value <- RoundBetaValue(configs$beta_value) %in% RoundBetaValue(beta_value_vec)
  out <- configs[keep_n & keep_tau & keep_type & keep_name & keep_value, , drop = FALSE]
  return(out)
}

#' Map a global configuration index to a one-based block identifier.
#'
#' @param config_index One-based global configuration index.
#' @param n_blocks Number of execution blocks.
#' @return Integer block id in \code{1:n_blocks}.
TerminalEventBlockId <- function(config_index, n_blocks) {
  out <- ((as.integer(config_index) - 1L) %% as.integer(n_blocks)) + 1L
  return(out)
}

#' Build the design-cell filename stem used for shards and merged files.
#'
#' @param config_index Global configuration index.
#' @param n Sample size per arm.
#' @param tau Truncation time.
#' @param type Pseudo-value type.
#' @param beta_name Varied coefficient name.
#' @param beta_value Varied coefficient value.
#' @param batch Batch label.
#' @return Character scalar.
TerminalEventConfigSlug <- function(
  config_index,
  n,
  tau,
  type,
  beta_name,
  beta_value,
  batch
) {
  out <- paste(
    sprintf("config_%03d", as.integer(config_index)),
    paste0("n_", as.integer(n)),
    paste0("tau_", NumericSlug(tau)),
    type,
    beta_name,
    paste0("value_", NumericSlug(beta_value)),
    batch,
    sep = "_"
  )
  return(out)
}

#' Split replicate IDs into contiguous shards.
#'
#' @param reps Total number of replicates.
#' @param shard_size Maximum replicates per shard.
#' @return Data frame with \code{shard}, \code{rep_start}, \code{rep_end},
#'   and \code{n_rep}.
TerminalEventShardPlan <- function(reps, shard_size) {
  reps <- as.integer(reps)
  shard_size <- as.integer(shard_size)
  if (reps < 1L || shard_size < 1L) {
    stop("reps and shard_size must be positive integers.", call. = FALSE)
  }
  n_shards <- as.integer(ceiling(reps / shard_size))
  starts <- as.integer((seq_len(n_shards) - 1L) * shard_size + 1L)
  ends <- as.integer(pmin(starts + shard_size - 1L, reps))
  out <- data.frame(
    shard = seq_len(n_shards),
    n_shards = n_shards,
    rep_start = starts,
    rep_end = ends,
    n_rep = ends - starts + 1L,
    stringsAsFactors = FALSE
  )
  return(out)
}

#' Convert a named list of scalars to a one-line JSON object.
#'
#' @param record Named list of length-one values.
#' @return Character scalar.
ScalarListToJson <- function(record) {
  keys <- names(record)
  if (is.null(keys) || any(!nzchar(keys))) {
    stop("JSON records must have named fields.", call. = FALSE)
  }
  EncodeJsonValue <- function(value) {
    if (length(value) != 1L) {
      stop("JSON record fields must be scalars.", call. = FALSE)
    }
    if (is.null(value) || (length(value) == 1L && is.na(value))) {
      return("null")
    }
    if (is.logical(value)) {
      if (isTRUE(value)) {
        return("true")
      }
      return("false")
    }
    if (is.numeric(value)) {
      return(format(value, scientific = TRUE, digits = 15, trim = TRUE))
    }
    text <- as.character(value)
    text <- gsub("\\", "\\\\", text, fixed = TRUE)
    text <- gsub("\"", "\\\"", text, fixed = TRUE)
    text <- gsub("\n", "\\n", text, fixed = TRUE)
    text <- gsub("\r", "\\r", text, fixed = TRUE)
    return(paste0("\"", text, "\""))
  }
  pieces <- vapply(
    keys,
    function(key) {
      paste0("\"", key, "\":", EncodeJsonValue(record[[key]]))
    },
    character(1L)
  )
  out <- paste0("{", paste(pieces, collapse = ","), "}")
  return(out)
}

#' Append one JSON object to a JSONL file.
#'
#' @param path Destination JSONL path.
#' @param record Named list of scalars.
#' @return Invisibly, \code{path}.
AppendJsonlRecord <- function(path, record) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  line <- ScalarListToJson(record)
  cat(line, "\n", file = path, append = TRUE, sep = "")
  return(invisible(path))
}

#' Write a finished-failed cell marker so resume will not recompute it.
#'
#' @param path Destination JSON path.
#' @param record Named list of scalars describing the failure.
#' @return Invisibly, \code{path}.
WriteFailedCellMarker <- function(path, record) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(ScalarListToJson(record), con = path)
  return(invisible(path))
}

#' Test whether a finished-failed marker exists for a cell.
#'
#' @param path Marker path.
#' @return Logical scalar.
CellFailedMarkerExists <- function(path) {
  out <- file.exists(path)
  return(out)
}

#' Required columns in terminal-event simulation output.
#'
#' @return Character vector.
TerminalEventRequiredColumns <- function() {
  out <- c(
    "rep",
    "seed",
    "attempt",
    "config_index",
    "block_id",
    "shard",
    "term",
    "estimate",
    "se",
    "lower",
    "upper",
    "pvalue",
    "variance_method",
    "n",
    "tau",
    "tau_gen",
    "type",
    "beta_name",
    "beta_value",
    "batch"
  )
  return(out)
}

#' Validate a terminal-event inference table.
#'
#' @param data Data frame of replicate-level output.
#' @param expected_rep_ids Integer vector of required replicate IDs.
#' @param config Optional named list of expected metadata scalars.
#' @return Invisibly, \code{TRUE}.
ValidateTerminalEventTable <- function(
  data,
  expected_rep_ids,
  config = NULL
) {
  if (!is.data.frame(data)) {
    stop("Terminal-event output must be a data frame.", call. = FALSE)
  }
  required <- TerminalEventRequiredColumns()
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "Terminal-event output missing columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  expected_rep_ids <- as.integer(expected_rep_ids)
  if (nrow(data) != length(expected_rep_ids)) {
    stop(
      "Expected ",
      length(expected_rep_ids),
      " rows but found ",
      nrow(data),
      ".",
      call. = FALSE
    )
  }
  if (!setequal(as.integer(data$rep), expected_rep_ids)) {
    stop("Replicate IDs do not match the assigned shard or cell.", call. = FALSE)
  }
  if (any(duplicated(as.integer(data$rep)))) {
    stop("Duplicate replicate IDs in terminal-event output.", call. = FALSE)
  }
  numeric_cols <- c("estimate", "se", "lower", "upper", "pvalue")
  values <- as.matrix(data[, numeric_cols, drop = FALSE])
  if (any(!is.finite(values))) {
    stop("Non-finite inference values in terminal-event output.", call. = FALSE)
  }
  if (any(data$se <= 0)) {
    stop("Standard errors must be positive.", call. = FALSE)
  }
  if (any(data$pvalue < 0 | data$pvalue > 1)) {
    stop("P-values must lie in [0, 1].", call. = FALSE)
  }
  if (any(data$lower > data$estimate | data$estimate > data$upper)) {
    stop("Interval bounds must satisfy lower <= estimate <= upper.", call. = FALSE)
  }
  if (!all(data$variance_method == "HC3 + second order")) {
    stop("variance_method must be 'HC3 + second order' for every row.", call. = FALSE)
  }
  if (any(as.integer(data$attempt) < 1L | as.integer(data$attempt) > 5L)) {
    stop("attempt must lie in 1:5.", call. = FALSE)
  }
  if (!is.null(config)) {
    CheckScalar <- function(column, expected) {
      observed <- unique(data[[column]])
      if (length(observed) != 1L || !isTRUE(all.equal(observed, expected))) {
        stop(
          "Metadata mismatch for ",
          column,
          ": expected ",
          expected,
          ".",
          call. = FALSE
        )
      }
      return(invisible(NULL))
    }
    if (!is.null(config$config_index)) {
      CheckScalar("config_index", as.integer(config$config_index))
    }
    if (!is.null(config$block_id)) {
      CheckScalar("block_id", as.integer(config$block_id))
    }
    if (!is.null(config$n)) {
      CheckScalar("n", as.integer(config$n))
    }
    if (!is.null(config$tau)) {
      CheckScalar("tau", as.numeric(config$tau))
    }
    if (!is.null(config$type)) {
      CheckScalar("type", as.character(config$type))
    }
    if (!is.null(config$beta_name)) {
      CheckScalar("beta_name", as.character(config$beta_name))
    }
    if (!is.null(config$beta_value)) {
      CheckScalar("beta_value", as.numeric(config$beta_value))
    }
    if (!is.null(config$tau_gen)) {
      CheckScalar("tau_gen", as.numeric(config$tau_gen))
    }
    if (!is.null(config$batch)) {
      CheckScalar("batch", as.character(config$batch))
    }
  }
  return(invisible(TRUE))
}

#' Validate a parquet file on disk, returning FALSE instead of throwing.
#'
#' @param path Parquet path.
#' @param expected_rep_ids Required replicate IDs.
#' @param config Optional expected metadata.
#' @return Logical scalar.
TryValidateTerminalEventParquet <- function(
  path,
  expected_rep_ids,
  config = NULL
) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  ok <- tryCatch(
    {
      data <- arrow::read_parquet(path)
      ValidateTerminalEventTable(
        data = data,
        expected_rep_ids = expected_rep_ids,
        config = config
      )
      TRUE
    },
    error = function(e) FALSE
  )
  return(isTRUE(ok))
}

#' Write a parquet file atomically and validate the result.
#'
#' Writes \code{path.tmp} in the destination directory, validates it, then
#' renames to \code{path}.
#'
#' @param data Data frame to write.
#' @param path Final parquet path.
#' @param expected_rep_ids Required replicate IDs.
#' @param config Optional expected metadata.
#' @return Invisibly, \code{path}.
AtomicWriteTerminalEventParquet <- function(
  data,
  path,
  expected_rep_ids,
  config = NULL
) {
  ValidateTerminalEventTable(
    data = data,
    expected_rep_ids = expected_rep_ids,
    config = config
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(path, ".tmp")
  if (file.exists(tmp_path)) {
    unlink(tmp_path)
  }
  arrow::write_parquet(data, tmp_path)
  checked <- arrow::read_parquet(tmp_path)
  ValidateTerminalEventTable(
    data = checked,
    expected_rep_ids = expected_rep_ids,
    config = config
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

#' Delete a path if it exists, including in-progress \code{.tmp} siblings.
#'
#' @param path File path.
#' @return Invisibly, \code{NULL}.
UnlinkTerminalEventPath <- function(path) {
  if (file.exists(path)) {
    unlink(path)
  }
  tmp_path <- paste0(path, ".tmp")
  if (file.exists(tmp_path)) {
    unlink(tmp_path)
  }
  return(invisible(NULL))
}

#' Classify errors that should stop the whole block rather than retry a cell.
#'
#' @param error An error condition.
#' @return Logical scalar.
IsSystemicPseudoRegError <- function(error) {
  message_text <- conditionMessage(error)
  pattern <- paste(
    "does not contain PseudoReg",
    "Unexpected PseudoReg variance",
    "there is no package called",
    "could not find function",
    sep = "|"
  )
  out <- grepl(pattern, message_text)
  return(out)
}

#' Compare cell-mean estimates to a previous simulation using Monte Carlo SEs.
#'
#' @param new_data New replicate-level data frame.
#' @param old_data Previous replicate-level data frame.
#' @param se_multiplier Flag threshold in Monte Carlo SEs of the difference.
#' @return Data frame of cell-level comparisons.
CompareTerminalEventMeans <- function(
  new_data,
  old_data,
  se_multiplier = 3
) {
  SummarizeMeans <- function(data, suffix) {
    data$beta_value_key <- RoundBetaValue(data$beta_value)
    out <- data %>%
      dplyr::group_by(
        .data$n,
        .data$tau,
        .data$type,
        .data$beta_name,
        .data$beta_value_key
      ) %>%
      dplyr::summarise(
        mean_estimate = mean(.data$estimate),
        sd_estimate = stats::sd(.data$estimate),
        n_rep = dplyr::n(),
        beta_value = mean(.data$beta_value),
        .groups = "drop"
      )
    names(out)[names(out) == "mean_estimate"] <- paste0("mean_", suffix)
    names(out)[names(out) == "sd_estimate"] <- paste0("sd_", suffix)
    names(out)[names(out) == "n_rep"] <- paste0("n_rep_", suffix)
    return(out)
  }
  new_sum <- SummarizeMeans(new_data, "new")
  old_sum <- SummarizeMeans(old_data, "old")
  out <- dplyr::inner_join(
    new_sum,
    old_sum,
    by = c("n", "tau", "type", "beta_name", "beta_value_key")
  )
  out$diff <- out$mean_new - out$mean_old
  out$se_diff <- sqrt(
    (out$sd_new^2) / out$n_rep_new + (out$sd_old^2) / out$n_rep_old
  )
  out$abs_z <- abs(out$diff / out$se_diff)
  out$flagged <- out$abs_z > se_multiplier
  return(out)
}
