#!/usr/bin/env Rscript
# Table 1: corrected normal-Wald rejection rates.
argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(dirname(normalizePath(script))), "utils", "bootstrap_regression.R"))
SetRegressionWd()
suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(ggplot2); library(optparse)
})
SourceRegressionUtil("ground_truth_io.R")
SourceRegressionUtil("t1e_sim.R")
SummarizeTypeIError <- function(dir_in, alpha = 0.05, expected_reps = 100000L) {
  if (!dir.exists(dir_in)) {
    stop("Type I error input directory does not exist: ", dir_in)
  }
  raw <- LoadCanonicalT1eData(dir_in)
  req <- T1eRequiredColumns()
  miss <- setdiff(req, names(raw))
  if (length(miss) > 0L) {
    stop("Canonical type I schema missing: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (any(raw$n_per_arm == 25L)) {
    stop("Canonical type I output must not contain 25-per-arm rows.", call. = FALSE)
  }
  if (!all(raw$n == raw$n_per_arm)) {
    stop("Type I n must equal n_per_arm.", call. = FALSE)
  }
  if (any(!is.finite(raw$pvalue_arm)) || any(raw$pvalue_arm < 0 | raw$pvalue_arm > 1)) {
    stop("Type I error p-values must be finite and lie in [0, 1].")
  }
  methods <- unique(raw$variance_method)
  if (!identical(methods, "HC3 + second order")) {
    stop(
      "Expected only variance_method = 'HC3 + second order'; found: ",
      paste(methods, collapse = ", ")
    )
  }
  expected_grid <- CanonicalT1eConfigs()
  cells <- raw %>%
    dplyr::group_by(.data$n_per_arm, .data$tau, .data$type, .data$source) %>%
    dplyr::summarise(
      replicates = dplyr::n(),
      n_reject = sum(.data$pvalue_arm < alpha),
      rejection_rate = mean(.data$pvalue_arm < alpha),
      n_retry = sum(.data$attempts > 1L, na.rm = TRUE),
      n_attempts_obs = sum(!is.na(.data$attempts)),
      .groups = "drop"
    )
  cells$mcse <- sqrt(cells$rejection_rate * (1 - cells$rejection_rate) / cells$replicates)
  wilson <- WilsonProportionCi(k = cells$n_reject, n = cells$replicates)
  cells$wilson_lower <- wilson$lower
  cells$wilson_upper <- wilson$upper
  cells$retry_rate <- ifelse(
    cells$n_attempts_obs > 0L,
    cells$n_retry / cells$n_attempts_obs,
    NA_real_
  )
  if (nrow(cells) != nrow(expected_grid)) {
    stop("Type I error output does not contain the expected 18 design cells.")
  }
  mapped <- expected_grid %>%
    dplyr::left_join(
      cells,
      by = c("n_per_arm", "tau", "type")
    )
  if (anyNA(mapped$replicates) || any(mapped$replicates != expected_reps)) {
    stop("Expected ", expected_reps, " replicates in every canonical cell.", call. = FALSE)
  }
  mcf <- cells %>%
    dplyr::filter(.data$type == "MCF") %>%
    dplyr::select(
      n = "n_per_arm",
      "tau",
      mcf_rejection_rate = "rejection_rate"
    )
  aumcf <- cells %>%
    dplyr::filter(.data$type == "AUC") %>%
    dplyr::select(
      n = "n_per_arm",
      "tau",
      aumcf_rejection_rate = "rejection_rate"
    )
  out <- mcf %>%
    dplyr::inner_join(aumcf, by = c("n", "tau")) %>%
    dplyr::arrange(.data$n, .data$tau)
  attr(out, "maximum_mcse") <- max(cells$mcse)
  attr(out, "alpha") <- alpha
  attr(out, "replicates") <- expected_reps
  attr(out, "cells") <- cells %>%
    dplyr::arrange(.data$n_per_arm, .data$tau, .data$type)
  return(out)
}

#' Independently re-aggregate canonical parquets and compare with a summary.
#'
#' @param dir_in Canonical directory.
#' @param summ Summary from \code{SummarizeTypeIError()}.
#' @param alpha Significance level.
#' @param expected_reps Replicates per cell.
#' @return Invisibly, \code{TRUE}.
ConfirmTypeIErrorSummary <- function(dir_in, summ, alpha, expected_reps) {
  again <- SummarizeTypeIError(
    dir_in = dir_in,
    alpha = alpha,
    expected_reps = expected_reps
  )
  cmp_cols <- c("n", "tau", "mcf_rejection_rate", "aumcf_rejection_rate")
  if (!isTRUE(all.equal(summ[, cmp_cols], again[, cmp_cols], tolerance = 0))) {
    stop("Independent type I aggregation does not match the written summary.", call. = FALSE)
  }
  cells <- attr(summ, "cells")
  cells2 <- attr(again, "cells")
  if (!isTRUE(all.equal(cells$n_reject, cells2$n_reject, tolerance = 0))) {
    stop("Independent type I rejection numerators do not match.", call. = FALSE)
  }
  if (!identical(nrow(cells), 18L) || sum(cells$replicates) != 18L * expected_reps) {
    stop("Independent type I row counts are incomplete.", call. = FALSE)
  }
  return(invisible(TRUE))
}

#' Build the manuscript LaTeX tabular for type I error.
#'
#' @param summ Output from code{SummarizeTypeIError()}.
#' @return Character scalar containing a LaTeX tabular environment.
BuildTypeIErrorLatex <- function(summ) {
  lines_out <- c(
    "\\begin{tabular}{rrrr}",
    "\\toprule",
    "\\makecell{Sample Size\\\\per Arm} &",
    "\\makecell{Truncation\\\\Time} &",
    "\\makecell{MCF Rejection\\\\Rate} &",
    "\\makecell{AUMCF Rejection\\\\Rate} \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(summ))) {
    row <- summ[i, ]
    if (i > 1L && row$n != summ$n[i - 1L]) {
      lines_out <- c(lines_out, "\\midrule")
    }
    lines_out <- c(
      lines_out,
      paste0(
        as.integer(row$n), " & ",
        as.integer(row$tau), " & ",
        sprintf("%.4f", row$mcf_rejection_rate), " & ",
        sprintf("%.4f", row$aumcf_rejection_rate), " \\\\"
      )
    )
  }
  lines_out <- c(lines_out, "\\bottomrule", "\\end{tabular}")
  return(paste(lines_out, collapse = "\n"))
}

#' Write the type I error table and an auditable CSV summary.
#'
#' @param summ Output from code{SummarizeTypeIError()}.
#' @param out_file Path to the LaTeX output.
#' @param overwrite Whether existing outputs may be replaced.
#' Write the type I error table and an auditable CSV summary.
#'
#' @param summ Output from code{SummarizeTypeIError()}.
#' @param out_file Path to the LaTeX output.
#' @param overwrite Whether existing outputs may be replaced.
#' @param dir_in Canonical parquet directory used to confirm the aggregation.
WriteTypeIErrorTable <- function(summ, out_file, overwrite, dir_in) {
  csv_file <- sub("[.]tex$", ".csv", out_file)
  if (identical(csv_file, out_file)) {
    csv_file <- paste0(out_file, ".csv")
  }
  cells_file <- sub("[.]csv$", "_cells.csv", csv_file)
  existing <- c(out_file, csv_file, cells_file)[file.exists(c(out_file, csv_file, cells_file))]
  if (!overwrite && length(existing) > 0L) {
    stop("Output file(s) already exist. Pass --overwrite or remove: ", paste(existing, collapse = "; "))
  }
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(BuildTypeIErrorLatex(summ), out_file)
  csv_summ <- data.frame(
    n = summ$n,
    tau = summ$tau,
    mcf_rejection_rate = summ$mcf_rejection_rate,
    aumcf_rejection_rate = summ$aumcf_rejection_rate
  )
  utils::write.csv(csv_summ, csv_file, row.names = FALSE)
  utils::write.csv(attr(summ, "cells"), cells_file, row.names = FALSE)
  ConfirmTypeIErrorSummary(
    dir_in = dir_in,
    summ = summ,
    alpha = attr(summ, "alpha"),
    expected_reps = attr(summ, "replicates")
  )
  reread <- utils::read.csv(csv_file)
  if (!isTRUE(all.equal(reread, csv_summ, check.attributes = FALSE, tolerance = 0))) {
    stop("Re-read type I CSV does not match the displayed summary.", call. = FALSE)
  }
  message("Wrote ", out_file)
  message("Wrote ", csv_file)
  message("Wrote ", cells_file)
}

opts <- parse_args(OptionParser(option_list=list(
  make_option("--dir",type="character",default="simulations/pseudoreg/t1e/null"),
  make_option("--out",type="character",default="tables/tab1.tex"),
  make_option("--expected_reps",type="integer",default=100000L),
  make_option("--alpha",type="double",default=0.05),
  make_option("--overwrite",action="store_true",default=FALSE)
)))
summ <- SummarizeTypeIError(opts$dir, opts$alpha, opts$expected_reps)
WriteTypeIErrorTable(summ, opts$out, opts$overwrite, opts$dir)
