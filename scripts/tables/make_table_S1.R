#!/usr/bin/env Rscript
# Table S1: corrected SE calibration and normal-Wald coverage.
argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(dirname(normalizePath(script))), "utils", "bootstrap_regression.R"))
SetRegressionWd()
suppressPackageStartupMessages({library(arrow); library(dplyr); library(optparse)})
SourceRegressionUtil("ground_truth_io.R")
SourceRegressionUtil("coverage_plot_helpers.R")
SourceRegressionUtil("se_calibration.R")

opts <- parse_args(OptionParser(option_list = list(
  make_option("--dir", type = "character", default = "simulations/pseudoreg/estimation"),
  make_option("--ground_truth", type = "character", default = "simulations/ground_truth"),
  make_option("--out", type = "character", default = "tables/tabS1.tex"),
  make_option("--manuscript_out", type = "character", default = NULL),
  make_option("--summary_out", type = "character", default = NULL),
  make_option("--overwrite", action = "store_true", default = FALSE)
)))
targets <- c(opts$out, opts$manuscript_out, opts$summary_out)
if (!opts$overwrite && any(file.exists(targets))) stop("Output exists; pass --overwrite.")
raw <- LoadCoverageData(opts$dir)
truth <- LoadGroundTruthRawAndAveraged(opts$ground_truth)$averaged
summary <- SummarizeSeCalibration(raw, truth)
if (nrow(summary) != 27L || any(summary$replicates != 2000L)) stop("Incomplete Table S1.")
for (path in targets) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
writeLines(BuildSeCalibrationTable(summary), opts$out)
if (!is.null(opts$manuscript_out)) writeLines(BuildSeCalibrationTable(summary, longtable = TRUE), opts$manuscript_out)
if (!is.null(opts$summary_out)) utils::write.csv(summary, opts$summary_out, row.names = FALSE)
message("Wrote corrected Table S1 from ", sum(summary$replicates), " replicates.")
