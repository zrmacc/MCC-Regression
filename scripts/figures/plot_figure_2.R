#!/usr/bin/env Rscript
# Figure 2: coverage; panel a MCF, panel b AUMCF.
argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(dirname(normalizePath(script))), "utils", "bootstrap_regression.R"))
SetRegressionWd()
suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(ggplot2); library(optparse)
})
SourceRegressionUtil("ground_truth_io.R")
SourceRegressionUtil("coverage_plot_helpers.R")
opts <- parse_args(OptionParser(option_list = list(
  make_option("--dir", type = "character", default = "simulations/pseudoreg/estimation"),
  make_option("--ground_truth", type = "character", default = "simulations/ground_truth"),
  make_option("--out_dir", type = "character", default = "figures")
)))
truth <- LoadGroundTruthRawAndAveraged(opts$ground_truth)$averaged
raw <- LoadCoverageData(opts$dir)
ValidateStudyInference(raw, "estimation")
raw <- JoinGroundTruthAndCoverage(raw, GroundTruthLong(truth))
plots <- AddPlotLabels(AggregateCoverage(raw))
if (anyNA(raw$truth) || any(!is.finite(raw$truth))) stop("Missing/nonfinite regression truth.")
for (typ in c("MCF", "AUC")) {
  if (!any(plots$type == typ)) stop("Missing pseudo-value type: ", typ)
  panel <- if (typ == "MCF") "a" else "b"
  p <- BuildCoveragePlot(plots[plots$type == typ, , drop = FALSE])
  SaveRegressionFigure(p, file.path(opts$out_dir, paste0("fig2", panel, ".pdf")), 10, 9)
}
