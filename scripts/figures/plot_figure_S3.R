#!/usr/bin/env Rscript
# Figure S3: terminal-event bias and coverage, MCF.
argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(dirname(normalizePath(script))), "utils", "bootstrap_regression.R"))
SetRegressionWd()
suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(ggplot2); library(optparse)
})
SourceRegressionUtil("ground_truth_io.R")
SourceRegressionUtil("terminal_plot_helpers.R")
PlotTerminalFigure("MCF", "S3")
