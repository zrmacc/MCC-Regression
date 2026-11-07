#!/usr/bin/env Rscript
# Verify active figure numbering, labels, filenames, and supplied asset copies.
argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(dirname(normalizePath(script))), "utils", "bootstrap_regression.R"))
suppressPackageStartupMessages(library(optparse))
SourceRegressionUtil("display_assets.R")
opts <- parse_args(OptionParser(option_list = list(
  make_option("--manuscript_dir", type = "character", default = NULL),
  make_option("--figure_dir", type = "character", default = "figures")
)))
if (is.null(opts$manuscript_dir)) stop("Supply --manuscript_dir with main.tex, supplement.tex, and figures/.", call. = FALSE)
assets <- ValidateManuscriptFigures(ResolveRegressionPath(opts$manuscript_dir), ResolveRegressionPath(opts$figure_dir))
message("PASS: ", sum(!is.na(assets$source)), " generated figure assets match the active main-text/supplement figures; the clinical figure is supplied separately.")
