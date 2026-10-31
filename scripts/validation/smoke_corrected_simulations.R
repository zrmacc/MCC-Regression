#!/usr/bin/env Rscript
# Tiny simulations in an explicitly isolated directory; no production reruns.
argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(dirname(normalizePath(script))), "utils", "bootstrap_regression.R"))
suppressPackageStartupMessages({library(arrow); library(optparse)})
opts <- parse_args(OptionParser(option_list = list(make_option("--out_dir", type = "character"))))
if (is.null(opts$out_dir)) stop("--out_dir is required and must be outside production results.")
out <- normalizePath(opts$out_dir, mustWork = FALSE)
if (dir.exists(out)) stop("Use a fresh smoke-test output directory.")
dir.create(out, recursive = TRUE)
run <- function(name, extra, resume = FALSE) {
  args <- c(file.path(RegressionDir(), "scripts", "sims", paste0("sim_", name, ".R")), "--n", "50", "--tau", "2", "--type", "AUC,MCF", "--reps", "2", "--batch", "validation", if (name == "power") "--out_dir" else "--outdir", file.path(out, name), extra, if (resume) "--resume")
  log <- file.path(out, paste0(name, if (resume) "_resume" else "", ".log"))
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(args), stdout = log, stderr = log)
  if (status != 0L) stop("Smoke failed: ", log)
}
jobs <- list(estimation = c("--beta_name", "beta_A", "--beta_value", "0"), frailty = c("--frailty_variance", "0,1"), power = paste0("--beta_A=", format(log(.8), digits = 15)), terminal_event = c("--beta_name", "beta_A", "--beta_value", "0", "--n_blocks", "1", "--block_id", "1", "--shard_size", "1"), type1 = character())
for (name in names(jobs)) {
  run(name, jobs[[name]])
  files <- list.files(file.path(out, name), "[.]parquet$", recursive = TRUE, full.names = TRUE)
  if (!length(files)) stop("Smoke run wrote no data: ", name)
  before <- tools::md5sum(files)
  times <- file.info(files)$mtime
  run(name, jobs[[name]], TRUE)
  stopifnot(identical(before, tools::md5sum(files)), identical(times, file.info(files)$mtime))
  raw <- dplyr::bind_rows(lapply(files, read_parquet))
  ValidateCorrectedInferenceData(raw, expected_reps = 2L)
  message(name, ": generation, normal-Wald schema, and unchanged resume: PASS")
}
message("All five corrected simulation drivers passed. Outputs are isolated under ", out)
