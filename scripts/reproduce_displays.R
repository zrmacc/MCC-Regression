#!/usr/bin/env Rscript
# Regenerate inferential displays and bias panels from validated saved runs.
argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(normalizePath(script)), "utils", "bootstrap_regression.R"))
SetRegressionWd()
suppressPackageStartupMessages({library(optparse); library(jsonlite)})
opts <- parse_args(OptionParser(option_list = list(
  make_option("--data_root", type = "character", default = "simulations"),
  make_option("--event_truth", type = "character", default = "simulations/ground_truth"),
  make_option("--terminal_truth", type = "character", default = "simulations/terminal_event_ground_truth"),
  make_option("--frailty_truth", type = "character", default = "simulations/frailty_ground_truth"),
  make_option("--out_dir", type = "character", default = "figures"),
  make_option("--table_dir", type = "character", default = "tables"),
  make_option("--manuscript_dir", type = "character", default = NULL),
  make_option("--manifest", type = "character", default = "provenance/display_manifest.json")
)))
for (k in names(opts)) if (!is.null(opts[[k]])) opts[[k]] <- ResolveRegressionPath(opts[[k]])
for (path in c(opts$out_dir, opts$table_dir, dirname(opts$manifest))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
jobs <- list()
run <- function(script, args) {
  log <- file.path(dirname(opts$manifest), paste0(sub("[.]R$", "", basename(script)), ".log"))
  command <- c(file.path(RegressionDir(), "scripts", script), args)
  tic <- Sys.time()
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(command), stdout = log, stderr = log)
  jobs[[length(jobs) + 1L]] <<- list(script = script, arguments = args, seconds = as.numeric(difftime(Sys.time(), tic, units = "secs")), status = status)
  if (status != 0L) stop("Display generation failed; see ", log, call. = FALSE)
  message("Completed ", script)
}
data_path <- function(study) file.path(opts$data_root, "pseudoreg", study)
for (figure in c("1", "2")) run(paste0("figures/plot_figure_", figure, ".R"), c("--dir", data_path("estimation"), "--ground_truth", opts$event_truth, "--out_dir", opts$out_dir))
run("figures/plot_figure_3.R", c("--dir", data_path("t1e/null"), "--out", file.path(opts$out_dir, "fig3.pdf")))
run("figures/plot_figure_4.R", c("--dir", data_path("power"), "--out_dir", opts$out_dir))
run("figures/plot_figure_S2.R", c("--dir", data_path("frailty"), "--ground_truth", opts$frailty_truth, "--out_dir", opts$out_dir))
for (figure in c("S3", "S4")) run(paste0("figures/plot_figure_", figure, ".R"), c("--dir", data_path("terminal_event"), "--ground_truth", opts$terminal_truth, "--out_dir", opts$out_dir))
run("tables/make_table_1.R", c("--dir", data_path("t1e/null"), "--out", file.path(opts$table_dir, "tab1.tex"), "--overwrite"))
table_args <- c("--dir", data_path("estimation"), "--ground_truth", opts$event_truth, "--out", file.path(opts$table_dir, "tabS1.tex"), "--summary_out", file.path(dirname(opts$manifest), "table_S1_summary.csv"), "--overwrite")
if (!is.null(opts$manuscript_dir)) table_args <- c(table_args, "--manuscript_out", file.path(opts$manuscript_dir, "tables/tabS1.tex"))
run("tables/make_table_S1.R", table_args)

mapping <- c("fig1a.pdf" = "fig2a.pdf", "fig1b.pdf" = "fig2b.pdf", "fig2a.pdf" = "fig3a.pdf", "fig2b.pdf" = "fig3b.pdf", "fig3.png" = "fig5.png", "fig4.pdf" = "fig6.pdf", "figS2a.pdf" = "fig4a.pdf", "figS2b.pdf" = "fig4b.pdf", "figS3.pdf" = "supp_fig_terminal_event_beta.pdf", "figS4.pdf" = "supp_fig_terminal_event_gamma.pdf")
if (!is.null(opts$manuscript_dir)) {
  for (file in names(mapping)) {
    target <- file.path(opts$manuscript_dir, "figures", mapping[[file]])
    if (!file.copy(file.path(opts$out_dir, file), target, overwrite = TRUE)) stop("Could not synchronize ", target)
  }
  if (!file.copy(file.path(opts$table_dir, "tab1.tex"), file.path(opts$manuscript_dir, "tables/tab1.tex"), overwrite = TRUE)) stop("Could not synchronize Table 1.")
}
sources <- c(unlist(lapply(c("estimation", "frailty", "power", "terminal_event", "t1e/null"), function(study) list.files(data_path(study), "[.]parquet$", recursive = TRUE, full.names = TRUE))), unlist(lapply(c(opts$event_truth, opts$terminal_truth, opts$frailty_truth), function(p) if (dir.exists(p)) list.files(p, "[.]parquet$", full.names = TRUE) else p)))
sources <- unique(sources[!grepl("/(shards|logs|failures)/|_shard_", sources)])
outputs <- c(file.path(opts$out_dir, names(mapping)), file.path(opts$out_dir, "fig3.pdf"), file.path(opts$table_dir, c("tab1.tex", "tabS1.tex")))
hashes <- function(paths) lapply(paths, function(path) list(path = path, md5 = unname(tools::md5sum(path))))
manifest <- list(variance_method = CorrectedVarianceMethod(), interval = "95% normal Wald", generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), R = R.version.string, jobs = jobs, inputs = hashes(sources), outputs = hashes(outputs), retained_display = "Web Figure 1 compares pseudo-values and contains no inference; its supplied artifacts are unchanged.")
jsonlite::write_json(manifest, opts$manifest, pretty = TRUE, auto_unbox = TRUE)
message("All displays generated; manifest: ", opts$manifest)
