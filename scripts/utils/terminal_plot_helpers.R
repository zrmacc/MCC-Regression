# Shared assembly for Figures S3 (MCF) and S4 (AUMCF).
PlotTerminalFigure <- function(type, number) {
  opts <- parse_args(OptionParser(option_list = list(
    make_option("--dir", type="character", default="simulations/pseudoreg/terminal_event"),
    make_option("--ground_truth", type="character", default="simulations/terminal_event_ground_truth"),
    make_option("--out_dir", type="character", default="figures")
  )))
  bias <- new.env(parent=.GlobalEnv)
  coverage <- new.env(parent=.GlobalEnv)
  sys.source(file.path(RegressionDir(), "scripts/utils/bias_plot_helpers.R"), bias)
  sys.source(file.path(RegressionDir(), "scripts/utils/coverage_plot_helpers.R"), coverage)
  raw <- bias$LoadBiasData(opts$dir)
  ValidateStudyInference(raw, "terminal_event")
  raw <- raw[toupper(raw$type)==type, , drop=FALSE]
  if (!nrow(raw)) stop("No data for ", type)
  truth <- LoadGroundTruthRawAndAveraged(opts$ground_truth)$averaged
  b <- bias$JoinGroundTruthAndBias(raw, bias$GroundTruthLong(truth))
  if (anyNA(b$truth) || any(!is.finite(b$truth))) stop("Missing/nonfinite terminal-event truth.")
  b <- bias$AddPlotLabels(bias$FilterMatchingTerm(b))
  c <- coverage$JoinGroundTruthAndCoverage(raw, coverage$GroundTruthLong(truth))
  c <- coverage$AddPlotLabels(coverage$AggregateCoverage(c))
  # Supply a null graphics device for layout calculations, avoiding Rplots.pdf.
  grDevices::pdf(file=NULL)
  on.exit(grDevices::dev.off(), add=TRUE)
  # Match main-figure typography and panel heights; retain every observation.
  bias_bound <- max(abs(pretty(c(-max(abs(b$bias)), max(abs(b$bias))), n = 2)))
  p_bias <- bias$BuildBiasPlot(b, base_size = MainFigureBaseSize(10), n_breaks = 3) +
    scale_y_continuous(name = "Bias", breaks = c(-bias_bound, 0, bias_bound)) +
    coord_cartesian(ylim = c(-bias_bound, bias_bound))
  p_coverage <- coverage$BuildCoveragePlot(c, base_size = MainFigureBaseSize(10), n_breaks = 3)
  SaveRegressionFigure(p_bias, file.path(opts$out_dir, paste0("fig",number,"a.pdf")),10,8)
  SaveRegressionFigure(p_coverage, file.path(opts$out_dir, paste0("fig",number,"b.pdf")),10,8)
  # Retain the combined export for readers who prefer a single source figure.
  p <- gridExtra::arrangeGrob(p_bias, p_coverage, ncol=1)
  SaveRegressionFigure(p, file.path(opts$out_dir, paste0("fig",number,".pdf")),10,16)
}
