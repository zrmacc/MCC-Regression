# Table S1 and its manuscript longtable body share one corrected summary.
SummarizeSeCalibration <- function(raw, truth, expected_reps = 2000L) {
  ValidateCorrectedInferenceData(raw, expected_reps)
  raw <- raw[toupper(raw$type) == "AUC" & raw$beta_name == "beta_A" & raw$term == "A", , drop = FALSE]
  expected <- CanonicalStudyConfigs("estimation")
  expected <- expected[expected$type == "AUC" & expected$beta_name == "beta_A", , drop = FALSE]
  if (!setequal(unique(InferenceKeys(raw, names(expected))), InferenceKeys(expected, names(expected)))) stop("Table S1 requires all 27 AUMCF treatment cells.", call. = FALSE)
  gt <- truth[truth$beta_name == "beta_A", , drop = FALSE]
  gt <- data.frame(tau = gt$tau, beta_value_key = GroundTruthJoinKey(gt$beta_value), truth = gt$beta_A_AUMCF)
  if (anyDuplicated(InferenceKeys(gt, c("tau", "beta_value_key")))) stop("Nonunique Table S1 truth keys.", call. = FALSE)
  raw$beta_value_key <- GroundTruthJoinKey(raw$beta_value)
  joined <- dplyr::left_join(raw, gt, by = c("tau", "beta_value_key"))
  if (any(!is.finite(joined$truth))) stop("Missing/nonfinite Table S1 truth.", call. = FALSE)
  joined$covers <- joined$lower <= joined$truth & joined$truth <= joined$upper
  joined |>
    dplyr::group_by(.data$n, .data$tau, .data$beta_value) |>
    dplyr::summarise(replicates = dplyr::n(), empirical_sd = stats::sd(.data$estimate), rms_corrected_se = sqrt(mean(.data$se^2)), coverage = mean(.data$covers), .groups = "drop") |>
    dplyr::arrange(.data$n, .data$tau, .data$beta_value)
}

# With longtable=TRUE, the supplement owns the environment, caption, label, and initial top rule.
BuildSeCalibrationTable <- function(summary, longtable = FALSE) {
  header <- c("\\makecell{N per\\\\Arm} &", "\\makecell{Truncation\\\\Time} &", "\\makecell{Rate\\\\Parameter} &", "\\makecell{Empirical\\\\SD} &", "\\makecell{Corrected\\\\SE} &", "\\makecell{Coverage} \\\\")
  start <- if (longtable) c(header, "\\midrule", "\\endfirsthead", "\\multicolumn{6}{l}{\\textit{Web Table \\thetable{} (continued)}}\\\\", "\\toprule", header, "\\midrule", "\\endhead", "\\bottomrule", "\\endfoot") else c("\\begin{tabular}{lrrccc}", "\\toprule", header, "\\midrule")
  rows <- vapply(seq_len(nrow(summary)), function(i) {
    r <- summary[i, ]
    rate <- if (abs(r$beta_value) < 1e-12) "ln(1.0)" else sprintf("ln(%.2f)", exp(r$beta_value))
    sprintf("%d & %d & %s & %.4f & %.4f & %.3f \\\\", r$n, r$tau, rate, r$empirical_sd, r$rms_corrected_se, r$coverage)
  }, character(1))
  c(start, rows, if (!longtable) c("\\bottomrule", "\\end{tabular}"))
}
