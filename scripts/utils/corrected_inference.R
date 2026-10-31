# Shared contract for all reported regression inference.
CorrectedVarianceMethod <- function() "HC3 + second order"

InferenceCellColumns <- function(data) {
  intersect(c("n", "tau", "type", "beta_name", "beta_value", "beta_A", "frailty_variance"), names(data))
}

InferenceKeys <- function(data, columns) {
  if (!length(columns) || !all(columns %in% names(data))) stop("Missing inference key columns.", call. = FALSE)
  values <- lapply(data[columns], function(x) if (is.numeric(x)) sprintf("%.12f", round(x, 12)) else as.character(x))
  do.call(paste, c(values, sep = "/"))
}

ValidateCorrectedInferenceData <- function(data, expected_reps = NULL, tolerance = 1e-10) {
  if (!nrow(data)) stop("No inference rows supplied.", call. = FALSE)
  if (!"variance_method" %in% names(data) || anyNA(data$variance_method) || any(data$variance_method != CorrectedVarianceMethod())) {
    stop("Every inference row must have variance_method = 'HC3 + second order'.", call. = FALSE)
  }
  cells <- InferenceCellColumns(data)
  keys <- c(cells, "rep")
  if (!all(keys %in% names(data)) || anyNA(data[keys])) stop("Missing design/replicate identifiers.", call. = FALSE)
  if (anyDuplicated(InferenceKeys(data, keys))) stop("Duplicate design/replicate identifiers.", call. = FALSE)
  if (any(!is.finite(data$rep)) || any(data$rep < 1 | data$rep != as.integer(data$rep))) stop("Invalid replicate identifiers.", call. = FALSE)
  p_name <- if ("pvalue" %in% names(data)) "pvalue" else "pvalue_arm"
  if (!p_name %in% names(data)) stop("Missing corrected p-values.", call. = FALSE)
  p <- data[[p_name]]
  if (any(!is.finite(p) | p < 0 | p > 1)) stop("Invalid corrected p-values.", call. = FALSE)
  coefficient_fields <- c("estimate", "se", "lower", "upper")
  if (any(coefficient_fields %in% names(data))) {
    if (!all(coefficient_fields %in% names(data))) stop("Incomplete coefficient inference fields.", call. = FALSE)
    if (any(!is.finite(as.matrix(data[coefficient_fields]))) || any(data$se <= 0)) stop("Nonfinite inference or nonpositive SE.", call. = FALSE)
    z <- stats::qnorm(.975)
    error <- max(abs(data$lower - (data$estimate - z * data$se)), abs(data$upper - (data$estimate + z * data$se)), abs(p - 2 * stats::pnorm(-abs(data$estimate / data$se))))
    if (error > tolerance) stop("Inference does not match corrected 95% normal Wald intervals/p-values.", call. = FALSE)
  }
  if (!is.null(expected_reps)) {
    groups <- split(data$rep, InferenceKeys(data, cells))
    if (any(!vapply(groups, function(x) identical(sort(as.integer(x)), seq_len(expected_reps)), logical(1)))) stop("Incomplete or unexpected replicate IDs within a design cell.", call. = FALSE)
  }
  invisible(data)
}

CanonicalStudyConfigs <- function(study) {
  switch(study,
    estimation = expand.grid(n = c(50L, 100L, 200L), tau = c(2, 3, 4), type = c("AUC", "MCF"), beta_name = c("beta_A", "beta_X", "beta_AX"), beta_value = signif(log(c(.75, 1, 1.25)), 15), stringsAsFactors = FALSE),
    terminal_event = CanonicalStudyConfigs("estimation"),
    frailty = expand.grid(n = c(50L, 100L, 200L), tau = c(2, 3, 4), type = c("AUC", "MCF"), frailty_variance = c(0, 1, 2, 3), stringsAsFactors = FALSE),
    power = expand.grid(n = seq(50L, 200L, 25L), tau = c(2, 3, 4), type = c("AUC", "MCF"), beta_A = signif(log(c(.6, .7, .8, .9)), 15), stringsAsFactors = FALSE),
    type1 = expand.grid(n = c(50L, 100L, 200L), tau = c(2, 3, 4), type = c("AUC", "MCF"), stringsAsFactors = FALSE),
    stop("Unknown study: ", study, call. = FALSE)
  )
}

AttachCanonicalConfigIndex <- function(configs, study) {
  canonical <- CanonicalStudyConfigs(study)
  columns <- names(canonical)
  index <- match(InferenceKeys(configs, columns), InferenceKeys(canonical, columns))
  if (anyNA(index)) stop("Requested configuration is outside the published ", study, " grid.", call. = FALSE)
  if (anyDuplicated(index)) stop("Duplicated requested configurations.", call. = FALSE)
  configs$config_index <- index
  configs
}

ValidateStudyInference <- function(data, study, expected_reps = NULL) {
  if (is.null(expected_reps)) expected_reps <- switch(study, power = 500L, type1 = 100000L, 2000L)
  ValidateCorrectedInferenceData(data, expected_reps)
  canonical <- CanonicalStudyConfigs(study)
  wanted <- InferenceKeys(canonical, names(canonical))
  observed <- unique(InferenceKeys(data, names(canonical)))
  if (!setequal(wanted, observed)) stop("Incomplete or unexpected ", study, " design grid.", call. = FALSE)
  invisible(data)
}

ValidateResumedInference <- function(data, config, expected_reps, base_seed = NULL) {
  ValidateCorrectedInferenceData(data, expected_reps)
  columns <- names(config)
  if (!all(columns %in% names(data))) stop("Resumed file lacks requested cell metadata.", call. = FALSE)
  if (any(!vapply(columns, function(k) all(InferenceKeys(data, k) == InferenceKeys(config, k)), logical(1)))) stop("Resumed file metadata do not match the requested cell.", call. = FALSE)
  if (!is.null(base_seed)) {
    wanted <- vapply(data$rep, function(i) RegressionRepSeed(base_seed, config$config_index, i), integer(1))
    if (!"seed" %in% names(data) || anyNA(data$seed) || any(data$seed != wanted)) stop("Resumed file seeds do not match the requested run.", call. = FALSE)
  }
  invisible(data)
}
