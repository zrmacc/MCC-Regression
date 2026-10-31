#!/usr/bin/env Rscript
# Fixed-seed checks of the covariance, provenance guards, and SE calibration.
argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(dirname(normalizePath(script))), "utils", "bootstrap_regression.R"))
suppressPackageStartupMessages({library(MCC); library(dplyr)})
SourceRegressionUtil("regression_dgp.R")
SourceRegressionUtil("ground_truth_io.R")
SourceRegressionUtil("se_calibration.R")
MustFail <- function(expr) {
  failed <- tryCatch({force(expr); FALSE}, error = function(e) TRUE)
  if (!failed) stop("Expected validation failure.", call. = FALSE)
}
ValidatePseudoReg()
set.seed(2026103002L)
w <- MakeDesignData(20L)
data <- MCC::GenData(covariates = cbind(1, w$A, w$X, w$AX), beta_event = c(0, -.1, .1, 0), beta_death = rep(0, 4), tau = 4)
for (type in c("MCF", "AUC")) {
  fit <- MCC::PseudoReg(~ A + X + A:X, data, tau = 2, type = type, covariates = w)
  pseudo <- MCC::GenPseudo(data, tau = 2, type = type)
  external <- stats::lm(pseudo ~ A + X + A:X, data = dplyr::inner_join(pseudo, w, by = "idx"))
  stopifnot(max(abs(coef(fit) - coef(external))) < 1e-12)
  bread <- solve(fit$Ahat)
  manual <- bread %*% crossprod(fit$h11_hc3 + fit$h12) %*% t(bread) / nrow(w)^2
  stopifnot(max(abs(manual - fit$vcov)) < 1e-12)
}

z <- qnorm(.975)
fixture <- data.frame(n = 50L, tau = 2, type = "AUC", beta_name = "beta_A", beta_value = 0, rep = 1:3, estimate = c(-.1, 0, .1), se = .2, variance_method = CorrectedVarianceMethod())
fixture$lower <- fixture$estimate - z * fixture$se
fixture$upper <- fixture$estimate + z * fixture$se
fixture$pvalue <- 2 * pnorm(-abs(fixture$estimate / fixture$se))
ValidateCorrectedInferenceData(fixture, 3L)
config <- fixture[1, InferenceCellColumns(fixture), drop = FALSE]
ValidateResumedInference(fixture, config, 3L)
bad_config <- config; bad_config$tau <- 3
MustFail(ValidateResumedInference(fixture, bad_config, 3L))
bad_config <- config; bad_config$tau_gen <- 4
MustFail(ValidateResumedInference(fixture, bad_config, 3L))
seeded <- fixture; seeded$config_index <- 1L; seeded$seed <- vapply(seeded$rep, function(i) RegressionRepSeed(1L, 1L, i), integer(1))
config$config_index <- 1L
ValidateResumedInference(seeded, config, 3L, 1L)
MustFail(ValidateResumedInference(seeded, config, 3L, 2L))
bad <- fixture; bad$variance_method <- NULL
MustFail(ValidateCorrectedInferenceData(bad))
bad <- fixture; bad$variance_method[1] <- "unsupported"
MustFail(ValidateCorrectedInferenceData(bad))
bad <- fixture; bad$variance_method[1] <- NA
MustFail(ValidateCorrectedInferenceData(bad))
MustFail(ValidateCorrectedInferenceData(rbind(fixture, fixture[1, ])))
MustFail(ValidateCorrectedInferenceData(fixture[-1, ], 3L))
bad <- fixture; bad$lower[1] <- bad$lower[1] - .01
MustFail(ValidateCorrectedInferenceData(bad))
bad <- fixture; bad$se[1] <- 0
MustFail(ValidateCorrectedInferenceData(bad))
bad <- fixture; bad$pvalue[1] <- NA
MustFail(ValidateCorrectedInferenceData(bad))
for (study in c("estimation", "frailty", "power")) {
  full <- AttachCanonicalConfigIndex(CanonicalStudyConfigs(study), study)
  subset <- full[c(1, 5, nrow(full)), setdiff(names(full), "config_index")]
  stopifnot(identical(AttachCanonicalConfigIndex(subset, study)$config_index, c(1L, 5L, nrow(full))))
}

grid <- CanonicalStudyConfigs("estimation")
grid <- grid[grid$type == "AUC" & grid$beta_name == "beta_A", ]
raw <- grid[rep(seq_len(nrow(grid)), each = 3), ]
raw$rep <- rep(1:3, nrow(grid))
raw$term <- "A"
raw$estimate <- rep(c(-.1, 0, .1), nrow(grid))
raw$se <- .2
raw$lower <- raw$estimate - z * raw$se
raw$upper <- raw$estimate + z * raw$se
raw$pvalue <- 2 * pnorm(-abs(raw$estimate / raw$se))
raw$variance_method <- CorrectedVarianceMethod()
truth <- unique(grid[c("beta_name", "beta_value", "tau")])
truth$beta_A_AUMCF <- 0
s <- SummarizeSeCalibration(raw, truth, expected_reps = 3L)
stopifnot(nrow(s) == 27L, all(abs(s$empirical_sd - .1) < 1e-12), all(abs(s$rms_corrected_se - .2) < 1e-12), all(s$coverage == 1))
bad_truth <- rbind(truth, truth[1, ])
MustFail(SummarizeSeCalibration(raw, bad_truth, 3L))
MustFail(SummarizeSeCalibration(raw, truth[-1, ], 3L))
rows <- function(x) x[grepl("^[0-9]+ &", x)]
stopifnot(identical(rows(BuildSeCalibrationTable(s)), rows(BuildSeCalibrationTable(s, TRUE))))
stopifnot(!any(grepl("\\caption", BuildSeCalibrationTable(s, TRUE), fixed = TRUE)), !any(grepl("\\label", BuildSeCalibrationTable(s, TRUE), fixed = TRUE)), !any(grepl("{longtable}", BuildSeCalibrationTable(s, TRUE), fixed = TRUE)))
message("Corrected covariance, rejection guards, stable indices, truth joins, and table formatting: PASS.")
