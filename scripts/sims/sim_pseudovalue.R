#!/usr/bin/env Rscript
# Supports Figure S1: exact-jackknife versus influence-function pseudo-value validation.
argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)][1L])
source(file.path(dirname(dirname(normalizePath(script))), "utils", "bootstrap_regression.R"))
SetRegressionWd()
suppressPackageStartupMessages({library(MCC); library(dplyr); library(optparse)})
opts <- parse_args(OptionParser(option_list = list(
  make_option("--n", type = "character", default = "50,100,200,400"),
  make_option("--tau", type = "character", default = "1,2,3,4"),
  make_option("--type", type = "character", default = "AUC,MCF"),
  make_option("--seed", type = "integer", default = 3141L),
  make_option("--out", type = "character", default = "simulations/pseudovalue/pseudovalue_calibration.rds")
)))
SourceRegressionUtil("regression_options.R")
#' Calculate Area
#'
#' @param data Data.frame.
#' @param tau Truncation time.
#' @return Scalar AUC.
.CalcArea <- function(data, tau) {
  auc <- MCC::SingleArmAUC(data, tau = tau)
  out <- auc@MargAreas$area[1]
  return(out)
}


#' Calculate MCF at truncation time tau
#'
#' @param data Data.frame with idx, status, time.
#' @param tau Truncation time.
#' @return Scalar MCF(tau).
.CalcMCFAtTau <- function(data, tau) {
  mcf <- MCC::CalcMCF(
    idx = data$idx,
    status = data$status,
    time = data$time,
    calc_var = FALSE
  )
  tau_eff <- max(mcf$time[mcf$time <= tau], na.rm = TRUE)
  if (!is.finite(tau_eff)) {
    tau_eff <- min(mcf$time)
  }
  out <- mcf$mcf[mcf$time == tau_eff][1]
  return(out)
}


#' Calculate Jackknife Pseudovalues (AUC)
#'
#' @param data Data.frame.
#' @param tau Truncation time.
#' @return Data.frame of jackknife pseudo-values.
CalcPseudoJKAUC <- function(data, tau) {

  n <- length(unique(data$idx))
  auc_full <- .CalcArea(data, tau = tau)

  out <- lapply(seq_len(n), function(i) {

    df_i <- data %>% dplyr::filter(idx != i)
    auc_i <- .CalcArea(df_i, tau = tau)
    pseudo_i <- n * auc_full - (n - 1) * auc_i
    out <- data.frame(idx = i, pseudo_jk = pseudo_i)
    return(out)
  })
  out <- do.call(rbind, out)
  return(out)
}


#' Calculate Jackknife Pseudovalues (MCF at tau)
#'
#' @param data Data.frame.
#' @param tau Truncation time.
#' @return Data.frame of jackknife pseudo-values.
CalcPseudoJKMCF <- function(data, tau) {

  n <- length(unique(data$idx))
  mcf_full <- .CalcMCFAtTau(data, tau = tau)

  out <- lapply(seq_len(n), function(i) {

    df_i <- data %>% dplyr::filter(idx != i)
    mcf_i <- .CalcMCFAtTau(df_i, tau = tau)
    pseudo_i <- n * mcf_full - (n - 1) * mcf_i
    out <- data.frame(idx = i, pseudo_jk = pseudo_i)
    return(out)
  })
  out <- do.call(rbind, out)
  return(out)
}


#' Simulation Instance
#'
#' @param n Sample size.
#' @param tau Truncation time.
#' @param type "AUC" or "MCF".
#' @return Data.frame with idx, pseudo_inf, pseudo_jk, and type.
SimInst <- function(n, tau, type = c("AUC", "MCF")) {

  type <- match.arg(type)
  data <- MCC::GenData(n = n, tau = 10)

  df1 <- MCC::GenPseudo(data, tau = tau, type = type)
  df1 <- df1 %>%
    dplyr::select(idx, pseudo) %>%
    dplyr::rename(pseudo_inf = pseudo)

  if (type == "AUC") {
    df2 <- CalcPseudoJKAUC(data, tau = tau)
  } else {
    df2 <- CalcPseudoJKMCF(data, tau = tau)
  }

  out <- merge(x = df1, y = df2, by = "idx")
  out$type <- type
  return(out)
}


# Retain the fixed grid order, including n=400 and tau=1, to reproduce
# the RNG state of the plotted n=50/100/200, tau=2/3/4 cells. Only S1 is plotted.
set.seed(opts$seed)
rows <- lapply(ParseType(opts$type), function(typ) {
  bind_rows(lapply(ParseNumericOption(opts$n), function(n) {
    bind_rows(lapply(ParseNumericOption(opts$tau), function(tau) {
      out <- SimInst(n=n, tau=tau, type=typ)
      out$n <- n; out$tau <- tau; out
    }))
  }))
})
dir.create(dirname(opts$out), recursive=TRUE, showWarnings=FALSE)
saveRDS(bind_rows(rows), opts$out)
message("Wrote ", opts$out)
