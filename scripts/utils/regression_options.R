# Shared command-line parsing helpers for Regression simulation scripts.

#' Parse a comma-separated pseudo-value type option.
#'
#' @param type_string Comma-separated types (e.g. \code{"AUC,MCF"}).
#' @return Character vector with values in \code{c("AUC", "MCF")}.
ParseType <- function(type_string) {
  out <- trimws(strsplit(type_string, ",")[[1L]])
  out <- toupper(out)
  out <- out[out %in% c("AUC", "MCF")]
  if (length(out) == 0L) {
    stop("--type must include at least one of: AUC, MCF", call. = FALSE)
  }
  return(out)
}

#' Parse a comma-separated beta_name option.
#'
#' @param beta_name_string Comma-separated names.
#' @return Character vector with values in \code{c("beta_A", "beta_X", "beta_AX")}.
ParseBetaName <- function(beta_name_string) {
  out <- trimws(strsplit(beta_name_string, ",")[[1L]])
  out <- out[out %in% c("beta_A", "beta_X", "beta_AX")]
  if (length(out) == 0L) {
    stop(
      "--beta_name must include at least one of: beta_A, beta_X, beta_AX",
      call. = FALSE
    )
  }
  return(out)
}

#' Parse a comma-separated numeric option.
#'
#' @param value_string Comma-separated numbers.
#' @return Numeric vector.
ParseNumericOption <- function(value_string) {
  out <- as.numeric(trimws(strsplit(value_string, ",")[[1L]]))
  if (anyNA(out)) {
    stop("Failed to parse a numeric option value.", call. = FALSE)
  }
  return(out)
}

#' Parse a comma-separated integer option.
#'
#' @param value_string Comma-separated integers.
#' @return Integer vector.
ParseIntegerOption <- function(value_string) {
  out <- as.integer(trimws(strsplit(value_string, ",")[[1L]]))
  if (anyNA(out)) {
    stop("Failed to parse an integer option value.", call. = FALSE)
  }
  return(out)
}
