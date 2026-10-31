# Path helpers for the Regression simulation subproject.
# Loaded automatically by R scripts via the bootstrap block at the top of each
# driver script (see README.md).

#' Path to the running Rscript file, if available.
#'
#' @return Character scalar or \code{NULL}.
ScriptPath <- function() {
  argv <- commandArgs(trailingOnly = FALSE)
  file_arg <- argv[grepl("^--file=", argv)]
  if (length(file_arg) > 0L) {
    out <- sub("^--file=", "", file_arg[1L])
  } else {
    out <- NULL
  }
  return(out)
}

#' Walk upward from \code{start_dir} to locate the Regression project root.
#'
#' The root is identified by \code{README.md}, \code{scripts/}, and
#' \code{scripts/utils/regression_paths.R} in the same directory.
#'
#' @param start_dir Directory to begin the search (default: script dir or cwd).
#' @return Normalized path to the Regression directory.
FindRegressionRoot <- function(start_dir = NULL) {
  if (is.null(start_dir)) {
    sp <- ScriptPath()
    if (is.null(sp)) {
      start_dir <- getwd()
    } else {
      start_dir <- dirname(normalizePath(sp, mustWork = FALSE))
    }
  }
  d <- normalizePath(start_dir, mustWork = FALSE)
  for (step in 0:12L) {
    cand <- if (step == 0L) {
      d
    } else {
      normalizePath(
        file.path(d, paste(rep("..", step), collapse = "/")),
        mustWork = FALSE
      )
    }
    if (
      file.exists(file.path(cand, "README.md")) &&
      dir.exists(file.path(cand, "scripts")) &&
      file.exists(file.path(cand, "scripts", "utils", "regression_paths.R"))
    ) {
      return(cand)
    }
  }
  stop(
    "Could not locate the Regression project root starting from: ",
    d,
    ". Run from the replication directory or use an absolute path to the script.",
    call. = FALSE
  )
}

#' Regression project root directory.
#'
#' Uses \code{.regression_root} when set by the bootstrap block; otherwise
#' calls \code{FindRegressionRoot()}.
#'
#' @return Character scalar path.
RegressionDir <- function() {
  if (exists(".regression_root", envir = .GlobalEnv, inherits = FALSE)) {
    out <- get(".regression_root", envir = .GlobalEnv)
  } else {
    out <- FindRegressionRoot()
  }
  return(out)
}

#' Set working directory to \code{RegressionDir()}.
#'
#' @return Invisibly, the Regression root path.
SetRegressionWd <- function() {
  reg <- RegressionDir()
  if (normalizePath(getwd(), mustWork = FALSE) != normalizePath(reg, mustWork = FALSE)) {
    setwd(reg)
    message("Working directory set to: ", reg)
  }
  return(invisible(reg))
}

#' Source a file under \code{scripts/utils/}.
#'
#' @param filename File name only (e.g. \code{"regression_dgp.R"}).
#' @return Invisibly, the sourced path.
SourceRegressionUtil <- function(filename) {
  path <- file.path(RegressionDir(), "scripts", "utils", filename)
  if (!file.exists(path)) {
    stop("Missing Regression utility: ", path, call. = FALSE)
  }
  source(path, local = FALSE)
  return(invisible(path))
}

# Resolve output paths relative to this replication checkout.
ResolveRegressionPath <- function(path) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) return(path)
  file.path(RegressionDir(), path)
}

# All included figure artifacts are kept directly under figures/.
FigureDir <- function(topic = NULL) {
  out <- file.path(RegressionDir(), "figures")
  return(out)
}

#' Save a ggplot as PDF, optionally also as PNG.
#'
#' By default only PDF is written. Set \code{also_png = TRUE} for dense plots
#' such as the type-I error QQ figure.
#'
#' @param plot ggplot object.
#' @param path_pdf Destination PDF path.
#' @param width Width in inches.
#' @param height Height in inches.
#' @param also_png If TRUE, also write a PNG sibling of \code{path_pdf}.
#' @param dpi PNG resolution when \code{also_png} is TRUE.
#' @return Invisibly, the PDF path.
SaveRegressionFigure <- function(
  plot,
  path_pdf,
  width,
  height,
  also_png = FALSE,
  dpi = 300
) {
  dir.create(
    dirname(path_pdf),
    recursive = TRUE,
    showWarnings = FALSE
  )
  ggplot2::ggsave(
    filename = path_pdf,
    plot = plot,
    width = width,
    height = height,
    device = "pdf"
  )
  message("Wrote ", path_pdf)
  if (isTRUE(also_png)) {
    path_png <- sub("\\.pdf$", ".png", path_pdf)
    ggplot2::ggsave(
      filename = path_png,
      plot = plot,
      width = width,
      height = height,
      device = "png",
      dpi = dpi
    )
    message("Wrote ", path_png)
  }
  return(invisible(path_pdf))
}
