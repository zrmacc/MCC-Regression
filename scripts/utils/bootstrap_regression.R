# Locate the Regression project root and load path utilities.
# Sourced by driver scripts via a short walk-up loop (see README.md).

if (!exists("RegressionDir", mode = "function")) {
  .__reg_argv <- commandArgs(trailingOnly = FALSE)
  .__reg_file <- .__reg_argv[grepl("^--file=", .__reg_argv)]
  .__reg_start <- if (length(.__reg_file) > 0L) {
    dirname(normalizePath(sub("^--file=", "", .__reg_file[1L]), mustWork = FALSE))
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
  .__reg_root <- NA_character_
  for (.__reg_step in 0:12L) {
    .__reg_cand <- if (.__reg_step == 0L) {
      .__reg_start
    } else {
      normalizePath(
        file.path(.__reg_start, paste(rep("..", .__reg_step), collapse = "/")),
        mustWork = FALSE
      )
    }
    if (file.exists(file.path(.__reg_cand, "scripts/utils/regression_paths.R"))) {
      .__reg_root <- .__reg_cand
      break
    }
  }
  if (is.na(.__reg_root)) {
    stop(
      "Could not find replication/scripts/utils/regression_paths.R. ",
      "Run scripts from the replication directory or pass an absolute script path.",
      call. = FALSE
    )
  }
  assign(".regression_root", .__reg_root, envir = .GlobalEnv)
  source(
    file.path(.__reg_root, "scripts/utils/regression_paths.R"),
    local = FALSE
  )
}

SourceRegressionUtil("corrected_inference.R")
