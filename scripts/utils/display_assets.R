# Active manuscript figures. Historical script/output names remain stable.
DisplayAssetMap <- function() {
  data.frame(
    document = c(rep("main", 7), rep("supplement", 6)),
    figure = c(1, 1, 2, 2, 3, 4, 5, 1, 1, 2, 2, 3, 3),
    label = c("fig:bias", "fig:bias", "fig:coverage", "fig:coverage", "fig:t1e", "fig:power", "fig:application", "supp:fig:frailty", "supp:fig:frailty", "supp:fig:terminal-event-beta", "supp:fig:terminal-event-beta", "supp:fig:terminal-event-gamma", "supp:fig:terminal-event-gamma"),
    script = c("plot_figure_1.R", "plot_figure_1.R", "plot_figure_2.R", "plot_figure_2.R", "plot_figure_3.R", "plot_figure_4.R", NA, "plot_figure_S2.R", "plot_figure_S2.R", "plot_figure_S3.R", "plot_figure_S3.R", "plot_figure_S4.R", "plot_figure_S4.R"),
    source = c("fig1a.pdf", "fig1b.pdf", "fig2a.pdf", "fig2b.pdf", "fig3.pdf", "fig4.pdf", NA, "figS2a.pdf", "figS2b.pdf", "figS3a.pdf", "figS3b.pdf", "figS4a.pdf", "figS4b.pdf"),
    manuscript = c("fig2a.pdf", "fig2b.pdf", "fig3a.pdf", "fig3b.pdf", "fig5.pdf", "fig6.pdf", "fig7.pdf", "fig4a.pdf", "fig4b.pdf", "supp_fig_terminal_event_beta_a.pdf", "supp_fig_terminal_event_beta_b.pdf", "supp_fig_terminal_event_gamma_a.pdf", "supp_fig_terminal_event_gamma_b.pdf"),
    stringsAsFactors = FALSE
  )
}

# Compare the active LaTeX references with the supplied/generated figure assets.
# Commented figures are excluded, and the collaborator-supplied clinical figure
# must be referenced but is never overwritten by the replication driver.
ValidateManuscriptFigures <- function(manuscript_dir, figure_dir, check_copies = TRUE) {
  assets <- DisplayAssetMap()
  figure_pattern <- "(?s)\\\\begin\\{figure\\*?\\}.*?\\\\end\\{figure\\*?\\}"
  graphics_pattern <- "\\\\includegraphics(?:\\[[^]]*\\])?\\{[^}]+\\}"
  extract_graphics <- function(text) {
    hits <- regmatches(text, gregexpr(graphics_pattern, text, perl = TRUE))[[1]]
    if (!length(hits)) return(character())
    basename(sub(".*\\{([^}]+)\\}$", "\\1", hits))
  }
  for (document in unique(assets$document)) {
    path <- file.path(manuscript_dir, paste0(document, ".tex"))
    if (!file.exists(path)) stop("Missing manuscript source: ", path, call. = FALSE)
    lines <- sub("(?<!\\\\)%.*$", "", readLines(path, warn = FALSE), perl = TRUE)
    text <- paste(lines, collapse = "\n")
    blocks <- regmatches(text, gregexpr(figure_pattern, text, perl = TRUE))[[1]]
    # Continued panels form one numbered figure even when they occupy separate pages.
    grouped <- character()
    for (block in blocks) {
      if (grepl("\\ContinuedFloat", block, fixed = TRUE)) {
        if (!length(grouped)) stop("A continued panel has no preceding figure.", call. = FALSE)
        grouped[length(grouped)] <- paste(grouped[length(grouped)], block, sep = "\n")
      } else grouped <- c(grouped, block)
    }
    blocks <- grouped
    expected <- assets[assets$document == document, , drop = FALSE]
    referenced <- extract_graphics(text)
    if (!setequal(referenced, expected$manuscript) || anyDuplicated(referenced)) {
      stop("Active figure references differ from DisplayAssetMap() in ", path,
           "; missing: ", paste(setdiff(expected$manuscript, referenced), collapse = ", "),
           "; unexpected: ", paste(setdiff(referenced, expected$manuscript), collapse = ", "), call. = FALSE)
    }
    for (i in seq_len(nrow(expected))) {
      row <- expected[i, ]
      matching <- which(vapply(blocks, function(b) grepl(paste0("\\label{", row$label, "}"), b, fixed = TRUE), logical(1)))
      if (length(matching) != 1L || matching != row$figure || !row$manuscript %in% extract_graphics(blocks[matching])) {
        stop("Figure number, label, or asset mismatch for ", document, " Figure ", row$figure, ": ", row$label, call. = FALSE)
      }
      target <- file.path(manuscript_dir, "figures", row$manuscript)
      if ((check_copies || is.na(row$source)) && !file.exists(target)) stop("Missing manuscript figure: ", target, call. = FALSE)
      if (check_copies && !is.na(row$source)) {
        source <- file.path(figure_dir, row$source)
        if (!file.exists(source)) stop("Missing replication figure: ", source, call. = FALSE)
        if (unname(tools::md5sum(source)) != unname(tools::md5sum(target))) {
          stop("Replication and manuscript figure files differ: ", row$source, " -> ", row$manuscript, call. = FALSE)
        }
      }
    }
  }
  invisible(assets)
}
