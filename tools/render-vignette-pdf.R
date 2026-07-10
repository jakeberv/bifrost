#!/usr/bin/env Rscript

find_repo_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(path, "_pkgdown.yml"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find repository root from ", getwd(), call. = FALSE)
    }
    path <- parent
  }
}

opt_value <- function(args, name, default = NULL) {
  hit <- which(args == name)
  if (length(hit) == 0L) {
    return(default)
  }
  if (hit[1L] == length(args)) {
    stop("Missing value for ", name, call. = FALSE)
  }
  args[[hit[1L] + 1L]]
}

main <- function(args) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required to render vignette PDFs.", call. = FALSE)
  }

  output_dir <- opt_value(args, "--output-dir", "docs/articles")
  skip <- logical(length(args))
  option_hits <- which(args == "--output-dir")
  skip[option_hits] <- TRUE
  skip[option_hits + 1L] <- TRUE
  slugs <- args[!skip]

  if (length(slugs) == 0L) {
    stop("Usage: render-vignette-pdf.R <slug> [<slug> ...] [--output-dir DIR]",
         call. = FALSE)
  }

  repo_root <- find_repo_root()
  setwd(repo_root)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (slug in slugs) {
    input <- file.path(repo_root, "vignettes", paste0(slug, ".Rmd"))
    if (!file.exists(input)) {
      stop("Missing vignette source for slug '", slug, "': ", input, call. = FALSE)
    }
    message("Rendering PDF for ", slug)
    rmarkdown::render(
      input = input,
      output_format = "rmarkdown::pdf_document",
      output_file = paste0(slug, ".pdf"),
      output_dir = output_dir,
      clean = TRUE,
      quiet = FALSE,
      envir = new.env(parent = globalenv())
    )
  }
}

main(commandArgs(trailingOnly = TRUE))
