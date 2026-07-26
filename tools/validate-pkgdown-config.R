#!/usr/bin/env Rscript

main <- function(args) {
  if (length(args) > 1L) {
    stop("Usage: validate-pkgdown-config.R [CONFIG]", call. = FALSE)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to validate pkgdown configuration.",
         call. = FALSE)
  }

  config_path <- if (length(args) == 1L) args[[1L]] else "_pkgdown.yml"
  config <- tryCatch(
    yaml::read_yaml(
      config_path,
      handlers = list(
        seq = function(value) structure(value, class = "yaml_sequence")
      )
    ),
    error = function(error) {
      stop("Could not parse pkgdown configuration: ", conditionMessage(error),
           call. = FALSE)
    }
  )
  if (!is.list(config)) {
    stop("pkgdown configuration must define a template mapping.", call. = FALSE)
  }
  template <- config[["template"]]
  if (!is.list(template)) {
    stop("pkgdown configuration must define a template mapping.", call. = FALSE)
  }
  if (is.null(template[["math-rendering"]])) {
    stop("pkgdown template must define math-rendering.", call. = FALSE)
  }
  if (!identical(template[["math-rendering"]], "katex")) {
    stop("pkgdown template math-rendering must be exactly 'katex'.", call. = FALSE)
  }
}

main(commandArgs(trailingOnly = TRUE))
