#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L || anyNA(args) || any(!nzchar(args))) {
  stop(
    "Usage: Rscript tools/checktor-audit.R <package-path> <output-dir>",
    call. = FALSE
  )
}

package_path <- args[[1L]]
output_dir <- args[[2L]]

if (!dir.exists(package_path)) {
  stop("Package path does not exist or is not a directory: ", package_path,
       call. = FALSE)
}

package_path <- normalizePath(package_path, winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(package_path, "DESCRIPTION"))) {
  stop("Package path has no DESCRIPTION file: ", package_path, call. = FALSE)
}

if (!dir.exists(output_dir)) {
  stop("Output path does not exist or is not a directory: ", output_dir,
       call. = FALSE)
}

output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
if (file.access(output_dir, mode = 2L) != 0L) {
  stop("Output directory is not writable: ", output_dir, call. = FALSE)
}

if (!requireNamespace("checktor", quietly = TRUE)) {
  stop("The checktor package is required to run this audit.", call. = FALSE)
}
if (!requireNamespace("generics", quietly = TRUE)) {
  stop("The generics package is required to write the tidy audit summary.",
       call. = FALSE)
}

results <- checktor::checktor(
  package_path,
  verbose = FALSE,
  progress = FALSE
)

expected_categories <- c(
  "code_issues",
  "description_issues",
  "documentation_issues",
  "general_issues",
  "policy_issues"
)
missing_categories <- setdiff(expected_categories, names(results))
if (length(missing_categories) > 0L) {
  stop(
    "checktor returned no result for: ",
    paste(missing_categories, collapse = ", "),
    call. = FALSE
  )
}

bare_failures <- vapply(expected_categories, function(category_name) {
  category <- results[[category_name]]
  is.list(category) &&
    length(category$passed) == 1L &&
    isFALSE(unname(category$passed)) &&
    is.null(names(category$passed))
}, logical(1))
if (any(bare_failures)) {
  failure_messages <- vapply(
    expected_categories[bare_failures],
    function(category_name) {
      message <- results[[category_name]]$message
      if (is.null(message) || !nzchar(message)) {
        message <- "category diagnostic failed without a message"
      }
      paste0(category_name, ": ", message)
    },
    character(1)
  )
  stop(
    "checktor could not complete one or more diagnostic categories: ",
    paste(failure_messages, collapse = "; "),
    call. = FALSE
  )
}

issue_results <- checktor::issues(results)
tidy_results <- generics::tidy(results)
errored_checks <- unique(c(
  tidy_results$check[
    !is.na(tidy_results$message) &
      grepl("\\(errored\\)$", tidy_results$message)
  ],
  issue_results$check[
    (!is.na(issue_results$message) &
       grepl("\\(errored\\)$", issue_results$message)) |
      (!is.na(issue_results$location) &
         startsWith(issue_results$location, "Diagnostic errored:"))
  ]
))
if (length(errored_checks) > 0L) {
  stop(
    "checktor diagnostic execution failed for: ",
    paste(errored_checks, collapse = ", "),
    call. = FALSE
  )
}

health_path <- file.path(output_dir, "checktor-health.md")
issues_path <- file.path(output_dir, "checktor-issues.csv")
summary_path <- file.path(output_dir, "checktor-summary.csv")

invisible(checktor::health_report(
  results,
  file = health_path,
  format = "markdown"
))
utils::write.csv(issue_results, issues_path,
                 row.names = FALSE, na = "")
utils::write.csv(tidy_results, summary_path,
                 row.names = FALSE, na = "")

report_paths <- c(health_path, issues_path, summary_path)
report_sizes <- file.info(report_paths)$size
if (anyNA(report_sizes) || any(report_sizes <= 0L)) {
  stop("One or more checktor reports were not written successfully.",
       call. = FALSE)
}

advisory_summary <- c(
  "## checktor advisory audit",
  "",
  paste0(
    "**Advisory only:** checktor findings do not fail CI; audit or report ",
    "generation errors do."
  ),
  "",
  paste0("- Package: `", basename(package_path), "`"),
  paste0("- checktor version: `", results$metadata$checktor_version, "`"),
  paste0("- Findings: **", results$metadata$total_issues, "**"),
  paste0("- Failed checks: **", results$metadata$failed_checks, "**"),
  paste0(
    "- The Markdown health report and CSV reports are available in the ",
    "workflow artifact."
  )
)

cat(paste(advisory_summary, collapse = "\n"), "\n")

step_summary <- Sys.getenv("GITHUB_STEP_SUMMARY", unset = "")
if (nzchar(step_summary)) {
  cat(
    paste(c("", advisory_summary, ""), collapse = "\n"),
    file = step_summary,
    append = TRUE
  )
}
