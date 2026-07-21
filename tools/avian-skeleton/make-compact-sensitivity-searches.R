#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

source_dir <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  Sys.getenv(
    "BIFROST_BODYPLAN_SEARCH_DIR",
    unset = file.path(
      "..",
      "passerine-bodyplan-evolution",
      "data",
      "temporal",
      "02_shift_search",
      "new_bifrost"
    )
  )
}

archive_search_dir <- file.path(
  "data",
  "temporal",
  "02_shift_search",
  "new_bifrost"
)

source_label <- paste(
  "Berv et al. (2026) Zenodo archive:",
  archive_search_dir
)

compact_provenance <- function(source_file = NULL) {
  provenance <- list(
    generated_on = as.character(Sys.Date()),
    source = source_label,
    source_dir = archive_search_dir,
    compact_script = "tools/avian-skeleton/make-compact-sensitivity-searches.R",
    notes = c(
      "Compact search objects keep final mapped trees, fitted BMM regime rates, shift nodes, IC summaries, compact IC history, VCVs, and IC weights.",
      "Individual candidate fits and proposal-level history are intentionally omitted.",
      "search$VCVs are proportional joint-model matrices; Part 5 post-hoc integration uses independent covariance refits instead."
    )
  )
  if (!is.null(source_file)) {
    provenance$source_file <- archive_source_file(source_file)
  }
  provenance
}

archive_source_file <- function(source_file) {
  file.path(
    "data",
    "temporal",
    "02_shift_search",
    "new_bifrost",
    basename(source_file)
  )
}

output_path <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  file.path(
    "inst",
    "extdata",
    "avian-skeleton",
    "passerine_bodyplan_search_sensitivity_compact.RDS"
  )
}

run_names <- c(
  "min30.ic40.gic", "min30.ic40.bic",
  "min20.ic40.gic", "min20.ic40.bic",
  "min10.ic40.gic", "min10.ic40.bic",
  "min30.ic20.gic", "min30.ic20.bic",
  "min20.ic20.gic", "min20.ic20.bic",
  "min10.ic20.gic", "min10.ic20.bic"
)

compact_model <- function(model) {
  out <- list(
    model = model$model,
    param = model$param
  )

  if (is.null(out$model) && !is.null(model$call$model)) {
    out$model <- as.character(model$call$model)
  }

  out
}

compact_search <- function(x, source_file) {
  out <- list(
    user_input = x$user_input,
    tree_no_uncertainty_untransformed = x$tree_no_uncertainty_untransformed,
    model_no_uncertainty = compact_model(x$model_no_uncertainty),
    shift_nodes_no_uncertainty = x$shift_nodes_no_uncertainty,
    optimal_ic = x$optimal_ic,
    baseline_ic = x$baseline_ic,
    IC_used = x$IC_used,
    num_candidates = x$num_candidates,
    model_fit_history = list(
      ic_acceptance_matrix = x$model_fit_history$ic_acceptance_matrix
    ),
    VCVs = x$VCVs,
    ic_weights = x$ic_weights
  )

  attr(out, "source_file") <- archive_source_file(source_file)
  attr(out, "provenance") <- compact_provenance(source_file)
  class(out) <- c("bifrost_search", "compact_bifrost_search", "list")
  out
}

if (!dir.exists(source_dir)) {
  stop("Source directory does not exist: ", source_dir, call. = FALSE)
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

searches <- vector("list", length(run_names))
names(searches) <- run_names

for (run_name in run_names) {
  source_file <- file.path(source_dir, paste0(run_name, ".RDS"))
  if (!file.exists(source_file)) {
    stop("Missing source run: ", source_file, call. = FALSE)
  }

  message("Compacting ", basename(source_file))
  full_search <- readRDS(source_file)
  searches[[run_name]] <- compact_search(full_search, source_file)
  rm(full_search)
  invisible(gc())
}

attr(searches, "source") <- source_label
attr(searches, "source_dir") <- archive_search_dir
attr(searches, "description") <- paste(
  "Compact Berv et al. (2026) sensitivity search outputs for",
  "rate-shift transition and magnitude-comparison examples."
)
attr(searches, "provenance") <- compact_provenance()
class(searches) <- c("bifrost_search_compact_bundle", "list")

saveRDS(searches, output_path, compress = "xz")

message("Wrote ", output_path)
message("Compressed size: ", file.info(output_path)$size, " bytes")
