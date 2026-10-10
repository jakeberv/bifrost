#!/usr/bin/env Rscript

# Run from the repository root; no searches or model refits are performed.
# Rscript tools/render-jaw-shape-figures.R SWEEP_CACHE HISTORY_RESULT
# The sweep cache lacks history. Supply a cached reproduction that includes it.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: render-jaw-shape-figures.R SWEEP_CACHE HISTORY_RESULT")
}
pkgload::load_all(".", quiet = TRUE)

cache_path <- args[[1L]]
history_path <- args[[2L]]
cache_key <- "GIC_threshold_20"
jaws_search <- readRDS(cache_path)[[cache_key]]
history_search <- readRDS(history_path)
trajectory <- icTrajectory(history_search)
accepted <- trajectory$status == "accepted"
near <- function(x, y) isTRUE(all.equal(x, y, tolerance = 1e-10))

# Refuse to combine the fitted model with an unrelated search trajectory.
stopifnot(
  length(jaws_search$shift_nodes_no_uncertainty) == 17L,
  identical(jaws_search$shift_nodes_no_uncertainty,
            history_search$shift_nodes_no_uncertainty),
  identical(as.integer(trajectory$candidate_node[accepted]),
            as.integer(jaws_search$shift_nodes_no_uncertainty)),
  nrow(trajectory) == 29L,
  sum(accepted) == 17L,
  sum(trajectory$status == "rejected") == 11L,
  all(trajectory$accepted[-1L] == (trajectory$delta_ic[-1L] >= 20)),
  abs(jaws_search$baseline_ic - history_search$baseline_ic) < 1e-7,
  abs(jaws_search$optimal_ic - history_search$optimal_ic) < 1e-7,
  abs(tail(trajectory$best_ic, 1L) - jaws_search$optimal_ic) < 1e-7,
  near(jaws_search$model_no_uncertainty$param,
       history_search$model_no_uncertainty$param),
  near(jaws_search$tree_no_uncertainty_untransformed,
       history_search$tree_no_uncertainty_untransformed),
  identical(names(jaws_search$VCVs), names(history_search$VCVs)),
  all(vapply(names(jaws_search$VCVs), function(nm) {
    near(unname(jaws_search$VCVs[[nm]]), unname(history_search$VCVs[[nm]]))
  }, logical(1L)))
)

# Also check the published numerical summary before replacing any figures.
sweep <- read.csv("vignettes/rate-map-jaw-shape/sweep_summary.csv")
reference <- sweep[sweep$threshold == 20, ]
stopifnot(
  nrow(reference) == 1L,
  reference$n_shifts == 17L,
  abs(reference$baseline_ic - jaws_search$baseline_ic) < 1e-7,
  abs(reference$optimal_ic - jaws_search$optimal_ic) < 1e-7
)

# Use the vignette's plotting code itself so the example and figures agree.
vignette_path <- "vignettes/jaw-shape-vignette.Rmd"
lines <- readLines(vignette_path)
chunk_code <- function(label) {
  start <- which(startsWith(lines, paste0("```{r ", label, ",")))
  stopifnot(length(start) == 1L)
  end <- which(seq_along(lines) > start & lines == "```")[1L]
  parse(text = lines[seq.int(start + 1L, end - 1L)])
}
jaws_search$model_fit_history <- history_search$model_fit_history
rm(history_search)
plot_env <- new.env(parent = globalenv())
plot_env$jaws_search <- jaws_search
plot_env$plotSimmap <- phytools::plotSimmap
out <- "vignettes/jaw-shape"
render_figure <- function(label, filename, width, height) {
  png(file.path(out, filename), width = width, height = height,
      res = 140, type = "cairo", bg = "white")
  on.exit(dev.off(), add = TRUE)
  eval(chunk_code(label), envir = plot_env)
}
render_figure("plot_cor", "cor_heatmap.png", 1100, 1100)
render_figure("plot_decay", "IC_decay.png", 1540, 951)
render_figure("plot_tree", "branch_rates.png", 1540, 951)

# Retain a small, inspectable numerical record without bundling large fits.
write.csv(as.data.frame(trajectory), file.path(out, "IC_trajectory.csv"),
          row.names = FALSE, na = "")
hash <- function(path) digest::digest(file = path, algo = "sha256")
artifacts <- file.path(out, c("cor_heatmap.png", "IC_decay.png",
                             "branch_rates.png", "IC_trajectory.csv"))
provenance <- list(
  regeneration_command = paste("Rscript tools/render-jaw-shape-figures.R",
                               shQuote(cache_path), shQuote(history_path)),
  cache = list(path = cache_path, sha256 = hash(cache_path), key = cache_key),
  history = list(path = history_path, sha256 = hash(history_path),
                 description = paste(
                   "Previously cached reproduction; the sweep cache omits history.",
                   "Accepted nodes, rates, tree, and VCV values checked against",
                   "the sweep cache; baseline and final GIC agree within 1e-7."
                 )),
  summary = list(n_proposals = 28L, n_accepted = 17L, n_rejected = 11L,
                 shift_nodes = jaws_search$shift_nodes_no_uncertainty,
                 baseline_ic = jaws_search$baseline_ic,
                 optimal_ic = jaws_search$optimal_ic,
                 delta_ic = jaws_search$baseline_ic - jaws_search$optimal_ic),
  plot_source = vignette_path,
  plot_chunks = c("plot_cor", "plot_decay", "plot_tree"),
  artifacts_sha256 = as.list(setNames(vapply(artifacts, hash, character(1L)),
                                     basename(artifacts)))
)
jsonlite::write_json(provenance, file.path(out, "figure-provenance.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = NA)
message("Regenerated Figures 3–5: 17 accepted / 28 proposals; final GIC ",
        format(jaws_search$optimal_ic, digits = 12))
