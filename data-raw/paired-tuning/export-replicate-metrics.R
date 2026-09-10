# Export replicate-level metrics from the recorded paired campaign.
# No simulations or model fits are run. Source this file to use the builder,
# or run Rscript export-replicate-metrics.R <campaign-full-dir> <output.rds>.
build_replicate_metrics <- function(root) {
  expected_commit <- "db18184ddb06a5019123647ce218f9b717f76e49"
  expected_design <- "c6327e5bbc59b6d5d948651d3142068ce3d732848e52851957d4fc0844fc518e"
  installed <- utils::packageDescription("bifrost")
  stopifnot(identical(installed$RemoteSha, expected_commit))
  files <- sort(list.files(file.path(root, "replicates"), "[.]rds$",
                          recursive = TRUE, full.names = TRUE))
  stopifnot(length(files) == 1500L)
  rows <- vector("list", 18000L)
  hashes <- character(length(files))
  j <- 0L
  for (i in seq_along(files)) {
    w <- readRDS(files[i])
    stopifnot(identical(w$design_fingerprint, expected_design),
              identical(w$provenance$package_remote_sha, expected_commit),
              length(w$results) == 12L,
              identical(names(w$results), w$search_configs$config_id))
    sim <- w$simulation
    hash_input <- sim
    hash_input$user_input <- NULL
    stopifnot(identical(digest::digest(hash_input, algo = "sha256"),
                        w$simulation_hash))
    hashes[i] <- digest::digest(file = files[i], algo = "sha256")
    if (identical(w$group$scenario, "null")) {
      # Null generators store baseline_tree, not the shifted paintedTree field.
      # Adapt only this in-memory evaluator input; preserve the saved object.
      sim$shiftNodes <- integer(0)
      sim$paintedTree <- sim$baseline_tree
    }
    sizes <- vapply(sim$shiftNodes, function(node) {
      sum(phytools::getDescendants(sim$paintedTree, node) <= ape::Ntip(sim$paintedTree))
    }, integer(1))
    for (k in seq_len(12L)) {
      config <- w$search_configs[k, ]
      result <- w$results[[k]]
      stopifnot(!any(!is.na(result$error) & nzchar(as.character(result$error))),
                length(result$candidate_nodes) == result$num_candidates,
                !anyDuplicated(result$candidate_nodes),
                all(result$shift_nodes_no_uncertainty %in% result$candidate_nodes))
      e <- bifrost::evaluateShiftRecovery(list(sim), list(result),
                                          fuzzy_distance = 2L, weighted = TRUE,
                                          verbose = FALSE)
      stopifnot(e$n_evaluable_replicates == 1L)
      row <- data.frame(
        scenario = w$group$scenario, replicate = w$job$replicate,
        dataset_id = w$job$job_id, simulation_hash = w$simulation_hash,
        simulation_seed = w$seeds$simulation, search_seed = w$seeds$search,
        config_id = config$config_id, IC = config$IC,
        threshold = config$shift_acceptance_threshold,
        min_descendant_tips = config$min_descendant_tips,
        status = "ok", n_candidates = result$num_candidates,
        n_true_shifts = length(sim$shiftNodes),
        n_true_shifts_10_19 = sum(sizes >= 10L & sizes <= 19L),
        n_true_shifts_20_40 = sum(sizes >= 20L & sizes <= 40L),
        n_inferred_shifts = length(result$shift_nodes_no_uncertainty),
        stringsAsFactors = FALSE)
      row$null_false_positive_rate <- if (row$n_true_shifts == 0L && row$n_candidates > 0L) {
        row$n_inferred_shifts / row$n_candidates
      } else NA_real_
      weights <- result$ic_weights
      total_weight <- if (is.null(weights)) 0 else sum(
        weights$ic_weight_withshift[weights$node %in% result$shift_nodes_no_uncertainty],
        na.rm = TRUE)
      for (mode in c("strict", "fuzzy")) {
        for (name in names(e$counts[[mode]])) {
          row[[paste0(mode, "_", name)]] <- unname(e$counts[[mode]][name])
        }
        for (name in names(e[[mode]])) row[[paste0(mode, "_", name)]] <- e[[mode]][[name]]
        for (name in names(e$weighted[[mode]])) {
          row[[paste0("weighted_", mode, "_", name)]] <- e$weighted[[mode]][[name]]
        }
        tp <- if (row$n_true_shifts == 0L) 0 else e$weighted[[mode]]$recall * row$n_true_shifts
        row[[paste0("weighted_", mode, "_TP")]] <- tp
        row[[paste0("weighted_", mode, "_FP")]] <- total_weight - tp
      }
      j <- j + 1L
      rows[[j]] <- row
    }
    if (i %% 100L == 0L) message("Evaluated saved outputs: ", i, "/1500 datasets")
  }
  rows <- do.call(rbind, rows)
  rownames(rows) <- NULL
  stopifnot(nrow(rows) == 18000L, length(unique(rows$dataset_id)) == 1500L,
            all(table(rows$dataset_id) == 12L),
            !anyDuplicated(rows[c("dataset_id", "config_id")]))
  equal <- function(a, b) stopifnot(isTRUE(all.equal(unname(a), unname(b), tolerance = 1e-12)))
  divide <- function(a, b) if (b == 0) NA_real_ else a / b
  f1 <- function(p, r) if (is.na(p + r) || p + r == 0) NA_real_ else 2 * p * r / (p + r)
  # Reproduce published pooled metrics from sums of counts, NOT means of
  # per-replicate ratios. Null node-level FPR is the mean of replicate ratios.
  for (id in unique(rows$config_id)) {
    d <- rows[rows$config_id == id, ]
    summary_dir <- if (dir.exists(file.path(root, "summaries-available"))) {
      "summaries-available"
    } else "summaries"
    reference <- readRDS(file.path(root, summary_dir, paste0(id, ".rds")))
    stopifnot(identical(reference$design_fingerprint, expected_design), nrow(d) == 500L)
    old <- reference$per_replicate
    idx <- match(d$simulation_hash, old$simulation_hash)
    stopifnot(!anyNA(idx), !anyDuplicated(idx))
    equal(d$n_candidates, old$n_candidates[idx])
    equal(d$n_inferred_shifts, old$n_inferred_shifts[idx])
    equal(d$n_true_shifts_10_19, old$n_true_shifts_10_19[idx])
    equal(d$n_true_shifts_20_40, old$n_true_shifts_20_40[idx])
    if (d$scenario[1] == "null") {
      equal(mean(d$null_false_positive_rate), reference$summary_row$mean_false_positive_rate)
      equal(mean(d$n_inferred_shifts > 0), reference$summary_row$fraction_any_false_positive)
    } else for (mode in c("strict", "fuzzy")) {
      counts <- colSums(d[paste0(mode, "_", c("TP", "FP", "FN", "TN"))])
      equal(counts, reference$evaluation$counts[[mode]])
      tp <- counts[1]; fp <- counts[2]; fn <- counts[3]; tn <- counts[4]
      p <- divide(tp, tp + fp); r <- divide(tp, tp + fn); s <- divide(tn, tn + fp)
      metrics <- list(precision = p, recall = r, f1 = f1(p, r), specificity = s,
                      fpr = divide(fp, fp + tn), balanced_accuracy = (r + s) / 2)
      for (name in names(metrics)) equal(metrics[[name]], reference$evaluation[[mode]][[name]])
      wtp <- sum(d[[paste0("weighted_", mode, "_TP")]])
      wfp <- sum(d[[paste0("weighted_", mode, "_FP")]])
      wp <- divide(wtp, wtp + wfp); wr <- divide(wtp, tp + fn)
      for (name in c("precision", "recall", "f1")) {
        equal(list(precision = wp, recall = wr, f1 = f1(wp, wr))[[name]],
              reference$evaluation$weighted[[mode]][[name]])
      }
    }
  }
  message("Validated all 36 scenario/settings summaries against saved campaign summaries.")
  list(schema_version = 1L, provenance = list(
    package_commit = expected_commit, design_fingerprint = expected_design,
    metric_accounting_version = "candidate-node-aware-v1", fuzzy_distance = 2L,
    replicate_files_sha256 = stats::setNames(hashes, basename(files)),
    aggregation = "Pool strict/fuzzy counts and weighted TP/FP before computing recovery metrics; average replicate FPR for null scenarios.",
    missing_metrics = "Undefined ratios retain NA, following evaluateShiftRecovery()."
  ), metrics = rows)
}


if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop("Usage: Rscript export-replicate-metrics.R <campaign-full-dir> <output.rds>")
  }
  if (file.exists(args[2])) stop("Output already exists: ", args[2])
  result <- build_replicate_metrics(args[1])
  dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("replicate-metrics-", tmpdir = dirname(args[2]))
  tryCatch({
    saveRDS(result, temporary, compress = "xz")
    stopifnot(identical(result, readRDS(temporary)))
    if (!file.rename(temporary, args[2])) stop("Could not publish output: ", args[2])
  }, finally = if (file.exists(temporary)) unlink(temporary))
  message("Wrote ", args[2], " (", file.info(args[2])$size, " bytes)")
}
