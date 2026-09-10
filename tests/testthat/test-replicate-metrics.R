test_that("replicate metrics preserve pairing and reconstruct the vignette cache", {
  root <- Sys.getenv("BIFROST_ARTIFACT_DIR", unset = "")
  skip_if(!nzchar(root) || !dir.exists(root), "repository artifacts unavailable")
  path <- file.path(root, "simulation-study-cache/passerine_replicate_metrics.rds")
  x <- readRDS(path)
  cache <- readRDS(file.path(root, "simulation-study-cache/passerine_preview_tables.rds"))
  d <- x$metrics
  expect_identical(x$schema_version, 1L)
  expect_identical(x$provenance$package_commit, cache$provenance$tuning$package_commit)
  expect_identical(x$provenance$design_fingerprint,
                   cache$provenance$tuning$paired_design$design_fingerprint)
  expect_identical(x$provenance$metric_accounting_version, "candidate-node-aware-v1")
  expect_identical(x$provenance$fuzzy_distance, 2L)
  expect_length(x$provenance$replicate_files_sha256, 1500L)
  expect_false(anyDuplicated(names(x$provenance$replicate_files_sha256)) > 0)
  expect_true(all(grepl("^[a-f0-9]{64}$", x$provenance$replicate_files_sha256)))
  expect_true(all(d$status == "ok"))
  expect_identical(ncol(d), 47L)
  expect_false(anyDuplicated(d[c("dataset_id", "config_id")]) > 0)
  stopifnot(nrow(d) == 18000L, length(unique(d$dataset_id)) == 1500L,
            all(table(d$scenario) == 6000L), all(table(d$config_id) == 500L))
  for (id in unique(d$dataset_id)) {
    z <- d[d$dataset_id == id, ]
    stopifnot(nrow(z) == 12L, length(unique(z$simulation_hash)) == 1L,
              length(unique(z$simulation_seed)) == 1L,
              length(unique(z$search_seed)) == 1L)
    stopifnot(setequal(paste(z$IC, z$threshold, z$min_descendant_tips),
                      with(expand.grid(IC = c("GIC", "BIC"), threshold = c(10, 20, 30),
                                       minimum = c(10, 20)), paste(IC, threshold, minimum))))
  }
  # Successful zero-shift searches remain present; misses are not dropped.
  zero_shifted <- d$n_true_shifts > 0 & d$n_inferred_shifts == 0
  expect_true(any(zero_shifted))
  expect_true(all(d$strict_FN[zero_shifted] == d$n_true_shifts[zero_shifted]))
  expect_true(all(d$fuzzy_FN[zero_shifted] == d$n_true_shifts[zero_shifted]))
  expect_true(all(d$strict_recall[zero_shifted] == 0))
  for (mode in c("strict", "fuzzy")) {
    counts <- d[paste0(mode, "_", c("TP", "FP", "FN", "TN"))]
    expect_true(all(as.matrix(counts) >= 0))
    expect_true(all(as.matrix(counts) == floor(as.matrix(counts))))
    expect_equal(counts[[1]] + counts[[2]], d$n_inferred_shifts)
    expect_equal(counts[[1]] + counts[[3]], d$n_true_shifts)
    expect_true(all(counts[[2]] + counts[[4]] <= d$n_candidates))
    expect_true(all(is.na(d[[paste0(mode, "_recall")]][d$n_true_shifts == 0])))
    divide <- function(a, b) ifelse(b == 0, NA_real_, a / b)
    harmonic <- function(p, r) ifelse(is.na(p + r) | p + r == 0, NA_real_, 2 * p * r / (p + r))
    p <- divide(counts[[1]], counts[[1]] + counts[[2]])
    r <- divide(counts[[1]], counts[[1]] + counts[[3]])
    s <- divide(counts[[4]], counts[[4]] + counts[[2]])
    expected <- list(precision = p, recall = r, f1 = harmonic(p, r),
                     specificity = s, fpr = divide(counts[[2]], counts[[2]] + counts[[4]]),
                     balanced_accuracy = (r + s) / 2)
    for (name in names(expected)) {
      expect_equal(d[[paste(mode, name, sep = "_")]], expected[[name]], tolerance = 1e-12)
    }
    tp <- d[[paste0("weighted_", mode, "_TP")]]
    fp <- d[[paste0("weighted_", mode, "_FP")]]
    # Subtraction can leave FP at machine-roundoff below zero.
    expect_true(all(tp >= -1e-12 & fp >= -1e-12))
    wp <- divide(tp, tp + fp); wr <- divide(tp, d$n_true_shifts)
    for (name in c("precision", "recall", "f1")) {
      expect_equal(d[[paste("weighted", mode, name, sep = "_")]],
                   list(precision = wp, recall = wr, f1 = harmonic(wp, wr))[[name]],
                   tolerance = 1e-12)
    }
  }
  equal <- function(a, b) {
    expect_equal(unname(a), unname(b), tolerance = 1e-12)
  }
  ratio <- function(a, b) if (b == 0) NA_real_ else a / b
  f1 <- function(p, r) if (is.na(p + r) || p + r == 0) NA_real_ else 2 * p * r / (p + r)
  for (ic in c("GIC", "BIC")) {
    grid <- cache$grid_summary[[tolower(ic)]]
    for (i in seq_len(nrow(grid))) {
      g <- grid[i, ]
      z <- d[d$IC == ic & d$threshold == g$shift_acceptance_threshold &
               d$min_descendant_tips == g$min_descendant_tips, ]
      null <- z[z$scenario == "null", ]
      equal(mean(null$null_false_positive_rate), g$null_mean_false_positive_rate)
      equal(mean(null$n_inferred_shifts > 0), g$null_fraction_any_false_positive)
      for (scenario in c("null", "proportional", "integration-rate")) {
        a <- z[z$scenario == scenario, ]
        prefix <- if (scenario == "integration-rate") "correlation" else scenario
        equal(mean(a$n_inferred_shifts), g[[paste0(prefix, "_mean_inferred_shifts")]])
        if (scenario == "null") next
        for (mode in c("strict", "fuzzy")) {
          counts <- colSums(a[paste0(mode, "_", c("TP", "FP", "FN", "TN"))])
          tp <- counts[1]; fp <- counts[2]; fn <- counts[3]; tn <- counts[4]
          p <- ratio(tp, tp + fp); r <- ratio(tp, tp + fn); s <- ratio(tn, tn + fp)
          metrics <- list(precision = p, recall = r, f1 = f1(p, r), specificity = s,
                          fpr = ratio(fp, fp + tn), balanced_accuracy = (r + s) / 2)
          for (name in names(metrics)) {
            equal(metrics[[name]], g[[paste(prefix, mode, name, sep = "_")]])
          }
        }
        tp <- sum(a$weighted_fuzzy_TP); fp <- sum(a$weighted_fuzzy_FP)
        equal(f1(ratio(tp, tp + fp), ratio(tp, sum(a$n_true_shifts))),
              g[[paste0(prefix, "_weighted_fuzzy_f1")]])
      }
    }
  }
})

test_that("the portable metrics exporter refuses invalid CLI calls without overwriting", {
  exporter <- test_path("../../data-raw/paired-tuning/export-replicate-metrics.R")
  skip_if_not(file.exists(exporter), "development exporter unavailable in package build")
  namespace <- new.env(parent = globalenv())
  sys.source(exporter, namespace)
  expect_true(is.function(namespace$build_replicate_metrics))
  rscript <- file.path(R.home("bin"), "Rscript")
  if (.Platform$OS.type == "windows") rscript <- paste0(rscript, ".exe")
  usage <- suppressWarnings(system2(rscript, shQuote(exporter), stdout = TRUE, stderr = TRUE))
  expect_identical(attr(usage, "status"), 1L)
  expect_match(paste(usage, collapse = "\n"), "Usage:")
  sentinel <- withr::local_tempfile(fileext = ".rds")
  saveRDS(list(preserve = TRUE), sentinel)
  before <- digest::digest(file = sentinel, algo = "sha256")
  refused <- suppressWarnings(system2(rscript,
    c(shQuote(exporter), shQuote("missing-campaign"), shQuote(sentinel)),
    stdout = TRUE, stderr = TRUE))
  expect_identical(attr(refused, "status"), 1L)
  expect_match(paste(refused, collapse = "\n"), "Output already exists")
  expect_identical(digest::digest(file = sentinel, algo = "sha256"), before)
})
