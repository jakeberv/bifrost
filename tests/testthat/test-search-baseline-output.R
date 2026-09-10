testthat::skip_on_cran()

baseline_output_inputs <- function() {
  tree <- ape::stree(12L, type = "left")
  tree$edge.length <- rep(1, nrow(tree$edge))
  traits <- withr::with_seed(221, matrix(rnorm(24), nrow = 12L))
  dimnames(traits) <- list(tree$tip.label, c("trait_a", "trait_b"))
  list(baseline_tree = tree, trait_data = traits)
}

expect_baseline_output <- function(result) {
  expect_length(result$shift_nodes_no_uncertainty, 0L)
  expect_identical(result$model_no_uncertainty$model, "BM")
  expect_true(is.na(result$model_no_uncertainty$param))
  expect_null(names(result$model_no_uncertainty$param))
  expect_identical(names(result$VCVs), "0")
  expect_identical(result$VCVs[["0"]], result$model_no_uncertainty$sigma$Pinv)
  expect_true(is.matrix(result$VCVs[["0"]]))
  trees <- list(
    result$tree_no_uncertainty_transformed,
    result$tree_no_uncertainty_untransformed,
    result$model_no_uncertainty$corrSt$phy,
    result$model_no_uncertainty$variables$tree
  )
  for (tree in trees) {
    expect_identical(colnames(tree$mapped.edge), "0")
    expect_true(all(unlist(lapply(tree$maps, names)) == "0"))
    expect_identical(unique(phytools::getStates(tree, type = "both")), "0")
  }
}

test_that("real baseline fits preserve numerics and expose the global covariance", {
  inputs <- baseline_output_inputs()
  baseline <- .bifrost_search_initialize_tree(inputs$baseline_tree)
  old_baseline <- generatePaintedTrees(baseline, min_tips = 12L)[[1L]]
  cases <- expand.grid(IC = c("GIC", "BIC"), cutoff = c(10L, 12L),
                       stringsAsFactors = FALSE)
  cases$method <- "LL"
  cases <- rbind(cases, data.frame(IC = "GIC", cutoff = 12L, method = "LOOCV"))

  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    reference <- withr::with_seed(221, .bifrost_search_fit_ic(
      case$IC, "trait_data ~ 1", old_baseline, inputs$trait_data,
      method = case$method, error = FALSE
    ))
    args <- c(inputs, list(
      min_descendant_tips = case$cutoff, shift_acceptance_threshold = 1e9,
      num_cores = 1L, IC = case$IC, method = case$method, error = FALSE,
      progress = FALSE, store_model_fit_history = FALSE
    ))
    if (case$cutoff == 12L) {
      expect_warning(
        result <- withr::with_seed(221, do.call(searchOptimalConfiguration, args)),
        "No non-root internal nodes"
      )
    } else {
      result <- withr::with_seed(221, do.call(searchOptimalConfiguration, args))
    }
    expect_baseline_output(result)
    expect_equal(result$baseline_ic, .bifrost_search_ic_value(reference, case$IC),
                 tolerance = 1e-8)
    expect_equal(result$optimal_ic, result$baseline_ic)
    expect_equal(result$model_no_uncertainty$logLik, reference$model$logLik,
                 tolerance = 1e-8)
    expect_equal(result$VCVs[["0"]], reference$model$sigma$Pinv, tolerance = 1e-8)
    expect_equal(result$tree_no_uncertainty_transformed$edge.length,
                 reference$model$corrSt$phy$edge.length, tolerance = 1e-8)
    expect_identical(result$tree_no_uncertainty_untransformed$edge.length,
                     inputs$baseline_tree$edge.length)
    expect_identical(result$num_candidates, if (case$cutoff == 12L) 0L else 2L)
    expect_identical(result$candidate_nodes,
                     if (case$cutoff == 12L) integer() else c(14L, 15L))
  }
})

baseline_output_mock_search <- function(IC, scenario, weights) {
  calls <- list()
  covariance <- matrix(c(2, 0.3, 0.3, 1), 2,
                       dimnames = list(c("trait_a", "trait_b"),
                                       c("trait_a", "trait_b")))
  testthat::local_mocked_bindings(
    .bifrost_search_fit_ic = function(IC, formula, tree, trait_data, ...) {
      states <- unique(phytools::getStates(tree, type = "both"))
      calls[[length(calls) + 1L]] <<- states
      proposal <- length(states) > 1L && !"shift" %in% states
      if (proposal && scenario == "fail") stop("Controlled proposal failure")
      if (proposal && scenario == "warn") warning("Controlled fitting warning")
      score <- if (scenario == "accept") 100 - 50 * (length(states) - 1L) else 100
      model <- list(
        model = if (length(states) == 1L) "BM" else "BMM",
        param = if (length(states) == 1L) NA else setNames(seq_along(states), states),
        sigma = list(Pinv = covariance), variables = list(tree = tree),
        corrSt = list(phy = tree), residuals = trait_data
      )
      list(model = model, GIC = list(GIC = score), BIC = list(BIC = score))
    }
  )
  warnings <- character()
  result <- withCallingHandlers(
    do.call(searchOptimalConfiguration, c(baseline_output_inputs(), list(
      IC = IC, min_descendant_tips = if (scenario == "none") 12L else 10L,
      shift_acceptance_threshold = 20, num_cores = 1L,
      progress = FALSE, store_model_fit_history = TRUE,
      uncertaintyweights = weights == "serial",
      uncertaintyweights_par = weights == "parallel"
    ))),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, warnings = warnings, calls = calls)
}

test_that("baseline fallback and accepted shifts preserve search and weight accounting", {
  for (IC in c("GIC", "BIC")) {
    for (scenario in c("none", "reject", "fail", "accept")) {
      for (weights in c("off", "serial", "parallel")) {
        run <- baseline_output_mock_search(IC, scenario, weights)
        result <- run$result
        expect_identical(run$calls[[1L]], "0")
        expect_identical(result$candidate_nodes,
                         if (scenario == "none") integer() else c(14L, 15L))
        expect_identical(result$num_candidates, length(result$candidate_nodes))
        expect_identical(length(run$calls), if (scenario == "none") 1L else
          if (scenario == "accept" && weights != "off") 7L else 5L)
        if (scenario == "accept") {
          expect_identical(result$model_no_uncertainty$model, "BMM")
          expect_identical(result$shift_nodes_no_uncertainty, c(14L, 15L))
          expect_equal(result$optimal_ic, 0)
          expect_identical(result$VCVs, extractRegimeVCVs(result$model_no_uncertainty))
        } else {
          expect_baseline_output(result)
          expect_equal(result$optimal_ic, 100)
        }
        if (weights == "off") {
          expect_false("ic_weights" %in% names(result))
        } else {
          expect_identical(nrow(result$ic_weights), if (scenario == "accept") 2L else 0L)
          if (scenario == "accept") {
            expect_equal(result$ic_weights$delta_ic, rep(-50, 2))
            expect_equal(result$ic_weights$ic_weight_withshift,
                         rep(1 / (1 + exp(-25)), 2))
          }
        }
        if (scenario == "fail") {
          expect_length(result$warnings, 2L)
          expect_true(all(grepl("Controlled proposal failure", run$warnings)))
        } else if (scenario == "none") {
          expect_length(run$warnings, 1L)
          expect_match(run$warnings, "No non-root internal nodes")
        } else {
          expect_length(run$warnings, 0L)
        }
        trajectory <- icTrajectory(result)
        expect_identical(trajectory$regime_id[1L], "0")
        expect_identical(trajectory$status,
                         c("baseline", rep(switch(scenario, accept = "accepted",
                           fail = "error", "rejected"), result$num_candidates)))
      }
    }
    run <- baseline_output_mock_search(IC, "warn", "off")
    expect_baseline_output(run$result)
    expect_length(run$result$warnings, 2L)
    expect_true(all(grepl("Controlled fitting warning", run$warnings)))
  }
})

test_that("baseline metadata does not change printing, recovery, or BMM-only contracts", {
  result <- baseline_output_mock_search("GIC", "reject", "off")$result
  legacy <- result
  legacy$tree_no_uncertainty_untransformed <- phytools::paintSubTree(
    result$tree_no_uncertainty_untransformed, 13L, state = "shift"
  )
  legacy$tree_no_uncertainty_transformed <- legacy$tree_no_uncertainty_untransformed
  legacy$model_no_uncertainty$corrSt$phy <- legacy$tree_no_uncertainty_untransformed
  legacy$model_no_uncertainty$variables$tree <- legacy$tree_no_uncertainty_untransformed
  legacy$VCVs <- list()
  expect_identical(capture.output(print(result)), capture.output(print(legacy)))
  for (truth in list(integer(), 14L)) {
    simdata <- list(list(shiftNodes = truth, paintedTree = result$tree_no_uncertainty_untransformed))
    capture.output(current <- evaluateShiftRecovery(simdata, list(result), weighted = FALSE))
    capture.output(previous <- evaluateShiftRecovery(simdata, list(legacy), weighted = FALSE))
    expect_identical(current, previous)
  }
  expect_error(rateMap(result, progress = FALSE), "named numeric parameter")
  expect_error(lineage_rates(result, progress = FALSE), "multi-regime BMM")
  expect_error(shift_transitions(result), "multi-regime BMM")
  expect_error(shift_waiting_times(result), "multi-regime BMM")
  expect_error(shift_magnitude_counts(result), "multi-regime BMM")
  expect_error(fit_rate_distribution(result), "multi-regime BMM")
})
