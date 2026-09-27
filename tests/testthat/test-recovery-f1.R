# Hand-counted fixtures catch zero F1 being mistaken for an undefined ratio.
f1_fixture <- function(truth = 6L, inferred = 7L, weights = rep(0.8, length(inferred))) {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  result <- list(shift_nodes_no_uncertainty = inferred, num_candidates = 2L,
                candidate_nodes = 6:7)
  if (!is.null(weights)) result$ic_weights <- data.frame(node = inferred, ic_weight_withshift = weights)
  list(sim = list(paintedTree = tree, shiftNodes = truth), result = result)
}

test_that("zero-hit recovery is zero for incorrect, absent, or spurious predictions", {
  cases <- list(f1_fixture(), f1_fixture(inferred = integer()),
                f1_fixture(truth = integer()), f1_fixture(inferred = NULL))
  for (case in cases) {
    out <- evaluateShiftRecovery(list(case$sim), list(case$result), fuzzy_distance = 0, verbose = FALSE)
    for (mode in c("strict", "fuzzy")) {
      expect_equal(out[[mode]]$f1, 0)
      expect_equal(out$weighted[[mode]]$f1, 0)
    }
    expect_equal(out$n_evaluable_replicates, 1L)
  }
})

test_that("empty or failed evaluations remain undefined and missing weights are not zero evidence", {
  empty <- f1_fixture(truth = integer(), inferred = integer())
  failed <- f1_fixture(); failed$result$error <- "fit failed"
  incomplete <- f1_fixture(); incomplete$result$shift_nodes_no_uncertainty <- NULL
  for (case in list(empty, failed, incomplete)) {
    out <- evaluateShiftRecovery(list(case$sim), list(case$result), fuzzy_distance = 0, verbose = FALSE)
    for (mode in c("strict", "fuzzy")) {
      expect_true(is.na(out[[mode]]$f1))
      expect_true(is.na(out$weighted[[mode]]$f1))
    }
  }
  missing_weights <- f1_fixture(weights = NULL)
  out <- evaluateShiftRecovery(list(missing_weights$sim), list(missing_weights$result), fuzzy_distance = 0, verbose = FALSE)
  expect_equal(out$strict$f1, 0)
  expect_true(is.na(out$weighted$strict$f1))
})

test_that("partial recovery pools missed shifts and respects IC weights", {
  hit <- f1_fixture(inferred = 6L, weights = 0.8)
  miss <- f1_fixture(inferred = integer())
  out <- evaluateShiftRecovery(list(hit$sim, miss$sim), list(hit$result, miss$result), fuzzy_distance = 0, verbose = FALSE)
  for (mode in c("strict", "fuzzy")) {
    expect_equal(out[[mode]]$f1, 2/3) # TP=1, FN=1, FP=0
    expect_equal(out$weighted[[mode]]$f1, 4/7) # 2*0.8 / (0.8 + 2 true)
    expect_equal(out[[mode]]$recall, 0.5)
  }
})

test_that("F1 tuning retains zero-hit settings as valid evidence", {
  case <- f1_fixture()
  out <- evaluateShiftRecovery(list(case$sim), list(case$result), fuzzy_distance = 0, verbose = FALSE)
  grid <- structure(list(IC = "GIC", base_search_options = list(), summary_table = data.frame(
    setting_id = 1:2, shift_acceptance_threshold = c(10,20), min_descendant_tips = c(2,2),
    null_mean_false_positive_rate = 0, null_fraction_any_false_positive = 0,
    null_evaluable_fraction = 1, proportional_evaluable_fraction = 1, correlation_evaluable_fraction = 1,
    proportional_strict_f1 = c(out$strict$f1, 0.1), correlation_strict_f1 = c(0.9,0.1)
  )), class = "bifrost_search_tuning_grid")
  expect_equal(selectTunedSearchParameters(grid, primary_metric = "strict_f1")$selected_row$setting_id, 1L)
})
