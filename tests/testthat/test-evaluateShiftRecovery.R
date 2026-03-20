testthat::skip_on_cran()

skip_if_eval_shift_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
}

test_that("evaluateShiftRecovery returns perfect strict metrics for exact matches", {
  skip_if_eval_shift_deps()

  set.seed(50)
  tr <- phytools::paintSubTree(ape::rtree(12), node = 13, state = "ancestral", anc.state = "ancestral")
  true_node <- 15L
  simdata <- list(list(paintedTree = tr, shiftNodes = true_node))
  simresults <- list(list(
    shift_nodes_no_uncertainty = true_node,
    num_candidates = 10L,
    ic_weights = data.frame(
      node = true_node,
      ic_with_shift = 1,
      ic_without_shift = 2,
      delta_ic = -1,
      ic_weight_withshift = 0.8,
      ic_weight_withoutshift = 0.2,
      evidence_ratio = 4
    )
  ))

  out <- evaluateShiftRecovery(simdata, simresults, fuzzy_distance = 2, weighted = TRUE, verbose = FALSE)

  testthat::expect_equal(out$strict$precision, 1)
  testthat::expect_equal(out$strict$recall, 1)
  testthat::expect_equal(out$fuzzy$precision, 1)
  testthat::expect_equal(out$fuzzy$recall, 1)
  testthat::expect_equal(out$weighted$strict$precision, 1)
  testthat::expect_equal(out$weighted$strict$recall, 0.8)
})

test_that("evaluateShiftRecovery distinguishes strict and fuzzy matches", {
  skip_if_eval_shift_deps()

  set.seed(51)
  tr0 <- ape::rtree(14)
  tr <- phytools::paintSubTree(tr0, node = ape::Ntip(tr0) + 1L, state = "ancestral", anc.state = "ancestral")
  internal_nodes <- setdiff(unique(tr$edge[, 1]), ape::Ntip(tr) + 1L)
  true_node <- internal_nodes[1]
  inferred_node <- tr$edge[tr$edge[, 1] == true_node, 2][1]

  simdata <- list(list(paintedTree = tr, shiftNodes = true_node))
  simresults <- list(list(
    shift_nodes_no_uncertainty = inferred_node,
    num_candidates = 12L,
    ic_weights = data.frame(
      node = inferred_node,
      ic_with_shift = 1,
      ic_without_shift = 2,
      delta_ic = -1,
      ic_weight_withshift = 0.6,
      ic_weight_withoutshift = 0.4,
      evidence_ratio = 1.5
    )
  ))

  out <- evaluateShiftRecovery(simdata, simresults, fuzzy_distance = 2, weighted = TRUE, verbose = FALSE)

  testthat::expect_equal(out$strict$recall, 0)
  testthat::expect_equal(out$fuzzy$recall, 1)
  testthat::expect_equal(out$weighted$fuzzy$recall, 0.6)
})

test_that("evaluateShiftRecovery validates inputs", {
  skip_if_eval_shift_deps()

  testthat::expect_error(
    evaluateShiftRecovery(simdata = 1, simresults = list()),
    "must both be lists"
  )
  testthat::expect_error(
    evaluateShiftRecovery(simdata = list(), simresults = list(list())),
    "same length"
  )
  testthat::expect_error(
    evaluateShiftRecovery(simdata = list(), simresults = list(), fuzzy_distance = -1),
    "single non-negative number"
  )
  testthat::expect_error(
    evaluateShiftRecovery(simdata = list(), simresults = list(), weighted = NA),
    "TRUE or FALSE"
  )
  testthat::expect_error(
    evaluateShiftRecovery(simdata = list(), simresults = list(), verbose = NA),
    "TRUE or FALSE"
  )
})

test_that("evaluateShiftRecovery handles edge cases and verbose output", {
  skip_if_eval_shift_deps()

  set.seed(52)
  tr <- phytools::paintSubTree(ape::rtree(10), node = 11, state = "ancestral", anc.state = "ancestral")
  simdata <- list(
    list(paintedTree = tr, shiftNodes = integer(0)),
    list(paintedTree = tr, shiftNodes = 11L),
    list()
  )
  simresults <- list(
    list(
      shift_nodes_no_uncertainty = integer(0),
      num_candidates = 6L,
      ic_weights = data.frame()
    ),
    list(
      shift_nodes_no_uncertainty = integer(0),
      num_candidates = 6L,
      ic_weights = data.frame()
    ),
    list()
  )

  out_text <- paste(
    capture.output(
      out <- evaluateShiftRecovery(
        simdata,
        simresults,
        fuzzy_distance = 1,
        weighted = FALSE,
        verbose = TRUE
      )
    ),
    collapse = "\n"
  )

  testthat::expect_match(out_text, "Strict Performance Metrics")
  testthat::expect_null(out$weighted)
  testthat::expect_equal(unname(out$counts$strict["FN"]), 1)
  testthat::expect_true(is.na(out$strict$precision) || out$strict$precision == 0)
})

test_that("evaluateShiftRecovery prints weighted verbose summaries", {
  skip_if_eval_shift_deps()

  set.seed(53)
  tr0 <- ape::rtree(12)
  tr <- phytools::paintSubTree(tr0, node = ape::Ntip(tr0) + 1L, state = "ancestral", anc.state = "ancestral")
  internal_nodes <- setdiff(unique(tr$edge[, 1]), ape::Ntip(tr) + 1L)
  true_node <- internal_nodes[1]
  inferred_node <- true_node

  simdata <- list(list(paintedTree = tr, shiftNodes = true_node))
  simresults <- list(list(
    shift_nodes_no_uncertainty = inferred_node,
    num_candidates = 10L,
    ic_weights = data.frame(
      node = inferred_node,
      ic_with_shift = 1,
      ic_without_shift = 2,
      delta_ic = -1,
      ic_weight_withshift = 0.75,
      ic_weight_withoutshift = 0.25,
      evidence_ratio = 3
    )
  ))

  out_text <- paste(
    capture.output(
      evaluateShiftRecovery(
        simdata,
        simresults,
        fuzzy_distance = 1,
        weighted = TRUE,
        verbose = TRUE
      )
    ),
    collapse = "\n"
  )

  testthat::expect_match(out_text, "Weighted Metrics")
  testthat::expect_match(out_text, "Weighted Precision \\(strict\\)")
  testthat::expect_match(out_text, "Weighted F1 Score  \\(fuzzy\\)")
})

test_that("evaluateShiftRecovery treats zero-candidate replicates as unevaluable for specificity", {
  skip_if_eval_shift_deps()

  set.seed(54)
  tr <- phytools::paintSubTree(ape::rtree(10), node = 11, state = "ancestral", anc.state = "ancestral")
  simdata <- list(list(paintedTree = tr, shiftNodes = 11L))
  simresults <- list(list(
    shift_nodes_no_uncertainty = integer(0),
    num_candidates = 0L,
    ic_weights = data.frame()
  ))

  out <- evaluateShiftRecovery(simdata, simresults, fuzzy_distance = 1, weighted = FALSE, verbose = FALSE)

  testthat::expect_equal(unname(out$counts$strict["TN"]), 0)
  testthat::expect_equal(unname(out$counts$fuzzy["TN"]), 0)
  testthat::expect_equal(out$strict$recall, 0)
  testthat::expect_true(is.na(out$strict$specificity))
  testthat::expect_true(is.na(out$strict$fpr))
  testthat::expect_true(is.na(out$strict$balanced_accuracy))
  testthat::expect_true(is.na(out$fuzzy$specificity))
  testthat::expect_true(is.na(out$fuzzy$balanced_accuracy))
})
