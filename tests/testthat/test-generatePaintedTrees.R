# tests/testthat/test-generatePaintedTrees.R

# Use testthat edition 3 if your project is set up for it
# testthat::local_edition(3)

# Guard tests if required packages aren't available
skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
}

# Optional: load deps explicitly (helpful in non-package projects)
# library(ape)
# library(phytools)

# Group: core generation and logging
# Test: Basic functionality works correctly
test_that("Basic functionality works correctly", {
  skip_if_missing_deps()

  set.seed(123)
  tree <- ape::rcoal(100)
  min_tips <- 10

  # Enable verbose messages for this test, then restore previous option
  old_verbose <- getOption("bifrost.verbose")
  options(bifrost.verbose = TRUE)
  on.exit(options(bifrost.verbose = old_verbose), add = TRUE)

  # Capture message() output
  out <- paste(
    testthat::capture_messages(
      painted_trees <- generatePaintedTrees(tree, min_tips)
    ),
    collapse = "\n"
  )

  # Structure checks
  expect_type(painted_trees, "list")
  expect_true(all(vapply(painted_trees, function(x) inherits(x, "phylo"), logical(1))))
  expect_true(all(grepl("^Node \\d+$", names(painted_trees))))

  # Output checks (because the function uses cat())
  expect_true(grepl("eligible nodes are detected", out))
  expect_true(grepl("sub-models generated", out))
})

# Group: tree preprocessing and parameter handling
# Test: Unrooted trees are properly rooted
test_that("Unrooted trees are properly rooted", {
  skip_if_missing_deps()

  set.seed(123)
  tree <- ape::rcoal(20)
  tree <- ape::unroot(tree)  # make unrooted
  expect_false(ape::is.rooted(tree))

  painted_trees <- generatePaintedTrees(tree, min_tips = 3)

  # All returned trees should be rooted
  expect_true(all(vapply(painted_trees, ape::is.rooted, logical(1))))
})

# Test: Minimum tips threshold is respected
test_that("Minimum tips threshold is respected", {
  skip_if_missing_deps()

  set.seed(123)
  tree <- ape::rcoal(50)
  min_tips_high <- 30
  min_tips_low  <- 5

  painted_trees_high <- generatePaintedTrees(tree, min_tips_high)
  painted_trees_low  <- generatePaintedTrees(tree, min_tips_low)

  expect_true(length(painted_trees_high) <= length(painted_trees_low))

  if (length(painted_trees_high) > 0) {
    for (i in seq_along(painted_trees_high)) {
      node_num <- as.numeric(sub("^Node\\s+", "", names(painted_trees_high)[i]))
      descendants <- phytools::getDescendants(tree, node_num)
      tip_desc    <- descendants[descendants <= ape::Ntip(tree)]
      expect_true(length(tip_desc) >= min_tips_high)
    }
  }
})

# Test: Edge cases are handled properly
test_that("Edge cases are handled properly", {
  skip_if_missing_deps()

  set.seed(123)
  small_tree <- ape::rcoal(5)

  result_small <- generatePaintedTrees(small_tree, min_tips = 2)
  expect_type(result_small, "list")

  result_impossible <- generatePaintedTrees(small_tree, min_tips = 10)
  expect_equal(length(result_impossible), 0)

  result_single <- generatePaintedTrees(small_tree, min_tips = 1)
  expect_type(result_single, "list")
})

# Test: Custom state parameter works
test_that("Custom state parameter works", {
  skip_if_missing_deps()

  set.seed(123)
  tree <- ape::rcoal(20)
  custom_state <- "custom_shift"

  painted_trees <- generatePaintedTrees(tree, min_tips = 3, state = custom_state)
  expect_type(painted_trees, "list")
  # Verifying the actual painting would require inspecting attributes,
  # which depends on paintSubTree's implementation.
})

# Test: Painted trees include the requested state in maps and mapped.edge
test_that("Painted trees include the requested state in maps and mapped.edge", {
  skip_if_missing_deps()

  set.seed(123)
  tree <- ape::rcoal(20)
  custom_state <- "shift_state"

  painted_trees <- generatePaintedTrees(tree, min_tips = 3, state = custom_state)
  expect_true(length(painted_trees) > 0)

  for (pt in painted_trees) {
    map_states <- unique(unlist(lapply(pt$maps, names)))
    expect_true(custom_state %in% map_states)
    expect_true(custom_state %in% colnames(pt$mapped.edge))
  }
})

# Group: stability and structure invariants
# Test: Results are consistent across runs with same input (smoke test with constructed inputs)
test_that("Results are consistent across runs with same input", {
  skip_if_missing_deps()

  tree <- ape::rcoal(30)
  min_tips <- 5

  result1 <- generatePaintedTrees(tree, min_tips)
  result2 <- generatePaintedTrees(tree, min_tips)

  expect_equal(length(result1), length(result2))
  expect_equal(names(result1), names(result2))
})

# Test: Original tree structure is preserved in painted trees
test_that("Original tree structure is preserved in painted trees", {
  skip_if_missing_deps()

  set.seed(123)
  tree <- ape::rcoal(25)
  min_tips <- 4

  painted_trees <- generatePaintedTrees(tree, min_tips)

  if (length(painted_trees) > 0) {
    for (painted_tree in painted_trees) {
      expect_equal(ape::Ntip(painted_tree),  ape::Ntip(tree))
      expect_equal(painted_tree$tip.label,    tree$tip.label)
      expect_equal(ape::Nnode(painted_tree), ape::Nnode(tree))
    }
  }
})

# Group: eligibility and reporting checks
# Test: Eligible nodes have at least min_tips descendants
test_that("Eligible nodes have at least min_tips descendants", {
  skip_if_missing_deps()

  set.seed(123)
  tree <- ape::rcoal(30)
  min_tips <- 8

  painted_trees <- generatePaintedTrees(tree, min_tips)

  for (nm in names(painted_trees)) {
    node_num <- as.numeric(sub("^Node\\s+", "", nm))
    descendants <- phytools::getDescendants(tree, node_num)
    tip_desc    <- descendants[descendants <= ape::Ntip(tree)]
    expect_true(length(tip_desc) >= min_tips)
  }
})

# Test: Console output reports matching counts
test_that("Console output reports matching counts", {
  skip_if_missing_deps()

  set.seed(123)
  tree <- ape::rcoal(40)
  min_tips <- 6

  old_verbose <- getOption("bifrost.verbose")
  options(bifrost.verbose = TRUE)
  on.exit(options(bifrost.verbose = old_verbose), add = TRUE)

  out <- paste(
    testthat::capture_messages(
      painted_trees <- generatePaintedTrees(tree, min_tips)
    ),
    collapse = "\n"
  )

  eligible_match  <- regmatches(out, regexpr("\\d+(?= eligible nodes)", out, perl = TRUE))
  generated_match <- regmatches(out, regexpr("\\d+(?= sub-models generated)", out, perl = TRUE))

  if (length(eligible_match) > 0 && length(generated_match) > 0) {
    eligible_count  <- as.numeric(eligible_match)
    generated_count <- as.numeric(generated_match)

    expect_equal(eligible_count, generated_count)
    expect_equal(generated_count, length(painted_trees))
  } else {
    testthat::skip("Output did not contain expected phrases; check generatePaintedTrees() cat() strings.")
  }
})
