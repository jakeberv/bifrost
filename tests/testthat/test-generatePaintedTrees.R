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

# -------------------------------------------------------------------
# Test 1: Basic functionality with a standard tree
# -------------------------------------------------------------------
test_that("Basic functionality works correctly", {
  skip_if_missing_deps()

  set.seed(123)
  tree <- ape::rcoal(100)
  min_tips <- 10

  # Capture stdout produced by cat()
  out <- testthat::capture_output(
    painted_trees <- generatePaintedTrees(tree, min_tips)
  )

  # Structure checks
  expect_type(painted_trees, "list")
  expect_true(all(vapply(painted_trees, function(x) inherits(x, "phylo"), logical(1))))
  expect_true(all(grepl("^Node \\d+$", names(painted_trees))))

  # Output checks (because the function uses cat())
  expect_true(grepl("eligible nodes are detected", out))
  expect_true(grepl("sub-models generated", out))
})

# -------------------------------------------------------------------
# Test 2: Tree rooting functionality
# -------------------------------------------------------------------
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

# -------------------------------------------------------------------
# Test 3: Minimum tips threshold
# -------------------------------------------------------------------
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

# -------------------------------------------------------------------
# Test 4: Edge cases
# -------------------------------------------------------------------
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

# -------------------------------------------------------------------
# Test 5: Custom state parameter
# -------------------------------------------------------------------
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

# -------------------------------------------------------------------
# Test 6: Consistency across runs with same input
# -------------------------------------------------------------------
test_that("Results are consistent across runs with same input", {
  skip_if_missing_deps()

  tree <- ape::rcoal(30)
  min_tips <- 5

  result1 <- generatePaintedTrees(tree, min_tips)
  result2 <- generatePaintedTrees(tree, min_tips)

  expect_equal(length(result1), length(result2))
  expect_equal(names(result1), names(result2))
})

# -------------------------------------------------------------------
# Test 7: Tree structure preservation
# -------------------------------------------------------------------
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

# -------------------------------------------------------------------
# Test 8: Indirect check of eligible node helper logic
# -------------------------------------------------------------------
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

# -------------------------------------------------------------------
# Test 9: Console output verification (counts)
# -------------------------------------------------------------------
test_that("Console output reports matching counts", {
  skip_if_missing_deps()

  set.seed(123)
  tree <- ape::rcoal(40)
  min_tips <- 6

  out <- testthat::capture_output(
    painted_trees <- generatePaintedTrees(tree, min_tips)
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
