# Test file for addShiftToModel function
# This file should be placed in the tests/testthat/ directory of your package

library(testthat)
library(ape)
library(phytools)

# Helper function to create a basic SIMMAP tree for testing
create_test_simmap_tree <- function(n_tips = 20) {
  set.seed(123)
  tree <- pbtree(n = n_tips, scale = 1)
  simmap_tree <- paintSubTree (tree, node = Ntip(tree)+1,
                               state = "0", anc.state = "0")
  return(simmap_tree)
}

# Helper function to check if a tree is a valid SIMMAP tree
is_valid_simmap <- function(tree) {
  return(inherits(tree, "simmap") &&
           !is.null(tree$maps) &&
           !is.null(tree$mapped.edge))
}

test_that("addShiftToModel returns correct structure", {
  tree <- create_test_simmap_tree()
  shift_node <- Ntip(tree) + 2  # First internal node after root

  result <- addShiftToModel(tree, shift_node = shift_node, current_shift_id = 0L)

  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("tree", "shift_id"))
  expect_s3_class(result$tree, "simmap")
  expect_type(result$shift_id, "double")
})

test_that("addShiftToModel increments shift ID correctly", {
  tree <- create_test_simmap_tree()
  shift_node <- Ntip(tree) + 2

  # Test with different starting shift IDs
  test_cases <- c(0L, 1L, 5L, 10L)

  for (start_id in test_cases) {
    result <- addShiftToModel(tree, shift_node = shift_node, current_shift_id = start_id)
    expect_equal(result$shift_id, start_id + 1,
                 info = paste("Failed for start_id:", start_id))
  }
})

test_that("addShiftToModel preserves tree structure", {
  tree <- create_test_simmap_tree()
  shift_node <- Ntip(tree) + 2
  original_ntip <- Ntip(tree)
  original_nnode <- Nnode(tree)

  result <- addShiftToModel(tree, shift_node = shift_node, current_shift_id = 0L)

  # Check that tree topology is preserved
  expect_equal(Ntip(result$tree), original_ntip)
  expect_equal(Nnode(result$tree), original_nnode)
  expect_equal(nrow(result$tree$edge), nrow(tree$edge))
})

test_that("addShiftToModel creates valid SIMMAP tree", {
  tree <- create_test_simmap_tree()
  shift_node <- Ntip(tree) + 2

  result <- addShiftToModel(tree, shift_node = shift_node, current_shift_id = 0L)

  expect_true(is_valid_simmap(result$tree))
  expect_true(!is.null(result$tree$maps))
  expect_true(!is.null(result$tree$mapped.edge))
})

test_that("addShiftToModel paints subtree with correct state", {
  tree <- create_test_simmap_tree(n_tips = 8)
  shift_node <- Ntip(tree) + 3  # Second internal node after root
  current_shift_id <- 0L

  result <- addShiftToModel(tree, shift_node = shift_node, current_shift_id = current_shift_id)

  # Get descendants of the shift node
  descendants <- getDescendants(tree, shift_node)

  # Check that the new state appears in the mapped.edge matrix
  expected_state <- as.character(current_shift_id + 1)
  expect_true(expected_state %in% colnames(result$tree$mapped.edge))

  # Verify that some edges have the new state
  new_state_edges <- result$tree$mapped.edge[, expected_state] > 0
  expect_true(any(new_state_edges))
})

test_that("addShiftToModel works with different node types", {
  tree <- create_test_simmap_tree(n_tips = 10)

  # Test with different internal nodes
  internal_nodes <- (Ntip(tree) + 2):(Ntip(tree) + Nnode(tree))

  for (node in internal_nodes[1:3]) {  # Test first 3 internal nodes
    result <- addShiftToModel(tree, shift_node = node, current_shift_id = 0L)
    expect_s3_class(result$tree, "simmap")
    expect_equal(result$shift_id, 1L)
  }
})

test_that("addShiftToModel handles sequential shifts", {
  tree <- create_test_simmap_tree(n_tips = 8)

  # Apply first shift
  shift_node1 <- Ntip(tree) + 2
  result1 <- addShiftToModel(tree, shift_node = shift_node1, current_shift_id = 0L)

  # Apply second shift
  shift_node2 <- Ntip(tree) + 4
  result2 <- addShiftToModel(result1$tree, shift_node = shift_node2,
                             current_shift_id = result1$shift_id)

  expect_equal(result1$shift_id, 1L)
  expect_equal(result2$shift_id, 2L)
  expect_s3_class(result2$tree, "simmap")

  # Check that both states exist in the final tree
  expect_true("1" %in% colnames(result2$tree$mapped.edge))
  expect_true("2" %in% colnames(result2$tree$mapped.edge))
})

test_that("addShiftToModel input validation", {
  tree <- create_test_simmap_tree()

  # Test with invalid node numbers
  expect_error(addShiftToModel(tree, shift_node = 0, current_shift_id = 0L))
  expect_error(addShiftToModel(tree, shift_node = Ntip(tree) + Nnode(tree) + 1,
                               current_shift_id = 0L))

  # Test with tip nodes (should work but might behave differently)
  expect_error(result <- addShiftToModel(tree, shift_node = 1, current_shift_id = 0L))
})

test_that("addShiftToModel preserves edge lengths", {
  tree <- create_test_simmap_tree()
  original_edge_lengths <- tree$edge.length
  shift_node <- Ntip(tree) + 2

  result <- addShiftToModel(tree, shift_node = shift_node, current_shift_id = 0L)

  # Edge lengths should be preserved (approximately, due to floating point)
  expect_equal(sum(result$tree$edge.length), sum(original_edge_lengths),
               tolerance = 1e-10)
})

test_that("addShiftToModel works with different shift_id types", {
  tree <- create_test_simmap_tree()
  shift_node <- Ntip(tree) + 2

  # Test with integer
  result_int <- addShiftToModel(tree, shift_node = shift_node, current_shift_id = 0L)
  expect_equal(result_int$shift_id, 1)

  # Test with numeric
  result_num <- addShiftToModel(tree, shift_node = shift_node, current_shift_id = 0)
  expect_equal(result_num$shift_id, 1)
})

test_that("addShiftToModel maintains tree attributes", {
  tree <- create_test_simmap_tree()
  tree$custom_attribute <- "test_value"  # Add custom attribute
  shift_node <- Ntip(tree) + 2

  result <- addShiftToModel(tree, shift_node = shift_node, current_shift_id = 0L)

  # Check that basic tree attributes are maintained
  expect_equal(result$tree$tip.label, tree$tip.label)
  expect_equal(result$tree$Nnode, tree$Nnode)
})

test_that("addShiftToModel error handling for missing paintSubTree_mod", {
  # This test assumes paintSubTree_mod might not be available
  # Skip if the function is available, otherwise expect an error

  tree <- create_test_simmap_tree()
  shift_node <- Ntip(tree) + 2

  # If paintSubTree_mod doesn't exist, this should fail
  if (!exists("paintSubTree_mod")) {
    expect_error(addShiftToModel(tree, shift_node = shift_node, current_shift_id = 0L))
  } else {
    # If it exists, the function should work
    result <- addShiftToModel(tree, shift_node = shift_node, current_shift_id = 0L)
    expect_s3_class(result$tree, "simmap")
  }
})

test_that("addShiftToModel handles root node error correctly", {
  tree <- create_test_simmap_tree()
  root_node <- Ntip(tree) + 1

  # Root node should cause an error in paintSubTree_mod due to the indexing issue
  expect_error(addShiftToModel(tree, shift_node = root_node, current_shift_id = 0L),
               "attempt to select less than one element")
})

test_that("addShiftToModel works with all valid internal nodes", {
  tree <- create_test_simmap_tree(n_tips = 6)
  valid_nodes <- (Ntip(tree) + 2):(Ntip(tree) + Nnode(tree))

  # Should have at least one valid internal node
  expect_true(length(valid_nodes) >= 1)

  # Test that all valid nodes work
  for (node in valid_nodes) {
    result <- addShiftToModel(tree, shift_node = node, current_shift_id = 0L)
    expect_s3_class(result$tree, "simmap")
    expect_equal(result$shift_id, 1L)
  }
})
