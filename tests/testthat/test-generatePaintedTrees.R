library(testthat)
library(ape)
library(phytools)  # Assuming paintSubTree comes from phytools

#' Unit Test Function for generatePaintedTrees
#'
#' This function tests various aspects of the generatePaintedTrees function
#' including input validation, node detection, tree painting, and output format.
test_generatePaintedTrees <- function() {
  
  # Test 1: Basic functionality with a standard tree
  test_that("Basic functionality works correctly", {
    set.seed(123)  # For reproducible results
    tree <- rcoal(100)
    min_tips <- 10
    
    # Capture console output
    result <- capture_output(
      painted_trees <- generatePaintedTrees(tree, min_tips)
    )
    
    # Check that function returns a list
    expect_type(painted_trees, "list")
    
    # Check that all elements in the list are phylo objects
    expect_true(all(sapply(painted_trees, function(x) inherits(x, "phylo"))))
    
    # Check that list names follow expected pattern
    expect_true(all(grepl("^Node \\d+$", names(painted_trees))))
    
    # Check that console output contains expected messages
    expect_true(grepl("eligible nodes are detected", result$output))
    expect_true(grepl("sub-models generated", result$output))
  })
  
  # Test 2: Tree rooting functionality
  test_that("Unrooted trees are properly rooted", {
    set.seed(123)
    tree <- rcoal(20)
    tree <- unroot(tree)  # Make tree unrooted
    
    expect_false(is.rooted(tree))
    
    painted_trees <- suppressMessages(generatePaintedTrees(tree, min_tips = 3))
    
    # All returned trees should be rooted
    expect_true(all(sapply(painted_trees, is.rooted)))
  })
  
  # Test 3: Minimum tips threshold
  test_that("Minimum tips threshold is respected", {
    set.seed(123)
    tree <- rcoal(50)
    min_tips_high <- 30  # Very restrictive threshold
    min_tips_low <- 5    # Permissive threshold
    
    painted_trees_high <- suppressMessages(generatePaintedTrees(tree, min_tips_high))
    painted_trees_low <- suppressMessages(generatePaintedTrees(tree, min_tips_low))
    
    # High threshold should produce fewer painted trees than low threshold
    expect_true(length(painted_trees_high) <= length(painted_trees_low))
    
    # With very high threshold, might get no eligible nodes
    if (length(painted_trees_high) > 0) {
      # Verify that selected nodes actually have enough descendants
      for (i in seq_along(painted_trees_high)) {
        node_num <- as.numeric(gsub("Node ", "", names(painted_trees_high)[i]))
        descendants <- getDescendants(tree, node_num)
        tip_descendants <- descendants[descendants <= Ntip(tree)]
        expect_true(length(tip_descendants) >= min_tips_high)
      }
    }
  })
  
  # Test 4: Edge cases
  test_that("Edge cases are handled properly", {
    set.seed(123)
    
    # Very small tree
    small_tree <- rcoal(5)
    result_small <- suppressMessages(generatePaintedTrees(small_tree, min_tips = 2))
    expect_type(result_small, "list")
    
    # Minimum tips larger than tree size
    result_impossible <- suppressMessages(generatePaintedTrees(small_tree, min_tips = 10))
    expect_equal(length(result_impossible), 0)
    
    # Single tip requirement
    result_single <- suppressMessages(generatePaintedTrees(small_tree, min_tips = 1))
    expect_type(result_single, "list")
  })
  
  # Test 5: Custom state parameter
  test_that("Custom state parameter works", {
    set.seed(123)
    tree <- rcoal(20)
    custom_state <- "custom_shift"
    
    painted_trees <- suppressMessages(generatePaintedTrees(tree, min_tips = 3, state = custom_state))
    
    expect_type(painted_trees, "list")
    # Note: Testing the actual painting would require inspecting tree attributes
    # which depends on how paintSubTree stores the state information
  })
  
  # Test 6: Consistency check
  test_that("Results are consistent across runs with same input", {
    tree <- rcoal(30)  # Don't set seed here to test with fixed tree
    min_tips <- 5
    
    result1 <- suppressMessages(generatePaintedTrees(tree, min_tips))
    result2 <- suppressMessages(generatePaintedTrees(tree, min_tips))
    
    # Should get same number of painted trees
    expect_equal(length(result1), length(result2))
    
    # Should get same node names
    expect_equal(names(result1), names(result2))
  })
  
  # Test 7: Tree structure preservation
  test_that("Original tree structure is preserved in painted trees", {
    set.seed(123)
    tree <- rcoal(25)
    min_tips <- 4
    
    painted_trees <- suppressMessages(generatePaintedTrees(tree, min_tips))
    
    if (length(painted_trees) > 0) {
      for (painted_tree in painted_trees) {
        # Same number of tips
        expect_equal(Ntip(painted_tree), Ntip(tree))
        
        # Same tip labels
        expect_equal(painted_tree$tip.label, tree$tip.label)
        
        # Same number of nodes
        expect_equal(Nnode(painted_tree), Nnode(tree))
      }
    }
  })
  
  # Test 8: Internal helper function
  test_that("getEligibleNodes helper function works correctly", {
    set.seed(123)
    tree <- rcoal(30)
    min_tips <- 8
    
    # We need to extract the helper function to test it
    # This is a bit tricky since it's defined inside the main function
    # We'll test indirectly through the main function's behavior
    
    painted_trees <- suppressMessages(generatePaintedTrees(tree, min_tips))
    
    # Each painted tree should correspond to a node with at least min_tips descendants
    for (tree_name in names(painted_trees)) {
      node_num <- as.numeric(gsub("Node ", "", tree_name))
      descendants <- getDescendants(tree, node_num)
      tip_descendants <- descendants[descendants <= Ntip(tree)]
      expect_true(length(tip_descendants) >= min_tips)
    }
  })
  
  # Test 9: Console output verification
  test_that("Console messages are informative", {
    set.seed(123)
    tree <- rcoal(40)
    min_tips <- 6
    
    output <- capture_output(
      painted_trees <- generatePaintedTrees(tree, min_tips)
    )
    
    # Extract numbers from output
    eligible_match <- regmatches(output$output, regexpr("\\d+(?= eligible nodes)", output$output, perl = TRUE))
    generated_match <- regmatches(output$output, regexpr("\\d+(?= sub-models generated)", output$output, perl = TRUE))
    
    if (length(eligible_match) > 0 && length(generated_match) > 0) {
      eligible_count <- as.numeric(eligible_match)
      generated_count <- as.numeric(generated_match)
      
      # Number of eligible nodes should equal number of generated painted trees
      expect_equal(eligible_count, generated_count)
      expect_equal(generated_count, length(painted_trees))
    }
  })
  
  cat("All tests completed successfully!\n")
}

# Helper function to run the tests
run_painted_trees_tests <- function() {
  cat("Running unit tests for generatePaintedTrees function...\n\n")
  test_generatePaintedTrees()
}

# Example usage:
# run_painted_trees_tests()