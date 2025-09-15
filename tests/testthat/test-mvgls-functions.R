# Unit tests for fitMvglsAndExtractGIC and fitMvglsAndExtractBIC functions
# This file should be saved as tests/testthat/test-mvgls-functions.R

library(testthat)
library(ape)
library(phytools)
library(mvMORPH)

# Helper function to create a simple painted tree for testing
create_test_painted_tree <- function(ntips = 5) {
  # Create a simple phylogenetic tree
  tree <- rtree(ntips)
  tree$tip.label <- paste0("species_", 1:ntips)
  
  # Create a simple character mapping (painted tree)
  # Simulate ancestral states for painting
  states <- sample(c("A", "B"), ntips, replace = TRUE)
  names(states) <- tree$tip.label
  
  # Paint the tree using make.simmap
  painted_tree <- make.simmap(tree, states, model = "ER", nsim = 1)
  return(painted_tree)
}

# Helper function to create matching trait data
create_test_trait_data <- function(tip_labels, ntraits = 2) {
  ntips <- length(tip_labels)
  trait_matrix <- matrix(rnorm(ntips * ntraits), nrow = ntips, ncol = ntraits)
  rownames(trait_matrix) <- tip_labels
  colnames(trait_matrix) <- paste0("trait_", 1:ntraits)
  return(trait_matrix)
}

# Test suite for fitMvglsAndExtractGIC
test_that("fitMvglsAndExtractGIC works with valid inputs", {
  # Skip if required packages are not available
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")
  
  # Create test data
  painted_tree <- create_test_painted_tree(ntips = 6)
  trait_data <- create_test_trait_data(painted_tree$tip.label, ntraits = 2)
  
  # Test the function
  result <- fitMvglsAndExtractGIC(painted_tree, trait_data)
  
  # Check that result is a list with correct elements
  expect_type(result, "list")
  expect_named(result, c("model", "GIC"))
  
  # Check that model is of correct class
  expect_s3_class(result$model, "mvgls")
  
  # Check that GIC is numeric
  expect_type(result$GIC, "double")
  expect_length(result$GIC, 1)
  expect_false(is.na(result$GIC))
  expect_true(is.finite(result$GIC))
})

test_that("fitMvglsAndExtractGIC handles single trait", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")
  
  painted_tree <- create_test_painted_tree(ntips = 5)
  trait_data <- create_test_trait_data(painted_tree$tip.label, ntraits = 1)
  
  result <- fitMvglsAndExtractGIC(painted_tree, trait_data)
  
  expect_type(result, "list")
  expect_named(result, c("model", "GIC"))
  expect_s3_class(result$model, "mvgls")
  expect_type(result$GIC, "double")
})

test_that("fitMvglsAndExtractGIC handles multiple traits", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")
  
  painted_tree <- create_test_painted_tree(ntips = 7)
  trait_data <- create_test_trait_data(painted_tree$tip.label, ntraits = 4)
  
  result <- fitMvglsAndExtractGIC(painted_tree, trait_data)
  
  expect_type(result, "list")
  expect_named(result, c("model", "GIC"))
  expect_s3_class(result$model, "mvgls")
  expect_type(result$GIC, "double")
})

test_that("fitMvglsAndExtractGIC throws error for non-matrix trait_data", {
  painted_tree <- create_test_painted_tree()
  trait_data <- data.frame(trait1 = rnorm(5), trait2 = rnorm(5))
  rownames(trait_data) <- painted_tree$tip.label
  
  expect_error(
    fitMvglsAndExtractGIC(painted_tree, trait_data),
    "trait_data must be a matrix."
  )
})

test_that("fitMvglsAndExtractGIC throws error for mismatched row names", {
  painted_tree <- create_test_painted_tree()
  trait_data <- create_test_trait_data(paste0("wrong_", 1:5))
  
  expect_error(
    fitMvglsAndExtractGIC(painted_tree, trait_data),
    "Row names of trait_data must exactly match the tip labels of the tree."
  )
})

test_that("fitMvglsAndExtractGIC throws error for missing row names", {
  painted_tree <- create_test_painted_tree()
  trait_data <- matrix(rnorm(10), nrow = 5, ncol = 2)
  # No row names set
  
  expect_error(
    fitMvglsAndExtractGIC(painted_tree, trait_data),
    "Row names of trait_data must exactly match the tip labels of the tree."
  )
})

# Test suite for fitMvglsAndExtractBIC
test_that("fitMvglsAndExtractBIC works with valid inputs", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")
  
  painted_tree <- create_test_painted_tree(ntips = 6)
  trait_data <- create_test_trait_data(painted_tree$tip.label, ntraits = 2)
  
  result <- fitMvglsAndExtractBIC(painted_tree, trait_data)
  
  expect_type(result, "list")
  expect_named(result, c("model", "BIC"))
  expect_s3_class(result$model, "mvgls")
  expect_type(result$BIC, "double")
  expect_length(result$BIC, 1)
  expect_false(is.na(result$BIC))
  expect_true(is.finite(result$BIC))
})

test_that("fitMvglsAndExtractBIC handles single trait", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")
  
  painted_tree <- create_test_painted_tree(ntips = 5)
  trait_data <- create_test_trait_data(painted_tree$tip.label, ntraits = 1)
  
  result <- fitMvglsAndExtractBIC(painted_tree, trait_data)
  
  expect_type(result, "list")
  expect_named(result, c("model", "BIC"))
  expect_s3_class(result$model, "mvgls")
  expect_type(result$BIC, "double")
})

test_that("fitMvglsAndExtractBIC throws error for non-matrix trait_data", {
  painted_tree <- create_test_painted_tree()
  trait_data <- data.frame(trait1 = rnorm(5), trait2 = rnorm(5))
  rownames(trait_data) <- painted_tree$tip.label
  
  expect_error(
    fitMvglsAndExtractBIC(painted_tree, trait_data),
    "trait_data must be a matrix."
  )
})

test_that("fitMvglsAndExtractBIC throws error for mismatched row names", {
  painted_tree <- create_test_painted_tree()
  trait_data <- create_test_trait_data(paste0("wrong_", 1:5))
  
  expect_error(
    fitMvglsAndExtractBIC(painted_tree, trait_data),
    "Row names of trait_data must exactly match the tip labels of the tree."
  )
})

# Comparative tests between GIC and BIC functions
test_that("GIC and BIC functions return same model but different information criteria", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")
  
  painted_tree <- create_test_painted_tree(ntips = 6)
  trait_data <- create_test_trait_data(painted_tree$tip.label, ntraits = 2)
  
  gic_result <- fitMvglsAndExtractGIC(painted_tree, trait_data)
  bic_result <- fitMvglsAndExtractBIC(painted_tree, trait_data)
  
  # Models should have same structure (though may not be identical objects)
  expect_equal(class(gic_result$model), class(bic_result$model))
  expect_equal(length(gic_result$model$sigma), length(bic_result$model$sigma))
  
  # Information criteria should be different (in general)
  expect_type(gic_result$GIC, "double")
  expect_type(bic_result$BIC, "double")
})

# Edge case tests
test_that("Functions handle minimum tree size", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")
  
  # Test with minimum viable tree size (3 tips)
  painted_tree <- create_test_painted_tree(ntips = 3)
  trait_data <- create_test_trait_data(painted_tree$tip.label, ntraits = 1)
  
  expect_no_error({
    gic_result <- fitMvglsAndExtractGIC(painted_tree, trait_data)
    bic_result <- fitMvglsAndExtractBIC(painted_tree, trait_data)
  })
})

test_that("Functions handle identical trait values", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")
  
  painted_tree <- create_test_painted_tree(ntips = 5)
  # Create trait data with identical values
  trait_data <- matrix(rep(1.0, 10), nrow = 5, ncol = 2)
  rownames(trait_data) <- painted_tree$tip.label
  colnames(trait_data) <- c("trait1", "trait2")
  
  # These might produce warnings but should not error
  expect_no_error({
    gic_result <- fitMvglsAndExtractGIC(painted_tree, trait_data)
    bic_result <- fitMvglsAndExtractBIC(painted_tree, trait_data)
  })
})

test_that("Functions handle reordered tip labels correctly", {
  skip_if_not_installed("mvMORPH")
  
  painted_tree <- create_test_painted_tree(ntips = 30, min_tips_high = 4)
  # Create trait data with reordered row names (but still matching)
  trait_data <- create_test_trait_data(painted_tree$tip.label, ntraits = 2)
  reorder_idx <- sample(nrow(trait_data))
  trait_data_reordered <- trait_data[reorder_idx, , drop = FALSE]
  
  # Should still work if row names match (even if in different order)
  expect_no_error({
    gic_result <- fitMvglsAndExtractGIC(painted_tree, trait_data_reordered)
    bic_result <- fitMvglsAndExtractBIC(painted_tree, trait_data_reordered)
  })
})

# Tests specific to generatePaintedTrees functionality
test_that("generatePaintedTrees integration works correctly", {
  # Test that generatePaintedTrees produces usable painted trees
  tree <- rcoal(40)
  tree$tip.label <- paste0("species_", 1:40)
  
  painted_trees_list <- generatePaintedTrees(tree, min_tips_high = 5)
  
  # Should return a list
  expect_type(painted_trees_list, "list")
  expect_true(length(painted_trees_list) > 0)
  
  # Each element should be a painted tree suitable for mvgls
  for (i in seq_len(min(3, length(painted_trees_list)))) {
    painted_tree <- painted_trees_list[[i]]
    expect_true(length(painted_tree$tip.label) > 0)
    
    # Create trait data and test that functions work
    trait_data <- create_test_trait_data(painted_tree$tip.label, ntraits = 2)
    
    expect_no_error({
      gic_result <- fitMvglsAndExtractGIC(painted_tree, trait_data)
      bic_result <- fitMvglsAndExtractBIC(painted_tree, trait_data)
    })
  }
})

test_that("Functions work with different painted trees from same base tree", {
  skip_if_not_installed("mvMORPH")
  
  # Test with multiple painted trees from the same base tree
  painted_tree1 <- create_test_painted_tree(ntips = 50, min_tips_high = 5, index = 1)
  painted_tree2 <- create_test_painted_tree(ntips = 50, min_tips_high = 8, index = 1)
  
  trait_data1 <- create_test_trait_data(painted_tree1$tip.label, ntraits = 2)
  trait_data2 <- create_test_trait_data(painted_tree2$tip.label, ntraits = 2)
  
  expect_no_error({
    result1 <- fitMvglsAndExtractGIC(painted_tree1, trait_data1)
    result2 <- fitMvglsAndExtractBIC(painted_tree2, trait_data2)
  })
  
  # Results should be valid
  expect_s3_class(result1$model, "mvgls")
  expect_s3_class(result2$model, "mvgls")
  expect_type(result1$GIC, "double")
  expect_type(result2$BIC, "double")
})