# Unit tests for fitMvglsAndExtractGIC and fitMvglsAndExtractBIC functions

library(testthat)
library(ape)
library(mvMORPH)

# Helper function to create test data
create_test_data <- function(n_tips = 20, n_traits = 2, min_tips_high = 5) {
  # Generate a random phylogenetic tree
  set.seed(123)  # For reproducibility
  base_tree <- rcoal(n_tips)
  
  # Generate painted trees using the custom function
  painted_trees <- generatePaintedTrees(base_tree, min_tips = min_tips_high)
  
  # Select the first painted tree for testing
  painted_tree <- painted_trees[[1]]
  
  # Generate simulated trait data
  trait_data <- mvSIM(painted_tree, nsim = n_traits, model = "BMM", param = list(sigma = 0.1, theta = 0))
  
  # Ensure trait_data is a matrix with proper row names
  if (!is.matrix(trait_data)) {
    trait_data <- as.matrix(trait_data)
  }
  rownames(trait_data) <- painted_tree$tip.label
  
  return(list(painted_tree = painted_tree, trait_data = trait_data))
}

# Test suite for fitMvglsAndExtractGIC
test_that("fitMvglsAndExtractGIC works with valid inputs", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("ape")
  
  # Create test data
  test_data <- create_test_data(n_tips = 20, n_traits = 2)
  painted_tree <- test_data$painted_tree
  trait_data <- test_data$trait_data
  
  # Test the function
  result <- fitMvglsAndExtractGIC(painted_tree, trait_data)
  
  # Check that result is a list with correct elements
  expect_type(result, "list")
  expect_named(result, c("model", "GIC"))
  
  # Check that model is of correct class (mvgls object)
  expect_s3_class(result$model, "mvgls")
  
  # Check that GIC is a numeric value
  expect_type(result$GIC$GIC, "matrix")
  expect_length(result$GIC, 7)
  expect_false(is.na(result$GIC))
})


test_that("fitMvglsAndExtractGIC handles multiple traits", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("ape")
  
  # Create test data with multiple traits
  test_data <- create_test_data(n_tips = 25, n_traits = 10)
  painted_tree <- test_data$painted_tree
  trait_data <- test_data$trait_data
  
  # Test the function
  result <- fitMvglsAndExtractGIC(painted_tree, trait_data)
  
  # Check results
  expect_type(result, "list")
  expect_named(result, c("model", "GIC"))
  expect_s3_class(result$model, "mvgls")
  expect_type(result$GIC$GIC, "matrix")
})

test_that("fitMvglsAndExtractGIC throws error for non-matrix trait_data", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("ape")
  
  # Create test data
  test_data <- create_test_data(n_tips = 20, n_traits = 2)
  painted_tree <- test_data$painted_tree
  trait_data <- as.data.frame(test_data$trait_data)  # Convert to data.frame
  
  # Test that error is thrown
  expect_error(
    fitMvglsAndExtractGIC(painted_tree, trait_data),
    "trait_data must be a matrix."
  )
})

test_that("fitMvglsAndExtractGIC does not work on univariate data", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("ape")
  
  # Create test data
  test_data <- create_test_data(n_tips = 20, n_traits = 1)
  painted_tree <- test_data$painted_tree
  trait_data <- test_data$trait_data

  # Test that error is thrown
  expect_error(
    fitMvglsAndExtractGIC(painted_tree, trait_data),
    "trait_data must be multivariate."
  )
})

test_that("fitMvglsAndExtractGIC throws error for mismatched row names", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("ape")
  
  # Create test data
  test_data <- create_test_data(n_tips = 20, n_traits = 2)
  painted_tree <- test_data$painted_tree
  trait_data <- test_data$trait_data
  
  # Modify row names to create mismatch
  rownames(trait_data)[1] <- "wrong_name"
  
  # Test that error is thrown
  expect_error(
    fitMvglsAndExtractGIC(painted_tree, trait_data),
    "Row names of trait_data must exactly match the tip labels of the tree."
  )
})

# Test suite for fitMvglsAndExtractBIC
test_that("fitMvglsAndExtractBIC works with valid inputs", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("ape")
  
  # Create test data
  test_data <- create_test_data(n_tips = 20, n_traits = 2)
  painted_tree <- test_data$painted_tree
  trait_data <- test_data$trait_data
  
  # Test the function
  result <- fitMvglsAndExtractBIC(painted_tree, trait_data)
  
  # Check that result is a list with correct elements
  expect_type(result, "list")
  expect_named(result, c("model", "BIC"))
  
  # Check that model is of correct class (mvgls object)
  expect_s3_class(result$model, "mvgls")
  
  # Check that BIC is a numeric value
  expect_type(result$BIC$BIC, "numeric")
  expect_length(result$BIC, 4)
  expect_false(is.na(result$BIC))
})


test_that("fitMvglsAndExtractBIC handles multiple traits", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("ape")
  
  # Create test data with multiple traits
  test_data <- create_test_data(n_tips = 25, n_traits = 4)
  painted_tree <- test_data$painted_tree
  trait_data <- test_data$trait_data
  
  # Test the function
  result <- fitMvglsAndExtractBIC(painted_tree, trait_data)
  
  # Check results
  expect_type(result, "list")
  expect_named(result, c("model", "BIC"))
  expect_s3_class(result$model, "mvgls")
  expect_type(result$BIC, "double")
})

test_that("fitMvglsAndExtractBIC throws error for non-matrix trait_data", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("ape")
  
  # Create test data
  test_data <- create_test_data(n_tips = 20, n_traits = 2)
  painted_tree <- test_data$painted_tree
  trait_data <- as.data.frame(test_data$trait_data)  # Convert to data.frame
  
  # Test that error is thrown
  expect_error(
    fitMvglsAndExtractBIC(painted_tree, trait_data),
    "trait_data must be a matrix."
  )
})

test_that("fitMvglsAndExtractBIC throws error for mismatched row names", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("ape")
  
  # Create test data
  test_data <- create_test_data(n_tips = 20, n_traits = 2)
  painted_tree <- test_data$painted_tree
  trait_data <- test_data$trait_data
  
  # Modify row names to create mismatch
  rownames(trait_data)[1] <- "wrong_name"
  
  # Test that error is thrown
  expect_error(
    fitMvglsAndExtractBIC(painted_tree, trait_data),
    "Row names of trait_data must exactly match the tip labels of the tree."
  )
})

# Comparative tests between GIC and BIC functions
test_that("GIC and BIC functions produce consistent models", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("ape")
  
  # Create test data
  test_data <- create_test_data(n_tips = 20, n_traits = 2)
  painted_tree <- test_data$painted_tree
  trait_data <- test_data$trait_data
  
  # Test both functions
  result_gic <- fitMvglsAndExtractGIC(painted_tree, trait_data)
  result_bic <- fitMvglsAndExtractBIC(painted_tree, trait_data)
  
  # Check that both models have same logLik (since they're fitting the same model)
  expect_equal(logLik(result_gic$model), logLik(result_bic$model))
  
  # Check that both produce valid information criteria
  expect_type(result_gic$GIC, "double")
  expect_type(result_bic$BIC, "double")
})

# Edge case tests
test_that("functions handle small trees", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("ape")
  
  # Create test data with minimum viable tree size
  test_data <- create_test_data(n_tips = 10, n_traits = 1, min_tips_high = 3)
  painted_tree <- test_data$painted_tree
  trait_data <- test_data$trait_data
  
  # Test both functions
  result_gic <- fitMvglsAndExtractGIC(painted_tree, trait_data)
  result_bic <- fitMvglsAndExtractBIC(painted_tree, trait_data)
  
  # Check results
  expect_type(result_gic, "list")
  expect_type(result_bic, "list")
  expect_s3_class(result_gic$model, "mvgls")
  expect_s3_class(result_bic$model, "mvgls")
})

test_that("functions handle trait data with missing row names", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("ape")
  
  # Create test data
  test_data <- create_test_data(n_tips = 20, n_traits = 2)
  painted_tree <- test_data$painted_tree
  trait_data <- test_data$trait_data
  
  # Remove row names
  rownames(trait_data) <- NULL
  
  # Test that error is thrown for both functions
  expect_error(
    fitMvglsAndExtractGIC(painted_tree, trait_data),
    "Row names of trait_data must exactly match the tip labels of the tree."
  )
  
  expect_error(
    fitMvglsAndExtractBIC(painted_tree, trait_data),
    "Row names of trait_data must exactly match the tip labels of the tree."
  )
})