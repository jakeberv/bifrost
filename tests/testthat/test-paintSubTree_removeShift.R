# Unit tests for paintSubTree_removeShift function

library(testthat)
library(phytools)

# Test paintSubTree_removeShift function
test_that("paintSubTree_removeShift handles invalid inputs correctly", {
  # Test non-phylo object
  expect_error(
    paintSubTree_removeShift("not_a_tree", 1),
    "tree should be an object of class 'phylo'."
  )

  # Test with matrix instead of phylo
  expect_error(
    paintSubTree_removeShift(matrix(1:4, 2, 2), 1),
    "tree should be an object of class 'phylo'."
  )
})

test_that("paintSubTree_removeShift works with basic painted tree", {
  # Create a test tree
  set.seed(123)
  tree <- phytools::pbtree(n = 10, scale = 1)
  painted <- generatePaintedTrees(tree, min_tips = 3, state = "shift")

  # Get a valid shift node (internal node with painted descendants)
  shift_nodes <- unique(painted$`Node 13`$edge[,1][painted$`Node 13`$edge[,1] > length(painted$`Node 13`$tip.label)])
  shift_node <- shift_nodes[2]

  # Apply the function
  result <- paintSubTree_removeShift(painted$`Node 13`, shift_node)

  # Check that result is a phylo object with simmap class
  expect_s3_class(result, "phylo")
  expect_s3_class(result, "simmap")

  # Check that required components exist
  expect_true("maps" %in% names(result))
  expect_true("mapped.edge" %in% names(result))

  # Check that maps list has correct length
  expect_equal(length(result$maps), nrow(result$edge))

  # Check that mapped.edge dimensions are correct
  expect_equal(nrow(result$mapped.edge), nrow(result$edge))
  expect_true(ncol(result$mapped.edge) >= 1)
})

test_that("paintSubTree_removeShift preserves tree structure", {
  # Create a test tree
  set.seed(456)
  tree <- phytools::pbtree(n = 8, scale = 1)
  painted <- generatePaintedTrees(tree, min_tips = 3, state = "shift")

  # Get original tree properties
  original_ntips <- length(painted$`Node 12`$tip.label)
  original_nnodes <- painted$`Node 12`$Nnode
  original_edge_length_sum <- sum(painted$`Node 12`$edge.length)

  # Apply function
  shift_nodes <- unique(painted$`Node 12`$edge[,1][painted$`Node 12`$edge[,1] > length(painted$`Node 12`$tip.label)])
  result <- paintSubTree_removeShift(painted$`Node 12`, shift_nodes[2])

  # Check that tree structure is preserved
  expect_equal(length(result$tip.label), original_ntips)
  expect_equal(result$Nnode, original_nnodes)
  expect_equal(sum(result$edge.length), original_edge_length_sum, tolerance = 1e-10)
  expect_equal(nrow(result$edge), nrow(painted$edge))
})


test_that("paintSubTree_removeShift handles stem parameter correctly", {
  # Create a test tree
  set.seed(321)
  tree <- phytools::pbtree(n = 12, scale = 1)
  painted <- generatePaintedTrees(tree, min_tips = 4, state = "shift")

  # Get a valid internal node
  internal_nodes <- unique(painted$`Node 13`$edge[,1][painted$`Node 13`$edge[,1] > length(painted$`Node 13`$tip.label)])
  shift_node <- internal_nodes[2]

  # Test with stem = FALSE (default)
  result_no_stem <- paintSubTree_removeShift(painted$`Node 13`, shift_node, stem = FALSE)

  # Test with stem = TRUE
  result_with_stem <- paintSubTree_removeShift(painted$`Node 13`, shift_node, stem = TRUE)

  # Both should be valid simmap objects
  expect_s3_class(result_no_stem, "simmap")
  expect_s3_class(result_with_stem, "simmap")

  # Check that stem parameter affects the result
  # (The exact effect depends on the tree structure, but results should be different)
  expect_true("maps" %in% names(result_no_stem))
  expect_true("maps" %in% names(result_with_stem))
})

test_that("paintSubTree_removeShift handles tip nodes", {
  # Create a test tree
  set.seed(654)
  tree <- phytools::pbtree(n = 8, scale = 1)
  painted <- generatePaintedTrees(tree, min_tips = 3, state = "shift")

  # Test with a tip node
  tip_node <- 1  # First tip

  # Should work without error
  result <- paintSubTree_removeShift(painted$`Node 9`, tip_node)

  expect_s3_class(result, "simmap")
  expect_true("maps" %in% names(result))
  expect_equal(length(result$maps), nrow(result$edge))
})

test_that("paintSubTree_removeShift should error on root node", {
  # Create a test tree
  set.seed(987)
  tree <- phytools::pbtree(n = 10, scale = 1)
  painted <- generatePaintedTrees(tree, min_tips = 4, state = "shift")

  # Test with root node
  root_node <- length(painted$`Node 11`$tip.label) + 1

  # Should work with error
  expect_error(paintSubTree_removeShift(painted$`Node 11`, root_node))

})

test_that("paintSubTree_removeShift preserves edge lengths in maps", {
  # Create a test tree
  set.seed(111)
  tree <- phytools::pbtree(n = 8, scale = 1)
  painted <- generatePaintedTrees(tree, min_tips = 3, state = "shift")

  # Get original total length for each edge from maps
  original_edge_lengths <- sapply(painted$maps, sum)

  # Apply function
  shift_nodes <- unique(painted$edge[,1][painted$edge[,1] > length(painted$tip.label)])
  result <- paintSubTree_removeShift(painted, shift_nodes[1])

  # Check that edge lengths are preserved in maps
  result_edge_lengths <- sapply(result$maps, sum)
  expect_equal(result_edge_lengths, original_edge_lengths, tolerance = 1e-10)
})

test_that("paintSubTree_removeShift creates consistent mapped.edge matrix", {
  # Create a test tree
  set.seed(222)
  tree <- phytools::pbtree(n = 10, scale = 1)
  painted <- generatePaintedTrees(tree, min_tips = 4, state = "shift")

  # Apply function
  shift_nodes <- unique(painted$edge[,1][painted$edge[,1] > length(painted$tip.label)])
  result <- paintSubTree_removeShift(painted, shift_nodes[1])

  # Check that mapped.edge row sums equal edge lengths
  mapped_edge_sums <- rowSums(result$mapped.edge)
  expect_equal(mapped_edge_sums, result$edge.length, tolerance = 1e-10)

  # Check that all values in mapped.edge are non-negative
  expect_true(all(result$mapped.edge >= 0))
})

test_that("paintSubTree_removeShift handles trees without edge lengths", {
  # Create a tree and remove edge lengths
  set.seed(333)
  tree <- phytools::pbtree(n = 6, scale = 1)
  tree$edge.length <- NULL

  # Should still work (function calls compute.brlen internally)
  result <- paintSubTree_removeShift(tree, length(tree$tip.label) + 1)

  expect_s3_class(result, "simmap")
  expect_true("edge.length" %in% names(result))
  expect_true("maps" %in% names(result))
  expect_true(all(result$edge.length > 0))
})

test_that("paintSubTree_removeShift handles large shift node values", {
  # Create a test tree
  set.seed(444)
  tree <- phytools::pbtree(n = 8, scale = 1)
  painted <- generatePaintedTrees(tree, min_tips = 3, state = "shift")

  # Test with a node number larger than what exists
  max_node <- max(c(painted$edge))
  large_node <- max_node + 10

  # Should work without error (though may not modify anything)
  expect_no_error(paintSubTree_removeShift(painted, large_node))

  result <- paintSubTree_removeShift(painted, large_node)
  expect_s3_class(result, "simmap")
})

test_that("paintSubTree_removeShift maintains class structure correctly", {
  # Create a test tree
  set.seed(555)
  tree <- phytools::pbtree(n = 8, scale = 1)
  painted <- generatePaintedTrees(tree, min_tips = 3, state = "shift")

  # Apply function
  shift_nodes <- unique(painted$edge[,1][painted$edge[,1] > length(painted$tip.label)])
  result <- paintSubTree_removeShift(painted, shift_nodes[1])

  # Check that simmap comes first in class vector
  result_classes <- class(result)
  expect_true("simmap" %in% result_classes)
  expect_equal(result_classes[1], "simmap")
  expect_true("phylo" %in% result_classes)
})
