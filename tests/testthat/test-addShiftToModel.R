# tests/testthat/test-addShiftToModel.R

skip_if_tree_editing_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
}

make_add_shift_tree <- function(n_tip = 8, seed = 123, baseline = "0") {
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tip, scale = 1)
  phytools::paintSubTree(
    tree,
    node = ape::Ntip(tree) + 1L,
    state = baseline,
    anc.state = baseline,
    stem = FALSE
  )
}

descendant_edges <- function(tree, node) {
  descendants <- phytools::getDescendants(tree, node)
  which(tree$edge[, 2] %in% descendants)
}

test_that("addShiftToModel increments the regime id and paints the target clade", {
  skip_if_tree_editing_deps()

  tree <- make_add_shift_tree()
  node <- ape::Ntip(tree) + 2L

  result <- addShiftToModel(tree, shift_node = node, current_shift_id = 0L)

  testthat::expect_named(result, c("tree", "shift_id"))
  testthat::expect_equal(result$shift_id, 1L)
  testthat::expect_s3_class(result$tree, "simmap")
  testthat::expect_s3_class(result$tree, "phylo")
  testthat::expect_equal(ape::Ntip(result$tree), ape::Ntip(tree))
  testthat::expect_equal(ape::Nnode(result$tree), ape::Nnode(tree))
  testthat::expect_equal(sum(result$tree$edge.length), sum(tree$edge.length))

  target_edges <- descendant_edges(result$tree, node)
  testthat::expect_gt(length(target_edges), 0L)
  testthat::expect_true(all(vapply(
    target_edges,
    function(edge) identical(names(result$tree$maps[[edge]]), "1"),
    logical(1)
  )))
  testthat::expect_equal(
    as.vector(rowSums(result$tree$mapped.edge)),
    result$tree$edge.length,
    tolerance = 1e-10
  )
})

test_that("addShiftToModel uses sequential ids without overwriting previous regimes", {
  skip_if_tree_editing_deps()

  tree <- make_add_shift_tree(n_tip = 10, seed = 55)
  first_node <- ape::Ntip(tree) + 2L
  second_node <- ape::Ntip(tree) + 4L

  first <- addShiftToModel(tree, shift_node = first_node, current_shift_id = 0L)
  second <- addShiftToModel(first$tree, shift_node = second_node, current_shift_id = first$shift_id)

  testthat::expect_equal(first$shift_id, 1L)
  testthat::expect_equal(second$shift_id, 2L)
  testthat::expect_true("1" %in% colnames(second$tree$mapped.edge))
  testthat::expect_true("2" %in% colnames(second$tree$mapped.edge))
})

test_that("addShiftToModel rejects nodes that cannot define a new internal shift", {
  skip_if_tree_editing_deps()

  tree <- make_add_shift_tree()
  root <- ape::Ntip(tree) + 1L

  testthat::expect_error(
    addShiftToModel(list(), shift_node = root, current_shift_id = 0L),
    "`tree` must be a phylo object"
  )
  testthat::expect_error(
    addShiftToModel(tree, shift_node = 1L, current_shift_id = 0L),
    "`shift_node` must be an internal node"
  )
  testthat::expect_error(
    addShiftToModel(tree, shift_node = NA_integer_, current_shift_id = 0L),
    "`shift_node` must be a single finite integer node id"
  )
  testthat::expect_error(
    addShiftToModel(tree, shift_node = c(root, root + 1L), current_shift_id = 0L),
    "`shift_node` must be a single finite integer node id"
  )
  testthat::expect_error(
    addShiftToModel(tree, shift_node = "node", current_shift_id = 0L),
    "`shift_node` must be a single finite integer node id"
  )
  testthat::expect_error(
    addShiftToModel(tree, shift_node = root, current_shift_id = 0L),
    "`shift_node` cannot be the root node"
  )
  testthat::expect_error(
    addShiftToModel(tree, shift_node = ape::Ntip(tree) + ape::Nnode(tree) + 1L, current_shift_id = 0L),
    "`shift_node` must identify a node in `tree`"
  )
})
