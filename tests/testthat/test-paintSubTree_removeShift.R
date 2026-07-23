# tests/testthat/test-paintSubTree_removeShift.R

skip_if_remove_shift_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
}

make_remove_shift_tree <- function() {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,(d:1,e:1):1):1);")
  baseline <- phytools::paintSubTree(
    tree,
    node = ape::Ntip(tree) + 1L,
    state = "0",
    anc.state = "0",
    stem = FALSE
  )
  shift_node <- ape::getMRCA(baseline, c("a", "b"))
  shifted <- paintSubTree_mod(
    baseline,
    node = shift_node,
    state = "1",
    anc.state = "0",
    stem = TRUE,
    overwrite = TRUE
  )
  list(
    baseline = baseline,
    shifted = shifted,
    shift_node = shift_node
  )
}

edge_indices_below <- function(tree, node) {
  descendants <- phytools::getDescendants(tree, node)
  which(tree$edge[, 2] %in% descendants)
}

test_that("paintSubTree_removeShift restores matching descendant edges to the parent state", {
  skip_if_remove_shift_deps()
  fixture <- make_remove_shift_tree()

  restored <- paintSubTree_removeShift(fixture$shifted, fixture$shift_node, stem = FALSE)
  descendant_edges <- edge_indices_below(restored, fixture$shift_node)

  testthat::expect_s3_class(restored, "simmap")
  testthat::expect_s3_class(restored, "phylo")
  testthat::expect_equal(ape::Ntip(restored), ape::Ntip(fixture$shifted))
  testthat::expect_equal(ape::Nnode(restored), ape::Nnode(fixture$shifted))
  testthat::expect_true(all(vapply(
    descendant_edges,
    function(edge) identical(names(restored$maps[[edge]]), "0"),
    logical(1)
  )))
  testthat::expect_equal(
    as.vector(rowSums(restored$mapped.edge)),
    restored$edge.length,
    tolerance = 1e-10
  )
})

test_that("paintSubTree_removeShift restores the incoming edge when stem is requested", {
  skip_if_remove_shift_deps()
  fixture <- make_remove_shift_tree()

  restored <- paintSubTree_removeShift(fixture$shifted, fixture$shift_node, stem = TRUE)
  stem_edge <- which(restored$edge[, 2] == fixture$shift_node)

  testthat::expect_length(stem_edge, 1L)
  testthat::expect_true(all(names(restored$maps[[stem_edge]]) == "0"))
  testthat::expect_equal(sum(restored$maps[[stem_edge]]), restored$edge.length[stem_edge])
})

test_that("paintSubTree_removeShift leaves tip requests as no-op tree updates", {
  skip_if_remove_shift_deps()
  fixture <- make_remove_shift_tree()

  restored <- paintSubTree_removeShift(fixture$shifted, shift_node = 1L)

  testthat::expect_s3_class(restored, "simmap")
  testthat::expect_equal(restored$maps, fixture$shifted$maps)
  testthat::expect_equal(restored$mapped.edge, fixture$shifted$mapped.edge)
})

test_that("paintSubTree_removeShift initializes plain phylo inputs before restoring states", {
  skip_if_remove_shift_deps()

  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  tree$edge.length <- NULL
  node <- ape::getMRCA(tree, c("a", "b"))

  restored <- paintSubTree_removeShift(tree, shift_node = node)

  testthat::expect_s3_class(restored, "simmap")
  testthat::expect_equal(length(restored$maps), nrow(restored$edge))
  testthat::expect_equal(
    as.vector(rowSums(restored$mapped.edge)),
    restored$edge.length,
    tolerance = 1e-10
  )
  testthat::expect_true(all(unique(unlist(lapply(restored$maps, names))) == "1"))
})

test_that("paintSubTree_removeShift validates tree and node inputs", {
  skip_if_remove_shift_deps()
  fixture <- make_remove_shift_tree()
  root <- ape::Ntip(fixture$shifted) + 1L

  testthat::expect_error(
    paintSubTree_removeShift("not-a-tree", fixture$shift_node),
    "phylo"
  )
  testthat::expect_error(
    paintSubTree_removeShift(fixture$shifted, root),
    "`shift_node` cannot be the root node"
  )
  testthat::expect_error(
    paintSubTree_removeShift(
      fixture$shifted,
      ape::Ntip(fixture$shifted) + ape::Nnode(fixture$shifted) + 1L
    ),
    "`shift_node` must identify a node in `tree`"
  )
})
