# tests/testthat/test-generatePaintedTrees.R

skip_if_painted_tree_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
}

make_painted_tree_fixture <- function() {
  ape::read.tree(text = "((a:1,b:1):1,(c:1,(d:1,e:1):1):1);")
}

eligible_tip_count <- function(tree, node) {
  descendants <- phytools::getDescendants(tree, node)
  sum(descendants <= ape::Ntip(tree))
}

node_id_from_name <- function(name) {
  as.integer(sub("^Node\\s+", "", name))
}

test_that("generatePaintedTrees returns one painted tree for each eligible clade including the root baseline", {
  skip_if_painted_tree_deps()

  tree <- make_painted_tree_fixture()
  root <- ape::Ntip(tree) + 1L
  painted <- generatePaintedTrees(tree, min_tips = 2, state = "shift_state")
  node_ids <- vapply(names(painted), node_id_from_name, integer(1))

  testthat::expect_type(painted, "list")
  testthat::expect_true(length(painted) > 0L)
  testthat::expect_true(root %in% node_ids)
  testthat::expect_true(all(vapply(node_ids, eligible_tip_count, integer(1), tree = tree) >= 2L))
  testthat::expect_true(all(vapply(painted, inherits, logical(1), what = "simmap")))
  testthat::expect_true(all(vapply(painted, function(x) {
    "shift_state" %in% unique(unlist(lapply(x$maps, names))) &&
      "shift_state" %in% colnames(x$mapped.edge) &&
      ape::Ntip(x) == ape::Ntip(tree) &&
      ape::Nnode(x) == ape::Nnode(tree) &&
      identical(x$tip.label, tree$tip.label)
  }, logical(1))))
})

test_that("generatePaintedTrees respects min_tips and returns empty lists when no clade qualifies", {
  skip_if_painted_tree_deps()

  tree <- make_painted_tree_fixture()
  low <- generatePaintedTrees(tree, min_tips = 2)
  high <- generatePaintedTrees(tree, min_tips = 3)
  impossible <- generatePaintedTrees(tree, min_tips = ape::Ntip(tree) + 1L)

  testthat::expect_gte(length(low), length(high))
  testthat::expect_length(impossible, 0L)
})

test_that("generatePaintedTrees roots unrooted input trees before painting", {
  skip_if_painted_tree_deps()

  tree <- ape::unroot(make_painted_tree_fixture())
  testthat::expect_false(ape::is.rooted(tree))

  painted <- generatePaintedTrees(tree, min_tips = 2)

  testthat::expect_true(all(vapply(painted, ape::is.rooted, logical(1))))
})

test_that("generatePaintedTrees validates user-facing arguments", {
  skip_if_painted_tree_deps()

  tree <- make_painted_tree_fixture()

  testthat::expect_error(
    generatePaintedTrees(list(), min_tips = 2),
    "`tree` must be a phylo object"
  )
  testthat::expect_error(
    generatePaintedTrees(tree, min_tips = 0),
    "`min_tips` must be a single finite integer"
  )
  testthat::expect_error(
    generatePaintedTrees(tree, min_tips = 2.5),
    "`min_tips` must be a single finite integer"
  )
  testthat::expect_error(
    generatePaintedTrees(tree, min_tips = 2, state = c("a", "b")),
    "`state` must be a single non-empty character string"
  )
})

test_that("bifrost integer validator enforces optional maximum bounds", {
  testthat::expect_equal(
    .bifrost_check_integer_scalar(3, "value", minimum = 1L, maximum = 5L),
    3L
  )
  testthat::expect_error(
    .bifrost_check_integer_scalar(6, "value", maximum = 5L),
    "`value` must be a single finite integer <= 5"
  )
})

test_that("generatePaintedTrees reports candidate and output counts only when verbose", {
  skip_if_painted_tree_deps()

  tree <- make_painted_tree_fixture()

  old_verbose <- getOption("bifrost.verbose")
  on.exit(options(bifrost.verbose = old_verbose), add = TRUE)

  options(bifrost.verbose = FALSE)
  quiet_messages <- testthat::capture_messages(
    quiet <- generatePaintedTrees(tree, min_tips = 2)
  )
  testthat::expect_length(quiet_messages, 0L)

  options(bifrost.verbose = TRUE)
  verbose_messages <- testthat::capture_messages(
    verbose <- generatePaintedTrees(tree, min_tips = 2)
  )
  testthat::expect_equal(verbose, quiet)
  testthat::expect_true(any(grepl("eligible nodes are detected", verbose_messages)))
  testthat::expect_true(any(grepl("sub-models generated", verbose_messages)))
})
