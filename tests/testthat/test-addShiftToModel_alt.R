# tests/testthat/test-addShiftToModel.R

# testthat::local_edition(3)

# ---- dependency guard --------------------------------------------------------
skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
}

# ---- helper: make baseline SIMMAP tree ---------------------------------------
make_simmap_tree <- function(n_tip = 12, seed = 101, baseline = "0") {
  set.seed(seed)
  tr <- phytools::pbtree(n = n_tip, scale = 1)
  root <- ape::Ntip(tr) + 1L
  phytools::paintSubTree(tr, node = root, state = baseline, anc.state = baseline, stem = FALSE)
}

# ---- Test: adding a shift increments ID and paints the correct clade ---------
test_that("addShiftToModel increments shift_id and paints a new regime", {
  skip_if_missing_deps()

  sim <- make_simmap_tree(12, seed = 42, baseline = "0")

  # Pick a non-root internal node
  internals <- (ape::Ntip(sim) + 2L):(ape::Ntip(sim) + ape::Nnode(sim))
  nd <- internals[1L]

  res <- addShiftToModel(sim, shift_node = nd, current_shift_id = 0L)

  # 1) shift_id increments
  expect_equal(res$shift_id, 1L)

  # 2) output is a list with correct structure
  expect_true(is.list(res))
  expect_true(all(c("tree", "shift_id") %in% names(res)))
  expect_true(inherits(res$tree, "phylo"))
  expect_true("maps" %in% names(res$tree))

  # 3) edges in that clade now have state "1"
  desc <- phytools::getDescendants(sim, nd)
  z <- which(res$tree$edge[, 2] %in% desc)
  expect_true(length(z) > 0)
  all_state_1 <- all(vapply(
    z,
    function(i) all(names(res$tree$maps[[i]]) == "1"),
    logical(1)
  ))
  expect_true(all_state_1)
})

# ---- Test: calling twice produces sequential shift_ids ----------------------
test_that("addShiftToModel produces sequential shift_ids when called twice", {
  skip_if_missing_deps()

  sim <- make_simmap_tree(10, seed = 55, baseline = "0")
  internals <- (ape::Ntip(sim) + 2L):(ape::Ntip(sim) + ape::Nnode(sim))

  # First shift
  res1 <- addShiftToModel(sim, shift_node = internals[1], current_shift_id = 0L)
  expect_equal(res1$shift_id, 1L)

  # Second shift on a different node
  res2 <- addShiftToModel(res1$tree, shift_node = internals[2], current_shift_id = res1$shift_id)
  expect_equal(res2$shift_id, 2L)

  # Ensure new clade is painted as "2"
  desc <- phytools::getDescendants(res1$tree, internals[2])
  z <- which(res2$tree$edge[, 2] %in% desc)
  expect_true(all(vapply(z, function(i) all(names(res2$tree$maps[[i]]) == "2"), logical(1))))
})
