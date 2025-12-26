# tests/testthat/test-removeShiftFromTree.R

# testthat::local_edition(3)

# ---- dependency guard --------------------------------------------------------
skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
}

# ---- helper: make a SIMMAP baseline (all edges "0") via paintSubTree --------
make_simmap_tree <- function(n_tip = 14, seed = 101, baseline = "0") {
  set.seed(seed)
  tr <- phytools::pbtree(n = n_tip, scale = 1)
  root <- ape::Ntip(tr) + 1L
  # Paint the whole tree from the root; stem = FALSE keeps single-segment edges
  phytools::paintSubTree(tr, node = root, state = baseline, anc.state = baseline, stem = FALSE)
}

# ---- helper: pick a non-root internal node with at least one descendant -----
pick_internal_node <- function(tree, offsets = c(2L, 3L, 4L)) {
  ntip <- ape::Ntip(tree)
  root <- ntip + 1L
  max_internal <- ntip + ape::Nnode(tree)

  for (off in offsets) {
    nd <- root + off
    if (nd < max_internal) {
      desc <- phytools::getDescendants(tree, nd)
      if (length(desc) > 0) return(nd)
    }
  }
  # fallback scan
  internals <- (root + 1L):max_internal
  for (nd in internals) {
    desc <- phytools::getDescendants(tree, nd)
    if (length(desc) > 0) return(nd)
  }
  NA_integer_
}

# Test: removeShiftFromTree restores parental state for a painted clade (non-root) (verifies state restoration)
test_that("removeShiftFromTree restores parental state for a painted clade (non-root)", {
  skip_if_missing_deps()

  # 1) SIMMAP baseline with global "0"
  sim0 <- make_simmap_tree(n_tip = 14, seed = 11, baseline = "0")

  # Choose a non-root internal node
  nd <- pick_internal_node(sim0, offsets = c(2L, 3L))
  if (is.na(nd)) testthat::skip("No suitable internal node found (non-root).")

  # Record parent and its state
  parent_nd    <- phytools::getParent(sim0, nd)
  node_states0 <- phytools::getStates(sim0, type = "nodes")
  expect_false(is.na(parent_nd))
  parent_state <- as.character(node_states0[as.character(parent_nd)])
  expect_false(is.na(parent_state))

  # 2) Paint the clade at nd as shift "1" (overwrite = TRUE, stem = TRUE)
  sim_shift <- paintSubTree_mod(
    tree      = sim0,
    node      = nd,
    state     = "1",
    anc.state = "0",
    stem      = TRUE,
    overwrite = TRUE
  )

  # Sanity: descendants are indeed "1"
  desc <- phytools::getDescendants(sim_shift, nd)
  z    <- which(sim_shift$edge[, 2] %in% desc)
  expect_true(length(z) > 0)
  all_one <- all(vapply(
    z,
    function(i) length(sim_shift$maps[[i]]) == 1 && names(sim_shift$maps[[i]])[1] == "1",
    logical(1)
  ))
  expect_true(all_one)

  # 3) Remove the shift and restore to parental state
  out <- removeShiftFromTree(sim_shift, shift_node = nd, stem = FALSE)

  # All descendant edges now match the parent's state (e.g., "0")
  z_out <- which(out$edge[, 2] %in% desc)
  expect_true(length(z_out) > 0)
  all_restored <- all(vapply(
    z_out,
    function(i) length(out$maps[[i]]) == 1 && names(out$maps[[i]])[1] == parent_state,
    logical(1)
  ))
  expect_true(all_restored)

  # Structure intact
  expect_true(inherits(out, "phylo"))
  expect_true("maps" %in% names(out))
  expect_true("mapped.edge" %in% names(out))
})

# Test: removeShiftFromTree is a no-op when the node is not a shift
test_that("removeShiftFromTree is a no-op when the node is not a shift", {
  skip_if_missing_deps()

  sim0 <- make_simmap_tree(n_tip = 10, seed = 99, baseline = "0")
  nd <- pick_internal_node(sim0, offsets = c(2L, 3L))
  if (is.na(nd)) testthat::skip("No suitable internal node found (non-root).")

  out <- removeShiftFromTree(sim0, shift_node = nd, stem = FALSE)

  expect_equal(out$maps, sim0$maps)
  expect_equal(out$mapped.edge, sim0$mapped.edge)
})
