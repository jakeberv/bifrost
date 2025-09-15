# tests/testthat/test-paintSubTree_mod.R

# testthat::local_edition(3)

# ---- dependency guard --------------------------------------------------------
skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
}

# ---- helper: make a SIMMAP baseline (all edges "0") via paintSubTree --------
make_simmap_tree <- function(n_tip = 16, seed = 101, baseline = "0") {
  set.seed(seed)
  tr <- phytools::pbtree(n = n_tip, scale = 1)
  root <- ape::Ntip(tr) + 1L
  # Paint whole tree from the root; stem=FALSE keeps single-segment edges
  phytools::paintSubTree(tr, node = root, state = baseline, anc.state = baseline, stem = FALSE)
}

# ---- helper: pick upstream (non-root) and nested internal nodes by offset ----
pick_up_and_down_by_offset <- function(tree, X_candidates = c(2L, 3L, 4L, 5L)) {
  ntip  <- ape::Ntip(tree)
  root  <- ntip + 1L
  max_i <- ntip + ape::Nnode(tree) # == ntip + (ntip - 1) for fully dichotomous trees

  for (X in X_candidates) {
    nd_up <- root + X
    if (nd_up < max_i) {
      desc <- phytools::getDescendants(tree, nd_up)
      int_desc <- desc[desc > ntip]
      if (length(int_desc) >= 1L) {
        return(list(nd_up = nd_up, nd_down = int_desc[1L]))
      }
    }
  }

  # fallback: scan other internals (excluding root)
  internals <- (root + 1L):max_i
  for (nd in internals) {
    desc <- phytools::getDescendants(tree, nd)
    int_desc <- desc[desc > ntip]
    if (length(int_desc) >= 1L) {
      return(list(nd_up = nd, nd_down = int_desc[1L]))
    }
  }
  NULL
}

# ---- Test 1: overwrite = TRUE (full repaint) --------------------------------
test_that("paintSubTree_mod overwrite=TRUE repaints the subtree to the new state", {
  skip_if_missing_deps()

  sim <- make_simmap_tree(n_tip = 16, seed = 11, baseline = "0")

  pick <- pick_up_and_down_by_offset(sim, X_candidates = c(2L, 3L))
  if (is.null(pick)) testthat::skip("No suitable upstream/downstream pair found.")
  nd_up <- pick$nd_up

  out <- paintSubTree_mod(sim, node = nd_up, state = "1",
                          anc.state = "0", stem = FALSE, overwrite = TRUE)

  desc_up <- phytools::getDescendants(sim, nd_up)
  z <- which(out$edge[, 2] %in% desc_up)

  expect_true(length(z) > 0)
  all_single_1 <- all(vapply(
    z,
    function(i) length(out$maps[[i]]) == 1 && names(out$maps[[i]])[1] == "1",
    logical(1)
  ))
  expect_true(all_single_1)
})

# ---- Test 2: overwrite = FALSE (nested clade preserved under upstream repaint)
test_that("paintSubTree_mod overwrite=FALSE preserves a prepainted nested clade", {
  skip_if_missing_deps()

  sim <- make_simmap_tree(n_tip = 16, seed = 22, baseline = "0")

  pick <- pick_up_and_down_by_offset(sim, X_candidates = c(2L, 3L))
  if (is.null(pick)) testthat::skip("No suitable upstream/downstream pair found.")
  nd_up   <- pick$nd_up
  nd_down <- pick$nd_down

  # Prepaint nested clade to "X" (single-segment edges; stem = FALSE)
  sim2 <- phytools::paintSubTree(sim, node = nd_down, state = "X", stem = FALSE)

  # Selective repaint upstream to "1": only baseline "0" edges flip; "X" must remain
  out <- paintSubTree_mod(sim2, node = nd_up, state = "1",
                          anc.state = "0", stem = FALSE, overwrite = FALSE)

  # Nested clade stays "X"
  desc_down <- phytools::getDescendants(sim2, nd_down)
  z_down <- which(out$edge[, 2] %in% desc_down)
  expect_true(length(z_down) > 0)
  all_down_X <- all(vapply(
    z_down,
    function(i) length(out$maps[[i]]) == 1 && names(out$maps[[i]])[1] == "X",
    logical(1)
  ))
  expect_true(all_down_X)

  # Some edges in the upstream clade but outside the nested clade flip to "1"
  desc_up <- phytools::getDescendants(sim2, nd_up)
  z_up_all  <- which(out$edge[, 2] %in% desc_up)
  z_up_only <- setdiff(z_up_all, z_down)
  expect_true(length(z_up_only) > 0)
  any_up_1 <- any(vapply(
    z_up_only,
    function(i) length(out$maps[[i]]) == 1 && names(out$maps[[i]])[1] == "1",
    logical(1)
  ))
  expect_true(any_up_1)
})



# ---- Test 4: stem = numeric (partial stem painting) -------------------------
test_that("paintSubTree_mod with numeric stem splits parent edge correctly", {
  skip_if_missing_deps()

  sim <- make_simmap_tree(n_tip = 16, seed = 44, baseline = "0")

  pick <- pick_up_and_down_by_offset(sim, X_candidates = c(2L, 3L))
  if (is.null(pick)) testthat::skip("No suitable upstream/downstream pair found.")
  nd_up <- pick$nd_up

  stem_fraction <- 0.3
  out <- paintSubTree_mod(sim, node = nd_up, state = "1",
                          anc.state = "0", stem = stem_fraction, overwrite = TRUE)

  # Check that parent edge has two segments with correct proportions
  stem_edge_idx <- which(out$edge[, 2] == nd_up)
  expect_equal(length(out$maps[[stem_edge_idx]]), 2)
  expect_equal(names(out$maps[[stem_edge_idx]]), c("0", "1"))

  # Check proportions (allowing for small numerical errors)
  total_length <- sum(out$maps[[stem_edge_idx]])
  anc_prop <- as.vector(out$maps[[stem_edge_idx]][1] / total_length)
  state_prop <- as.vector(out$maps[[stem_edge_idx]][2] / total_length)
  expect_equal(anc_prop, 1 - stem_fraction, tolerance = 1e-10)
  expect_equal(state_prop, stem_fraction, tolerance = 1e-10)
})

# ---- Test 6: Error handling for invalid stem with tips ---------------------
test_that("paintSubTree_mod throws error when stem=FALSE for tip nodes", {
  skip_if_missing_deps()

  sim <- make_simmap_tree(n_tip = 8, seed = 66, baseline = "0")

  tip_node <- 2L
  expect_error(
    paintSubTree_mod(sim, node = tip_node, state = "1",
                     anc.state = "0", stem = FALSE, overwrite = TRUE),
    "stem must be TRUE for node <= N"
  )
})

# ---- Test 7: Tree without edge lengths (should be computed) ----------------
test_that("paintSubTree_mod handles trees without edge.length", {
  skip_if_missing_deps()

  set.seed(77)
  tr <- phytools::pbtree(n = 8, scale = 1)
  tr$edge.length <- NULL  # Remove edge lengths

  root <- ape::Ntip(tr) + 1L
  expect_no_error({
    out <- paintSubTree_mod(tr, node = root, state = "1",
                            anc.state = "0", stem = FALSE, overwrite = TRUE)
  })

  # Should have edge lengths after processing
  expect_true(!is.null(out$edge.length))
  expect_true(all(out$edge.length > 0))
})

# ---- Test 8: Non-phylo object error ----------------------------------------
test_that("paintSubTree_mod throws error for non-phylo objects", {
  skip_if_missing_deps()

  not_a_tree <- list(tip.label = c("A", "B"), edge = matrix(c(3, 3, 1, 2), 2, 2))

  expect_error(
    paintSubTree_mod(not_a_tree, node = 1, state = "1"),
    "tree should be an object of class \"phylo\"."
  )
})

# ---- Test 9: Complex multi-state preservation ------------------------------
test_that("paintSubTree_mod preserves complex multi-state mappings with overwrite=FALSE", {
  skip_if_missing_deps()

  sim <- make_simmap_tree(n_tip = 16, seed = 88, baseline = "0")

  pick <- pick_up_and_down_by_offset(sim, X_candidates = c(2L, 3L))
  if (is.null(pick)) testthat::skip("No suitable upstream/downstream pair found.")
  nd_up <- pick$nd_up
  nd_down <- pick$nd_down

  # Create a more complex scenario with multiple states
  sim2 <- phytools::paintSubTree(sim, node = nd_down, state = "X", stem = FALSE)

  # Add another state to a different part of the tree
  other_internals <- (ape::Ntip(sim2) + 2L):(ape::Ntip(sim2) + ape::Nnode(sim2))
  other_internals <- setdiff(other_internals, c(nd_up, nd_down))
  if (length(other_internals) > 0) {
    sim3 <- phytools::paintSubTree(sim2, node = other_internals[1], state = "Y", stem = FALSE)
  } else {
    sim3 <- sim2
  }

  # Now selectively paint - should preserve both X and Y states
  out <- paintSubTree_mod(sim3, node = nd_up, state = "1",
                          anc.state = "0", stem = FALSE, overwrite = FALSE)

  # Check that we have all expected states
  all_states <- unique(unlist(lapply(out$maps, names)))
  expect_true("X" %in% all_states)  # Should be preserved
  expect_true("1" %in% all_states)  # Should be added
  expect_true("0" %in% all_states)  # Should remain in unmodified parts
})

# ---- Test 10: mapped.edge consistency --------------------------------------
test_that("paintSubTree_mod maintains consistency between maps and mapped.edge", {
  skip_if_missing_deps()

  sim <- make_simmap_tree(n_tip = 16, seed = 99, baseline = "0")

  pick <- pick_up_and_down_by_offset(sim, X_candidates = c(2L, 3L))
  if (is.null(pick)) testthat::skip("No suitable upstream/downstream pair found.")
  nd_up <- pick$nd_up

  out <- paintSubTree_mod(sim, node = nd_up, state = "1",
                          anc.state = "0", stem = 0.4, overwrite = TRUE)

  # Check that mapped.edge sums match edge lengths
  for (i in 1:nrow(out$edge)) {
    maps_sum <- sum(out$maps[[i]])
    mapped_edge_sum <- sum(out$mapped.edge[i, ])
    expect_equal(maps_sum, mapped_edge_sum, tolerance = 1e-10)
    expect_equal(maps_sum, out$edge.length[i], tolerance = 1e-10)
  }
})

# ---- Test 11: Root node painting -------------------------------------------
test_that("paintSubTree_mod can paint from root node", {
  skip_if_missing_deps()

  sim <- make_simmap_tree(n_tip = 8, seed = 111, baseline = "0")
  root <- ape::Ntip(sim) + 1L

  out <- paintSubTree_mod(sim, node = root, state = "1",
                          anc.state = "0", stem = FALSE, overwrite = TRUE)

  # All edges should be painted as "1"
  all_single_1 <- all(vapply(
    seq_along(out$maps),
    function(i) length(out$maps[[i]]) == 1 && names(out$maps[[i]])[1] == "1",
    logical(1)
  ))
  expect_true(all_single_1)
})

# ---- Test 12: Empty subtree edge case --------------------------------------
test_that("paintSubTree_mod handles single-tip subtrees", {
  skip_if_missing_deps()

  # Create a simple tree where we can isolate a single tip
  set.seed(122)
  tr <- phytools::pbtree(n = 4, scale = 1)
  root <- ape::Ntip(tr) + 1L
  sim <- phytools::paintSubTree(tr, node = root, state = "0",
                                anc.state = "0", stem = FALSE)

  # Paint just one tip
  out <- paintSubTree_mod(sim, node = 1L, state = "1",
                          anc.state = "0", stem = TRUE, overwrite = TRUE)

  # Only the edge leading to tip 1 should change
  tip1_edge_idx <- which(out$edge[, 2] == 1L)
  expect_equal(names(out$maps[[tip1_edge_idx]])[1], "1")

  # Other edges should remain "0"
  other_edges <- setdiff(seq_along(out$maps), tip1_edge_idx)
  all_others_zero <- all(vapply(
    other_edges,
    function(i) all(names(out$maps[[i]]) == "0"),
    logical(1)
  ))
  expect_true(all_others_zero)
})

# ---- Test 13: Numeric vs character states ----------------------------------
test_that("paintSubTree_mod works with numeric states", {
  skip_if_missing_deps()

  sim <- make_simmap_tree(n_tip = 8, seed = 133, baseline = 0)  # numeric baseline

  pick <- pick_up_and_down_by_offset(sim, X_candidates = c(2L))
  if (is.null(pick)) testthat::skip("No suitable upstream/downstream pair found.")
  nd_up <- pick$nd_up

  out <- paintSubTree_mod(sim, node = nd_up, state = 1,
                          anc.state = 0, stem = FALSE, overwrite = TRUE)

  # Check that numeric states work
  desc_up <- phytools::getDescendants(sim, nd_up)
  z <- which(out$edge[, 2] %in% desc_up)

  expect_true(length(z) > 0)
  all_numeric_1 <- all(vapply(
    z,
    function(i) length(out$maps[[i]]) == 1 && as.numeric(names(out$maps[[i]])[1]) == 1,
    logical(1)
  ))
  expect_true(all_numeric_1)
})
