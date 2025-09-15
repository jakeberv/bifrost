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
