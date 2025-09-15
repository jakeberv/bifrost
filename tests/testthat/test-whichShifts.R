# tests/testthat/test-whichShifts.R

# testthat::local_edition(3)

# ---- dependency guard --------------------------------------------------------
skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
}

# ---- helper: make a SIMMAP baseline (all edges "0") via paintSubTree --------
make_simmap_tree <- function(n_tip = 16, seed = 123, baseline = "0") {
  set.seed(seed)
  tr <- phytools::pbtree(n = n_tip, scale = 1)
  root <- ape::Ntip(tr) + 1L
  # Paint the whole tree from the root; stem = FALSE avoids splitting root edge
  phytools::paintSubTree(tr, node = root, state = baseline, anc.state = baseline, stem = FALSE)
}

# ---- Test 1: single global state -> MRCA for baseline state -----------------
test_that("whichShifts returns MRCA for the baseline state when it covers all tips", {
  skip_if_missing_deps()

  sim <- make_simmap_tree(n_tip = 16, seed = 11, baseline = "0")

  out <- whichShifts(sim)

  # MRCA of *all* tips with state "0" should be the root node
  mrca0 <- ape::getMRCA(sim, sim$tip.label)

  expect_type(out, "integer")
  expect_length(out, 1)
  expect_identical(out, mrca0)
})

# ---- Test 2: two disjoint painted clades -> returns both + baseline MRCA ----
test_that("whichShifts returns MRCA nodes for all states with >=2 tips (including baseline)", {
  skip_if_missing_deps()

  sim <- make_simmap_tree(n_tip = 18, seed = 33, baseline = "0")

  ntip <- ape::Ntip(sim)
  internals <- (ntip + 2L):(ntip + ape::Nnode(sim))  # exclude root (ntip+1)
  # Pick two reasonably separated internal nodes
  nd1 <- internals[2L]
  nd2 <- internals[length(internals) - 1L]

  sim2 <- phytools::paintSubTree(sim, node = nd1, state = "1", stem = FALSE)
  sim3 <- phytools::paintSubTree(sim2, node = nd2, state = "2", stem = FALSE)

  # Expected: MRCA of all "1" tips, MRCA of all "2" tips, and MRCA of remaining "0" tips
  tip_states <- phytools::getStates(sim3, type = "tips")

  tips1 <- names(tip_states[tip_states == "1"])
  tips2 <- names(tip_states[tip_states == "2"])
  tips0 <- names(tip_states[tip_states == "0"])

  exp_nodes <- integer(0)
  if (length(tips1) >= 2) exp_nodes <- c(exp_nodes, ape::getMRCA(sim3, tips1))
  if (length(tips2) >= 2) exp_nodes <- c(exp_nodes, ape::getMRCA(sim3, tips2))
  if (length(tips0) >= 2) exp_nodes <- c(exp_nodes, ape::getMRCA(sim3, tips0))
  exp_nodes <- unique(exp_nodes)

  shifts <- whichShifts(sim3)

  expect_type(shifts, "integer")
  expect_setequal(shifts, exp_nodes)
})

# ---- Test 3: singleton state (one tip only) is ignored; baseline MRCA present
test_that("whichShifts ignores singleton states but still returns baseline MRCA", {
  skip_if_missing_deps()

  sim <- make_simmap_tree(n_tip = 14, seed = 44, baseline = "0")

  # Paint a single tip to a unique state "X"
  tip_id <- 1L
  simX <- phytools::paintSubTree(sim, node = tip_id, state = "X", stem = TRUE)

  # Baseline "0" is still present on many tips -> include its MRCA; ignore "X" singleton
  tip_states <- phytools::getStates(simX, type = "tips")
  tips0 <- names(tip_states[tip_states == "0"])
  mrca0 <- if (length(tips0) >= 2) ape::getMRCA(simX, tips0) else integer(0)

  shifts <- whichShifts(simX)

  expect_type(shifts, "integer")
  # Should be exactly the baseline MRCA (typically the root)
  expect_identical(shifts, mrca0)
})
