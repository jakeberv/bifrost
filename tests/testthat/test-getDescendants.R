# tests/testthat/test-getDescendants.R

# testthat::local_edition(3)

# ---- dependency guard --------------------------------------------------------
skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("ape")
}

# ---- helper: a fixed, known tree --------------------------------------------
# Newick: ((t1,t2),(t3,(t4,t5)));
# Tip labels are in order, so tips are 1..5 as "t1","t2","t3","t4","t5".
make_known_tree <- function() {
  tr <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")
  # Ensure standard postorder (usually default, but be explicit)
  ape::reorder.phylo(tr, order = "postorder")
}

# Test: Tip node: no descendants; include.node=TRUE returns the tip (checks return value)
test_that("Tip node: no descendants; include.node=TRUE returns the tip", {
  skip_if_missing_deps()
  tr <- make_known_tree()
  # Tip "t1" is tip index which() returns (should be 1)
  tip1 <- which(tr$tip.label == "t1")
  out1 <- getDescendants(tr, tip1, include.node = FALSE)
  expect_length(out1, 0)

  out2 <- getDescendants(tr, tip1, include.node = TRUE)
  expect_setequal(out2, tip1)
})

# Test: MRCA(t4,t5) descendants are exactly {t4, t5}
test_that("MRCA(t4,t5) descendants are exactly {t4, t5}", {
  skip_if_missing_deps()
  tr <- make_known_tree()
  t4 <- which(tr$tip.label == "t4")
  t5 <- which(tr$tip.label == "t5")

  n45 <- ape::getMRCA(tr, c("t4", "t5"))
  expect_true(is.finite(n45))
  out <- getDescendants(tr, n45, include.node = FALSE)

  expect_setequal(out, c(t4, t5))
})

# Test: MRCA(t3, MRCA(t4,t5)) descendants include {t3, t4, t5, n45}
test_that("MRCA(t3, MRCA(t4,t5)) descendants include {t3, t4, t5, n45}", {
  skip_if_missing_deps()
  tr <- make_known_tree()
  t3  <- which(tr$tip.label == "t3")
  t4  <- which(tr$tip.label == "t4")
  t5  <- which(tr$tip.label == "t5")

  n45 <- ape::getMRCA(tr, c("t4", "t5"))
  n345 <- ape::getMRCA(tr, c("t3", tr$tip.label[t4], tr$tip.label[t5])) # == MRCA(t3, n45)
  expect_true(is.finite(n345))

  out <- getDescendants(tr, n345, include.node = FALSE)

  # Should include the two tips {t4, t5}, their MRCA n45, and tip t3
  expect_true(all(c(t3, t4, t5, n45) %in% out))
})

# Test: Root descendants are all nodes except the root itself
test_that("Root descendants are all nodes except the root itself", {
  skip_if_missing_deps()
  tr <- make_known_tree()
  ntip <- ape::Ntip(tr)
  nnode <- ape::Nnode(tr)
  root <- ntip + 1L

  out <- getDescendants(tr, root, include.node = FALSE)

  all_nodes <- seq_len(ntip + nnode)
  expect_setequal(out, setdiff(all_nodes, root))
})

# Test: Out-of-range node returns empty unless include.node=TRUE (checks return value)
test_that("Out-of-range node returns empty unless include.node=TRUE", {
  skip_if_missing_deps()
  tr <- make_known_tree()
  bogus <- ape::Ntip(tr) + ape::Nnode(tr) + 100L

  out1 <- getDescendants(tr, bogus, include.node = FALSE)
  expect_length(out1, 0)

  # NOTE: current implementation will return the 'bogus' id if include.node=TRUE,
  # even though it isn't a real node id in the tree. We assert that behavior here.
  out2 <- getDescendants(tr, bogus, include.node = TRUE)
  expect_setequal(out2, bogus)
})
