# tests/testthat/test-rate-map.R

.rate_map_skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  testthat::skip_if_not_installed("progressr")
}

.rate_map_fixture <- function() {
  set.seed(20260521)
  tr <- ape::rtree(8)
  base <- phytools::paintSubTree(
    tree = ape::as.phylo(tr),
    node = ape::Ntip(tr) + 1L,
    state = "0"
  )
  shifted_a <- phytools::paintSubTree(
    tree = base,
    node = ape::Ntip(base) + 2L,
    state = "1"
  )
  shifted_b <- phytools::paintSubTree(
    tree = base,
    node = ape::Ntip(base) + 3L,
    state = "1"
  )
  list(base = base, shifted_a = shifted_a, shifted_b = shifted_b)
}

.rate_map_fit <- function(tree, param = c("0" = 0.1, "1" = 0.5)) {
  list(
    variables = list(tree = tree),
    param = param
  )
}

test_that("rateMap summarizes small synthetic stochastic-map fits", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  fits <- list(
    .rate_map_fit(fx$shifted_a, c("0" = 0.10, "1" = 0.50)),
    .rate_map_fit(fx$shifted_b, c("0" = 0.20, "1" = 0.80))
  )

  out <- rateMap(
    fits,
    res = 8,
    plot = FALSE,
    progress = FALSE,
    ncolors = 16,
    palette = "Viridis"
  )

  testthat::expect_s3_class(out, "rateMap")
  testthat::expect_s3_class(out$tree, "simmap")
  testthat::expect_length(out$cols, 16)
  testthat::expect_length(out$lims, 2)
  testthat::expect_true(is.data.frame(out$intervals))
  testthat::expect_true(all(c(
    "edge",
    "parent",
    "child",
    "depth_start",
    "depth_end",
    "interval_length",
    "value",
    "color_bin"
  ) %in% names(out$intervals)))
  testthat::expect_equal(out$n_fits, 2L)
})

test_that("rateMap auto-detects bifrost_search-style objects", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  fit <- structure(
    list(
      tree_no_uncertainty_untransformed = fx$shifted_a,
      model_no_uncertainty = list(param = c("0" = 0.25, "1" = 0.75))
    ),
    class = c("bifrost_search", "list")
  )

  out <- rateMap(
    list(fit, fit),
    res = 4,
    plot = FALSE,
    progress = FALSE,
    ncolors = 8
  )

  testthat::expect_s3_class(out, "rateMap")
  testthat::expect_equal(out$n_fits, 2L)
  testthat::expect_equal(range(out$intervals$value), out$lims)
})

test_that("rateMap validates mapped state and parameter compatibility", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$shifted_a, c("0" = 0.1))),
      plot = FALSE,
      progress = FALSE
    ),
    "missing rate parameters"
  )

  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$shifted_a, c("0" = 0.1, "1" = 0))),
      log = TRUE,
      plot = FALSE,
      progress = FALSE
    ),
    "non-positive"
  )
})

test_that("rateMap errors on no-shift NA params by default and can omit them", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  invalid <- .rate_map_fit(fx$base, NA_real_)
  valid <- .rate_map_fit(fx$shifted_a, c("0" = 0.1, "1" = 0.5))

  testthat::expect_error(
    rateMap(list(invalid), plot = FALSE, progress = FALSE),
    "named numeric parameter"
  )

  out <- rateMap(
    list(invalid, valid),
    res = 4,
    plot = FALSE,
    progress = FALSE,
    na_action = "omit"
  )

  testthat::expect_equal(out$n_fits, 1L)
  testthat::expect_equal(out$omitted, 1L)
})

test_that("rateMap mapped.edge row sums match edge lengths", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()
  out <- rateMap(
    list(.rate_map_fit(fx$shifted_a), .rate_map_fit(fx$shifted_a)),
    res = 10,
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(
    as.numeric(rowSums(out$tree$mapped.edge)),
    out$tree$edge.length,
    tolerance = 1e-10
  )
})

test_that("rateMap palette resolution handles names, vectors, and reversal", {
  .rate_map_skip_if_missing_deps()

  vir <- .rateMap_colors(5, "Viridis")
  vir_rev <- .rateMap_colors(5, "Viridis", reverse_palette = TRUE)
  testthat::expect_equal(vir_rev, rev(vir))

  ramp <- .rateMap_colors(6, c("black", "white"))
  ramp_rev <- .rateMap_colors(6, c("black", "white"), reverse_palette = TRUE)
  testthat::expect_equal(ramp_rev, rev(ramp))

  fun_cols <- .rateMap_colors(4, function(n) grDevices::hcl.colors(n, "Inferno"))
  testthat::expect_length(fun_cols, 4)
  testthat::expect_type(fun_cols, "character")
})

test_that("plotRateMap draws rate maps on a graphics device", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()
  out <- rateMap(
    list(.rate_map_fit(fx$shifted_a), .rate_map_fit(fx$shifted_b)),
    res = 4,
    plot = FALSE,
    progress = FALSE
  )

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  testthat::expect_error(
    plotRateMap(out, type = "arc", show_tip_labels = FALSE, legend = FALSE, hold = FALSE),
    NA
  )
})
