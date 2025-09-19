# tests/testthat/test-generateViridisColorScale.R

# testthat::local_edition(3)

# ---- dependency guard --------------------------------------------------------
skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("viridis")
}

# ---- helper: a fixed, known parameter vector --------------------------------
make_named_params <- function() {
  c(slow = 0.10, medium = 0.50, fast = 0.90)
}

# ---- helper: hex validator (accepts #RRGGBB or #RRGGBBAA) -------------------
.is_hex_color <- function(x) {
  grepl("^#(?:[0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$", x)
}

# ---- Test 1: basic structure and lengths ------------------------------------
test_that("Returns expected list structure with correct lengths", {
  skip_if_missing_deps()

  params <- make_named_params()
  out <- generateViridisColorScale(params)

  expect_type(out, "list")
  expect_named(out, c("NamedColors", "ParamColorMapping"))

  expect_length(out$NamedColors, length(params))
  expect_length(out$ParamColorMapping, length(params))
})

# ---- Test 2: names are preserved and aligned with sorted order --------------
test_that("Names are preserved and correspond to increasing parameter values", {
  skip_if_missing_deps()

  params <- make_named_params()
  ord <- order(params)
  expected_names_sorted <- names(params)[ord]
  expected_values_sorted <- as.numeric(params[ord])

  out <- generateViridisColorScale(params)

  # Check that ParamColorMapping values are sorted ascending
  expect_equal(unname(out$ParamColorMapping), sort(as.numeric(params)))

  # Check that names match the sorted order
  expect_equal(names(out$ParamColorMapping), expected_names_sorted)
  expect_equal(names(out$NamedColors), expected_names_sorted)
})

# ---- Test 3: colors look like valid hex strings and count matches -----------
test_that("Colors are valid hex codes and count matches input", {
  skip_if_missing_deps()

  params <- make_named_params()
  out <- generateViridisColorScale(params)

  expect_true(all(.is_hex_color(out$NamedColors)))
  expect_length(unique(out$NamedColors), length(params))
})

# ---- Test 4: order invariance (permuted input yields same sorted mapping) ---
test_that("Permutation of input still yields a correctly sorted mapping", {
  skip_if_missing_deps()

  params <- make_named_params()
  perm <- sample(seq_along(params))
  params_perm <- params[perm]

  out_perm <- generateViridisColorScale(params_perm)

  expect_equal(unname(out_perm$ParamColorMapping), sort(as.numeric(params)))
  expect_equal(names(out_perm$ParamColorMapping), names(params)[order(params)])
})

# ---- Test 5: duplicate/tied values are handled (non-decreasing order) -------
test_that("Duplicate values are allowed; mapping is non-decreasing and names are a permutation", {
  skip_if_missing_deps()

  params <- c(a = 1, b = 2, c = 2, d = 0.5)
  out <- generateViridisColorScale(params)

  vals <- unname(out$ParamColorMapping)
  expect_true(all(diff(vals) >= 0))   # non-decreasing

  # Names in output should be a permutation of input names
  expect_setequal(names(out$ParamColorMapping), names(params))
  expect_setequal(names(out$NamedColors), names(params))
})

# ---- Test 6: unnamed input still returns colors (names may be NULL) ---------
test_that("Unnamed input returns colors; output names may be NULL", {
  skip_if_missing_deps()

  params <- c(0.2, 0.8, 0.5)
  out <- generateViridisColorScale(params)

  expect_length(out$NamedColors, length(params))
  expect_length(out$ParamColorMapping, length(params))

  # If input had no names, outputs typically have NULL names
  expect_true(is.null(names(out$NamedColors)) || all(nchar(names(out$NamedColors)) > 0))
  expect_true(is.null(names(out$ParamColorMapping)) || all(nchar(names(out$ParamColorMapping)) > 0))
})

# ---- Test 7: constant vector does not error and returns colors --------------
test_that("Constant-valued input does not error and returns valid colors", {
  skip_if_missing_deps()

  # Note: normalization would produce NaN for zero range, but those normalized
  # values aren't used to compute colors in the current implementation.
  params <- c(a = 5, b = 5, c = 5)
  expect_silent({
    out <- generateViridisColorScale(params)
  })

  expect_length(out$NamedColors, length(params))
  expect_length(out$ParamColorMapping, length(params))
  expect_true(all(unname(out$ParamColorMapping) == 5))

  expect_true(all(.is_hex_color(out$NamedColors)))
})

# ---- Test 8: deterministic palette length and endpoints ---------------------
test_that("Color vector matches viridis(n) for endpoints", {
  skip_if_missing_deps()

  params <- c(a = 0.3, b = 0.1, c = 1.0, d = 0.7)
  out <- generateViridisColorScale(params)

  n <- length(params)
  ref <- viridis::viridis(n)

  # Because function uses viridis(n) directly, the set of colors should match,
  # though names/orders are by sorted params.
  expect_setequal(unname(out$NamedColors), ref)

  # Check endpoints are present
  expect_true(ref[1] %in% out$NamedColors)
  expect_true(ref[n] %in% out$NamedColors)
})
