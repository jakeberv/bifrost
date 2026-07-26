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

# Group: basic mapping and ordering
# Test: Returns expected list structure with correct lengths (checks return value)
test_that("Returns expected list structure with correct lengths", {
  skip_if_missing_deps()

  params <- make_named_params()
  out <- generateViridisColorScale(params)

  expect_type(out, "list")
  expect_named(out, c("NamedColors", "ParamColorMapping"))

  expect_length(out$NamedColors, length(params))
  expect_length(out$ParamColorMapping, length(params))
})

# Test: Names are preserved and correspond to increasing parameter values
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

# Test: Colors are valid hex codes and count matches input
test_that("Colors are valid hex codes and count matches input", {
  skip_if_missing_deps()

  params <- make_named_params()
  out <- generateViridisColorScale(params)

  expect_true(all(.is_hex_color(out$NamedColors)))
  expect_length(unique(out$NamedColors), length(params))
})

# Test: Permutation of input still yields a correctly sorted mapping
test_that("Permutation of input still yields a correctly sorted mapping", {
  skip_if_missing_deps()

  params <- make_named_params()
  perm <- sample(seq_along(params))
  params_perm <- params[perm]

  out_perm <- generateViridisColorScale(params_perm)

  expect_equal(unname(out_perm$ParamColorMapping), sort(as.numeric(params)))
  expect_equal(names(out_perm$ParamColorMapping), names(params)[order(params)])
})

# Test: Color assignments depend on sorted rank, not numeric spacing
test_that("Color assignments are identical for inputs with the same rank ordering", {
  skip_if_missing_deps()

  evenly_spaced <- c(slow = 1, medium = 2, fast = 3)
  radically_spaced <- c(slow = -1000000, medium = -1, fast = 10000000)

  evenly_spaced_out <- generateViridisColorScale(evenly_spaced)
  radically_spaced_out <- generateViridisColorScale(radically_spaced)

  expect_equal(
    evenly_spaced_out$NamedColors,
    radically_spaced_out$NamedColors
  )
  expect_equal(names(evenly_spaced_out$NamedColors), c("slow", "medium", "fast"))
  expect_equal(names(radically_spaced_out$NamedColors), c("slow", "medium", "fast"))
})

# Test: Only numeric parameter vectors are accepted
test_that("Non-numeric parameter vectors produce a clear error", {
  skip_if_missing_deps()

  expect_error(
    generateViridisColorScale(c(slow = "low", fast = "high")),
    "params must be a numeric vector"
  )
  expect_error(
    generateViridisColorScale(factor(c("slow", "fast"))),
    "params must be a numeric vector"
  )
})

# Group: edge cases and input variants
# Test: Duplicate values are allowed; mapping is non-decreasing and names are a permutation
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

# Test: Unnamed input returns colors; output names may be NULL (checks return value)
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

# Test: Constant-valued input does not error and returns valid colors (expect_silent on constant params)
test_that("Constant-valued input does not error and returns valid colors", {
  skip_if_missing_deps()

  # Tied values still receive one color per sorted input position.
  params <- c(a = 5, b = 5, c = 5)
  expect_silent({
    out <- generateViridisColorScale(params)
  })

  expect_length(out$NamedColors, length(params))
  expect_length(out$ParamColorMapping, length(params))
  expect_true(all(unname(out$ParamColorMapping) == 5))

  expect_true(all(.is_hex_color(out$NamedColors)))
})

# Test: Color vector matches viridis(n) for endpoints
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
