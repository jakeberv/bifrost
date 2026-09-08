# Exercise real dataset generation and recovery evaluation. Only the expensive
# model search is replaced; its RNG consumption must not change later datasets.
paired_tuning_fixture <- function() {
  set.seed(762)
  tree <- ape::rtree(20)
  traits <- matrix(rnorm(40), 20, 2, dimnames = list(tree$tip.label, c("x", "y")))
  createSimulationTemplate(tree, traits, formula = "trait_data ~ 1", method = "LL", error = FALSE)
}

paired_tuning_args <- function(template) list(
  template = template, IC = "GIC", shift_acceptance_thresholds = c(10, 20),
  min_descendant_tips_values = c(3L, 4L), tree_tip_count = 16L,
  null_replicates = 2L, recovery_replicates = 2L,
  null_simulation_options = list(simulation_generator = "empirical"),
  proportional_simulation_options = list(simulation_generator = "empirical",
    num_shifts = 1L, min_shift_tips = 4L, max_shift_tips = 8L,
    scale_factor_range = c(2, 4), exclude_range = c(2.5, 3), buffer = 0L),
  base_search_options = list(formula = "trait_data ~ 1", method = "LL", error = FALSE),
  weighted = FALSE, num_cores = 1L, seed = 763L, store_studies = TRUE
)

paired_tuning_search <- function(baseline_tree, trait_data, min_descendant_tips, ...) {
  painted <- generatePaintedTrees(baseline_tree, min_tips = min_descendant_tips)[-1L]
  candidates <- as.integer(sub("^Node ", "", names(painted)))
  # Include a random search result to detect differing worker RNG kinds, too.
  list(shift_nodes_no_uncertainty = NULL, candidate_nodes = candidates,
       num_candidates = length(candidates), optimal_ic = sum(trait_data) + sum(runif(25)))
}

paired_tuning_contents <- function(grid) lapply(grid$studies, function(studies) {
  lapply(studies, function(study) list(
    datasets = lapply(study$simdata, function(sim) { sim$user_input <- NULL; sim }),
    results = study$results
  ))
})

test_that("tuning pairs actual datasets across settings and IC without reusing replicates", {
  args <- paired_tuning_args(paired_tuning_fixture())
  local_rebind("searchOptimalConfiguration", paired_tuning_search, environment(runSearchTuningGrid))
  grid <- do.call(runSearchTuningGrid, args)
  expect_true(grid$paired_settings)
  contents <- paired_tuning_contents(grid)
  for (scenario in c("null", "proportional", "correlation")) {
    reference <- contents[[1L]][[scenario]]$datasets
    expect_false(identical(reference[[1L]], reference[[2L]]))
    for (row in contents[-1L]) expect_identical(row[[scenario]]$datasets, reference)
    expect_length(unique(grid$summary_table[[paste0(scenario, "_seed")]]), 1L)
    expect_identical(grid$study_seeds[[scenario]], grid$summary_table[[paste0(scenario, "_seed")]][[1L]])
  }
  args$IC <- "BIC"
  bic <- do.call(runSearchTuningGrid, args)
  expect_identical(paired_tuning_contents(bic), contents)
  args$shift_acceptance_thresholds <- rev(args$shift_acceptance_thresholds)
  reordered <- do.call(runSearchTuningGrid, args)
  expect_identical(paired_tuning_contents(reordered)[[1L]], contents[[1L]])
  # A fresh run without an explicit seed still pairs rows, but not successive runs.
  args$seed <- NULL
  unseeded <- do.call(runSearchTuningGrid, args)
  unseeded_contents <- paired_tuning_contents(unseeded)
  expect_false(anyNA(unseeded$study_seeds))
  expect_false(identical(unseeded_contents[[1L]], contents[[1L]]))
  for (row in unseeded_contents[-1L]) {
    for (scenario in names(row)) expect_identical(row[[scenario]]$datasets,
                                               unseeded_contents[[1L]][[scenario]]$datasets)
  }
})

test_that("paired tuning matches serial and parallel datasets, results and summaries", {
  skip_on_cran()
  skip_on_covr()
  skip_on_os("windows")
  skip_if_not(future::supportsMulticore())
  withr::local_envvar(c(RSTUDIO = NA, RSTUDIO_SESSION_INITIALIZED = NA))
  args <- paired_tuning_args(paired_tuning_fixture())
  local_rebind("searchOptimalConfiguration", paired_tuning_search, environment(runSearchTuningGrid))
  serial <- do.call(runSearchTuningGrid, args)
  args$num_cores <- 2L
  parallel <- do.call(runSearchTuningGrid, args)
  expect_identical(paired_tuning_contents(parallel), paired_tuning_contents(serial))
  expect_identical(parallel$summary_table, serial$summary_table)
})
