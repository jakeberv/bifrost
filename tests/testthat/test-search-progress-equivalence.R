test_that("progress preserves complete searches with serial and parallel IC weights", {
  # Forking also keeps the development namespace available without installing it.
  skip_if_not(future::supportsMulticore(), "forked workers unavailable")
  case <- readRDS(test_path("fixtures", "mvgls-starts-historical.rds"))$cases$bm_hl
  Y <- case$args$data$Y[, 1:2]
  Y[1:8, ] <- Y[1:8, ] * 20
  withr::local_seed(91)
  for (cores in 1:2) {
    for (parallel_weights in c(FALSE, TRUE)) {
      run <- function(progress) suppressMessages(suppressWarnings(
        searchOptimalConfiguration(
          case$args$tree, Y, min_descendant_tips = 4, num_cores = cores,
          shift_acceptance_threshold = 0, method = "LL", verbose = FALSE,
          uncertaintyweights = !parallel_weights,
          uncertaintyweights_par = parallel_weights, progress = progress
        )
      ))
      set.seed(91)
      plain <- run(FALSE)
      expected_seed <- .Random.seed
      set.seed(91)
      displayed <- run(TRUE)
      expect_identical(.Random.seed, expected_seed)
      expect_gt(length(plain$shift_nodes_no_uncertainty), 0L)
      for (field in c("shift_nodes_no_uncertainty", "optimal_ic", "ic_weights")) {
        expect_identical(displayed[[field]], plain[[field]], info = field)
      }
      expect_identical(displayed$model_no_uncertainty$start_values,
                       plain$model_no_uncertainty$start_values)
      expect_identical(displayed$model_no_uncertainty$opt,
                       plain$model_no_uncertainty$opt)
      fits <- function(result) lapply(result$model_fit_history$fits, function(x) {
        list(node = x$candidate_node, accepted = x$accepted,
             starts = x$model$model$start_values, opt = x$model$model$opt)
      })
      expect_identical(fits(displayed), fits(plain))
    }
  }
})
