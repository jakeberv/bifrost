test_that("worker thread limits never change the calling process", {
  # Native thread counts cannot be queried portably. Catch either setter at
  # the library boundary without changing the test runner's thread settings.
  testthat::local_mocked_bindings(
    blas_set_num_threads = function(...) {
      stop("Attempted to change the caller's BLAS threads")
    },
    omp_set_num_threads = function(...) {
      stop("Attempted to change the caller's OpenMP threads")
    },
    .package = "RhpcBLASctl"
  )

  expect_invisible(.bifrost_search_limit_worker_threads(Sys.getpid()))
})
