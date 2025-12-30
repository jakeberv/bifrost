# bifrost 0.1.2

* Addressed CRAN reviewer feedback following review of 0.1.1:
  - `plot_ic_acceptance_matrix()` now saves and restores the user’s graphical parameters via an immediate `on.exit()` (prevents leaking `par()` settings across calls).

* Plotting:
  - Added `rate_limits` argument to `plot_ic_acceptance_matrix()` (default `c(-400, 150)`) to control the secondary y-axis limits for the rate-of-improvement overlay (validated numeric length-2, finite).

* Search results output:
  - Added a `bifrost_search` S3 class and `print.bifrost_search()` method for `searchOptimalConfiguration()` results (compact console summary; optional ASCII IC-history plot via `txtplot` when `store_model_fit_history = TRUE`; prints IC weights when present).
  - Print output includes a citation hint (`citation("bifrost")`); package citation metadata updated in `inst/CITATION`.

* IC weights / no-shift behavior:
  - Standardized `ic_weights` output across serial and parallel uncertainty-weight modes; always returns a `data.frame` with consistent columns, and returns an empty `data.frame` with the same schema when no shifts are detected.
  - When no shifts are detected, `model_no_uncertainty` now returns the baseline `mvgls` model (instead of `NULL`).

* Documentation / vignettes / tests:
  - Updated jaw-shape vignette chunk printing of `ic_weights` to avoid RStudio paged/Unicode rendering issues.
  - Expanded and stabilized unit tests and CI configuration (including `Config/testthat/parallel: false`).

# bifrost 0.1.1

* Addressed CRAN reviewer feedback following review of 0.1.0:
  - Replaced all uses of shorthand `T`/`F` with `TRUE`/`FALSE`.
  - Ensured all informational output is suppressible via `message()`/`warning()` and controlled by a `verbose` flag.
  - Redirected all on-disk output generated during model fitting to `tempdir()` to comply with CRAN file system policies and avoid writing to the user’s working directory.
  - Ensured graphical parameters and global options are restored using immediate `on.exit()` calls.
  - Refined parallelization behavior to be CRAN-safe and cross-platform:
    - Parallel candidate evaluation uses `future` with `multicore` on Unix outside RStudio and `multisession` otherwise.
    - BLAS/OpenMP threads are capped to one per worker during parallel execution to avoid CPU oversubscription.
    - Sequential execution remains the default when `num_cores = 1`.
  - Improved documentation clarity around parallel execution, verbosity, and model fit history storage.

# bifrost 0.1.0

* Initial CRAN submission.
