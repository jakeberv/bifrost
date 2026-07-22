# bifrost (development version)

* Search progress:
  - `searchOptimalConfiguration()` now displays persistent, Future-compatible CLI progress for candidate scoring, greedy shift evaluation, and IC-weight re-estimation by default.
  - Reached stage rows remain stacked at the bottom while `verbose = TRUE` output streams above them.
  - Stage spinners now redraw continuously during long model fits while completion counts, percentages, and ETA advance only after a fit finishes.
  - Added `progress = FALSE` as an explicit opt-out independent of detailed `verbose` messages.
  - Captured expression-valued search inputs now print safely in `print.bifrost_search()` output.

* Search inputs and diagnostics:
  - Formula-based searches now accept formula objects as well as character strings, numeric response-only data frames for intercept-only searches, and named-column data-frame formulas for pGLS-style workflows.
  - Added `icTrajectory()` and its `plot()` method for inspecting stored search histories; `plot_ic_acceptance_matrix()` remains available as a compatibility wrapper.
  - Stored model-fit histories now retain richer accepted, rejected, and errored candidate records.

* Branch-rate summaries:
  - Added the `rateMap()` workflow and supporting view, control, flagging, print, and plot methods for summarizing branch-rate patterns across completed searches.
  - Improved category legends for uneven rate breaks and strengthened validation of category colors.

* Empirically calibrated simulations:
  - Added empirical null, proportional-shift, and integration-rate robustness workflows centered on the fitted residual covariance.
  - Added `simulation_generator = c("original", "empirical")` for explicit generator selection.
  - Simulation generators now default to `"original"` for exact reproduction of the published operations; the full-covariance Wishart/spectral generator remains available explicitly as `simulation_generator = "empirical"`.
  - Reduced multisession transfer size by using compact namespace-level workers and by avoiding a complete calibration template inside every replicate's call record.

* Documentation / vignettes:
  - Added two rate-map jaw-shape workflows and refreshed the existing jaw-shape vignette.
  - Added a two-part empirically calibrated simulation guide covering performance assessment, search tuning, and empirical application.
  - Added tooling and CI workflows for generated vignette PDFs and executable Colab notebooks.
  - Updated Berv et al. (2026) citation metadata and avian skeleton references for the published *Nature Ecology & Evolution* article DOI.

* Maintenance:
  - Added an automated CRAN downloads tracker and generated chart for the README and development website.

# bifrost 0.1.4

* Documentation / vignettes:
  - Added a new "Quick Start with bifrost" vignette with a minimal end-to-end simulated example.
  - Clarified `searchOptimalConfiguration()` documentation around acceptable tree inputs, recommended `mvgls()` methods (`"H&L"` vs `"LL"`), and the role of `error = TRUE`.
  - Reworked the README to foreground installation, documentation, and citation guidance.
  - Added two pkgdown-only background articles on multivariate Brownian motion / shifts and on whole-tree PCA / model-selection issues.

* Citation / metadata:
  - Updated package authorship metadata to reflect the current author list.
  - Updated `citation("bifrost")` for the live bioRxiv preprint and the application paper.
  - Added the foundational `mvMORPH` citations to the package citation metadata.
  - Added a formatted citation section and dynamic bioRxiv badge to the README.

* Maintenance:
  - Disabled a deprecated vignette-preview step in GitHub Actions CI.

# bifrost 0.1.3

* Addressed CRAN reviewer feedback following review of 0.1.2:
  - Added explicit return-value documentation (`@return` / `\value{}`) for the exported
    `print.bifrost_search()` method, clarifying that the function returns the input object
    invisibly and is called for its printing side effects.

* Plotting:
  - `plot_ic_acceptance_matrix()` gains an optional `baseline_ic` argument to plot and compute
    `diff(IC)` relative to the true no-shift baseline (useful when `matrix_data` begins at the
    first evaluated shift model rather than the true baseline).

* Documentation / vignettes:
  - Updated the jaw-shape vignette with additional static figures (evolutionary correlation heatmap,
    IC-trajectory plot, and branch-rate visualization) and improved plotting annotations.

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
