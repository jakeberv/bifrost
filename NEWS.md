# bifrost 0.2.0

## Changes for existing users

* Increased the minimum supported R version from 4.1 to 4.2.
* Search defaults are now `min_descendant_tips = 10` and
  `shift_acceptance_threshold = 20`, matching the focal settings of Berv et al.
  (2026). These defaults can change search results when arguments are omitted.
* Searches with `IC = "BIC"` now default to `method = "LL"`. An explicitly
  supplied method takes precedence; GIC searches retain the `mvgls()` default.
* Removed `plot_ic_acceptance_matrix()`. For a `bifrost_search` or compatible
  search-result list, use `plot(icTrajectory(x))`. Legacy callers using a raw
  two-column `matrix_data` object can migrate with:

  ```r
  legacy <- list(
    baseline_ic = baseline_ic,
    IC_used = "GIC",
    model_fit_history = list(ic_acceptance_matrix = matrix_data)
  )
  plot(icTrajectory(legacy))
  ```

  Plotting arguments map from `plot_title` to `main`,
  `plot_rate_of_improvement` to `show_delta`, and `rate_limits` to
  `delta_limits`. Supply `baseline_ic` to `icTrajectory()`.
* Vignettes and empirical datasets are now distributed through the package
  website rather than the CRAN package. Replace former
  `system.file("extdata", ...)` paths with `bifrost_example_file()`. The first
  uncached request downloads the checksum-verified artifact tracked on GitHub
  `main`; subsequent calls use the verified cache unless `refresh = TRUE`.
  Installation, attachment, and package examples do not require these downloads.
* `runSearchTuningGrid()` now pairs simulated datasets across settings in both
  serial and parallel runs. Separate GIC and BIC calls with matching simulation
  inputs and seeds are also paired. Previously seeded grid results will change;
  existing cached results are not replaced. Returned objects record
  `paired_settings` and `study_seeds`.

## Search and result inspection

* Fixed a regression in formula normalization that added an intercept to
  formulas explicitly using `0 +` or `- 1`. Searches now preserve the requested
  intercept setting for numeric and factor predictors.
* Expanded formula support to accept formula objects, numeric response-only
  data frames for intercept-only searches, and named-column data-frame formulas
  for pGLS-style workflows.
* Added `icTrajectory()` and its plot method for inspecting search histories.
  Results now record candidate nodes, resolved progress settings, and detailed
  proposal histories, including accepted, rejected, and errored fits.
* Added persistent, Future-compatible progress displays for candidate scoring,
  shift evaluation, and IC-weight re-estimation. Use `progress = FALSE` to
  disable them independently of `verbose` output.
* Strengthened validation of descendant-tip cutoffs, shift-acceptance thresholds,
  and conflicting uncertainty-weight options. Diagnostics flag small clade
  sizes, low acceptance thresholds, and searches with no eligible candidates.
  Simulation runs muffle repeated settings advisories while preserving fitting
  and optimizer warnings.
* Baseline-only results consistently retain regime label `"0"` and the fitted
  global BM covariance in `VCVs[["0"]]`. Covariance summaries no longer issue
  the proportional multi-regime warning for this single covariance.
* Fixed printing of expression-valued search inputs.

## Rates, shifts, and regime covariance

* Added `rateMap()` and supporting methods for inspecting branch-rate patterns.
  Legends respect uneven category breaks. `generateViridisColorScale()` requires
  numeric input and uses sorted rank rather than numeric distance.
* Added `lineage_rates()` and tools for summarizing shift nodes, transitions,
  waiting times, and magnitudes; fitting and bootstrapping rate distributions;
  and comparing shift magnitudes.
* Added `fit_regime_covariances()`, `fit_regime_covariance_runs()`, module
  diagnostics, correlation-matrix PCA, integration summaries, and
  `regime_integration_pgls()` for post-hoc analyses of fitted regimes.

## Simulation and tuning

* Fixed F1 scores incorrectly reported as `NA` when recovery is zero, including
  searches that miss every true shift. Undefined cases retain `NA`. Corrected
  the supplementary replicate metrics and their export pipeline; pooled
  vignette summaries and selected settings are unchanged.

* Tuning recommendations now retain `method` and `error` settings inherited
  from the simulation template. Explicit search overrides remain authoritative;
  simulation fits, scores, and the selection rule are unchanged.
* Added reproducible simulation templates, null and shifted datasets,
  false-positive and shift-recovery studies, recovery evaluation, and fixed-IC
  tuning grids.
* Added empirical null, proportional-shift, and integration-rate robustness
  workflows based on fitted residual covariance. Generators default to
  `simulation_generator = "original"` to reproduce the published operations;
  the full-covariance Wishart/spectral generator is available with `"empirical"`.
* Added `selectTunedSearchParameters()` to filter settings using false-positive
  and evaluability safeguards and rank feasible settings by fuzzy balanced
  accuracy by default.
* Fixed recovery evaluation for successful searches with a `NULL` shift-node
  vector. Zero-shift results now contribute missed shifts to strict, fuzzy,
  and weighted summaries. Saved results can be reassessed without refitting;
  failed or incomplete records remain excluded.
* Reduced data transfer to parallel workers. Parallel search and simulation
  preserve the caller's Future plan and reproducible RNG state while avoiding
  nested worker oversubscription.

## Documentation and maintenance

* Added runnable help examples and website guides for rate maps, avian skeleton
  analyses, simulation, and tuning, with downloadable PDFs and Colab notebooks.
  Longer help examples use bounded `\donttest{}` blocks exercised in CI.
* Updated citation metadata for the published *Nature Ecology & Evolution*
  article and added a CRAN downloads chart to the README and website.
* Added minimum dependency versions `future (>= 1.49.0)` and
  `phytools (>= 2.0-3)`, added `plotrix`, corrected dependency declarations,
  and moved website-only dependencies to `Config/Needs/website`.
* Improved validation and error handling across lineage-rate, regime-integration,
  simulation, and tuning workflows, and reorganized search internals while
  preserving existing positional arguments.
* Replaced the fragile method-forwarding test affected by mvMORPH 1.2.2 with
  balanced examples that check the requested method and a finite GIC. This
  changes the tests, not the BMM starting-value calculation in mvMORPH.

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
