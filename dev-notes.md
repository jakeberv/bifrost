# Development Notes

### 2025-12-30 — Commit 95afbf8: Add coverage test for `rate_limits` validation (Codecov patch fix)

**Files touched**
- tests/testthat/test-plot_ic_acceptance_matrix.R
- dev-notes.md

**Summary**
- Added a targeted unit test exercising `plot_ic_acceptance_matrix()` input validation for `rate_limits` when `plot_rate_of_improvement = TRUE`.
- Covers the new error path (`stop(...)`) introduced with `rate_limits`, restoring Codecov patch coverage for the plotting changes.

**Intent**
- Keep coverage stable while preserving strict argument validation for the new `rate_limits` option.


### 2025-12-30 — Commit 96028bb: Restore plotting par; add rate-limits control; vignette print fix

**Files touched**
- R/plot_ic_acceptance_matrix.R
- man/plot_ic_acceptance_matrix.Rd
- vignettes/jaw-shape-vignette.Rmd
- dev-notes.md

**Summary**
- `plot_ic_acceptance_matrix()` plotting robustness + configurability:
  - Restores user graphical parameters on exit (`par(no.readonly=TRUE)` + `on.exit(par(oldpar))`) to avoid leaking `par()` settings across calls (addresses CRAN-style note).
  - Added `rate_limits` argument (default `c(-400, 150)`) to control the secondary y-axis limits for the rate-of-improvement overlay.
  - Added input validation for `rate_limits` (numeric length-2, finite; sorted internally).
- Documentation:
  - Regenerated Rd; updated details to reflect restored `par()` behavior and the new `rate_limits` argument.
- Vignette:
  - Print `ic_weights` as a matrix inside the chunk to avoid RStudio rendering issues (paged/Unicode table display).

**Intent**
- Make the IC acceptance plot function safe to call inside packages/vignettes (no persistent `par()` side effects) and provide user control over the diff(IC) axis scale.


### 2025-12-27 — Commit 380134c: Linting cleanup + CRAN compliance fixes (no behavior change)

**Files touched**
- R/plot_ic_acceptance_matrix.R
- R/searchOptimalConfiguration.R
- R/utils.R
- tests/testthat/test-addShiftToModel.R
- tests/testthat/test-fitMvglsAndExtractBIC.formula.R
- tests/testthat/test-fitMvglsAndExtractGIC.formula.r
- tests/testthat/test-plot_ic_acceptance_matrix.R
- .Rbuildignore
- .gitignore

**Summary**
- Performed a focused linting and style cleanup across core R files and tests:
  - Standardized spacing, argument formatting, and quoting (no functional changes).
  - Improved readability and consistency with tidyverse / CRAN style expectations.
- Minor test refactors to make `on.exit()` usage explicit and robust in plotting tests.
- No changes to model logic, algorithms, or numerical behavior.
- This commit intentionally does *not* address the outstanding CRAN note about restoring graphical parameters; that fix will be handled in a subsequent commit.

**Intent**
- Reduce diff noise and improve maintainability before implementing CRAN-requested fixes.
- Ensure future functional changes are isolated from formatting-only edits.

### 2025-12-26 — Commit 2ca1922: Update print output, citation hint, and vignette demo
**Files touched**
- R/bifrost_search-methods.R
- inst/CITATION
- vignettes/jaw-shape-vignette.Rmd
- dev-notes.md

**Summary**
- Print method: reorganized output blocks for readability, added versioned header underline, and added a one-line citation hint (`citation("bifrost")`) at the end of the printout.
- Citation: updated `inst/CITATION` entries to reflect current package metadata.
- Vignette: added a small snippet demonstrating that `bifrost_search` has a custom print method.

### 2025-12-26 — Commit 96c9941: Refactor print helpers; expand/clean test suite; minor CI + vignette updates
**Files touched**
- .Rbuildignore
- .gitignore
- .github/workflows/R-CMD-check.yaml
- R/bifrost_search-methods.R
- dev-notes.md
- tests/testthat/test-addShiftToModel.R
- tests/testthat/test-addShiftToModel_alt.R
- tests/testthat/test-extractRegimeVCVs.R
- tests/testthat/test-fitMvglsAndExtractBIC.formula.R
- tests/testthat/test-fitMvglsAndExtractGIC.formula.r
- tests/testthat/test-generatePaintedTrees.R
- tests/testthat/test-generateViridisColorScale.R
- tests/testthat/test-getDescendants.R
- tests/testthat/test-mvgls-functions.R
- tests/testthat/test-paintSubTree_mod.R
- tests/testthat/test-paintSubTree_removeShift.R
- tests/testthat/test-plot_ic_acceptance_matrix.R
- tests/testthat/test-print-bifrost_search.R
- tests/testthat/test-removeShiftFromTree.R
- tests/testthat/test-searchOptimalConfiguration.R
- tests/testthat/test-whichShifts.R
- vignettes/jaw-shape-vignette.Rmd

**Summary**
- Refactored `print.bifrost_search` into internal helper functions (behavior/output preserved) to reduce complexity and improve maintainability.
- Continued improving test coverage and stability:
  - Added/updated targeted tests across core helpers (`addShiftToModel*`, `fitMvglsAndExtract*`, `generatePaintedTrees`, plotting helpers, shift utilities, and `whichShifts`).
  - Consolidated/cleaned print-method tests in `test-print-bifrost_search.R`.
  - Added additional coverage-oriented assertions while keeping tests deterministic.
- Minor repository/CI housekeeping:
  - Updated ignore rules and the R-CMD-check workflow.
  - Updated jaw-shape vignette content as part of ongoing documentation refinement.

### 2025-12-20 — Commit 108b19e: Add `bifrost_search` print method with IC-history plot; reorganize tests
**Files touched**
- DESCRIPTION
- NAMESPACE
- R/searchOptimalConfiguration.R
- R/bifrost_search-methods.R
- man/print.bifrost_search.Rd
- R/utils.R
- tests/testthat/test-searchOptimalConfiguration.R
- tests/testthat/test-print-bifrost_search.R
- tests/testthat/test-extractRegimeVCVs.R

**Summary**
- Added an S3 class (`bifrost_search`) for `searchOptimalConfiguration()` results and a new `print.bifrost_search()` method:
  - Prints a compact console summary (IC baseline/optimal/ΔIC, search settings, mvgls fit info, shift nodes).
  - Optionally prints an ASCII IC-history plot using `txtplot` when `store_model_fit_history=TRUE`.
  - Prints IC weights (support table) when present.
- Added `txtplot` to `Imports` to support the IC-history console plot.
- Reorganized tests to keep core search tests focused:
  - Moved print-method tests into a dedicated file (`test-print-bifrost_search.R`) and expanded branch coverage.
- Cleaned up `extractRegimeVCVs` by removing/commenting unreachable error code and added tests for missing-component behavior and scaling logic.
- Added `Config/testthat/parallel: false` to stabilize test execution behavior in CI.

### 2025-12-19 — Commit 1e3af6d: Improve IC-weights output consistency and expand tests
**Files touched**
- R/searchOptimalConfiguration.R
- README.md
- tests/testthat/test-searchOptimalConfiguration.R
- vignettes/jaw-shape-vignette.Rmd

**Summary**
- Standardized `ic_weights` output across serial (`uncertaintyweights`) and parallel (`uncertaintyweights_par`) modes:
  - Always returns a `data.frame` with columns:
    `node`, `ic_with_shift`, `ic_without_shift`, `delta_ic`,
    `ic_weight_withshift`, `ic_weight_withoutshift`, `evidence_ratio`
  - When no shifts are detected, returns an empty `data.frame` with the same schema (instead of `NA`).
- Improved no-shifts behavior: `model_no_uncertainty` now returns the baseline `mvgls` model rather than `NULL`.
