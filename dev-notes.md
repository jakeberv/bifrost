# Development Notes

## Unreleased

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
