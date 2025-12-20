# Development Notes

## Unreleased

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
