## Resubmission

This is a new submission of **bifrost 0.1.4**. The current CRAN version is **bifrost 0.1.3**.

### Changes since 0.1.3

- Added a new package vignette, **"Quick Start with bifrost"**, with a minimal end-to-end simulated example.
- Clarified `searchOptimalConfiguration()` documentation around tree input requirements, recommended `mvgls()` methods, and the role of `error = TRUE`.
- Reworked the README to foreground installation, documentation, and citation guidance.
- Added two pkgdown-only background articles on multivariate Brownian motion / shifts and whole-tree PCA / model-selection issues. These articles are excluded from the CRAN tarball via `.Rbuildignore`.
- Updated package authorship metadata and refreshed `citation("bifrost")` for the live bioRxiv preprint, the in-press application paper, and the foundational `mvMORPH` references.
- Disabled a deprecated vignette-preview step in GitHub Actions CI. This does not affect package runtime or CRAN checks.

## Test environments

- macOS Sequoia 15.x, R 4.4.2 (local, aarch64)

## R CMD check results

0 errors | 0 warnings | 0 notes

Checked locally with `R CMD build .` and `R CMD check --as-cran bifrost_0.1.4.tar.gz`.

## Downstream dependencies

- None.
