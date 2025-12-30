## Resubmission

This is a resubmission following CRAN feedback on **bifrost 0.1.1**. The previous submission was **bifrost 0.1.1**. The current submission is **bifrost 0.1.2**.

### Changes made in response to CRAN comments

- `plot_ic_acceptance_matrix()` now saves and restores the user’s graphical parameters via an immediate `on.exit()`:
  `oldpar <- par(no.readonly = TRUE)` followed immediately by `on.exit(par(oldpar), add = TRUE)`.
  This ensures `par()` settings are reset even if the function exits early or errors.

### Other updates since 0.1.1 (unrelated to the CRAN note)

- `plot_ic_acceptance_matrix()` gained a user-facing `rate_limits` argument (default `c(-400, 150)`) to control the secondary y-axis limits for the rate-of-improvement overlay, with input validation.
- `bifrost_search` print output was refactored for readability/maintainability (behavior/output preserved), with additional tests.
- Vignette: `ic_weights` is printed as a matrix in the chunk to avoid RStudio paged/Unicode table rendering issues.
- Added a unit test exercising the `rate_limits` validation error path.
- Various maintenance/refactoring and test-suite improvements; no changes to core model algorithms or numerical behavior.

## Test environments

- macOS Sequoia 15.x, R 4.4.2 (local, aarch64)
- GitHub Actions CI: macOS-latest, ubuntu-latest, windows-latest
- R-hub (run via GitHub)

## R CMD check results

0 errors | 0 warnings | 0 notes

Checked locally with `R CMD check --as-cran` / `devtools::check(args = "--as-cran")`.
GitHub Actions CI checks also pass on macOS, Linux, and Windows.

## Downstream dependencies

- None.
