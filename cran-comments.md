## Resubmission

This is a resubmission following CRAN feedback on **bifrost 0.1.2**. The previous submission was **bifrost 0.1.2**. The current submission is **bifrost 0.1.3**.

### Changes made in response to CRAN comments

- Added explicit return-value documentation (`@return` / `\value{}`) for the exported
  `print.bifrost_search()` method, documenting that the function returns the input object
  invisibly and is called for its printing side effects.

### Other updates since 0.1.2 (unrelated to the CRAN note)

- `plot_ic_acceptance_matrix()` gained an optional `baseline_ic` argument to use the true
  no-shift baseline IC for the baseline annotation and for computing the rate-of-improvement
  series (`diff(IC)`), which is useful when `matrix_data` begins at the first evaluated shift model.
- Jaw-shape vignette: added additional static figures (evolutionary correlation heatmap, IC
  trajectory plot, and branch-rate visualization) and updated plotting/legend annotations for
  clarity. Vignette continues to build under `--as-cran`.
- Expanded unit tests to cover the new `baseline_ic` behavior and input validation.

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
