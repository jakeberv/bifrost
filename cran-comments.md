## Resubmission

This is a resubmission following CRAN feedback on **bifrost 0.1.0**. The current submission is **bifrost 0.1.1**. Changes made:

- Replaced all uses of `T`/`F` with `TRUE`/`FALSE` (including documentation), and avoided using `T`/`F` as names.
- Made console output suppressible and quiet by default: informational output now uses `message()`/`warning()` and is controlled by `verbose`. When `plot = TRUE` in an interactive RStudio session, some progress output is written via `cat()` so it remains visible while plots are updating.
- Ensured the package does not write to the user's home filespace / working directory: when `store_model_fit_history = TRUE`, per-iteration results are written under a `tempdir()` subdirectory and read back in at the end of the run.
- Ensured any temporary changes to user settings (options / graphical parameters) are restored via an *immediate* `on.exit()` (including `par()` in `plot_ic_acceptance_matrix()`).
- Refined parallelization behavior: parallel candidate scoring uses the `future` framework with platform-appropriate backends, and BLAS/OpenMP threads are temporarily capped to 1 per worker to avoid CPU oversubscription; thread limits are restored immediately after parallel sections and the `future` plan is reset to sequential.

## Test environments

- macOS Sequoia 15.x, R 4.4.2 (local, aarch64)
- macOS-latest (GitHub Actions CI)
- ubuntu-latest (GitHub Actions CI)
- windows-latest (GitHub Actions CI)
- Windows (win-builder devel/release)
- Linux (rhub)

## R CMD check results

0 errors | 0 warnings | 2 notes

Notes observed locally:
- checking CRAN incoming feasibility ... NOTE (resubmission metadata)
- checking for future file timestamps ... NOTE (unable to verify current time)

## Downstream dependencies

- None.

## Additional details

- Vignettes build successfully and documentation renders correctly.
- Unit tests pass across all tested platforms.
- Parallel code paths are supported via the `future` ecosystem; parallel workers cap BLAS/OpenMP threads to 1 to avoid CPU oversubscription.
