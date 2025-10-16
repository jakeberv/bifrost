## Test environments
- macOS Sequoia 15.6.1, R 4.4.2 (local, aarch64)  
- macOS-latest (GitHub Actions CI)  
- ubuntu-latest (GitHub Actions CI)  
- windows-latest (GitHub Actions CI)  
- Windows (win-builder devel/release)  
- Linux (rhub)

## R CMD check results

0 errors | 0 warnings | 0 notes

- All checks pass cleanly on all tested platforms.
- Vignettes build successfully and documentation renders correctly.
- Continuous integration and test coverage are verified via GitHub Actions and Codecov.

## Downstream dependencies

- None (new package).

## Additional details

- This is the first CRAN submission of **bifrost**.  
- All unit tests pass across all major platforms.  
- The package builds and checks without NOTES on macOS, Linux, and Windows.  
- Parallel tests and vignette builds run successfully under `--as-cran`.
