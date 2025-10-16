## Test environments
- macOS Sequoia 15.6.1, R 4.4.2 (local, aarch64)
- macOS-latest (GitHub Actions CI)
- ubuntu-latest (GitHub Actions CI)
- windows-latest (GitHub Actions CI)
- Windows (win-builder devel/release)
- Linux (rhub)

## R CMD check results
0 errors | 0 warnings | 2 notes

* NOTE: New submission  
  – This is the first CRAN submission of **bifrost**.

* NOTE: "checking for future file timestamps ... unable to verify current time"  
  – This appears on some macOS/APFS installations and in CI due to filesystem
    timestamp resolution and clock checks. No code relies on file modification
    times, and this NOTE is not reproducible on other platforms. Nothing in the
    package writes or modifies files during checks.

## Downstream dependencies
- None (new package).

## Additional details
- All unit tests pass across all major platforms.
- Continuous integration and test coverage are verified via GitHub Actions and Codecov.
