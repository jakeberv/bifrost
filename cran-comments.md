## Test environments
- macOS Sequoia 15.6.1, R 4.4.2 (local, aarch64)
- MacOS CI (if applicable)
- Windows (win-builder devel/release) [optional: paste results]
- Linux (rhub) [optional: paste results]

## R CMD check results
0 errors | 0 warnings | 2 notes

* NOTE: New submission
  – This is the first CRAN submission of **bifrost**.

* NOTE: "checking for future file timestamps ... unable to verify current time"
  – This appears on some macOS/APFS installations and CI due to filesystem
    timestamp resolution/clock checks. No code relies on file mtimes, and this
    NOTE is not reproducible on other platforms. Nothing in the package writes
    or modifies files during checks.

## Downstream dependencies
- None (new package).
