# bifrost 0.2.3

This is a maintenance update to CRAN version 0.2.2 addressing the CRAN
notification of November 3, 2026.

## Changes

Replace the deprecated `future::supportsMulticore()` re-export with
`parallelly::supportsMulticore()` and declare `parallelly` in Imports. The
deprecated call produces warnings with parallelly 1.50.0, causing test failures
in CRAN's R-devel checks. Tests have also been updated to use `parallelly`
directly. Thanks to Henrik Bengtsson for identifying and fixing the issue.

The release also clarifies search and rate-map documentation and worked
examples; it does not change the search algorithm.

## Checks

`R CMD check --as-cran --run-donttest` on the 0.2.3 source archive completed
on macOS arm64 with R 4.6.1 and parallelly 1.50.0:

**0 errors | 0 warnings | 0 notes**

The check included PDF and HTML manuals, examples, and tests. The test suite
reported 2,980 passing expectations, no failures or warnings, and 138 skips.

GitHub Actions checked the merged fix on Ubuntu R-devel 4.7.0 with parallelly
1.50.0: `R CMD check --as-cran --run-donttest --no-manual` reported
**Status: OK**, with 2,980 passing expectations, no test failures or warnings,
and 138 skips. That check used the PR source before the 0.2.3 release metadata,
not the final 0.2.3 source archive. Windows, macOS, Ubuntu release, and Ubuntu
oldrel-1 checks for the PR also passed. The CRAN package index listed no
reverse dependencies when checked.
