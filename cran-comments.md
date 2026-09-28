# CRAN submission

This submission updates bifrost from CRAN version 0.1.4 to 0.2.0.

## CRAN check failure

This release addresses the test failure reported with mvMORPH 1.2.2, with a
correction deadline of 2026-10-16. The old method-forwarding test used a BMM
fit with only two tips in one regime and produced a nonfinite GIC. It has
been replaced with balanced examples that verify the requested fitting method
and a finite GIC. These tests pass with CRAN mvMORPH 1.2.2.

The redirecting SharedIt URL flagged by incoming pretests has also been
replaced in README.md with the publisher's direct SharedIt ePDF URL.

## Release changes

This release adds search-history inspection, branch- and lineage-rate
summaries, post-hoc regime analyses, and simulation and tuning workflows.
It also fixes formula-intercept handling, inherited tuning settings, F1
calculations, and covariance validation. NEWS.md describes the changes,
including migration from `plot_ic_acceptance_matrix()` to
`plot(icTrajectory(x))`.

## Package contents and network access

Worked articles and empirical data are distributed through the package website
and excluded from the source package. `bifrost_example_file()` downloads data
only on explicit request, verifies SHA-256 and byte size, and reports download
failures. Package installation, attachment, and examples do not require these
downloads.

## Test environments and results

- Local: macOS Sequoia 15.7.9, aarch64, R 4.4.2. Compatibility was checked
  with CRAN mvMORPH 1.2.2.
- GitHub Actions: macOS and Windows R release; Ubuntu R release, devel, and
  oldrel-1. All five package-check jobs passed for the merged release changes.

Additional CI checks passed for test coverage (100%, zero uncovered lines),
parallel-worker smoke tests on Linux and Windows, vignette artifacts, and the
pkgdown website build. The advisory checktor audit also completed; it is a
supplementary diagnostic rather than a substitute for R CMD check.
CI results: https://github.com/jakeberv/bifrost/pull/232/checks

The final local `devtools::check()` with CRAN incoming checks, the manual, and
`run_dont_test = TRUE` completed on 2026-09-27 with:

0 errors | 0 warnings | 1 note

The sole note is the local environment's inability to verify the current time.
Examples (including `\donttest`), tests, and the manual passed.
