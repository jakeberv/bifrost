# Resubmission

This is a resubmission of the `bifrost` 0.2.0 release, updating CRAN version
0.1.4.

## Changes in this resubmission

This release addresses the test failure reported for CRAN bifrost 0.1.4
with mvMORPH 1.2.2, with a correction deadline of 2026-10-16. The old
method-forwarding test used a BMM fit with only two tips in one regime and
produced a nonfinite GIC. The mvMORPH maintainer has confirmed an unintended
change in its starting-value calculation and is preparing an upstream fix.

The development version of bifrost had already replaced that test. Its
current method-forwarding tests now use balanced, better-supported examples
and check both the requested fitting method and a finite GIC. No workaround
for mvMORPH's starting-value calculation is included in bifrost.

The permanently redirecting `rdcu.be` SharedIt URL reported by the incoming
pretests has been replaced in `README.md` by the publisher's canonical,
tokenized SharedIt ePDF URL. This preserves free access to the article while
avoiding the permanent redirect.

## Major changes

This release adds search-history inspection, branch- and lineage-rate
summaries, post-hoc regime analyses, and simulation and tuning workflows. It
removes `plot_ic_acceptance_matrix()`; the supported migration is
`plot(icTrajectory(x))`. There are no reverse dependencies on CRAN.

## Package contents and network access

Worked articles and empirical data are website-only and excluded from the
source package. `bifrost_example_file()` accesses the network only after an
explicit user request, verifies artifact SHA-256 and byte size, and fails
gracefully. Installation, attachment, examples, and checks remain network-free.

## Test environments

- macOS Sequoia 15.7.9, aarch64, R 4.4.2, checked on 2026-09-26
  using mvMORPH 1.2.2 installed from the official CRAN source tarball.
- GitHub Actions for main before the final CRAN-preparation fixes: macOS and Windows R release;
  Ubuntu R release, devel, and oldrel-1. All five jobs passed. These CI runs
  preceded the CRAN publication of mvMORPH 1.2.2; the local check above used
  that published version explicitly.

## Baseline R CMD check results

These results precede the final formula-intercept, tuning-recommendation, and
F1 calculation and covariance-validation fixes. Rebuild and check the final release tarball before submission.

`R CMD check --as-cran --run-donttest`:

0 errors | 0 warnings | 1 note

The note is the local environment's inability to verify the current time.
Examples (including `\donttest`), tests, and the PDF manual all passed.
The CRAN-mode test suite reported 2393 passes, no failures or warnings, and
141 skips for CRAN-excluded tests and unavailable development-only, optional,
or website-only inputs. The revised method-forwarding tests passed.

## Checks after the final preparation fixes

On 2026-09-27, the development test suite with CRAN mvMORPH 1.2.2 and
local website artifacts reports 4798 passing expectations across 539 tests,
with no failures, warnings, or skips. The
`runSearchTuningGrid()` and `selectTunedSearchParameters()` help examples pass.
The formula-intercept, tuning-recommendation, F1, and covariance-validation
fixes have regression tests.
The supplementary replicate metrics were corrected from saved counts; all
non-F1 data and the pooled vignette cache are unchanged. Other findings from
the release review remain to be resolved before submission.
