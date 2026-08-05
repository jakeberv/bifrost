# Resubmission

This is a resubmission of the `bifrost` 0.2.0 release, updating CRAN version
0.1.4.

## Changes in this resubmission

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

- macOS Sequoia 15.7.8, aarch64, R 4.4.2

## R CMD check results

`R CMD check --as-cran --run-donttest`:

0 errors | 0 warnings | 1 note

The note is the local environment's inability to verify the current time.
