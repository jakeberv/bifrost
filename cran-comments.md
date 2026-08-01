# Submission

This is the `bifrost` 0.2.0 release, updating CRAN version 0.1.4.

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

## URL note

The publisher's free-access SharedIt URL <https://rdcu.be/fuxB1> redirects over
HTTPS to the corresponding article PDF and returns HTTP status 200.

## Test environments

- macOS Sequoia 15.7.8, aarch64, R 4.4.2

## R CMD check results

`R CMD check --as-cran --run-donttest`:

0 errors | 0 warnings | 2 notes

The notes are the SharedIt redirect described above and the local environment's
inability to verify the current time.
