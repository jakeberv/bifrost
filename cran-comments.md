# Submission

This is the `bifrost` 0.2.0 release, updating CRAN version 0.1.4.

## Major changes

- Adds search-history inspection, branch-rate maps, lineage-rate and
  shift-distribution summaries, post-hoc regime covariance and integration
  analyses, and simulation/tuning workflows.
- Removes the superseded `plot_ic_acceptance_matrix()` function. The supported
  migration is `plot(icTrajectory(x))`.
- Raises the minimum R version from 4.1 to 4.2 and requires
  `phytools (>= 2.0-3)`. Rate-map plotting now declares `plotrix` directly and
  does not write compatibility functions into the user's global environment.

## Package contents and example data

The worked vignettes are website-only. Their sources, generated files, and
empirical payloads are excluded from the CRAN source package; package
installation, attachment, examples, and routine checks remain network-free.

`bifrost_example_file()` performs network access only when a user explicitly
requests an empirical artifact. The first uncached request downloads the
manifest and artifact currently tracked on GitHub `main`, verifies the recorded
SHA-256 checksum and byte size, and caches the verified result. Later calls use
that cache unless `refresh = TRUE` checks `main` again. A local verified mirror
can be selected with `BIFROST_ARTIFACT_DIR`, and network failures produce an
actionable error without corrupting a prior verified cache entry.

## URL note

The README intentionally uses <https://rdcu.be/fuxB1>. This is the publisher's
free-access SharedIt URL for the application paper. It redirects over HTTPS to
the corresponding Nature article PDF and returns HTTP status 200.

## Test environments

- Local: macOS Sequoia 15.7.8, aarch64, R 4.4.2
- GitHub Actions is configured for R-devel, release, and oldrel on Ubuntu, plus
  release R on macOS and Windows. The final remote results will be recorded
  after the local release-candidate review.

## R CMD check results

The local release-candidate check was run with
`R CMD check --as-cran --run-donttest`.

0 errors | 0 warnings | 2 notes

The notes are the intentional SharedIt redirect described above and an
environmental inability to verify the current time.

## Downstream dependencies

There are no reverse dependencies on CRAN.
