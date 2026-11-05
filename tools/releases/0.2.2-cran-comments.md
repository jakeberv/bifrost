# bifrost 0.2.2

This updates CRAN version 0.2.0; version 0.2.1 was released on GitHub only.

## Changes

This release restores the data-dependent starting-value strategy used in
earlier analyses, while retaining an option for native mvMORPH initialization.
It also fixes progress-related random-state and numerical-thread handling.
The adapted starting-value code credits Julien Clavel. New imports are
`corpcor` for initialization and `RhpcBLASctl` for worker thread control.
Further changes and reproducibility limitations are described in NEWS.md.

## Checks

`R CMD check --as-cran --run-donttest`, including PDF and HTML manuals,
completed on macOS arm64 (R 4.4.2) on 2026-10-17 UTC with:

**0 errors | 0 warnings | 1 note**

The sole note was `unable to verify current time`, an environment-related
timestamp check. Installation, examples and tests passed, with 2,980 passing
test assertions and no test warnings or failures. The RNG chunk-invariance
test now exercises two and three chunks using at most two worker processes.

[CI](https://github.com/jakeberv/bifrost/pull/255/checks)
also checks macOS and Windows R release, and Ubuntu R release, devel and
oldrel-1 (`--as-cran --run-donttest --no-manual`). A Windows temporary-file
note was traced to the worker-cancellation test; its worker files now use a
test-owned temporary directory with scoped cleanup.

## Additional notes

- Worked articles and empirical data are hosted separately. Data downloads
  occur only on explicit request, not during installation, attachment or checks.
- High-dimensional fits remain sensitive to initialization and random state;
  some intermediate fits can fail to converge (issue #242).
- The existing Gumbel KL-divergence bootstrap can understate uncertainty because
  of overflow/underflow. Its correction remains deferred to a later release:
  https://github.com/jakeberv/bifrost/issues/231.
