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
completed on macOS arm64 (R 4.4.2) on 2026-10-15 with:

**0 errors | 0 warnings | 1 note**

The sole note was `unable to verify current time`, an environment-related
timestamp check. Installation, examples and tests passed.

[CI checks](https://github.com/jakeberv/bifrost/actions/runs/38823234675)
also passed on macOS and Windows R release, and Ubuntu R release, devel and
oldrel-1 (`--as-cran --run-donttest --no-manual`).

## Additional notes

- Worked articles and empirical data are hosted separately. Data downloads
  occur only on explicit request, not during installation, attachment or checks.
- High-dimensional fits remain sensitive to initialization and random state;
  some intermediate fits can fail to converge (issue #242).
- The existing Gumbel KL-divergence bootstrap can understate uncertainty because
  of overflow/underflow. Its correction remains deferred to a later release:
  https://github.com/jakeberv/bifrost/issues/231.
