# CRAN submission preparation: bifrost 0.2.2

Draft for an update from CRAN version 0.2.0 to 0.2.2. Version 0.2.1 was
published on GitHub. Version 0.2.2 has not yet been submitted or published.

## Release changes

Searches now default to data-dependent BM/BMM starting values adapted from
the mvMORPH 1.2.1 development initializer used in earlier analyses. The
compatibility helper passes these values through `start` to the installed
mvMORPH; it does not pin or modify that package. Native initialization remains
available through `start_strategy = "native"` and is required for EmpBayes.
The adapted code credits Julien Clavel in its source and in `Authors@R`.
`corpcor` and `RhpcBLASctl` are now explicit imports.

This update also corrects progress-related random-stream and numerical-thread
handling and refines the console displays. It includes the citation and source
archive updates from the GitHub-only 0.2.1 release. NEWS.md describes the
changed defaults and reproducibility limitations.

## Package contents and network access

Worked articles and empirical data remain on the package website and are
excluded from the source package. `bifrost_example_file()` downloads data
only on explicit request, verifies SHA-256 and byte size, and reports download
failures. Installation, attachment, and examples do not require these downloads.

## Verified CI results

The following results apply to main commit
`204c3c5a7e71c6cb30a4414ef07e36fbd58561c7` (2026-10-15 UTC):

- [R CMD check](https://github.com/jakeberv/bifrost/actions/runs/38679236397)
  passed on macOS and Windows R release, and Ubuntu R release, devel, and
  oldrel-1, using `--as-cran --run-donttest --no-manual`.
- [Coverage and parallel smoke tests](https://github.com/jakeberv/bifrost/actions/runs/38679236413)
  passed. Whole-package R coverage was exactly 100%, with zero uncovered rows.
  The coverage run recorded 4,947 passing assertions, zero failures or warnings,
  and 43 expected skips for coverage-specific cases or development files absent
  from the built package. Separate Linux and Windows parallel smoke jobs passed.
- [The pkgdown build and deployment](https://github.com/jakeberv/bifrost/actions/runs/38679236325)
  passed.
- [The advisory audit](https://github.com/jakeberv/bifrost/actions/runs/38679236377)
  completed. It still reports print/cat usage and documentation/style findings;
  a successful advisory job does not mean every advisory finding is cleared.
  URL-liveness and advisory spelling checks were skipped.

## Known limitations retained in this update

- High-dimensional searches remain sensitive to initialization and random
  state, and some intermediate fits can fail to converge. The compatibility
  policy restores a reproducible starting strategy, not guaranteed optimality.
  Follow-up: https://github.com/jakeberv/bifrost/issues/242.
- The existing Gumbel KL-divergence bootstrap can discard valid large draws
  because of overflow/underflow and understate uncertainty. Its correction is
  deferred and is not part of 0.2.2 preparation.
  Follow-up: https://github.com/jakeberv/bifrost/issues/231.

## Pending before submission

The final 0.2.2 source archive has not yet received the local check including
its PDF manual. Do not reuse the September 0.2.0 check as evidence for this
submission. After the final preparation PR is merged, record the final commit,
archive checksum, actual check results and any explanations for notes here;
then remove this draft/pending wording before submission.
