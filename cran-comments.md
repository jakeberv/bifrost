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
`30e58e82af46f6f0ef35ebf69f3a563f8b459774` (2026-10-15 UTC):

- [R CMD check](https://github.com/jakeberv/bifrost/actions/runs/38688925126)
  passed on macOS and Windows R release, and Ubuntu R release, devel, and
  oldrel-1, using `--as-cran --run-donttest --no-manual`.
- [Coverage and parallel smoke tests](https://github.com/jakeberv/bifrost/actions/runs/38688925143)
  passed. Whole-package R coverage was exactly 100%, with zero uncovered rows.
  The coverage run recorded 4,947 passing assertions, zero failures or warnings,
  and 43 expected skips for coverage-specific cases or development files absent
  from the built package. Separate Linux and Windows parallel smoke jobs passed.
- [The pkgdown build and deployment](https://github.com/jakeberv/bifrost/actions/runs/38688925130)
  passed.
- [The advisory audit](https://github.com/jakeberv/bifrost/actions/runs/38688925163)
  completed. It still reports print/cat usage and documentation/style findings;
  a successful advisory job does not mean every advisory finding is cleared.
  URL-liveness and advisory spelling checks were skipped.

## Final source-package check

On 2026-10-15, `bifrost_0.2.2.tar.gz` built from the same commit passed
`R CMD check --as-cran --run-donttest --library=<isolated library>` with:

**0 errors | 0 warnings | 1 note**

The sole note was `unable to verify current time`, from the local environment's
timestamp check. Installation, loading, examples (including `\donttest`),
tests, and both PDF and HTML manual checks passed. The test suite recorded
2,969 passing assertions, zero failures or warnings, and 138 expected skips
(including CRAN-mode skips and development files absent from the tarball).
The 101-page PDF manual also passed a separate `qpdf --check`.

Environment: macOS Sequoia 15.7.9, arm64, R 4.4.2, mvMORPH 1.2.3, and OpenBLAS
0.3.34. Bifrost was installed in an isolated library; other dependencies were
accessed read-only. Numerical libraries used one thread per worker, with at
most eight workers. Archive contents were inspected to confirm that local
outputs, caches, and development tooling were excluded.

Archive SHA-256:

```text
153dbcf17195ad2c0e727fda3acad1e47c91a4ae67f32f5a609b525b0dbfe92b
```

The first attempt passed tests and the PDF manual but crashed inside the
existing V8 8.2.0 installation during the HTML manual check. A standalone V8
call reproduced the crash without loading bifrost. A private build of V8 8.2.0
with bundled engine 14.6.202.26 resolved it. The complete check then passed on
the identical archive, without disabling checks or changing bifrost source.
Existing installations were unchanged during this validation.

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

Merge the check-evidence PR through normal review and confirm CI on its merge
commit. Select the release commit and verify that its package contents still
match the checked source; rerun the final check if package contents change.
These maintainer notes and `tools/releases/0.2.2.md` are excluded from source
archives. Keep the checksum above tied to the actual checked archive, not an
unverified rebuild. Remove draft/pending wording only when submission is
authorized. No tag, GitHub release, or CRAN submission has been made for 0.2.2.
