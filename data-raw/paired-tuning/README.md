# Paired tuning-grid reproduction scripts

This directory publishes the execution sources used for the 12-dimensional,
10--40-tip planted-clade paired tuning grid in the Part 2 simulation-study vignette. They are
reproduction scripts, not a new `bifrost` package API. The public API tutorial
shows the study design, but it does not itself implement checkpoint recovery
after a simulation generator failure.

## Recorded source identity

The three files here are byte-for-byte copies of the original local campaign
sources. The campaign used `bifrost` package commit
`db18184ddb06a5019123647ce218f9b717f76e49` (package version 0.2.0).

| Published file | SHA-256 |
|---|---|
| `run-vignette-tuning-grid.R` | `35c6a555b65476ea7890f3b1a9be82b55c9b4b409fdbcaf9a8973d5b37f1591a` |
| `run-manuscript-simulation-suite.R` | `65a2de0f0d63363e9879f17cad97eeb9c4384dc9b8145ce3be8d879660313f2d` |
| `summarize-available-tuning-grid.R` | `580e65d28b20736926c431a7883f2637e04fb5cca13e63ab590925e43b37465f` |

`run-manuscript-simulation-suite.R` is a shared implementation dependency and
must remain beside `run-vignette-tuning-grid.R`. The node-specific launcher and
the original campaign log are intentionally not published.

## Package prerequisite

Install the recorded `bifrost` commit into the R library that will be active
when the scripts run. For example, with `remotes` already installed:

```r
remotes::install_github(
  "jakeberv/bifrost@db18184ddb06a5019123647ce218f9b717f76e49"
)
```

The runners load `bifrost` from the active R library. `--repo=PATH` records
checkout provenance and controls default output placement; it does not load or
install the package from that checkout. Confirm the installed source identity
before a reproduction run, for example with
`packageDescription("bifrost")$RemoteSha`.

## Full design

The full grid uses a 12-trait empirical covariance template and three 250-tip
scenarios: null, five proportional-rate shifts, and five integration-rate
(correlation) shifts. Planted clades span 10--40 descendant tips. Each scenario
has 500 attempted replicates. Within a
replicate, the same generated dataset is searched under all 12 combinations of
information criterion (GIC or BIC), shift-acceptance threshold (10, 20, or 30),
and minimum descendant-tip count (10 or 20). Thus the complete design contains
1,500 replicate jobs and 18,000 searches. Searches are response-only fits with
`error = FALSE`; the empirical calibration template uses `error = TRUE`.

Inspect the design without loading the package, generating data, or fitting a
model:

```sh
Rscript data-raw/paired-tuning/run-vignette-tuning-grid.R --dry-run --mode=full --repo=.
```

The runner uses the installed `bifrost` package from the active R library. By
default it obtains the example tree and traits through
`bifrost_example_file()`. For an offline copy, add `--data-dir=/path/to/data`;
that directory must contain the expected example tree and trait files.

## Run, resume, and checkpoints

Execution is opt-in. From a package checkout, a portable full-run command is:

```sh
Rscript data-raw/paired-tuning/run-vignette-tuning-grid.R \
  --execute --mode=full --repo=. --cores=8
```

`--cores=N` optionally runs multiple replicate jobs on one machine. Each worker
is single-core and processes the 12 searches for its assigned replicate
serially. Choose a worker count appropriate for available memory and physical
cores. No scheduler is required. `--only=GROUP_ID` and `--job-index=N` can
restrict execution; `--list` and `--list-jobs` show valid selections.

The default checkout-relative output root is
`local-cache/vignette-tuning-grid-12d-clades10-40`. Each replicate has a saved simulation
checkpoint and one checkpoint per completed search. Re-running the same command
validates and resumes compatible results. Existing final replicate files are
validated and skipped. There is deliberately no overwrite option: use a new
`--output-dir` if package, data, design, or runner identity changes. A single-replicate
smoke mode exists, writes separately, and is not evidence for the
full study:

```sh
Rscript data-raw/paired-tuning/run-vignette-tuning-grid.R \
  --execute --mode=smoke --repo=. --cores=1
```

For a fully completed campaign, rebuild summaries without simulation or model
fitting with:

```sh
Rscript data-raw/paired-tuning/run-vignette-tuning-grid.R \
  --summarize-only --mode=full --repo=.
```

## Summarize retained results after generator failures

The separate summarizer is only for a campaign situation in which every
missing replicate is accounted for by a known shift-placement generator
failure. Supply the output root and that run's original stderr log:

```sh
Rscript data-raw/paired-tuning/summarize-available-tuning-grid.R \
  /path/to/vignette-tuning-grid-12d-clades10-40 /path/to/original-stderr.log
```

Unlike the execution runner, this summarizer has no `--data-dir` parameter. It
re-resolves and verifies the original example inputs through the installed
package's `bifrost_example_file()` resolver. Matching original inputs must be
available either through a valid `BIFROST_ARTIFACT_DIR` artifact directory or
through the verified downloader cache before offline summarization. Running the
campaign with only `--data-dir` does not by itself prepare the retained-results
summarizer for offline use.

It accepts only the expected “failed to place all requested shifts after 100
attempts” records, rejects missing, duplicate, contradictory, or unknown job
records, validates saved provenance and checkpoints, and writes a new
`full/summaries-available` directory without modifying raw campaign files. Its
search-performance denominators are successfully generated datasets;
generation success is reported against all attempted datasets. Because the
campaign log is machine/run-specific, users reproduce this route with the log
from their own run rather than a repository copy.
