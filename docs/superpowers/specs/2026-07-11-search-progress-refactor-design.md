# Search Progress Refactor Design

## Goal

Reduce the size and maintenance burden of the progress-bar upgrade while retaining its public API, terminal experience, parallel behavior, error handling, and 100% changed-executable-line coverage.

The refactor starts from commit `61618d`, which is the reviewed, passing implementation and rollback point.

## Non-negotiable behavior

- `progress = TRUE` remains default-on and keyword-only after `...`.
- Reached stage rows appear incrementally: candidate scoring first, greedy search second, and IC-weight re-estimation third.
- Completed rows remain bottom-pinned until `searchOptimalConfiguration()` exits.
- With `verbose = TRUE`, search messages stream above all active rows.
- Dynamic terminals retain continuously animated spinners during long fits.
- Non-dynamic or redirected output remains bounded rather than logging every heartbeat.
- Candidate scoring and parallel IC-weight re-estimation retain their existing Future parallelism; greedy search remains sequential.
- Future plans, thread variables, RNG state, and unresolved work are restored or cancelled on success, error, and interrupt.
- `progress = FALSE` retains the old low-overhead execution paths and quiet/verbose semantics.
- Captured expressions such as `shift_acceptance_threshold = -Inf` continue to print safely.
- Dual uncertainty-weight flags continue to fail before fitting or progress-row creation.

## Selected refactor

### Remove dead compatibility code

Delete `.bifrost_search_legacy_progress_handler()`. It has no production caller: the no-argument branch exists only for internal compatibility tests, while standalone stages already create an owned shared session.

Make `.bifrost_search_progress_handler()` require a session and label, or move its reporter construction directly into the session's `handler()` method. Remove tests and temporary coverage paths that exist only for the deleted legacy implementation.

Delete `.bifrost_search_future_value()`, which has no production caller. Keep `.bifrost_search_await_work()` because greedy search and serial IC-weight refits use it.

### Consolidate the progress UI state machine

Keep the renderer/session boundary because it isolates CLI side effects and enables deterministic behavioral tests. Reduce duplication inside it by:

- extracting one row-creation/update path shared by active and skipped stages;
- storing row state in one consistent record;
- constructing the stage-bound `progressr` reporter inside the session;
- keeping one idempotent finalizer;
- retaining explicit multi-row redraw after semantic CLI output because released CLI restores only its most recent active row.

Do not replace the session with manual ANSI cursor control or an unreleased CLI multiline feature.

### Simplify stage integration where it remains clear

Keep the three explicit stage bodies because candidate scoring, sequential greedy decisions, and weight refits have different outputs and failure semantics. Extract only small repeated mechanics—heartbeat message construction, completion summaries, or session delegation—when the resulting helper has one clear responsibility.

Do not hide the search algorithm inside a generic callback framework merely to reduce line count.

### Reduce test implementation coupling

Retain behavioral coverage for:

- progressive one-to-three row creation and pinned completion;
- verbose output above multiple active rows;
- dynamic heartbeat animation and bounded non-dynamic output;
- skip, failure, interrupt, and idempotent cleanup;
- disabled progress and positional `...` compatibility;
- serial/parallel tick counts, RNG stability, cancellation, and thread/Future restoration;
- public print behavior and all three public stage paths.

Remove tests that only inspect deleted private closure bodies or repeat the same lifecycle branch through multiple renderers. Prefer table-driven state transitions and public/injected-renderer assertions over deparsed function-body assertions.

Changed executable lines must still measure 100%; reduced production branching should reduce the test burden naturally.

### Remove process-only artifacts from the final branch diff

Delete the tracked files under `docs/superpowers/specs/` and `docs/superpowers/plans/` after implementation and verification. They are development-process records rather than package documentation and currently account for hundreds of branch-added lines. User-facing documentation remains in roxygen/Rd, `NEWS.md`, and the quick-start vignette.

## Size targets

These are cleanup targets, not reasons to obscure code:

- remove at least 250 net runtime R lines relative to commit `61618d`;
- remove at least 1,100 total branch-added lines relative to commit `61618d`;
- reduce `tests/testthat/test-search-progress.R` materially while retaining all behavioral and coverage gates;
- introduce no new dependency or public argument.

If a target conflicts with clarity or correctness, preserve clarity and document the shortfall.

## Verification

- Demonstrate RED before each behavior-affecting production change and GREEN afterward.
- Run focused progress, print, and public-search tests.
- Run the full package suite with zero failures, warnings, and skips.
- Measure current `main...HEAD` changed executable R lines with `covr`; require 100% and zero uncovered lines.
- Run `R CMD build` and `R CMD check --no-manual` from a clean direct source tree; require `Status: OK`.
- Run the external PTY example with two cores, `verbose = TRUE`, `progress = TRUE`, parallel weights, and `-Inf`; verify the progressive pinned stack and successful `print(result)`.
- Compare final diff statistics with commit `61618d` and report actual production/test/process-line reductions.

## Alternatives rejected

- Stock `progressr::handler_cli(clear = FALSE)` is much smaller but cannot retain and redraw all completed stage rows beneath later verbose output.
- One active current-stage bar would be smaller but violates the approved pinned-stack behavior.
- Dropping continuous animation would remove Future polling code but reverses an explicitly requested feature.
- Requiring an unreleased CLI multiline implementation is not suitable for a package dependency contract.
