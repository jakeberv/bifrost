# Bottom-Pinned Search Progress Design

## Goal

Keep every search progress row that has been reached in a bottom-pinned terminal status area while `searchOptimalConfiguration()` is running. Verbose output must stream above that status area. Also make `print.bifrost_search()` robust to expression-valued arguments captured by `match.call()`, including `shift_acceptance_threshold = -Inf`.

## User experience

Progress rows appear incrementally rather than all at once:

1. Candidate scoring creates the first row when that stage is reached.
2. The completed candidate row remains pinned, and greedy search creates a second row beneath it.
3. The completed candidate and greedy rows remain pinned, and IC-weight re-estimation creates a third row beneath them.

At any point, verbose output is inserted above all currently visible progress rows. A completed row shows a checkmark, a full bar, its final count, elapsed time, and summary while remaining part of the active status area. A skipped stage creates its row only when the stage is reached and immediately changes it to a durable skipped state.

The rows remain terminal-managed only for the duration of `searchOptimalConfiguration()`. On normal return they are finalized in creation order as three persistent console lines. Later output from the caller, such as `print(result)`, is ordinary output and may follow those finalized lines.

## Selected architecture

Use CLI's native active-progress output contract instead of manually moving the cursor. CLI documents that semantic output emitted while progress bars are active temporarily removes the bars, writes the output, and restores the bars at the bottom:

<https://cli.r-lib.org/reference/cli_progress_bar.html#cli-output-while-the-progress-bar-is-active>

A private, call-scoped search progress session owns every CLI bar created during one search. The session exposes operations to:

- create a stage row when that stage is reached;
- update completion and heartbeat state;
- mark a stage completed, skipped, or failed without terminating its CLI bar;
- emit verbose text through the active CLI progress handler;
- finalize all created bars at search exit.

Each stage keeps its existing `progressr` progressor. Future workers still emit progress conditions only; they never own terminal output. A stage-specific `progressr` handler forwards lifecycle events into the shared session instead of creating and terminating an independent CLI bar.

## Row state and rendering

The session stores rows in stage order and assigns one explicit CLI progress-bar id and environment to each row. Bars use `current = FALSE`, `auto_terminate = FALSE`, and explicit ids so multiple rows may remain active simultaneously.

One custom format selects its leading symbol and status from row state:

- active: animated spinner, percentage, current/total, ETA, and current status;
- complete: checkmark, 100%, final current/total, elapsed time, and completion summary;
- skipped: information symbol and the documented skip reason;
- failed: cross, partial percentage/count, elapsed time, and failure summary.

Completion changes row state and forces a redraw but does not call `cli_progress_done()`. This is what keeps earlier rows pinned while later rows run. The session finalizer calls `cli_progress_done()` for all created rows with persistent output enabled. Cleanup is idempotent so normal completion, errors, and interrupts cannot finalize a row twice.

## Verbose and plotting output

When progress is enabled and at least one stage row exists, the internal verbose logger sends text through `cli_progress_output()` for an active row. CLI then hides the complete status area, emits the text above it, and restores every active row. Before the first stage row exists, or when progress is disabled, the existing `message()` behavior remains.

The special interactive RStudio plotting path remains supported. Search-owned verbose text uses the session output path whenever a progress row is active. Direct plotting output is not redirected; the session's next forced redraw restores the status area. Non-dynamic terminals retain CLI's documented line-oriented fallback rather than receiving raw cursor-control sequences.

## Skips, failures, and restoration

- A zero-work stage creates a row only when reached and immediately marks it skipped.
- A recoverable greedy proposal error remains a normal stage update and increments the errored-proposal count.
- An uncaught stage error marks the current row failed, leaves previously completed rows intact, finalizes all rows that were actually created, restores Future and thread settings, and propagates the original condition.
- Stages not reached because an error aborts the search do not appear.
- `progress = FALSE` creates no session rows and preserves the current quiet and verbose-only behavior.

## Print-method fix

`user_input` remains the captured `match.call()` list for reproducibility. Consequently, an argument such as `-Inf` is stored as a language expression rather than a numeric scalar.

The print method must not evaluate or numerically coerce captured expressions. The search-block threshold uses the existing safe value formatter, which formats numeric scalars to the current precision and deparses language objects such as `-Inf`. This fixes the observed error without changing the result schema or evaluating arbitrary caller expressions.

## Compatibility and scope

- The public API is unchanged.
- Candidate scoring and parallel IC-weight re-estimation retain their existing parallel Future infrastructure.
- Greedy search remains sequential.
- Heartbeats continue to update the spinner without advancing completion.
- Completed rows remain colored according to CLI and terminal capabilities.
- No global CLI handler or option is permanently changed.
- `rateMap()` and its progress handling remain out of scope.

## Testing

Implementation follows red-green TDD.

- A print regression test constructs or runs a search with `shift_acceptance_threshold = -Inf`, verifies that printing does not error, and verifies that `Threshold: -Inf` is displayed.
- An injected renderer test proves rows are created only when their stage is reached and that completion does not terminate them.
- A state-order test proves the visible row sequence grows from one to two to three.
- A verbose-output test proves output is sent through the active CLI progress handler while prior completed rows remain registered.
- Skip and failure tests verify durable row state, idempotent cleanup, error propagation, and restoration of Future/thread settings.
- A real CLI integration test verifies verbose lines render above the active stack in a dynamic terminal and that final rows persist.
- Existing serial/parallel tick-count, heartbeat, RNG, plotting, quiet-mode, full-suite, coverage, build, and `R CMD check` verification remain required.

## Alternatives rejected

- Manual ANSI cursor management is too fragile across terminal widths, RStudio, redirected output, and interrupts.
- Creating all three rows at search start violates the desired incremental appearance.
- Finalizing each row as soon as its stage completes turns it into ordinary transcript output, allowing later verbose text to appear below it.
- A single multiline composite bar fights CLI's per-bar lifecycle and makes stage-specific ETA and failure handling harder.
