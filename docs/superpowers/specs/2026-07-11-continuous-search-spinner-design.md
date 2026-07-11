# Continuous Search Spinner Design

## Goal

Keep the spinner in each `searchOptimalConfiguration()` progress line moving while a long model fit is running, without advancing the progress count until that fit actually completes.

## Selected approach

Expensive fits run in Future workers while the main R process polls their state every 100 milliseconds. Each poll emits a zero-increment `progressr` event. The CLI handler redraws the spinner when its own display timer permits, while ordinary completion events remain the only events that advance `current/total`, percentage, and ETA.

Candidate scoring and parallel IC-weight re-estimation use a small number of chunk Futures so the large tree and trait objects are serialized once per worker rather than once per candidate. Greedy search remains algorithmically sequential: it launches exactly one proposal fit, waits for it with heartbeats, applies the accept/reject decision in the main process, and only then launches the next proposal.

## Concurrency and compatibility

- `num_cores` continues to cap simultaneous expensive fits. Candidate scoring and parallel IC weights may run up to `num_cores` fits concurrently; greedy search and serial IC weights run one fit at a time.
- With progress enabled and `num_cores = 1`, one model fit runs in a background Future so the main process remains available to animate the terminal. The second process is the main-process renderer, not an additional model fit.
- With progress disabled, existing serial and `future.apply` paths remain unchanged and no heartbeat polling overhead is introduced.
- Unix terminals outside RStudio retain the multicore backend. RStudio and non-Unix platforms retain multisession.
- Workers emit conditions only. All CLI rendering, plotting, verbose output, and search-state mutation remain in the main process.

## Internal interfaces

The search helpers gain three private units:

1. A Future-plan scope that caps BLAS/OpenMP threads, selects the existing backend, and restores Future and thread state in `finally`.
2. A Future wait loop that polls without blocking, invokes a supplied heartbeat callback at a 100 ms cadence, returns the Future value, and cancels outstanding work on errors or interrupts.
3. A chunked Future map that preserves input/result order and relays worker progress conditions while the main process emits heartbeat events.

Stage work derives a heartbeat callback from the existing progressor by calling it with `amount = 0`. No public heartbeat frequency or style argument is added.

## Failure behavior

An error returned by a worker propagates through the same existing stage error path. Recoverable greedy fit errors still count as one errored proposal. Uncaught errors and interrupts cancel unresolved Futures, restore Future and thread state, leave the partial failed progress line visible, and propagate unchanged.

## Testing

- A focused asynchronous wait test proves that several zero-increment heartbeats occur during a deliberately slow Future and that the returned value is unchanged.
- Progress collector tests distinguish heartbeat amounts from completion amounts and prove that heartbeats do not advance a stage.
- Candidate, greedy, and IC-weight tests prove that completed-fit counts and ordering remain unchanged.
- Error/interrupt tests prove outstanding work is cancelled and Future/thread settings are restored.
- A CLI integration test confirms the spinner receives repeated redraw opportunities without creating additional persistent lines.
- The full test suite, package build, and `R CMD check --no-manual` remain required.

## Alternatives rejected

- A worker-owned terminal spinner risks interleaved or corrupted output and violates the existing main-process rendering contract.
- A dedicated heartbeat worker would consume one of the requested compute workers and reduce model-fitting parallelism.
- `later` callbacks cannot execute while the main R process is blocked inside a synchronous model fit.
- A subprocess that writes ANSI sequences directly would be fragile across terminals and platforms.
