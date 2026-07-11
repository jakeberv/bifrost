# Bottom-Pinned Search Progress Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Keep each reached search-stage row pinned at the bottom while verbose output streams above the progressively growing stack, and make `print.bifrost_search()` handle captured expressions such as `-Inf`.

**Architecture:** A call-scoped progress session owns all CLI bars created during one search. Stage-specific `progressr` handlers update rows in that shared session but do not terminate completed rows; the session sends verbose text through `cli_progress_output()` and finalizes all rows when the public search exits. The print fix formats captured call expressions without evaluating or coercing them.

**Tech Stack:** R 4.1+, `cli` 3.6+, `progressr`, `future`, `future.apply`, `testthat`, `covr`.

## Global Constraints

- Stage rows appear only when their stage is reached; the visible stack grows from one row to two rows to three rows.
- Completed rows remain active and pinned until `searchOptimalConfiguration()` exits.
- Future workers continue to signal conditions only; all terminal output remains in the main R process.
- `progress = FALSE` retains the current silent and verbose-only behavior without CLI-session overhead.
- The public API and `bifrost_search` result schema remain unchanged.
- Uncaught errors and interrupts propagate unchanged after progress, Future, RNG, and thread cleanup.
- Changed executable R lines must retain 100% line coverage.

---

### Task 1: Format captured threshold expressions safely

**Files:**
- Modify: `tests/testthat/test-search-progress.R`
- Modify: `R/bifrost_search-methods.R:213-233`

**Interfaces:**
- Consumes: `.bifrost_fmt_val(v, digits = 2, cutoff = 80)`, which formats numeric values and deparses language objects without evaluation.
- Produces: `print.bifrost_search()` output containing `Threshold: -Inf` for a result whose captured call contains `shift_acceptance_threshold = -Inf`.

- [ ] **Step 1: Write the failing regression test**

Add a dependency-free test that stores the same language object produced by `match.call()`:

```r
test_that("print.bifrost_search formats captured threshold expressions", {
  object <- structure(
    list(
      user_input = list(
        shift_acceptance_threshold = quote(-Inf),
        progress = TRUE
      ),
      shift_nodes_no_uncertainty = integer(),
      warnings = character()
    ),
    class = c("bifrost_search", "list")
  )

  testthat::expect_no_error(
    output <- paste(testthat::capture_output(print(object)), collapse = "\n")
  )
  testthat::expect_match(output, "Threshold: -Inf", fixed = TRUE)
})
```

- [ ] **Step 2: Run the regression test and verify RED**

Run:

```bash
Rscript -e 'devtools::test(filter = "search-progress")'
```

Expected: FAIL with `"'language' object cannot be coerced to type 'double'"` from `.bifrost_print_search_block()`.

- [ ] **Step 3: Implement the minimal safe formatter change**

Replace numeric coercion of `d$thr` with the existing value formatter:

```r
.bifrost_kv(
  "Threshold",
  if (!is.null(d$thr)) .bifrost_fmt_val(d$thr, digits = 2) else NULL
)
```

Do not evaluate the captured expression and do not change `user_input`.

- [ ] **Step 4: Run the regression and print suites and verify GREEN**

Run:

```bash
Rscript -e 'devtools::test(filter = "print-bifrost_search|search-progress")'
```

Expected: all assertions pass with zero failures and warnings.

- [ ] **Step 5: Commit the print fix**

```bash
git add R/bifrost_search-methods.R tests/testthat/test-search-progress.R
git commit -m "Fix printing captured search thresholds"
```

---

### Task 2: Add the shared progressive CLI status session

**Files:**
- Modify: `R/searchOptimalConfiguration-helpers.R:1-230`
- Modify: `tests/testthat/test-search-progress.R`

**Interfaces:**
- Produces: `.bifrost_search_cli_renderer()` with `create()`, `update()`, `output()`, and `done()` operations.
- Produces: `.bifrost_search_progress_session(enabled, renderer)` with `handler(label)`, `update_stage(label, config, state, progression, row_state, status = NULL)`, `skip(label, reason)`, `output(text)`, `has_rows()`, `rows()`, and `finalize()` operations.
- Produces: `.bifrost_search_progress_handler(session, label)`, a `progressr` handler bound to one stage row.
- Consumes: existing progression state (`config$max_steps`, `state$step`, `state$message`, and `progression$amount`).

- [ ] **Step 1: Add an injectable renderer and failing progressive-stack test**

Define this test-only renderer in `test-search-progress.R`:

```r
.recording_search_renderer <- function(events) {
  list(
    create = function(label, total) {
      id <- paste0("row-", length(events$created) + 1L)
      events$created <- c(events$created, label)
      list(id = id, label = label, total = total)
    },
    update = function(row, state, current, total, status, force = TRUE) {
      events$updates <- c(events$updates, list(list(
        id = row$id,
        state = state,
        current = current,
        total = total,
        status = status
      )))
      invisible(NULL)
    },
    output = function(row, text) {
      events$output <- c(events$output, text)
      invisible(NULL)
    },
    done = function(row, result) {
      events$done <- c(events$done, row$id)
      invisible(NULL)
    }
  )
}

test_that("search status rows appear progressively and finalize together", {
  events <- new.env(parent = emptyenv())
  events$created <- character()
  events$updates <- list()
  events$output <- character()
  events$done <- character()
  session <- .bifrost_search_progress_session(
    TRUE,
    .recording_search_renderer(events)
  )

  run_stage <- function(label) {
    .bifrost_search_run_stage(
      enabled = TRUE,
      steps = 1L,
      initial_message = label,
      work = function(tick) {
        tick(message = paste(label, "complete"))
        list(value = label, done = paste(label, "done"))
      },
      handler = session$handler(label)
    )
  }

  run_stage("[1/3] Candidate scoring")
  testthat::expect_length(events$created, 1L)
  testthat::expect_length(events$done, 0L)
  session$output("verbose between stages")

  run_stage("[2/3] Greedy shift search")
  testthat::expect_length(events$created, 2L)
  testthat::expect_length(events$done, 0L)

  session$skip("[3/3] IC-weight re-estimation", "not requested")
  testthat::expect_length(events$created, 3L)
  testthat::expect_equal(tail(events$updates, 1L)[[1L]]$state, "skipped")

  session$finalize()
  session$finalize()
  testthat::expect_identical(events$done, paste0("row-", 1:3))
  testthat::expect_identical(events$output, "verbose between stages")
})
```

- [ ] **Step 2: Run the progressive-stack test and verify RED**

Run:

```bash
Rscript -e 'devtools::test(filter = "search-progress")'
```

Expected: FAIL because `.bifrost_search_progress_session()` does not exist.

- [ ] **Step 3: Implement the CLI renderer**

Create a renderer whose row environment stores `stage_state`, `status`, and `output_text`. Its `create()` operation calls:

```r
id <- cli::cli_progress_bar(
  total = total,
  format = format,
  format_done = format,
  format_failed = format,
  clear = FALSE,
  current = FALSE,
  auto_terminate = FALSE,
  .auto_close = FALSE,
  .envir = row_env
)
```

Use one conditional custom format:

```r
format <- paste0(
  "{if (stage_state == 'active') cli::pb_spin else ",
  "if (stage_state == 'failed') cli::symbol$cross else ",
  "if (stage_state == 'skipped') cli::symbol$info else cli::symbol$tick} ",
  "{cli::pb_bar} {cli::pb_percent} ",
  "[{cli::pb_current}/{cli::pb_total}] ",
  "{if (stage_state == 'active') paste0('(ETA: ', cli::pb_eta, ')') ",
  "else paste0('[', cli::pb_elapsed, ']')} {status}"
)
```

`update()` mutates the row environment and calls `cli_progress_update(id = row$id, set = current, total = total, status = status, force = force, .envir = row$envir)`. `output()` stores `output_text` and calls `cli_progress_output("{output_text}", id = row$id, .envir = row$envir)`. `done()` calls `cli_progress_done()` with `result = "failed"` only for failed rows and `"done"` otherwise.

- [ ] **Step 4: Implement the session state machine and bound handler**

Implement the session as a closure with an ordered `rows` list and an idempotent `finalized` flag. The stage handler must:

```r
reporter <- list(
  reset = function(...) invisible(NULL),
  hide = function(...) invisible(NULL),
  unhide = function(...) invisible(NULL),
  initiate = function(config, state, progression, ...) {
    session$update_stage(label, config, state, progression, "active")
  },
  update = function(config, state, progression, ...) {
    session$update_stage(label, config, state, progression, "active")
  },
  finish = function(config, state, progression, ...) {
    session$update_stage(label, config, state, progression, "complete")
  },
  interrupt = function(config, state, progression, ...) {
    session$update_stage(
      label,
      config,
      state,
      progression,
      "failed",
      status = conditionMessage(progression)
    )
  }
)
```

Retain the sentinel adjustment `display_total <- config$max_steps - 1L` and clamp `display_step <- min(state$step, display_total)`. Heartbeats with `amount = 0` force redraws without changing the displayed step. `finish` changes row state but never calls `renderer$done()`.

- [ ] **Step 5: Add failing skip, failure, output, and finalization tests**

Add focused assertions that:

```r
testthat::expect_identical(session$rows(), c(
  "[1/3] Candidate scoring",
  "[2/3] Greedy shift search",
  "[3/3] IC-weight re-estimation"
))
testthat::expect_length(events$done, 0L)      # before finalize
testthat::expect_length(events$done, 3L)      # after finalize twice
testthat::expect_true(any(vapply(
  events$updates,
  function(event) identical(event$state, "failed"),
  logical(1)
)))
```

Run the failure case through `.bifrost_search_run_stage()` and verify the original synthetic error is rethrown after the row is marked failed.

- [ ] **Step 6: Run session tests and verify GREEN**

Run:

```bash
Rscript -e 'devtools::test(filter = "search-progress")'
```

Expected: all session, heartbeat, Future, RNG, lifecycle, and weight assertions pass.

- [ ] **Step 7: Commit the shared status session**

```bash
git add R/searchOptimalConfiguration-helpers.R tests/testthat/test-search-progress.R
git commit -m "Add shared search progress status session"
```

---

### Task 3: Integrate the session with the public search and verbose logger

**Files:**
- Modify: `R/searchOptimalConfiguration.R:333-661`
- Modify: `R/searchOptimalConfiguration-helpers.R:185-230`
- Modify: `tests/testthat/test-search-progress.R`
- Modify: `tests/testthat/test-search-progress-silent.R`

**Interfaces:**
- Consumes: `.bifrost_search_progress_session(enabled)` from Task 2.
- Changes: `.bifrost_search_run_stage(..., session = NULL)` selects `session$handler(initial_message)` when a session is supplied.
- Changes: `.bifrost_search_report_skipped_stage(..., session = NULL)` delegates skipped-row rendering to the session.
- Produces: one call-scoped session used by all three public search stages and the verbose logger.

- [ ] **Step 1: Write failing public integration tests**

Inject a recording session constructor with `testthat::local_mocked_bindings()` and run the existing small deterministic public-search fixtures. Assert:

```r
testthat::expect_identical(events$created, c(
  "[1/3] Candidate scoring",
  "[2/3] Greedy shift search",
  "[3/3] IC-weight re-estimation"
))
testthat::expect_true(any(grepl(
  "Sorting and evaluating shifts",
  events$output,
  fixed = TRUE
)))
testthat::expect_identical(events$done, paste0("row-", 1:3))
```

Retain the quiet test proving `progress = FALSE, verbose = FALSE` emits no messages and creates no rows.

- [ ] **Step 2: Run public integration tests and verify RED**

Run:

```bash
Rscript -e 'devtools::test(filter = "search-progress")'
```

Expected: FAIL because the public search does not create or pass a shared session and verbose output is not routed through it.

- [ ] **Step 3: Create and finalize the session in `searchOptimalConfiguration()`**

Immediately after verbose option setup, add:

```r
progress_session <- .bifrost_search_progress_session(isTRUE(progress))
on.exit(progress_session$finalize(), add = TRUE)
```

Pass `session = progress_session` to all three `.bifrost_search_run_stage()` calls and all `.bifrost_search_report_skipped_stage()` calls.

- [ ] **Step 4: Route verbose output through CLI while rows are active**

Change the internal logger ordering to:

```r
if (progress_session$has_rows()) {
  progress_session$output(txt)
} else if (isTRUE(plot) && interactive() && identical(Sys.getenv("RSTUDIO"), "1")) {
  cat(txt, "\n", sep = "")
  if (sink.number(type = "output") == 0) utils::flush.console()
} else {
  message(txt)
}
```

This preserves pre-stage and progress-disabled message behavior while using CLI's active-bar output contract once a row exists.

- [ ] **Step 5: Update stage and skipped-stage helpers**

Extend `.bifrost_search_run_stage()` with `session = NULL` and choose its handler lazily:

```r
owns_session <- is.null(session) && is.null(handler)
if (owns_session) {
  session <- .bifrost_search_progress_session(TRUE)
  on.exit(session$finalize(), add = TRUE)
}
if (is.null(handler)) {
  handler <- session$handler(initial_message)
}
```

Use a `handler = NULL` default so `progress = FALSE` never constructs a handler. Extend `.bifrost_search_report_skipped_stage()` so a supplied session calls `session$skip(label, reason)`; retain the standalone CLI alert fallback for existing internal tests.

- [ ] **Step 6: Add a real CLI redraw integration test**

With `options(cli.dynamic = TRUE)` and the real renderer, capture message output from two one-step stages separated by `session$output("verbose detail")`, then finalize. Strip ANSI styling and assert the final occurrences of both stage labels appear after `verbose detail`, proving CLI redrew the active stack below the output:

```r
plain <- cli::ansi_strip(paste(rendered, collapse = "\n"))
verbose_pos <- max(gregexpr("verbose detail", plain, fixed = TRUE)[[1L]])
stage1_pos <- max(gregexpr("[1/3] Candidate scoring", plain, fixed = TRUE)[[1L]])
stage2_pos <- max(gregexpr("[2/3] Greedy shift search", plain, fixed = TRUE)[[1L]])
testthat::expect_gt(stage1_pos, verbose_pos)
testthat::expect_gt(stage2_pos, verbose_pos)
```

- [ ] **Step 7: Run focused tests and verify GREEN**

Run:

```bash
Rscript -e 'devtools::test(filter = "search-progress|print-bifrost_search")'
```

Expected: all assertions pass with zero failures and warnings.

- [ ] **Step 8: Commit public integration**

```bash
git add R/searchOptimalConfiguration.R R/searchOptimalConfiguration-helpers.R tests/testthat/test-search-progress.R tests/testthat/test-search-progress-silent.R
git commit -m "Pin search progress below verbose output"
```

---

### Task 4: Document and verify the completed behavior

**Files:**
- Modify: `NEWS.md:3-8`
- Modify: `vignettes/quick-start-vignette.Rmd:116-124`
- Test: `tests/testthat/test-search-progress.R`
- Test: `tests/testthat/test-print-bifrost_search.R`

**Interfaces:**
- Documents the unchanged public arguments `progress` and `verbose`.
- Verifies the complete package and 100% branch-changed executable-line coverage.

- [ ] **Step 1: Update user-facing documentation**

Add to `NEWS.md` that reached stage rows remain bottom-pinned while verbose output streams above them, and that captured expression-valued inputs print safely. Replace the quick-start progress description with:

```markdown
- `progress`: display search-stage rows as each stage is reached. Completed rows
  remain stacked at the bottom while `verbose = TRUE` output streams above them;
  spinners animate during long fits, while counts and ETA advance only when fits
  complete. Set it to `FALSE` for quiet batch or document builds.
```

- [ ] **Step 2: Run formatting and focused verification**

Run:

```bash
git diff --check
Rscript -e 'devtools::test(filter = "search-progress|print-bifrost_search")'
```

Expected: no whitespace errors; all focused assertions pass.

- [ ] **Step 3: Run the full package suite**

Run:

```bash
Rscript -e 'devtools::test()'
```

Expected: zero failures, warnings, and skips.

- [ ] **Step 4: Verify changed-line coverage**

Run the existing `covr` installed-package trace plus focused progress/print traces, intersect `covr:::tally_coverage(..., by = "line")` with `git diff --unified=0 main...HEAD -- R`, and require:

```text
Changed executable R coverage: 100.00%
Uncovered changed executable R lines: 0
```

- [ ] **Step 5: Build and check the source package**

Run from a temporary directory:

```bash
R CMD build /Users/cotinga/jacob.berv@gmail.com/Code/bifrost
R CMD check --no-manual bifrost_0.1.4.tar.gz
```

Expected: `Status: OK`.

- [ ] **Step 6: Run an external-terminal demonstration**

Run the 14-tip, two-core example with `verbose = TRUE`, `progress = TRUE`, `shift_acceptance_threshold = -Inf`, and `uncertaintyweights_par = TRUE`. Confirm visually that rows appear one at a time, the stack grows to three, verbose text stays above it, all bars retain color, and `print(result)` completes with `Threshold: -Inf`.

- [ ] **Step 7: Commit documentation and verification-facing changes**

```bash
git add NEWS.md vignettes/quick-start-vignette.Rmd
git commit -m "Document bottom-pinned search progress"
```
