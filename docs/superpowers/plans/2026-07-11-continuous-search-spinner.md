# Continuous Search Spinner Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Continuously animate the search spinner during long fits while preserving truthful completion counts, existing parallelism, and failure cleanup.

**Architecture:** Private Future helpers keep expensive fits off the main R process and poll them every 100 ms. Polls emit zero-increment `progressr` events; fit completion remains the only event that advances a stage. Candidate and parallel weight work is chunked, while greedy work remains one-fit-at-a-time.

**Tech Stack:** R, future, future.apply, progressr, cli, testthat

## Global Constraints

- Do not add or change public arguments.
- Keep `num_cores` as the maximum number of concurrent expensive fits.
- Keep worker processes free of terminal rendering and plotting.
- Keep the `progress = FALSE` serial/parallel behavior unchanged.
- Restore Future plans and BLAS/OpenMP environment variables after success, errors, and interrupts.
- Keep partial failed progress lines persistent.

---

### Task 1: Add asynchronous Future polling primitives

**Files:**
- Modify: `R/searchOptimalConfiguration-helpers.R`
- Test: `tests/testthat/test-search-progress.R`

**Interfaces:**
- Produces: `.bifrost_search_with_future_plan(workers, is_rstudio_flag, work, animate)`
- Produces: `.bifrost_search_await_futures(futures, heartbeat, interval = 0.1)`
- Produces: `.bifrost_search_future_lapply(X, FUN, workers, is_rstudio_flag, heartbeat)`
- Produces: `.bifrost_search_future_value(work, is_rstudio_flag, heartbeat)`

- [ ] **Step 1: Write a failing heartbeat test**

Add a test that starts a 350 ms Future, calls `.bifrost_search_await_futures()`, records heartbeat timestamps, and expects the unchanged value plus at least two heartbeat calls.

```r
beats <- 0L
value <- .bifrost_search_future_value(
  work = function() { Sys.sleep(0.35); 42L },
  is_rstudio_flag = TRUE,
  heartbeat = function() beats <<- beats + 1L
)
expect_identical(value, 42L)
expect_gte(beats, 2L)
```

- [ ] **Step 2: Run the focused test and confirm the missing-helper failure**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-search-progress.R", filter = "heartbeat")'`

Expected: FAIL because `.bifrost_search_future_value()` does not exist.

- [ ] **Step 3: Implement scoped plan setup and polling**

Refactor the current thread-variable setup into `.bifrost_search_with_future_plan()`. When animation needs one logical worker, create a two-slot backend but launch only one compute Future. Implement polling with `future::resolved(..., timeout = 0)`, `Sys.sleep(interval)`, `future::value()`, and cancellation in `on.exit`.

```r
while (!all(done)) {
  done[!done] <- vapply(futures[!done], future::resolved, logical(1), timeout = 0)
  if (!all(done)) {
    heartbeat()
    Sys.sleep(interval)
  }
}
values <- lapply(futures, future::value)
```

- [ ] **Step 4: Run the heartbeat test and existing progress lifecycle tests**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-search-progress.R")'`

Expected: PASS with no warnings.

### Task 2: Animate candidate scoring and greedy fits

**Files:**
- Modify: `R/searchOptimalConfiguration-helpers.R`
- Modify: `R/searchOptimalConfiguration.R`
- Test: `tests/testthat/test-search-progress.R`

**Interfaces:**
- Consumes: Future helpers from Task 1.
- Changes: `.bifrost_search_lapply(..., heartbeat = NULL)` uses the old path when `heartbeat` is `NULL` and the chunked Future path otherwise.
- Changes: `.bifrost_search_forward(..., heartbeat = NULL, is_rstudio = FALSE)` evaluates greedy fits asynchronously only when a heartbeat is supplied.

- [ ] **Step 1: Extend the collector and write failing amount assertions**

Record `progression$amount` in `.collect_progress_handler()`. Add slow candidate and greedy fits and assert at least one zero amount plus exactly one unit amount per completed fit.

```r
completion_amounts <- events$amounts[events$types == "update" & events$amounts > 0]
heartbeat_amounts <- events$amounts[events$types == "update" & events$amounts == 0]
expect_length(completion_amounts, candidate_count)
expect_gt(length(heartbeat_amounts), 0L)
```

- [ ] **Step 2: Run focused tests and confirm that no heartbeat events exist**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-search-progress.R", filter = "candidate|greedy")'`

Expected: FAIL because the current fits block the main process.

- [ ] **Step 3: Wire candidate heartbeat callbacks**

Inside the candidate stage, pass a callback that calls the progressor without advancing it.

```r
heartbeat <- function() {
  tick(amount = 0, message = "[1/3] Candidate scoring - fitting candidates")
}
```

Use chunk Futures for animated candidate scoring and reassemble results in input order.

- [ ] **Step 4: Wire greedy heartbeat callbacks**

Run the greedy loop in the main process, but evaluate each `fit()` via `.bifrost_search_future_value()`. Keep warning handling, search-state mutation, verbose messages, plotting, history writes, and completion ticks in the main process.

- [ ] **Step 5: Run candidate and greedy tests**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-search-progress.R", filter = "candidate|greedy")'`

Expected: PASS, including accepted/rejected/recoverable-error counts.

### Task 3: Animate serial and parallel IC-weight fits

**Files:**
- Modify: `R/searchOptimalConfiguration-helpers.R`
- Modify: `R/searchOptimalConfiguration.R`
- Test: `tests/testthat/test-search-progress.R`

**Interfaces:**
- Consumes: animated `.bifrost_search_lapply()` from Task 2.
- Changes: `.bifrost_search_calculate_ic_weights(..., heartbeat = NULL)` uses one logical worker for serial weights and `num_cores` workers for parallel weights.

- [ ] **Step 1: Write failing serial and parallel heartbeat assertions**

Use slow fake weight fits. Assert zero-increment updates occur in both modes and exactly one positive completion amount occurs per accepted shift.

- [ ] **Step 2: Run the focused weight test and confirm failure**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-search-progress.R", filter = "IC-weight")'`

Expected: FAIL because weight fits do not emit heartbeats.

- [ ] **Step 3: Consolidate weight fitting through the animated map**

Build leave-one-shift-out trees in the main process, evaluate them through `.bifrost_search_lapply()`, calculate rows in the main process, and keep the existing completion message for every shift.

- [ ] **Step 4: Run the focused weight and skipped-stage tests**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-search-progress.R", filter = "IC-weight|skipped")'`

Expected: PASS with unchanged skipped-stage output.

### Task 4: Prove cleanup and document the behavior

**Files:**
- Modify: `tests/testthat/test-search-progress.R`
- Modify: `R/searchOptimalConfiguration.R`
- Modify: `man/searchOptimalConfiguration.Rd`
- Modify: `NEWS.md`
- Modify: `vignettes/quick-start-vignette.Rmd`

**Interfaces:**
- No public API changes.
- Documentation promises continuous animation only when `progress = TRUE`.

- [ ] **Step 1: Write an error cleanup regression test**

Start a slow Future alongside one that throws, assert the original error propagates, and verify `future::plan("list")` plus all five thread environment variables are restored.

- [ ] **Step 2: Run the cleanup test and verify it fails before final cleanup wiring**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-search-progress.R", filter = "restores")'`

Expected: FAIL until outstanding Futures are cancelled and state restoration is complete.

- [ ] **Step 3: Finish cancellation and cleanup paths**

Cancel unresolved Futures from `on.exit`, preserve the original worker error/interrupt, and restore state in `finally`.

- [ ] **Step 4: Update documentation**

Explain that enabled progress uses main-process heartbeat redraws during Future-hosted fits, that completion counts remain fit-based, and that `progress = FALSE` retains the quiet path. Regenerate Rd with roxygen2.

- [ ] **Step 5: Run focused and full verification**

Run:

```sh
Rscript -e 'testthat::test_file("tests/testthat/test-search-progress.R")'
Rscript -e 'testthat::test_file("tests/testthat/test-searchOptimalConfiguration.R")'
Rscript -e 'roxygen2::roxygenise(roclets = c("rd", "namespace"))'
Rscript -e 'devtools::test()'
R CMD build .
R CMD check --no-manual bifrost_0.1.4.tar.gz
git diff --check
```

Expected: zero test failures, successful build, `Status: OK`, and no whitespace errors.
