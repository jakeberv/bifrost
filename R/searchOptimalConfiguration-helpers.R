.bifrost_search_cli_renderer <- function() {
  rows <- list()
  visible <- 0L
  cursor_hidden <- FALSE
  wrap_disabled <- FALSE
  dynamic <- cli::is_dynamic_tty("stderr")
  spinner <- cli::get_spinner()
  now <- function() proc.time()[["elapsed"]]

  duration <- function(seconds) sprintf("%.1fs", seconds)

  format_row <- function(row) {
    end <- if (is.null(row$finished)) now() else row$finished
    elapsed <- end - row$started
    ratio <- if (row$total > 0L) min(row$current / row$total, 1) else 1
    filled <- round(30 * ratio)
    bar <- paste0(
      strrep(if (cli::is_utf8_output()) "\u2588" else "=", filled),
      strrep(if (cli::is_utf8_output()) " " else "-", 30L - filled)
    )
    icon <- switch(
      row$state,
      active = spinner$frames[[floor(1000 * elapsed / spinner$interval) %%
        length(spinner$frames) + 1L]],
      failed = cli::symbol$cross,
      skipped = cli::symbol$info,
      cli::symbol$tick
    )
    timing <- if (row$state == "active") {
      eta <- if (row$current == 0L) "?" else duration(
        elapsed * (row$total - row$current) / row$current
      )
      paste0("(ETA: ", eta, ")")
    } else paste0("[", duration(elapsed), "]")
    cli::ansi_strtrim(sprintf(
      "%s %s %3d%% [%d/%d] %s %s",
      icon, cli::col_green(bar), round(100 * ratio), row$current,
      row$total, timing, row$status
    ), max(cli::console_width() - 1L, 1L))
  }

  clear <- function() {
    if (!dynamic || visible == 0L) return(invisible(NULL))
    cat("\r\033[2K", file = stderr())
    for (i in seq_len(visible - 1L)) cat("\033[1A\r\033[2K", file = stderr())
    visible <<- 0L
  }

  set_wrap <- function(enabled) {
    cat(if (enabled) "\033[?7h" else "\033[?7l", file = stderr())
    wrap_disabled <<- !enabled
  }

  draw <- function(final = FALSE) {
    if (!final && !cursor_hidden && isTRUE(getOption("cli.hide_cursor", TRUE))) {
      cli::ansi_hide_cursor("stderr")
      cursor_hidden <<- TRUE
    }
    lines <- vapply(rows, format_row, character(1))
    clear()
    if (final && wrap_disabled) set_wrap(TRUE)
    if (!final && !wrap_disabled) set_wrap(FALSE)
    cat(paste(lines, collapse = "\n"), if (final) "\n" else "\r",
        sep = "", file = stderr())
    if (!final) visible <<- length(lines)
    if (final && cursor_hidden) {
      cli::ansi_show_cursor("stderr")
      cursor_hidden <<- FALSE
    }
  }

  list(
    create = function(label, total) {
      id <- paste0("stage-", length(rows) + 1L)
      row <- list2env(list(
        state = "active", current = 0L, total = total, status = label,
        started = now(), finished = NULL
      ), parent = baseenv())
      rows[[id]] <<- row
      row
    },
    update = function(row, state, current, total, status, force = TRUE) {
      row$state <- state
      row$current <- current
      row$total <- total
      row$status <- status
      if (state != "active" && is.null(row$finished)) row$finished <- now()
      if (dynamic) draw() else if (force) cat(
        format_row(row), "\n", sep = "", file = stderr()
      )
    },
    output = function(text) {
      clear()
      if (wrap_disabled) set_wrap(TRUE)
      cat(text, "\n", sep = "", file = stderr())
      if (dynamic) draw()
    },
    done = function() if (dynamic) draw(TRUE)
  )
}

.bifrost_search_progress_session <- function(
    enabled,
    renderer = .bifrost_search_cli_renderer()) {
  stage_rows <- list()
  finalized <- FALSE

  inactive <- function() !isTRUE(enabled) || finalized

  add_row <- function(label,
                      total,
                      current = 0L,
                      state = "active",
                      status = label) {
    row <- list2env(list(
      renderer_row = renderer$create(label, total),
      current = current,
      total = total,
      state = state,
      status = status
    ), parent = emptyenv())
    stage_rows[[label]] <<- row
    row
  }

  render <- function(row, force = TRUE) {
    renderer$update(
      row$renderer_row,
      state = row$state,
      current = row$current,
      total = row$total,
      status = row$status,
      force = force
    )
  }

  update_stage <- function(label,
                           config,
                           state,
                           progression,
                           row_state,
                           status = NULL) {
    if (inactive()) {
      return(invisible(NULL))
    }

    display_total <- config$max_steps - 1L
    display_step <- min(state$step, display_total)
    row <- stage_rows[[label]]
    if (is.null(row)) row <- add_row(label, display_total)
    heartbeat <- isTRUE(progression$amount == 0)
    if (!heartbeat) row$current <- display_step
    row$total <- display_total
    row$state <- row_state
    row$status <- if (is.null(status)) paste(state$message, collapse = "") else status
    render(
      row,
      force = !heartbeat || cli::is_dynamic_tty("stderr")
    )
    invisible(NULL)
  }

  handler <- function(label) {
    report <- function(config,
                       state,
                       progression,
                       row_state = "active",
                       status = NULL,
                       ...) {
      update_stage(label, config, state, progression, row_state, status)
    }
    reporter <- list(
      reset = function(...) invisible(NULL),
      hide = function(...) invisible(NULL),
      unhide = function(...) invisible(NULL),
      initiate = report,
      update = report,
      finish = function(config, state, progression, ...) {
        report(config, state, progression, "complete")
      },
      interrupt = function(config, state, progression, ...) {
        report(
          config,
          state,
          progression,
          "failed",
          conditionMessage(progression)
        )
      }
    )
    progressr::make_progression_handler(
      "cli",
      reporter = reporter,
      enable = TRUE,
      enable_after = 0,
      interval = 0,
      clear = FALSE,
      target = "terminal"
    )
  }

  list(
    handler = handler,
    skip = function(label, reason) {
      if (inactive()) return(invisible(NULL))
      row <- stage_rows[[label]]
      if (is.null(row)) {
        row <- add_row(label, 1L, current = 1L, state = "skipped")
      }
      row$state <- "skipped"
      row$status <- paste(sub(" .*", "", label), reason)
      render(row)
      invisible(NULL)
    },
    output = function(text) {
      if (inactive() || length(stage_rows) == 0L) {
        return(invisible(NULL))
      }
      renderer$output(text)
      if (!cli::is_dynamic_tty("stderr")) {
        for (row in stage_rows) render(row)
      }
      invisible(NULL)
    },
    has_rows = function() !inactive() && length(stage_rows) > 0L,
    finalize = function() {
      if (finalized) return(invisible(NULL))
      finalized <<- TRUE
      if (length(stage_rows)) renderer$done()
      invisible(NULL)
    }
  )
}

.bifrost_search_run_stage <- function(enabled,
                                      steps,
                                      initial_message,
                                      work,
                                      done_message = initial_message,
                                      handler = NULL,
                                      session = NULL,
                                      skipped = NULL) {
  if (!is.null(skipped)) {
    .bifrost_search_report_skipped_stage(
      enabled, initial_message, skipped, session
    )
    enabled <- FALSE
  }
  if (!isTRUE(enabled)) {
    no_op_tick <- function(...) invisible(NULL)
    return(work(no_op_tick))
  }

  owns_session <- is.null(session) && is.null(handler)
  if (owns_session) {
    session <- .bifrost_search_progress_session(TRUE)
    on.exit(session$finalize(), add = TRUE)
  }
  if (is.null(handler)) {
    handler <- session$handler(initial_message)
  }

  progressr::with_progress(
    {
      tick <- progressr::progressor(
        # Reserve one internal milestone so a final heartbeat can arrive after
        # the last worker completion but before the Future value is collected.
        # The CLI reporter subtracts this sentinel from its displayed total.
        steps = steps + 1L,
        message = initial_message,
        auto_finish = FALSE,
        on_exit = FALSE
      )
      value <- work(tick)
      tick(type = "finish", message = done_message)
      value
    },
    handlers = handler,
    cleanup = FALSE,
    delay_conditions = c("message", "warning"),
    enable = TRUE
  )
}

.bifrost_search_report_skipped_stage <- function(enabled,
                                                 label,
                                                 reason,
                                                 session = NULL) {
  if (!isTRUE(enabled)) return(invisible(NULL))
  if (is.null(session)) cli::cli_alert_info("{sub(' .*', '', label)} {reason}")
  else session$skip(label, reason)

  invisible(NULL)
}

.bifrost_search_heartbeat <- function(enabled, tick, message) {
  if (!isTRUE(enabled)) return(NULL)
  force(tick)
  force(message)
  function() tick(amount = 0, message = message)
}

.bifrost_search_with_future_plan <- function(workers,
                                            is_rstudio_flag,
                                            work,
                                            ensure_async = FALSE) {
  thread_vars <- c(
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS"
  )
  old_threads <- Sys.getenv(thread_vars, unset = NA_character_)
  old_plan <- future::plan("list")

  backend_workers <- if (isTRUE(ensure_async) && workers <= 1L) I(1L) else workers

  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )

  restore_threads <- function() {
    for (nm in thread_vars) {
      val <- old_threads[[nm]]
      if (is.na(val)) {
        Sys.unsetenv(nm)
      } else {
        do.call(Sys.setenv, setNames(list(val), nm))
      }
    }
  }

  tryCatch(
    {
      if (.Platform$OS.type == "unix" &&
          !identical(Sys.info()[["sysname"]], "SunOS") &&
          !is_rstudio_flag) {
        plan(multicore, workers = backend_workers)
      } else {
        plan(multisession, workers = backend_workers)
      }

      work()
    },
    finally = {
      plan(old_plan)
      restore_threads()
    }
  )
}

.bifrost_search_await_futures <- function(futures,
                                          heartbeat = function() invisible(NULL),
                                          interval = 0.1) {
  if (length(futures) == 0L) {
    return(list())
  }

  resolved <- rep(FALSE, length(futures))
  values <- vector("list", length(futures))
  finished <- FALSE

  on.exit({
    if (!finished) {
      for (i in which(!resolved)) {
        try(future::cancel(futures[[i]]), silent = TRUE)
      }
    }
  }, add = TRUE)

  while (!all(resolved)) {
    for (i in which(!resolved)) {
      if (future::resolved(futures[[i]], timeout = 0)) {
        values[i] <- list(future::value(futures[[i]]))
        resolved[i] <- TRUE
      }
    }

    if (!all(resolved)) {
      heartbeat()
      Sys.sleep(interval)
    }
  }

  finished <- TRUE
  values
}

.bifrost_search_rng_seeds <- function(count) {
  if (count == 0L) {
    return(list())
  }

  random_env <- globalenv()
  had_seed <- exists(".Random.seed", envir = random_env, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = random_env) else NULL
  old_kind <- RNGkind()

  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = random_env)
    } else if (exists(".Random.seed", envir = random_env, inherits = FALSE)) {
      rm(".Random.seed", envir = random_env)
    }
  }, add = TRUE)

  if (!had_seed) {
    stats::runif(1L)
  }

  seed <- get(".Random.seed", envir = random_env)
  is_lecuyer_seed <- length(seed) == 7L && seed[[1L]] %% 10000L == 407L
  if (!is_lecuyer_seed) {
    RNGkind("L'Ecuyer-CMRG")
    seed <- get(".Random.seed", envir = random_env)
  }

  seeds <- vector("list", count)
  for (i in seq_len(count)) {
    seeds[[i]] <- parallel::nextRNGSubStream(seed)
    seed <- parallel::nextRNGStream(seed)
  }
  seeds
}

.bifrost_search_with_rng_seed <- function(seed, work) {
  random_env <- globalenv()
  had_seed <- exists(".Random.seed", envir = random_env, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = random_env) else NULL

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = random_env)
    } else if (exists(".Random.seed", envir = random_env, inherits = FALSE)) {
      rm(".Random.seed", envir = random_env)
    }
  }, add = TRUE)

  assign(".Random.seed", seed, envir = random_env)
  work()
}

.bifrost_search_await_work <- function(work,
                                       heartbeat = function() invisible(NULL),
                                       seed = NULL,
                                       interval = 0.1) {
  if (is.null(seed)) {
    seed <- .bifrost_search_rng_seeds(1L)[[1L]]
  }
  fit_future <- future::future(work(), seed = seed)
  .bifrost_search_await_futures(
    list(fit_future),
    heartbeat = heartbeat,
    interval = interval
  )[[1L]]
}

.bifrost_search_future_lapply <- function(X,
                                          FUN,
                                          workers,
                                          is_rstudio_flag,
                                          heartbeat,
                                          interval = 0.1) {
  if (length(X) == 0L) {
    return(list())
  }
  if (length(workers) != 1L ||
      !is.numeric(workers) ||
      !isTRUE(is.finite(workers))) {
    stop("`workers` must be a single finite number.", call. = FALSE)
  }
  workers <- max(1L, as.integer(workers))

  chunk_count <- min(workers, length(X))
  chunks <- split(
    seq_along(X),
    rep(seq_len(chunk_count), length.out = length(X))
  )
  item_seeds <- .bifrost_search_rng_seeds(length(X))

  .bifrost_search_with_future_plan(
    workers = workers,
    is_rstudio_flag = is_rstudio_flag,
    ensure_async = TRUE,
    work = function() {
      seeded_eval <- .bifrost_search_with_rng_seed
      chunk_futures <- lapply(chunks, function(indices) {
        future::future({
          lapply(indices, function(i) {
            seeded_eval(
              item_seeds[[i]],
              function() FUN(X[[i]])
            )
          })
        }, seed = item_seeds[[indices[[1L]]]])
      })

      chunk_values <- .bifrost_search_await_futures(
        chunk_futures,
        heartbeat = heartbeat,
        interval = interval
      )

      values <- vector("list", length(X))
      for (i in seq_along(chunks)) {
        values[chunks[[i]]] <- chunk_values[[i]]
      }
      values
    }
  )
}

.bifrost_run_future_lapply_safe <- function(X, FUN, workers, is_rstudio_flag) {
  .bifrost_search_with_future_plan(
    workers = workers,
    is_rstudio_flag = is_rstudio_flag,
    work = function() {
      future.apply::future_lapply(
        X,
        FUN,
        future.seed = TRUE,
        future.scheduling = TRUE
      )
    }
  )
}

.bifrost_search_lapply <- function(X,
                                  FUN,
                                  num_cores,
                                  is_rstudio,
                                  heartbeat = NULL) {
  if (!is.null(heartbeat)) {
    return(.bifrost_search_future_lapply(
      X,
      FUN,
      workers = num_cores,
      is_rstudio_flag = is_rstudio,
      heartbeat = heartbeat
    ))
  }

  if (num_cores <= 1L) {
    return(lapply(X, FUN))
  }

  .bifrost_run_future_lapply_safe(
    X,
    FUN,
    workers = num_cores,
    is_rstudio_flag = is_rstudio
  )
}

.bifrost_search_validate_ic <- function(IC) {
  if (IC != "GIC" && IC != "BIC") {
    stop("IC must be GIC or BIC")
  }

  invisible(IC)
}

.bifrost_search_model_fun <- function(IC) {
  .bifrost_search_validate_ic(IC)

  if (IC == "GIC") {
    fitMvglsAndExtractGIC.formula
  } else {
    fitMvglsAndExtractBIC.formula
  }
}

.bifrost_search_fit_ic <- function(IC, formula, tree, trait_data, ...) {
  model_fun <- .bifrost_search_model_fun(IC)
  model_fun(formula, tree, trait_data, ...)
}

.bifrost_search_ic_value <- function(model_result, IC) {
  .bifrost_search_validate_ic(IC)

  if (IC == "GIC") {
    model_result$GIC$GIC
  } else {
    model_result$BIC$BIC
  }
}

.bifrost_search_score_candidates <- function(candidate_trees_shifts,
                                           baseline_ic,
                                           IC,
                                           formula,
                                           trait_data,
                                           args_list,
                                           num_cores,
                                           is_rstudio,
                                           tick = function(...) invisible(NULL),
                                           heartbeat = NULL,
                                           fit = .bifrost_search_fit_ic) {
  candidate_results <- .bifrost_search_lapply(
    seq_along(candidate_trees_shifts),
    function(i) {
      tree <- candidate_trees_shifts[[i]]
      result <- do.call(fit, c(list(IC, formula, tree, trait_data), args_list))
      tick(message = paste(
        "[1/3]", names(candidate_trees_shifts)[i], "scored"
      ))
      result
    },
    num_cores = num_cores,
    is_rstudio = is_rstudio,
    heartbeat = heartbeat
  )

  delta_ic_list <- sapply(candidate_results, function(res) {
    baseline_ic - .bifrost_search_ic_value(res, IC)
  })

  list(
    candidate_results = candidate_results,
    delta_ic_list = delta_ic_list,
    sorted_candidates = candidate_trees_shifts[order(delta_ic_list, decreasing = TRUE)]
  )
}

.bifrost_search_history_dir <- function(store_model_fit_history) {
  if (!isTRUE(store_model_fit_history)) {
    return(NULL)
  }

  base_dir <- file.path(tempdir(), "bifrost_fit_history")
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

  date_str <- format(Sys.Date(), "%Y-%m-%d")
  sub_dir <- file.path(base_dir, date_str)
  counter <- 1L
  while (dir.exists(sub_dir)) {
    counter <- counter + 1L
    sub_dir <- file.path(base_dir, paste0(date_str, "_", counter))
  }
  dir.create(sub_dir, recursive = TRUE)

  sub_dir
}

.bifrost_search_save_history <- function(model_fit_history, sub_dir, iteration_num) {
  saveRDS(
    model_fit_history,
    file = file.path(sub_dir, paste0("iteration_", iteration_num, ".rds"))
  )

  invisible(NULL)
}

.bifrost_search_load_history <- function(sub_dir, IC) {
  rds_files <- list.files(
    sub_dir,
    pattern = "^iteration_\\d+\\.rds$",
    full.names = TRUE
  )
  rds_files <- rds_files[order(as.numeric(gsub("\\D", "", basename(rds_files))))]

  model_fit_history <- lapply(rds_files, readRDS)

  ic_acceptance_matrix <- do.call(rbind, lapply(model_fit_history, function(x) {
    stored_ic <- if (is.null(x$ic)) {
      NA_real_
    } else {
      suppressWarnings(as.numeric(x$ic[1L]))
    }
    if (is.finite(stored_ic)) {
      c(stored_ic, x$accepted)
    } else if (is.null(x$model)) {
      c(NA_real_, x$accepted)
    } else {
      c(.bifrost_search_ic_value(x$model, IC), x$accepted)
    }
  }))

  list(
    fits = model_fit_history,
    ic_acceptance_matrix = ic_acceptance_matrix
  )
}

.bifrost_search_forward <- function(sorted_candidates,
                                   current_best_tree,
                                   current_best_ic,
                                   shift_id,
                                   IC,
                                   formula,
                                   trait_data,
                                   shift_acceptance_threshold,
                                   store_model_fit_history,
                                   sub_dir,
                                   plot,
                                   verbose_log,
                                   tick = function(...) invisible(NULL),
                                   heartbeat = NULL,
                                   is_rstudio = FALSE,
                                   fit = .bifrost_search_fit_ic,
                                   .future_plan_active = FALSE,
                                   ...) {
  if (length(sorted_candidates) > 0L &&
      !is.null(heartbeat) &&
      !isTRUE(.future_plan_active)) {
    dots <- list(...)
    forward_args <- c(list(
      sorted_candidates = sorted_candidates,
      current_best_tree = current_best_tree,
      current_best_ic = current_best_ic,
      shift_id = shift_id,
      IC = IC,
      formula = formula,
      trait_data = trait_data,
      shift_acceptance_threshold = shift_acceptance_threshold,
      store_model_fit_history = store_model_fit_history,
      sub_dir = sub_dir,
      plot = plot,
      verbose_log = verbose_log,
      tick = tick,
      heartbeat = heartbeat,
      is_rstudio = is_rstudio,
      fit = fit,
      .future_plan_active = TRUE
    ), dots)

    return(.bifrost_search_with_future_plan(
      workers = 1L,
      is_rstudio_flag = is_rstudio,
      ensure_async = TRUE,
      work = function() do.call(.bifrost_search_forward, forward_args)
    ))
  }

  shift_vec <- list() #initialize shift_vec
  model_with_shift_no_uncertainty <- NULL #initialize output
  best_tree_no_uncertainty <- NULL #initialize output
  # Initialize the list to collect warning messages
  warnings_list <- list()
  fit_seeds <- if (is.null(heartbeat)) {
    NULL
  } else {
    .bifrost_search_rng_seeds(length(sorted_candidates))
  }

  # In-memory accumulator (lightweight; actual fits stored on disk)
  model_fit_history <- list()

  #Run the primary shift configuration search

  for (i in seq_along(sorted_candidates)) {
    shift_node_name <- names(sorted_candidates)[i]
    shift_node_number <- as.integer(sub("Node ", "", shift_node_name))
    percent_complete <- round((i / length(sorted_candidates)) * 100, 2)
    verbose_log("Evaluating shift at node %d (%.2f%% complete)", shift_node_number, percent_complete)
    outcome <- "error"

    add_shift_result <- addShiftToModel(current_best_tree, shift_node_number, shift_id)
    shifted_tree <- add_shift_result$tree
    shift_id <- add_shift_result$shift_id

    if(plot == TRUE){
      nodelabels(text = shift_id, node = shift_node_number)
    }

    tryCatch({
      fit_args <- c(list(IC, formula, shifted_tree, trait_data), list(...))
      fit_work <- function() do.call(fit, fit_args)
      model_with_shift <- withCallingHandlers(
        if (is.null(heartbeat)) {
          fit_work()
        } else {
          .bifrost_search_await_work(
            fit_work,
            heartbeat = heartbeat,
            seed = fit_seeds[[i]]
          )
        },
        warning = function(w) {
          warning_message <- paste("Warning in evaluating shift at node", shift_node_number, ":", w$message)
          warning(warning_message)
          warnings_list[[length(warnings_list) + 1]] <<- warning_message
          invokeRestart("muffleWarning")
        }
      )
      new_ic <- .bifrost_search_ic_value(model_with_shift, IC)

      # Calculate delta IC
      delta_ic <- current_best_ic - new_ic

      # Store model fit and acceptance status (including delta_ic)
      if (store_model_fit_history) {
        model_fit_history <- list(
          step = i,
          candidate_node = shift_node_number,
          regime_id = as.character(shift_id),
          model = model_with_shift,
          ic = new_ic,
          accepted = delta_ic >= shift_acceptance_threshold,
          delta_ic = delta_ic,
          status = if (delta_ic >= shift_acceptance_threshold) "accepted" else "rejected"
        )
      }

      # Decision logic (unchanged)
      if (delta_ic >= shift_acceptance_threshold) {
        current_best_tree <- shifted_tree
        current_best_ic <- new_ic
        outcome <- "accepted"
        verbose_log(
          "Shift at node %d accepted. Updated %s: %.2f; Delta %s: %.2f",
          shift_node_number, IC, current_best_ic, IC, delta_ic
        )

        shift_vec[[length(shift_vec) + 1]] <- shift_node_number

        best_tree_no_uncertainty <- current_best_tree
        model_with_shift_no_uncertainty <- model_with_shift
      } else {
        outcome <- "rejected"
        verbose_log(
          "Shift at node %d rejected. Delta %s: %.2f < threshold: %.2f",
          shift_node_number, IC, delta_ic, shift_acceptance_threshold
        )
      }
    }, error = function(e) {
      # Handle errors (unchanged)
      warning_message <- paste("Error in evaluating shift at node", shift_node_number, ":", e$message)
      warning(warning_message)
      warnings_list[[length(warnings_list) + 1]] <<- warning_message

      # Also store the error in the model fit history
      if (store_model_fit_history) {
        model_fit_history <<- list(
          step = i,
          candidate_node = shift_node_number,
          regime_id = as.character(shift_id),
          model = NULL,
          ic = NA_real_,
          accepted = FALSE,
          delta_ic = NA,
          status = "error",
          error = e$message
        )
      }

    })
    tick(message = sprintf(
      "[2/3] Node %d %s",
      shift_node_number,
      if (outcome == "error") "failed" else outcome
    ))
    if (isTRUE(store_model_fit_history) && !is.null(sub_dir)) {
      iteration_num <- i + 1L
      .bifrost_search_save_history(model_fit_history, sub_dir, iteration_num)
    }
    if(plot == TRUE){
      colorvec <- setNames(object = c("black", rainbow(length(unique(getStates(shifted_tree, type = "both"))) - 1)),
                           nm = sort(as.numeric(unique(getStates(shifted_tree, type = "both")))))
      plotSimmap(current_best_tree, colors = colorvec, ftype = "off")
    }
  }

  list(
    current_best_tree = current_best_tree,
    current_best_ic = current_best_ic,
    shift_id = shift_id,
    shift_vec = shift_vec,
    shifts_no_uncertainty = unlist(shift_vec),
    model_with_shift_no_uncertainty = model_with_shift_no_uncertainty,
    best_tree_no_uncertainty = best_tree_no_uncertainty,
    warnings_list = warnings_list,
    model_fit_history = model_fit_history
  )
}

.bifrost_search_empty_ic_weights_df <- function() {
  data.frame(
    node = integer(),
    ic_with_shift = numeric(),
    ic_without_shift = numeric(),
    delta_ic = numeric(),
    ic_weight_withshift = numeric(),
    ic_weight_withoutshift = numeric(),
    evidence_ratio = numeric()
  )
}

.bifrost_search_ic_weights_row <- function(shift_node_number, original_ic, ic_without_shift) {
  delta_ic <- original_ic - ic_without_shift

  icw <- aicw(c(original_ic, ic_without_shift))$aicweights
  w_with <- icw[1]
  w_without <- icw[2]
  er <- w_with / w_without

  data.frame(
    node = shift_node_number,
    ic_with_shift = original_ic,
    ic_without_shift = ic_without_shift,
    delta_ic = delta_ic,
    ic_weight_withshift = w_with,
    ic_weight_withoutshift = w_without,
    evidence_ratio = er
  )
}

.bifrost_search_calculate_ic_weights <- function(uncertaintyweights,
                                                uncertaintyweights_par,
                                                shift_vec,
                                                best_tree_no_uncertainty,
                                                model_with_shift_no_uncertainty,
                                                IC,
                                                formula,
                                                trait_data,
                                                args_list,
                                                num_cores,
                                                is_rstudio,
                                                verbose_log,
                                                tick = function(...) invisible(NULL),
                                                heartbeat = NULL,
                                                fit = .bifrost_search_fit_ic) {
  ic_weights_df <- .bifrost_search_empty_ic_weights_df()

  if (xor(uncertaintyweights, uncertaintyweights_par)) {

    # If no shifts, return empty df (consistent in both modes)
    if (length(unlist(shift_vec)) == 0) {
      verbose_log("%s", "No shifts were detected in the initial search; skipping IC weights calculation.")
      ic_weights_df <- .bifrost_search_empty_ic_weights_df()

    } else {

      # Retrieve the IC of the optimized model before uncertainty analysis
      original_ic <- .bifrost_search_ic_value(model_with_shift_no_uncertainty, IC)

      verbose_log("Considering %d shifts in the candidate set", length(shift_vec))
      verbose_log(
        "There are %d shifts in the mapped tree",
        length(unique(getStates(best_tree_no_uncertainty, type = "both"))) - 1
      )

      if (uncertaintyweights) {
        verbose_log("%s", "Calculating IC weights for initially identified shifts...")

        calculate_serial_weights <- function() {
          serial_weights <- .bifrost_search_empty_ic_weights_df()
          shift_nodes <- unlist(shift_vec)
          fit_seeds <- if (is.null(heartbeat)) {
            NULL
          } else {
            .bifrost_search_rng_seeds(length(shift_nodes))
          }

          for (i in seq_along(shift_nodes)) {
            shift_node_number <- shift_nodes[[i]]
            verbose_log("Re-estimating model without shift at node %d", shift_node_number)

            tree_without_current_shift <- removeShiftFromTree(
              best_tree_no_uncertainty,
              shift_node_number
            )
            fit_args <- c(
              list(IC, formula, tree_without_current_shift, trait_data),
              args_list
            )
            fit_work <- function() do.call(fit, fit_args)
            model_without_current_shift <- if (is.null(heartbeat)) {
              fit_work()
            } else {
              .bifrost_search_await_work(
                fit_work,
                heartbeat = heartbeat,
                seed = fit_seeds[[i]]
              )
            }
            ic_without_current_shift <- .bifrost_search_ic_value(
              model_without_current_shift,
              IC
            )

            ic_weight_row <- .bifrost_search_ic_weights_row(
              shift_node_number,
              original_ic,
              ic_without_current_shift
            )

            verbose_log("IC weight for the shift is %.2f", ic_weight_row$ic_weight_withshift)

            serial_weights <- rbind(
              serial_weights,
              ic_weight_row
            )
            tick(message = sprintf(
              "[3/3] Weighted node %d", shift_node_number
            ))
          }

          serial_weights
        }

        if (is.null(heartbeat)) {
          ic_weights_df <- calculate_serial_weights()
        } else {
          ic_weights_df <- .bifrost_search_with_future_plan(
            workers = 1L,
            is_rstudio_flag = is_rstudio,
            ensure_async = TRUE,
            work = calculate_serial_weights
          )
        }
      }

      if (uncertaintyweights_par) {
        verbose_log("%s", "Calculating IC weights for initially identified shifts...")

        ic_weights_df <- .bifrost_search_empty_ic_weights_df()

        shift_removed_trees <- lapply(unlist(shift_vec), function(shift_node_number) {
          removeShiftFromTree(best_tree_no_uncertainty, shift_node_number)
        })

        ic_results <- .bifrost_search_lapply(
          seq_along(shift_removed_trees),
          function(i) {
            tree <- shift_removed_trees[[i]]
            model_without_shift <- do.call(
              fit,
              c(list(IC, formula, tree, trait_data), args_list)
            )
            ic_without_shift <- .bifrost_search_ic_value(model_without_shift, IC)
            delta_ic <- original_ic - ic_without_shift

            icw <- aicw(c(original_ic, ic_without_shift))$aicweights

            tick(message = sprintf(
              "[3/3] Weighted node %d", unlist(shift_vec)[i]
            ))

            c(
              ic_without_shift = ic_without_shift,
              delta_ic = delta_ic,
              ic_weight_withshift = icw[1],
              ic_weight_withoutshift = icw[2]
            )
          },
          num_cores = num_cores,
          is_rstudio = is_rstudio,
          heartbeat = heartbeat
        )

        for (i in seq_along(shift_removed_trees)) {
          shift_node_number <- unlist(shift_vec)[i]
          ic_res <- ic_results[[i]]

          # scalar extraction (avoids named-vector quirks)
          ic_without <- as.numeric(ic_res[["ic_without_shift"]])

          ic_weights_df <- rbind(
            ic_weights_df,
            .bifrost_search_ic_weights_row(shift_node_number, original_ic, ic_without)
          )
        }
      }
    }

  } else {
    if (isTRUE(uncertaintyweights) && isTRUE(uncertaintyweights_par)) {
      stop("Exactly one of uncertaintyweights or uncertaintyweights_par must be TRUE.")
    }
    # If both are FALSE, do nothing (IC weights not requested).
  }

  ic_weights_df
}

.bifrost_search_no_uncertainty_components <- function(shifts_no_uncertainty,
                                                     model_with_shift_no_uncertainty,
                                                     best_tree_no_uncertainty,
                                                     baseline_model,
                                                     baseline_candidate_tree) {
  if (length(shifts_no_uncertainty) > 0L) {
    # use the accepted-shift model/tree
    transformed <- model_with_shift_no_uncertainty$model$corrSt$phy
    untransformed <- model_with_shift_no_uncertainty$model$corrSt$phy
    untransformed$edge.length <- best_tree_no_uncertainty$edge.length
    model <- model_with_shift_no_uncertainty$model
  } else {
    # fallback to baseline model/tree when no shifts were accepted
    transformed <- baseline_model$model$corrSt$phy
    untransformed <- baseline_model$model$corrSt$phy
    untransformed$edge.length <- baseline_candidate_tree$edge.length
    model <- baseline_model$model
  }

  list(
    tree_transformed = transformed,
    tree_untransformed = untransformed,
    model = model
  )
}
