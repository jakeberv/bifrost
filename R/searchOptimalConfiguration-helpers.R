.bifrost_search_progress_handler <- function() {
  progressr::handler_cli(
    enable = TRUE,
    show_after = 0,
    clear = FALSE,
    format = paste0(
      "{cli::pb_spin} {cli::pb_bar} {cli::pb_percent} ",
      "[{cli::pb_current}/{cli::pb_total}] ",
      "(ETA: {cli::pb_eta}) {cli::pb_status}"
    ),
    format_done = paste0(
      "{cli::symbol$tick} {cli::pb_bar} {cli::pb_percent} ",
      "[{cli::pb_current}/{cli::pb_total}] ",
      "[{cli::pb_elapsed}] {cli::pb_status}"
    ),
    format_failed = paste0(
      "{cli::symbol$cross} {cli::pb_bar} {cli::pb_percent} ",
      "[{cli::pb_current}/{cli::pb_total}] ",
      "[{cli::pb_elapsed}] {cli::pb_status}"
    )
  )
}

.bifrost_search_run_stage <- function(enabled,
                                      steps,
                                      initial_message,
                                      work,
                                      handler = .bifrost_search_progress_handler()) {
  if (!isTRUE(enabled)) {
    no_op_tick <- function(...) invisible(NULL)
    stage_result <- work(no_op_tick)
    return(stage_result$value)
  }

  progressr::with_progress(
    {
      tick <- progressr::progressor(
        steps = steps,
        message = initial_message,
        auto_finish = FALSE,
        on_exit = FALSE
      )
      stage_result <- work(tick)
      tick(type = "finish", message = stage_result$done)
      stage_result$value
    },
    handlers = handler,
    cleanup = FALSE,
    delay_conditions = c("message", "warning"),
    enable = TRUE
  )
}

.bifrost_search_report_skipped_stage <- function(enabled, label, reason) {
  if (isTRUE(enabled)) {
    cli::cli_alert_info("{label} - skipped: {reason}")
  }

  invisible(NULL)
}

.bifrost_run_future_lapply_safe <- function(X, FUN, workers, is_rstudio_flag) {
  thread_vars <- c(
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS"
  )
  old_threads <- Sys.getenv(thread_vars, unset = NA_character_)

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
        plan(multicore, workers = workers)
      } else {
        plan(multisession, workers = workers)
      }

      future.apply::future_lapply(
        X,
        FUN,
        future.seed = TRUE,
        future.scheduling = TRUE
      )
    },
    finally = {
      plan(sequential)
      restore_threads()
    }
  )
}

.bifrost_search_lapply <- function(X, FUN, num_cores, is_rstudio) {
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
                                           fit = .bifrost_search_fit_ic) {
  candidate_results <- .bifrost_search_lapply(
    seq_along(candidate_trees_shifts),
    function(i) {
      tree <- candidate_trees_shifts[[i]]
      result <- do.call(fit, c(list(IC, formula, tree, trait_data), args_list))
      tick(message = paste(
        "[1/3] Candidate scoring - completed",
        names(candidate_trees_shifts)[i]
      ))
      result
    },
    num_cores = num_cores,
    is_rstudio = is_rstudio
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
                                   fit = .bifrost_search_fit_ic,
                                   ...) {
  shift_vec <- list() #initialize shift_vec
  model_with_shift_no_uncertainty <- NULL #initialize output
  best_tree_no_uncertainty <- NULL #initialize output
  # Initialize the list to collect warning messages
  warnings_list <- list()
  outcome_counts <- c(accepted = 0L, rejected = 0L, error = 0L)

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
      model_with_shift <- withCallingHandlers(
        fit(IC, formula, shifted_tree, trait_data, ...),
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
    outcome_counts[[outcome]] <- outcome_counts[[outcome]] + 1L
    tick(message = sprintf(
      "[2/3] Greedy shift search - node %d %s",
      shift_node_number,
      outcome
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
    outcome_counts = outcome_counts,
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

        ic_weights_df <- .bifrost_search_empty_ic_weights_df()

        for (shift_node_number in unlist(shift_vec)) {
          verbose_log("Re-estimating model without shift at node %d", shift_node_number)

          tree_without_current_shift <- removeShiftFromTree(best_tree_no_uncertainty, shift_node_number)
          model_without_current_shift <- do.call(
            fit,
            c(list(IC, formula, tree_without_current_shift, trait_data), args_list)
          )
          ic_without_current_shift <- .bifrost_search_ic_value(model_without_current_shift, IC)

          ic_weight_row <- .bifrost_search_ic_weights_row(
            shift_node_number,
            original_ic,
            ic_without_current_shift
          )

          verbose_log("IC weight for the shift is %.2f", ic_weight_row$ic_weight_withshift)

          ic_weights_df <- rbind(
            ic_weights_df,
            ic_weight_row
          )
          tick(message = sprintf(
            "[3/3] IC-weight re-estimation - completed node %d",
            shift_node_number
          ))
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
              "[3/3] IC-weight re-estimation - completed node %d",
              unlist(shift_vec)[i]
            ))

            c(
              ic_without_shift = ic_without_shift,
              delta_ic = delta_ic,
              ic_weight_withshift = icw[1],
              ic_weight_withoutshift = icw[2]
            )
          },
          num_cores = num_cores,
          is_rstudio = is_rstudio
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
