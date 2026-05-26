bifrost_run_future_lapply_safe <- function(X, FUN, workers, is_rstudio_flag) {
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

bifrost_search_validate_ic <- function(IC) {
  if (IC != "GIC" && IC != "BIC") {
    stop("IC must be GIC or BIC")
  }

  invisible(IC)
}

bifrost_search_model_fun <- function(IC) {
  bifrost_search_validate_ic(IC)

  if (IC == "GIC") {
    fitMvglsAndExtractGIC.formula
  } else {
    fitMvglsAndExtractBIC.formula
  }
}

bifrost_search_fit_ic <- function(IC, formula, tree, trait_data, ...) {
  model_fun <- bifrost_search_model_fun(IC)
  model_fun(formula, tree, trait_data, ...)
}

bifrost_search_ic_value <- function(model_result, IC) {
  bifrost_search_validate_ic(IC)

  if (IC == "GIC") {
    model_result$GIC$GIC
  } else {
    model_result$BIC$BIC
  }
}

bifrost_search_score_candidates <- function(candidate_trees_shifts,
                                           baseline_ic,
                                           IC,
                                           formula,
                                           trait_data,
                                           args_list,
                                           num_cores,
                                           is_rstudio) {
  candidate_results <- bifrost_run_future_lapply_safe(
    candidate_trees_shifts,
    function(tree) {
      do.call(bifrost_search_fit_ic, c(list(IC, formula, tree, trait_data), args_list))
    },
    workers = num_cores,
    is_rstudio_flag = is_rstudio
  )

  delta_ic_list <- sapply(candidate_results, function(res) {
    baseline_ic - bifrost_search_ic_value(res, IC)
  })

  list(
    candidate_results = candidate_results,
    delta_ic_list = delta_ic_list,
    sorted_candidates = candidate_trees_shifts[order(delta_ic_list, decreasing = TRUE)]
  )
}

bifrost_search_history_dir <- function(store_model_fit_history) {
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

bifrost_search_save_history <- function(model_fit_history, sub_dir, iteration_num) {
  saveRDS(
    model_fit_history,
    file = file.path(sub_dir, paste0("iteration_", iteration_num, ".rds"))
  )

  invisible(NULL)
}

bifrost_search_load_history <- function(sub_dir, IC) {
  rds_files <- list.files(
    sub_dir,
    pattern = "^iteration_\\d+\\.rds$",
    full.names = TRUE
  )
  rds_files <- rds_files[order(as.numeric(gsub("\\D", "", basename(rds_files))))]

  model_fit_history <- lapply(rds_files, readRDS)

  ic_acceptance_matrix <- do.call(rbind, lapply(model_fit_history, function(x) {
    if (is.null(x$model)) {
      c(NA_real_, x$accepted)
    } else {
      c(bifrost_search_ic_value(x$model, IC), x$accepted)
    }
  }))

  list(
    fits = model_fit_history,
    ic_acceptance_matrix = ic_acceptance_matrix
  )
}

bifrost_search_empty_ic_weights_df <- function() {
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

bifrost_search_ic_weights_row <- function(shift_node_number, original_ic, ic_without_shift) {
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

bifrost_search_calculate_ic_weights <- function(uncertaintyweights,
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
                                                progress) {
  ic_weights_df <- bifrost_search_empty_ic_weights_df()

  if (xor(uncertaintyweights, uncertaintyweights_par)) {

    # If no shifts, return empty df (consistent in both modes)
    if (length(unlist(shift_vec)) == 0) {
      progress("%s", "No shifts were detected in the initial search; skipping IC weights calculation.")
      ic_weights_df <- bifrost_search_empty_ic_weights_df()

    } else {

      # Retrieve the IC of the optimized model before uncertainty analysis
      original_ic <- bifrost_search_ic_value(model_with_shift_no_uncertainty, IC)

      progress("Considering %d shifts in the candidate set", length(shift_vec))
      progress(
        "There are %d shifts in the mapped tree",
        length(unique(getStates(best_tree_no_uncertainty, type = "both"))) - 1
      )

      if (uncertaintyweights) {
        progress("%s", "Calculating IC weights for initially identified shifts...")

        ic_weights_df <- bifrost_search_empty_ic_weights_df()

        for (shift_node_number in unlist(shift_vec)) {
          progress("Re-estimating model without shift at node %d", shift_node_number)

          tree_without_current_shift <- removeShiftFromTree(best_tree_no_uncertainty, shift_node_number)
          model_without_current_shift <- do.call(
            bifrost_search_fit_ic,
            c(list(IC, formula, tree_without_current_shift, trait_data), args_list)
          )
          ic_without_current_shift <- bifrost_search_ic_value(model_without_current_shift, IC)

          ic_weight_row <- bifrost_search_ic_weights_row(
            shift_node_number,
            original_ic,
            ic_without_current_shift
          )

          progress("IC weight for the shift is %.2f", ic_weight_row$ic_weight_withshift)

          ic_weights_df <- rbind(
            ic_weights_df,
            ic_weight_row
          )
        }
      }

      if (uncertaintyweights_par) {
        progress("%s", "Calculating IC weights for initially identified shifts in parallel...")

        ic_weights_df <- bifrost_search_empty_ic_weights_df()

        shift_removed_trees <- lapply(unlist(shift_vec), function(shift_node_number) {
          removeShiftFromTree(best_tree_no_uncertainty, shift_node_number)
        })

        ic_results <- bifrost_run_future_lapply_safe(
          shift_removed_trees,
          function(tree) {
            model_without_shift <- do.call(
              bifrost_search_fit_ic,
              c(list(IC, formula, tree, trait_data), args_list)
            )
            ic_without_shift <- bifrost_search_ic_value(model_without_shift, IC)
            delta_ic <- original_ic - ic_without_shift

            icw <- aicw(c(original_ic, ic_without_shift))$aicweights

            c(
              ic_without_shift = ic_without_shift,
              delta_ic = delta_ic,
              ic_weight_withshift = icw[1],
              ic_weight_withoutshift = icw[2]
            )
          },
          workers = num_cores,
          is_rstudio_flag = is_rstudio
        )

        for (i in seq_along(shift_removed_trees)) {
          shift_node_number <- unlist(shift_vec)[i]
          ic_res <- ic_results[[i]]

          # scalar extraction (avoids named-vector quirks)
          ic_without <- as.numeric(ic_res[["ic_without_shift"]])

          ic_weights_df <- rbind(
            ic_weights_df,
            bifrost_search_ic_weights_row(shift_node_number, original_ic, ic_without)
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

bifrost_search_no_uncertainty_components <- function(shifts_no_uncertainty,
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
