#' Search for an Optimal Multi-Regime (Shift) Configuration on a Phylogeny
#'
#' @description
#' Greedy, stepwise search for evolutionary regime shifts on a SIMMAP-style phylogeny
#' using multivariate \code{mvgls} fits from \pkg{mvMORPH}. The routine:
#' \enumerate{
#'   \item builds one-shift candidate trees for all internal nodes meeting a tip-size threshold
#'         (via \code{\link{generatePaintedTrees}}),
#'   \item fits each candidate in parallel and ranks them by improvement in the chosen
#'         information criterion (IC; \code{GIC} or \code{BIC}),
#'   \item iteratively adds shifts that pass a user-defined acceptance threshold,
#'   \item optionally revisits accepted shifts to prune overfitting using a small IC tolerance window,
#'   \item optionally computes per-shift IC weights by refitting the model with each shift removed.
#' }
#'
#' Models are fitted directly in multivariate trait space (no PCA), assuming a multi-rate
#' Brownian Motion with proportional VCV scaling across regimes. Extra arguments in \code{...}
#' are forwarded to \code{\link[mvMORPH]{mvgls}} (e.g., \code{method = "LL"} or
#' \code{method = "PL-LOOCV"}, \code{penalty}, \code{error = TRUE}, etc.).
#'
#' @param baseline_tree A rooted SIMMAP/\code{phylo} object representing the baseline
#'   (single-regime) tree. If not SIMMAP-initialized, it should already be painted to a
#'   single baseline state and have tip order matching \code{trait_data}.
#' @param trait_data A \code{matrix} or \code{data.frame} of trait values with row names
#'   matching \code{baseline_tree$tip.label} (same order).
#' @param formula Character formula passed to \code{mvgls}; defaults to \code{"trait_data ~ 1"}
#'   (intercept-only). Use \code{cbind()} for multivariate responses
#'   (e.g., \code{"cbind(t1, t2, ...) ~ 1"}).
#' @param min_descendant_tips Integer (≥1). Minimum number of tips required for an internal node
#'   to be considered as a candidate shift (forwarded to \code{\link{generatePaintedTrees}}).
#' @param num_cores Integer. Number of workers for parallel candidate scoring. Uses
#'   \code{future::plan(multicore)} on Unix and \code{future::plan(multisession)} on Windows.
#' @param ic_uncertainty_threshold Numeric (≥0). IC tolerance used in the optional post-search
#'   pruning step (\code{uncertainty = TRUE}). Shifts whose removal changes the current best IC
#'   by \eqn{\le} this value are pruned.
#' @param shift_acceptance_threshold Numeric (≥0). Minimum IC improvement
#'   (baseline − new) required to accept a candidate shift during the forward search.
#'   Larger values yield more conservative models.
#' @param uncertaintyweights Logical. If \code{TRUE}, compute per-shift IC weights serially by
#'   refitting the optimized model with each shift removed in turn. Exactly one of
#'   \code{uncertaintyweights} or \code{uncertaintyweights_par} must be \code{TRUE} to trigger
#'   IC-weight calculations.
#' @param uncertaintyweights_par Logical. As above, but compute per-shift IC weights in parallel
#'   using \pkg{future.apply}.
#' @param plot Logical. If \code{TRUE}, draw/update a SIMMAP plot as the search proceeds
#'   (requires \pkg{phytools}).
#' @param IC Character. Which information criterion to use, one of \code{"GIC"} or \code{"BIC"}
#'   (case-sensitive).
#' @param store_model_fit_history Logical. If \code{TRUE}, store a per-iteration record of fitted
#'   models, acceptance decisions, and IC values; warnings/errors are also kept.
#' @param ... Additional arguments passed to \code{\link[mvMORPH]{mvgls}} (e.g., \code{method},
#'   \code{penalty}, \code{target}, \code{error}, etc.).
#'
#' @details
#' \strong{Search outline.}
#' \enumerate{
#'   \item \emph{Baseline:} Fit \code{mvgls} on the baseline tree (single regime) to obtain the baseline IC.
#'   \item \emph{Candidates:} Build one-shift trees for eligible internal nodes
#'         (\code{\link{generatePaintedTrees}}); fit each with
#'         \code{\link{fitMvglsAndExtractGIC.formula}} or \code{\link{fitMvglsAndExtractBIC.formula}}
#'         and rank by ΔIC.
#'   \item \emph{Greedy add:} Add the top candidate, refit, and accept if
#'         ΔIC \eqn{\ge} \code{shift_acceptance_threshold}; continue down the ranked list.
#'   \item \emph{Optional IC weights:} If \code{uncertaintyweights} (or \code{uncertaintyweights_par})
#'         is \code{TRUE}, compute an IC weight for each accepted shift by refitting the final model with that
#'         shift removed and comparing the two ICs via \code{\link[mvMORPH]{aicw}}.
#' }
#'
#' \strong{Parallelization.} Candidate sub-model fits are distributed with \pkg{future} + \pkg{future.apply}.
#' On Unix, \code{multicore} is used; on Windows, \code{multisession}. A sequential plan is restored afterward.
#'
#' \strong{Plotting.} If \code{plot = TRUE}, trees are rendered with
#' \code{\link[phytools]{plotSimmap}()}; shift IDs are labeled with \code{\link[ape]{nodelabels}()}.
#'
#' \strong{Regime VCVs.} The returned \code{$VCVs} are extracted from the fitted multi-regime model via
#' \code{\link{extractRegimeVCVs}} and reflect regime-specific covariance
#' estimates (when \code{mvgls} is fitted under a PL/ML method).
#'
#' @return A named \code{list} with (at minimum):
#' \itemize{
#'   \item \code{user_input}: captured call (as a list) for reproducibility.
#'   \item \code{tree_no_uncertainty_transformed}: SIMMAP tree from the optimal (no-uncertainty) model
#'         on the transformed scale used internally by \code{mvgls}.
#'   \item \code{tree_no_uncertainty_untransformed}: same topology with original edge lengths restored.
#'   \item \code{model_no_uncertainty}: the final \code{mvgls} model object.
#'   \item \code{shift_nodes_no_uncertainty}: integer vector of accepted shift nodes.
#'   \item \code{optimal_ic}: final IC value; \code{baseline_ic}: baseline IC.
#'   \item \code{IC_used}: \code{"GIC"} or \code{"BIC"}; \code{num_candidates}: count of candidate one-shift models evaluated.
#'   \item \code{model_fit_history}: if \code{store_model_fit_history = TRUE}, a list of per-iteration fits
#'         and an \code{ic_acceptance_matrix} (IC value and acceptance flag per step).
#'   \item \code{VCVs}: named list of regime-specific VCV matrices extracted from the final model
#'         (penalized-likelihood estimates if PL was used).
#' }
#' Additional components appear conditionally:
#' \itemize{
#'   \item \code{ic_weights}: a \code{data.frame} of per-shift IC weights and evidence ratios when
#'         \code{uncertaintyweights} or \code{uncertaintyweights_par} is \code{TRUE}.
#'   \item \code{tree_uncertainty_transformed}, \code{tree_uncertainty_untransformed},
# #'         \code{model_uncertainty}, \code{shift_nodes_uncertainty}: returned when
# #'         \code{uncertainty = TRUE} and pruning removed at least one shift.
#'   \item \code{warnings}: character vector of warnings/errors encountered during fitting (if any).
#' }
#'
#' @section Convergence and robustness:
#' The search is greedy and may converge to a local optimum. Use a stricter
#' \code{shift_acceptance_threshold} and/or enable \code{uncertainty = TRUE} to reduce overfitting.
#' Re-run with different \code{min_descendant_tips} and IC choices (GIC vs BIC) to assess stability.
#'
#' @seealso
#' \code{\link{generatePaintedTrees}},
#' \code{\link{fitMvglsAndExtractGIC.formula}},
#' \code{\link{fitMvglsAndExtractBIC.formula}},
#' \code{\link{addShiftToModel}},
#' \code{\link{removeShiftFromTree}},
#' \code{\link{extractRegimeVCVs}};
#' packages: \pkg{mvMORPH}, \pkg{phytools}, \pkg{ape}, \pkg{future}, \pkg{future.apply}.
#'
#' @examples
#' \dontrun{
#' library(ape)
#' library(phytools)
#' set.seed(1)
#' tr <- pbtree(n = 80, scale = 1)
#' # Paint a single global baseline state "0"
#' base <- phytools::paintBranches(tr, edge = unique(tr$edge[,2]),
#'                                 state = "0", anc.state = "0")
#'
#' # Fake multivariate data (2 traits)
#' X <- cbind(trait1 = rnorm(Ntip(base)), trait2 = rnorm(Ntip(base)))
#' rownames(X) <- base$tip.label
#'
#' res <- searchOptimalConfiguration(
#'   baseline_tree = base,
#'   trait_data    = as.data.frame(X),
#'   formula       = "cbind(trait1, trait2) ~ 1",
#'   min_descendant_tips = 10,
#'   num_cores = 2,
#'   shift_acceptance_threshold = 10,  # conservative
#'   IC = "GIC",
#'   plot = FALSE
#' )
#'
#' res$shift_nodes_no_uncertainty
#' res$optimal_ic - res$baseline_ic
#' str(res$VCVs)
#' }
#' @param baseline_tree A \code{phylo} baseline tree.
#' @param trait_data Data frame of trait values with tip rownames.
#' @param formula Model formula, e.g. \code{cbind(y1, y2) ~ x}.
#' @param min_descendant_tips Integer minimum tips per candidate shift.
#' @param num_cores Integer, cores for parallel search.
#' @param ic_uncertainty_threshold Numeric, IC change threshold.
#' @param shift_acceptance_threshold Numeric, acceptance cutoff.
# #' @param uncertainty Character, how to treat uncertainty. #commented out for now
#' @param uncertaintyweights Numeric vector of weights.
#' @param uncertaintyweights_par List of parameters for weights.
# #' @param postorder_traversal Logical, traverse postorder. commented out for now
#' @param plot Logical, produce plots during search.
#' @param IC Character, information criterion, e.g. \code{"BIC"}.
#' @param store_model_fit_history Logical, keep fit history.
#' @param ... Passed to lower-level fitting functions.
#' @importFrom future plan multicore multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom mvMORPH mvgls GIC aicw
#' @importFrom phytools plotSimmap getStates
#' @importFrom ape nodelabels Ntip
#' @importFrom stats setNames BIC
#' @importFrom grDevices rainbow
#' @export
searchOptimalConfiguration <-
  function(baseline_tree,
           trait_data,
           formula = 'trait_data~1',
           min_descendant_tips,
           num_cores = 2,
           ic_uncertainty_threshold = 1.0,
           shift_acceptance_threshold = 1.0,
           #uncertainty = F,
           uncertaintyweights = F,
           uncertaintyweights_par = F,
           #postorder_traversal = F,
           plot = T,
           IC = 'GIC',
           store_model_fit_history = TRUE, ...) {

    # Capture user input
    user_input <- as.list(match.call())

    #generate initial set of painted candidate trees with shifts at each sub-node
    cat('Generating candidate shift models...\n')
    candidate_trees <- generatePaintedTrees(baseline_tree, min_descendant_tips)
    candidate_trees_shifts <- candidate_trees[-1]

    #fit the initial baseline model to the baseline tree with a global regime (state=0)
    cat('Fitting baseline model...\n')

    #select which information criterion to use
    if (IC != "GIC" && IC != "BIC") {
      stop("IC must be GIC or BIC")
    }
    if(IC=="GIC"){
      baseline_model <- fitMvglsAndExtractGIC.formula(formula, candidate_trees[[1]], trait_data, ...)
      baseline_ic <- baseline_model$GIC$GIC
    }
    if(IC=="BIC"){
      baseline_model <- fitMvglsAndExtractBIC.formula(formula, candidate_trees[[1]], trait_data, ...)
      baseline_ic <- baseline_model$BIC$BIC
    }
    cat(paste('Baseline ', IC, ':', round(baseline_ic, digits = 2), '\n'))

    #evaluate all of the candidate trees under GIC or BIC
    # Capture additional arguments into a list
    args_list <- list(...)

    if (.Platform$OS.type == "unix") {
      plan(multicore, workers = num_cores)
    } else {
      plan(multisession, workers = num_cores)
    }
    cat('Fitting sub-models in parallel...\n')
    candidate_results <-
      future.apply::future_lapply(candidate_trees_shifts, function(tree) {
        if (IC == "GIC") {
          do.call(fitMvglsAndExtractGIC.formula, c(list(formula, tree, trait_data), args_list))
        } else if (IC == "BIC") {
          do.call(fitMvglsAndExtractBIC.formula, c(list(formula, tree, trait_data), args_list))
        }
      }, future.seed = TRUE, future.scheduling = TRUE)
    plan(sequential)

    #generate the delta IC lists
    if (IC == "GIC"){
      delta_ic_list <- sapply(candidate_results, function(res) baseline_ic - res$GIC$GIC)
    } else if (IC == "BIC") {
      delta_ic_list <- sapply(candidate_results, function(res) baseline_ic - res$BIC$BIC)
    }

    cat('Sorting and evaluating shifts...\n')
    sorted_candidates <- candidate_trees_shifts[order(delta_ic_list, decreasing=T)]
    current_best_tree <- baseline_tree
    current_best_ic <- baseline_ic
    shift_id <- 0

    #plot current shift tree
    if(plot==T){
      plotSimmap(current_best_tree, ftype='off')
    }

    #in case the user wants postorder traversal for shift searching (default off)
    # if(postorder_traversal){
    #   print('Candidate shifts are sorted in Postorder')
    #   sorted_candidates <- (candidate_trees_shifts)
    #   current_best_tree <- baseline_tree
    #   current_best_ic <- baseline_ic
    #   shift_id <- 0
    # }

    shift_vec <- list() #initialize shift_vec
    model_with_shift_no_uncertainty<-NULL #initialize output
    best_tree_no_uncertainty<-NULL #initialize output
    # Initialize the list to collect warning messages
    warnings_list <- list()
    model_fit_history<-list() #new object to capture the full history of the search

    #Run the primary shift configuration search

    for (i in seq_along(sorted_candidates)) {
      shift_node_name <- names(sorted_candidates)[i]
      shift_node_number <- as.integer(sub("Node ", "", shift_node_name))
      percent_complete <- round((i / length(sorted_candidates)) * 100, 2)
      cat(paste('Evaluating shift at node:', shift_node_number, '-', percent_complete, '% complete', '\n'))

      add_shift_result <- addShiftToModel(current_best_tree, shift_node_number, shift_id)
      shifted_tree <- add_shift_result$tree
      shift_id <- add_shift_result$shift_id

      if(plot==T){
        nodelabels(text = shift_id, node=shift_node_number)
      }

      tryCatch({
        if (IC == "GIC") {
          model_with_shift <- withCallingHandlers(
            fitMvglsAndExtractGIC.formula(formula, shifted_tree, trait_data, ...),
            warning = function(w) {
              warning_message <- paste("Warning in evaluating shift at node", shift_node_number, ":", w$message)
              warning(warning_message)
              warnings_list[[length(warnings_list) + 1]] <<- warning_message
              invokeRestart("muffleWarning")
            }
          )
          new_ic <- model_with_shift$GIC$GIC
        } else if (IC == "BIC") {
          model_with_shift <- withCallingHandlers(
            fitMvglsAndExtractBIC.formula(formula, shifted_tree, trait_data, ...),
            warning = function(w) {
              warning_message <- paste("Warning in evaluating shift at node", shift_node_number, ":", w$message)
              warning(warning_message)
              warnings_list[[length(warnings_list) + 1]] <<- warning_message
              invokeRestart("muffleWarning")
            }
          )
          new_ic <- model_with_shift$BIC$BIC
        }

        # Calculate delta IC
        delta_ic <- current_best_ic - new_ic

        # Store model fit and acceptance status (including delta_ic)
        if (store_model_fit_history) {
          model_fit_history[[length(model_fit_history) + 1]] <- list(
            model = model_with_shift,
            accepted = delta_ic >= shift_acceptance_threshold,
            delta_ic = delta_ic
          )
        }

        # Decision logic (unchanged)
        if (delta_ic >= shift_acceptance_threshold) {
          current_best_tree <- shifted_tree
          current_best_ic <- new_ic
          cat(paste('Shift at node', shift_node_number, 'accepted. Updated', IC, ':', round(current_best_ic, digits = 2), 'Delta', IC, ':', round(delta_ic, digits = 2), '\n'))

          shift_vec[[length(shift_vec) + 1]] <- shift_node_number

          best_tree_no_uncertainty <- current_best_tree
          model_with_shift_no_uncertainty <- model_with_shift
        } else {
          cat(paste('Shift at node', shift_node_number, 'rejected. Delta', IC, ':', round(delta_ic, digits = 2), 'is less than threshold:', shift_acceptance_threshold, '\n'))
        }
      }, error = function(e) {
        # Handle errors (unchanged)
        warning_message <- paste("Error in evaluating shift at node", shift_node_number, ":", e$message)
        warning(warning_message)
        warnings_list[[length(warnings_list) + 1]] <<- warning_message

        # Also store the error in the model fit history
        if (store_model_fit_history) {
          model_fit_history[[length(model_fit_history) + 1]] <- list(
            model = NULL,
            accepted = FALSE,
            delta_ic = NA,
            error = e$message
          )
        }

      })

      if(plot==T){
        colorvec <- setNames(object = c('black', rainbow(length(unique(getStates(shifted_tree, type = 'both')))-1)),
                             nm = sort(as.numeric(unique(getStates(shifted_tree, type = 'both')))))
        plotSimmap(current_best_tree, colors=colorvec, ftype='off')
      }
    }

    #print(paste(shift_vec))
    shifts_no_uncertainty<-unlist(shift_vec)
    cat(paste("Shifts detected at nodes:", paste(shift_vec, collapse = ", "), '\n'))

    # If activated, this section removes shifts after re-evaluating
    model_without_shift <- NULL # Initialize as NULL
    # if (uncertainty) {
    #   cat('Post-search re-evaluation to reduce overfitting...')
    #   shift_nodes <- unlist(shift_vec)
    #   print(paste("Re-evaluating nodes", shift_nodes))
    #   shift_vec_uncertainty <- shift_nodes # Tracker (to remove shifts from)
    #
    #   root_node <- Ntip(baseline_tree) + 1
    #   for (shift_node_number in shift_nodes) {
    #     if (shift_node_number != root_node) {
    #       cat(paste('Re-evaluating shift at node:', shift_node_number))
    #       tree_without_shift <- removeShiftFromTree(current_best_tree, shift_node_number)
    #
    #       tryCatch({
    #         # Evaluate the model without the current shift using the selected IC
    #         if (IC == "GIC") {
    #           temp_model <- fitMvglsAndExtractGIC.formula(formula, tree_without_shift, trait_data, ...)
    #           ic_without_shift <- temp_model$GIC$GIC
    #         } else if (IC == "BIC") {
    #           temp_model <- fitMvglsAndExtractBIC.formula(formula, tree_without_shift, trait_data, ...)
    #           ic_without_shift <- temp_model$BIC$BIC
    #         }
    #
    #         if (abs(current_best_ic - ic_without_shift) <= ic_uncertainty_threshold) {
    #           current_best_tree <- tree_without_shift
    #           current_best_ic <- ic_without_shift
    #           cat(paste('Shift at node', shift_node_number, 'removed. Updated', IC, ':', round(current_best_ic, digits=2)))
    #           shift_vec_uncertainty <- shift_vec_uncertainty[!shift_vec_uncertainty == shift_node_number]
    #           model_without_shift <- temp_model # Update only if the condition is met
    #         }
    #       }, error = function(e) {
    #         warning(paste("Error in re-evaluating shift at node", shift_node_number, ":", e$message))
    #       })
    #     } else {
    #       cat(paste('Skipping re-evaluation for the root node:', shift_node_number))
    #     }
    #   }
    # }

    # New Section for Calculating Information Criterion Weights Post Optimization
    ic_weights_df <- NA  # Initialize the results vector
    if (xor(uncertaintyweights, uncertaintyweights_par)) {
      if (uncertaintyweights) {
        if (length(unlist(shift_vec)) > 0) {
          cat('Calculating IC weights for initially identified shifts...\n')
          ic_weights_df <- data.frame(node = integer(), ic_with_shift = numeric(), ic_without_shift = numeric(), delta_ic = numeric(), ic_weight = numeric())

          # Retrieve the IC of the optimized model before uncertainty analysis
          original_ic <- if (IC == "GIC") model_with_shift_no_uncertainty$GIC$GIC else model_with_shift_no_uncertainty$BIC$BIC

          cat(paste("Considering", length(shift_vec), "shifts in the candidate set",'\n'))
          cat(paste("There are", length(unique(getStates(best_tree_no_uncertainty, type = 'both'))) - 1, "shifts in the mapped tree",'\n'))

          for (shift_node_number in unlist(shift_vec)) {
            cat(paste('Re-estimating model without shift at node:', shift_node_number, '\n'))

            # Remove the shift temporarily from best_tree_no_uncertainty and re-estimate the IC
            tree_without_current_shift <- removeShiftFromTree(best_tree_no_uncertainty, shift_node_number)
            model_without_current_shift_function <- if (IC == "GIC") fitMvglsAndExtractGIC.formula else fitMvglsAndExtractBIC.formula
            model_without_current_shift <- model_without_current_shift_function(formula, tree_without_current_shift, trait_data, ...)
            ic_without_current_shift <- if (IC == "GIC") model_without_current_shift$GIC$GIC else model_without_current_shift$BIC$BIC

            # Calculate the difference in IC
            delta_ic <- original_ic - ic_without_current_shift

            # Use aicw function to calculate the IC weight
            ic_weights <- aicw(c(original_ic, ic_without_current_shift))$aicweights
            ic_weight <- ic_weights[1]  # First element is the weight for the model with the shift

            cat(paste("IC weight for the shift is", round(ic_weight, digits=2)))
            ic_weights_df <- rbind(
              ic_weights_df,
              data.frame(
                node = shift_node_number,
                ic_with_shift = original_ic,
                ic_without_shift = ic_without_current_shift,
                delta_ic = delta_ic,
                ic_weight = ic_weight
              )
            )
          }
        } else {
          cat("No shifts were detected in the initial search, skipping IC weights calculation.\n")
          ic_weights_df <- NA
        }
      }

      if (uncertaintyweights_par) {
        if (length(unlist(shift_vec)) > 0) {
          cat('Calculating IC weights for initially identified shifts in parallel...\n')
          ic_weights_df <- data.frame(node = integer(), ic_with_shift = numeric(), ic_without_shift = numeric(), delta_ic = numeric(), ic_weight = numeric())

          # Retrieve the IC of the optimized model before uncertainty analysis
          original_ic <- if (IC == "GIC") model_with_shift_no_uncertainty$GIC$GIC else model_with_shift_no_uncertainty$BIC$BIC

          # Prepare a list of trees with shifts removed
          shift_removed_trees <- lapply(unlist(shift_vec), function(shift_node_number) {
            return(removeShiftFromTree(best_tree_no_uncertainty, shift_node_number))
          })

          # Enable parallel processing
          if (.Platform$OS.type == "unix") {
            plan(multicore, workers = num_cores)
          } else {
            plan(multisession, workers = num_cores)
          }
          # Capture additional arguments into a list
          args_list <- list(...)
          # Calculate IC weights in parallel
          ic_results <- future.apply::future_lapply(shift_removed_trees, function(tree) {
            model_function <- if (IC == "GIC") fitMvglsAndExtractGIC.formula else fitMvglsAndExtractBIC.formula
            model_without_shift <- do.call(model_function, c(list(formula, tree, trait_data), args_list))
            ic_without_shift <- if (IC == "GIC") model_without_shift$GIC$GIC else model_without_shift$BIC$BIC
            delta_ic <- original_ic - ic_without_shift
            ic_weights <- aicw(c(original_ic, ic_without_shift))$aicweights
            return(c(ic_weight_withshift = ic_weights[1], ic_weight_withoutshift = ic_weights[2], delta_ic = delta_ic))
          }, future.seed = TRUE, future.scheduling = T)

          # Reset to sequential plan
          future::plan(future::sequential)

          # Add results to the dataframe
          for (i in seq_along(shift_removed_trees)) {
            shift_node_number <- unlist(shift_vec)[i]
            ic_res <- ic_results[[i]]
            ic_weights_df <- rbind(ic_weights_df, data.frame(
              node = shift_node_number,
              ic_with_shift = original_ic,
              ic_without_shift = original_ic - ic_res['delta_ic'],
              delta_ic = ic_res['delta_ic'],
              ic_weight_withshift = ic_res['ic_weight_withshift'],
              ic_weight_withoutshift = ic_res['ic_weight_withoutshift'],
              evidence_ratio = ic_res['ic_weight_withshift'] / ic_res['ic_weight_withoutshift']
            ))
          }
        } else {
          cat("No shifts were detected in the initial search, skipping IC weights calculation.\n")
          ic_weights_df <- NA
        }
      }
    } else {
      print("Only one of uncertaintyweights or uncertaintyweights_par can be set to TRUE")
    }

    # Print statements for the optimal configuration and delta GIC/BIC
    if (IC == "GIC") {
      cat(paste('Optimal configuration found with GIC:', round(current_best_ic, digits=2), '\n'))
      cat(paste('Global Delta GIC:', round(baseline_ic, digits=2) - round(current_best_ic, digits=2), '\n'))
    } else if (IC == "BIC") {
      cat(paste('Optimal configuration found with BIC:', round(current_best_ic, digits=2), '\n'))
      cat(paste('Global Delta BIC:', round(baseline_ic, digits=2) - round(current_best_ic, digits=2), '\n'))
    }

    # Assembling the results
    {
      # Taking directly from the model fit to ensure the correct tree is transferred to the results
      # if(uncertainty) {
      #   opt_uncertainty_transformed <- model_without_shift$model$corrSt$phy
      #   opt_uncertainty_untransformed <- model_without_shift$model$corrSt$phy
      #   opt_uncertainty_untransformed$edge.length <- best_tree_no_uncertainty$edge.length
      # }

      if (length(shifts_no_uncertainty) > 0L) {
        # use the accepted-shift model/tree
        opt_nouncertainty_transformed   <- model_with_shift_no_uncertainty$model$corrSt$phy
        opt_nouncertainty_untransformed <- model_with_shift_no_uncertainty$model$corrSt$phy
        opt_nouncertainty_untransformed$edge.length <- best_tree_no_uncertainty$edge.length
      } else {
        # fallback to baseline model/tree when no shifts were accepted
        opt_nouncertainty_transformed   <- baseline_model$model$corrSt$phy
        opt_nouncertainty_untransformed <- baseline_model$model$corrSt$phy
        opt_nouncertainty_untransformed$edge.length <- candidate_trees[[1]]$edge.length
      }

      # Create the main list that will always be returned
      result_list <- list(
        user_input = user_input,
        tree_no_uncertainty_transformed = opt_nouncertainty_transformed,
        tree_no_uncertainty_untransformed = opt_nouncertainty_untransformed,
        model_no_uncertainty = model_with_shift_no_uncertainty$model,
        shift_nodes_no_uncertainty = shifts_no_uncertainty,
        optimal_ic = current_best_ic,
        baseline_ic = baseline_ic,
        IC_used = IC,
        num_candidates = length(sorted_candidates),
        model_fit_history = model_fit_history
      )

      # Create the IC and acceptance matrix from the model fit history
      if (store_model_fit_history) {
        ic_acceptance_matrix <- do.call(rbind, lapply(model_fit_history, function(x) {
          if (is.null(x$model)) {
            c(NA, x$accepted)
          } else {
            if (IC == "GIC") c(x$model$GIC$GIC, x$accepted) else c(x$model$BIC$BIC, x$accepted)
          }
        }))

        result_list$model_fit_history <- list(
          fits = model_fit_history,
          ic_acceptance_matrix = ic_acceptance_matrix
        )
      }

      # Generate the VCVs per regime from the overall model fit
      model_output <- result_list$model_no_uncertainty
      result_list$VCVs <- extractRegimeVCVs(model_output)

      # Add the ic_weights to the list conditionally
      if (uncertaintyweights | uncertaintyweights_par) {
        result_list$ic_weights <- ic_weights_df
      }

      # Add the uncertainty outputs to the list conditionally
      # if (uncertainty) {
      #   result_list$tree_uncertainty_transformed <- opt_uncertainty_transformed
      #   result_list$tree_uncertainty_untransformed <- opt_uncertainty_untransformed
      #   result_list$model_uncertainty <- model_without_shift$model
      #   result_list$shift_nodes_uncertainty <- shift_vec_uncertainty
      # }

      # Add the warnings to the output list conditionally
      if(length(warnings_list) > 0) {
        result_list$warnings <- warnings_list
      }

    }

    return(result_list)

  }
#... is optional arguments to be passed to mvgls

