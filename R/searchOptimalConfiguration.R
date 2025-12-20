#' Search for an Optimal Multi-Regime (Shift) Configuration on a Phylogeny
#'
#' @description
#' Greedy, stepwise search for evolutionary regime shifts on a SIMMAP-style phylogeny
#' using multivariate \code{mvgls} fits from \pkg{mvMORPH}. The routine:
#' \enumerate{
#'   \item builds one-shift candidate trees for all internal nodes meeting a tip-size threshold
#'         (via \code{generatePaintedTrees}),
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
#' @param trait_data A \code{matrix} or \code{data.frame} of continuous trait values with row
#'   names matching \code{baseline_tree$tip.label} (same order). For the default
#'   \code{formula = "trait_data ~ 1"}, \code{trait_data} is typically supplied as a numeric
#'   matrix so that the multivariate response is interpreted correctly by \code{mvgls()}.
#'   When using more general formulas (e.g., pGLS-style models), a \code{data.frame} with
#'   named columns can be used instead.
#' @param formula Character formula passed to \code{mvgls}. Defaults to
#'   \code{"trait_data ~ 1"}, which fits an intercept-only model treating the supplied
#'   multivariate trait matrix as the response. This is the appropriate choice for most
#'   morphometric data where there are no predictor variables. For more general models,
#'   \code{formula} can reference subsets of \code{trait_data} explicitly, for example
#'   \code{"trait_data[, 1:5] ~ 1"} to treat columns 1–5 as a multivariate response, or
#'   \code{"trait_data[, 1:5] ~ trait_data[, 6]"} to fit a multivariate pGLS with column 6
#'   as a predictor.
#' @param min_descendant_tips Integer (\eqn{\ge}1). Minimum number of tips required for an internal node
#'   to be considered as a candidate shift (forwarded to \code{generatePaintedTrees}). Larger values
#'   reduce the number of candidate shifts by excluding very small clades. For empirical datasets,
#'   values around \code{10} are a reasonable starting choice and can be tuned in sensitivity analyses.
#' @param num_cores Integer. Number of workers for parallel candidate scoring. Uses
#'   \code{future::plan(multicore)} on Unix outside \code{RStudio}; otherwise uses
#'   \code{future::plan(multisession)}. During the parallel candidate-scoring blocks, BLAS/OpenMP
#'   threads are capped to 1 (per worker) to avoid CPU oversubscription.
#' @param ic_uncertainty_threshold Numeric (\eqn{\ge}0). Reserved for future development
#'   in post-search pruning and uncertainty analysis; currently not used by
#'   \code{searchOptimalConfiguration()}.
#' @param shift_acceptance_threshold Numeric (\eqn{\ge}0). Minimum IC improvement
#'   (baseline - new) required to accept a candidate shift during the forward search.
#'   Larger values yield more conservative models. For analyses based on the Generalized
#'   Information Criterion (\code{"GIC"}), a threshold on the order of \code{20} units is a
#'   conservative choice that tends to admit only strongly supported shifts. Simulation
#'   studies (Berv et al., in preparation) suggest that this choice yields good balanced
#'   accuracy between detecting true shifts and avoiding false positives, but users should
#'   explore alternative thresholds in sensitivity analyses for their own datasets.
#' @param uncertaintyweights Logical. If \code{TRUE}, compute per-shift IC weights serially by
#'   refitting the optimized model with each shift removed in turn. Exactly one of
#'   \code{uncertaintyweights} or \code{uncertaintyweights_par} must be \code{TRUE} to trigger
#'   IC-weight calculations; setting both to \code{TRUE} will result in an error. When enabled,
#'   the per-shift weights are returned in the \code{$ic_weights} component of the result.
#' @param uncertaintyweights_par Logical. As above, but compute per-shift IC weights in parallel
#'   using \pkg{future.apply}. Exactly one of \code{uncertaintyweights} or
#'   \code{uncertaintyweights_par} must be \code{TRUE} to trigger IC-weight calculations.
#' @param plot Logical. If \code{TRUE}, draw/update a SIMMAP plot as the search proceeds
#'   (requires \pkg{phytools}).
#' @param IC Character. Which information criterion to use, one of \code{"GIC"} or \code{"BIC"}
#'   (case-sensitive).
#' @param store_model_fit_history Logical. If \code{TRUE}, store a per-iteration record of fitted
#'   models, acceptance decisions, and IC values. To keep memory usage low during the search,
#'   per-iteration results are written to a temporary directory (\code{tempdir()}) and read back
#'   into memory at the end of the run.
#' @param verbose Logical. If \code{TRUE}, report progress during candidate generation and model
#'   fitting. By default, progress is emitted via \code{message()}. When \code{plot = TRUE} in an
#'   interactive \code{RStudio} session, progress is written via \code{cat()} so it remains visible
#'   while plots are updating. Set to \code{FALSE} to run quietly (default). Use
#'   \code{suppressMessages()} (and \code{capture.output()} if needed) to silence or capture output.
#' @param ... Additional arguments passed to \code{\link[mvMORPH]{mvgls}} (e.g., \code{method},
#'   \code{penalty}, \code{target}, \code{error}, etc.).
#'
#' @details
#' \strong{Input requirements.}
#' \itemize{
#'   \item \emph{Tree:} \code{baseline_tree} should be a rooted \code{phylo} (or SIMMAP-style) tree
#'         with branch lengths interpreted in units of time. An ultrametric tree is not required.
#'   \item \emph{Trait data alignment:} \code{rownames(trait_data)} must match
#'         \code{baseline_tree$tip.label} in both names and order; any tips without data should be
#'         pruned beforehand.
#'   \item \emph{Data type:} \code{trait_data} is typically a numeric matrix of continuous traits;
#'         high-dimensional settings (p \eqn{\ge} n) are supported via penalized-likelihood
#'         \code{mvgls()} fits.
#' }
#'
#' \strong{Search outline.}
#' \enumerate{
#'   \item \emph{Baseline:} Fit \code{mvgls} on the baseline tree (single regime) to obtain the baseline IC.
#'   \item \emph{Candidates:} Build one-shift trees for eligible internal nodes
#'         (\code{generatePaintedTrees}); fit each with
#'         \code{fitMvglsAndExtractGIC.formula} or \code{fitMvglsAndExtractBIC.formula}
#'         (internal helpers; not exported) and rank by \eqn{\Delta}IC.
#'   \item \emph{Greedy add:} Add the top candidate, refit, and accept if
#'         \eqn{\Delta}IC \eqn{\ge} \code{shift_acceptance_threshold}; continue down the ranked list.
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
#' \code{extractRegimeVCVs} and reflect regime-specific covariance
#' estimates (when \code{mvgls} is fitted under a PL/ML method).
#'
#' For high-dimensional trait datasets (p \eqn{\ge} n), penalized-likelihood settings in
#' \code{mvgls()} are often required for stable estimation. In practice, methods such as
#' \code{method = "LL"} or \code{method = "H&L"} combined with appropriate penalties (e.g.,
#' ridge-type penalties) have proven effective for intercept-only multivariate Brownian
#' motion models, as illustrated in the package vignettes. Users should consult the
#' \pkg{mvMORPH} documentation for details on available methods and penalties and
#' tune these choices to the structure of their data.
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
#'         (loaded from temporary files written during the run) and an \code{ic_acceptance_matrix}
#'         (IC value and acceptance flag per step).
#'   \item \code{VCVs}: named list of regime-specific VCV matrices extracted from the final model
#'         (penalized-likelihood estimates if PL was used).
#' }
#' Additional components appear conditionally:
#' \itemize{
#'   \item \code{ic_weights}: a \code{data.frame} of per-shift IC weights and evidence ratios when
#'         \code{uncertaintyweights} or \code{uncertaintyweights_par} is \code{TRUE}.
#'   \item \code{warnings}: character vector of warnings/errors encountered during fitting (if any).
#' }
#'
#' @section Convergence and robustness:
#' The search is greedy and may converge to a local optimum. Use a stricter
#' \code{shift_acceptance_threshold} to reduce overfitting, and re-run the search
#' with different \code{min_descendant_tips} and IC choices (\code{"GIC"} vs \code{"BIC"})
#' to assess stability of the inferred shifts. For a given run, the optional IC-weight
#' calculations (\code{uncertaintyweights} or \code{uncertaintyweights_par}) can be used
#' to quantify support for individual shifts. It is often helpful to repeat the analysis
#' under slightly different settings (e.g., thresholds or candidate-size constraints) and
#' compare the resulting sets of inferred shifts.
#'
#' @seealso
#' \code{\link[mvMORPH]{mvgls}}, \code{\link[mvMORPH]{GIC}}, \code{\link[stats]{BIC}},
#' \code{\link{plot_ic_acceptance_matrix}} for visualizing IC trajectories and shift
#' acceptance decisions, and \code{\link{generateViridisColorScale}} for mapping
#' regime-specific rates or parameters to a viridis color scale when plotting trees;
#' packages: \pkg{mvMORPH}, \pkg{future}, \pkg{future.apply}, \pkg{phytools}, \pkg{ape}.
#'
#' @note
#' Internally, this routine coordinates multiple unexported helper functions:
#' \code{generatePaintedTrees}, \code{fitMvglsAndExtractGIC.formula},
#' \code{fitMvglsAndExtractBIC.formula}, \code{addShiftToModel},
#' \code{removeShiftFromTree}, and \code{extractRegimeVCVs}. Through these,
#' it may also invoke lower-level utilities such as \code{paintSubTree_mod}
#' and \code{paintSubTree_removeShift}. These helpers are internal
#' implementation details and are not part of the public API.
#'
#' @examples
#' library(ape)
#' library(phytools)
#' library(mvMORPH)
#' set.seed(1)
#'
#' # Simulate a tree
#' tr <- pbtree(n = 50, scale = 1)
#'
#' # Define two regimes: "0" (baseline) and "1" (high-rate) on a subset of tips
#' states <- setNames(rep("0", Ntip(tr)), tr$tip.label)
#' high_clade_tips <- tr$tip.label[1:20]
#' states[high_clade_tips] <- "1"
#'
#' # Make a SIMMAP tree for the BMM simulation
#' simmap <- phytools::make.simmap(tr, states, model = "ER", nsim = 1)
#'
#' # Simulate traits under a BMM model with ~10x higher rate in regime "1"
#' sigma <- list(
#'   "0" = diag(0.1, 2),
#'   "1" = diag(1.0, 2)
#' )
#' theta <- c(0, 0)
#'
#' sim <- mvMORPH::mvSIM(
#'   tree  = simmap,
#'   nsim  = 1,
#'   model = "BMM",
#'   param = list(
#'     ntraits = 2,
#'     sigma   = sigma,
#'     theta   = theta
#'   )
#' )
#'
#' # mvSIM returns either a matrix or a list of matrices depending on mvMORPH version
#' X <- if (is.list(sim)) sim[[1]] else sim
#' rownames(X) <- simmap$tip.label
#'
#' # Run the search on the unpainted tree (single baseline regime)
#' res <- searchOptimalConfiguration(
#'   baseline_tree              = as.phylo(simmap),
#'   trait_data                 = X,
#'   formula                    = "trait_data ~ 1",
#'   min_descendant_tips        = 10,
#'   num_cores                  = 1,   # keep it simple / CRAN-safe
#'   shift_acceptance_threshold = 20,  # conservative GIC threshold
#'   IC                         = "GIC",
#'   plot                       = FALSE,
#'   store_model_fit_history    = FALSE,
#'   verbose                    = FALSE
#' )
#'
#' res$shift_nodes_no_uncertainty
#' res$optimal_ic - res$baseline_ic
#' str(res$VCVs)
#' @importFrom future plan multicore multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom mvMORPH mvgls GIC aicw
#' @importFrom phytools plotSimmap getStates paintSubTree
#' @importFrom ape nodelabels Ntip as.phylo
#' @importFrom stats setNames BIC
#' @importFrom grDevices rainbow
#' @export
searchOptimalConfiguration <-
  function(baseline_tree,
           trait_data,
           formula = "trait_data ~ 1",
           min_descendant_tips,
           num_cores = 2,
           ic_uncertainty_threshold = 1.0,
           shift_acceptance_threshold = 1.0,
           #uncertainty = FALSE,
           uncertaintyweights = FALSE,
           uncertaintyweights_par = FALSE,
           #postorder_traversal = FALSE,
           plot = FALSE,
           IC = 'GIC',
           store_model_fit_history = TRUE,
           verbose = FALSE,
           ...) {

    #capturing global option for verbose to pass to internal helpers
    old_verbose_opt <- getOption("bifrost.verbose")
    on.exit(options(bifrost.verbose = old_verbose_opt), add = TRUE)
    options(bifrost.verbose = isTRUE(verbose))

    #internal helper for capturing verbose output
    .progress <- function(...) {
      if (!isTRUE(verbose)) return(invisible(NULL))

      txt <- sprintf(...)

      # In RStudio interactive plotting, use stdout (cat) because messages can get swallowed.
      if (isTRUE(plot) && interactive() && identical(Sys.getenv("RSTUDIO"), "1")) {
        cat(txt, "\n", sep = "")
        if (sink.number(type = "output") == 0) utils::flush.console()
      } else {
        message(txt)
      }

      invisible(NULL)
    }

    # Capture user input
    user_input <- as.list(match.call())

    # Input tree should be painted SIMMAP tree with global state zero
    baseline_tree <- paintSubTree(((as.phylo(baseline_tree))),
                                  node = length(baseline_tree$tip.label) + 1,
                                  state = 0)

    #generate initial set of painted candidate trees with shifts at each sub-node
    .progress("%s", "Generating candidate shift models...")
    candidate_trees <- generatePaintedTrees(baseline_tree, min_descendant_tips)
    candidate_trees_shifts <- candidate_trees[-1]

    #fit the initial baseline model to the baseline tree with a global regime (state=0)
    .progress("%s", "Fitting baseline model...")

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
    .progress("Baseline %s: %.2f", IC, baseline_ic)

    #evaluate all of the candidate trees under GIC or BIC
    # Capture additional arguments into a list
    args_list <- list(...)

    is_rstudio <- identical(Sys.getenv("RSTUDIO"), "1")

    # --- local helper: run a future_lapply with safe plan + capped BLAS/OpenMP threads ---
    .run_future_lapply_safe <- function(X, FUN, workers, is_rstudio_flag) {
      .thread_vars <- c(
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS"
      )
      .old_threads <- Sys.getenv(.thread_vars, unset = NA_character_)

      Sys.setenv(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        NUMEXPR_NUM_THREADS = "1"
      )

      .restore_threads <- function() {
        for (nm in .thread_vars) {
          val <- .old_threads[[nm]]
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
          .restore_threads()
        }
      )
    }

    .progress("%s", "Fitting sub-models in parallel...")

    candidate_results <- .run_future_lapply_safe(
      candidate_trees_shifts,
      function(tree) {
        if (IC == "GIC") {
          do.call(fitMvglsAndExtractGIC.formula, c(list(formula, tree, trait_data), args_list))
        } else if (IC == "BIC") {
          do.call(fitMvglsAndExtractBIC.formula, c(list(formula, tree, trait_data), args_list))
        }
      },
      workers = num_cores,
      is_rstudio_flag = is_rstudio
    )

    #generate the delta IC lists
    if (IC == "GIC"){
      delta_ic_list <- sapply(candidate_results, function(res) baseline_ic - res$GIC$GIC)
    } else if (IC == "BIC") {
      delta_ic_list <- sapply(candidate_results, function(res) baseline_ic - res$BIC$BIC)
    }

    .progress("%s", "Sorting and evaluating shifts...")
    sorted_candidates <- candidate_trees_shifts[order(delta_ic_list, decreasing = TRUE)]
    current_best_tree <- baseline_tree
    current_best_ic <- baseline_ic
    shift_id <- 0

    #plot current shift tree
    if (isTRUE(plot)) {
      .progress("%s", "Plotting initial tree...")
      plotSimmap(current_best_tree, ftype = "off")
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

    # Where to store on-disk history (CRAN-safe temp location)
    sub_dir <- NULL

    if (isTRUE(store_model_fit_history)) {
      base_dir <- file.path(tempdir(), "bifrost_fit_history")
      if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

      # Dated (and uniquified) subdirectory per run
      date_str <- format(Sys.Date(), "%Y-%m-%d")
      sub_dir <- file.path(base_dir, date_str)
      counter <- 1L
      while (dir.exists(sub_dir)) {
        counter <- counter + 1L
        sub_dir <- file.path(base_dir, paste0(date_str, "_", counter))
      }
      dir.create(sub_dir, recursive = TRUE)
    }

    # In-memory accumulator (lightweight; actual fits stored on disk)
    model_fit_history <- list()

    #Run the primary shift configuration search

    for (i in seq_along(sorted_candidates)) {
      shift_node_name <- names(sorted_candidates)[i]
      shift_node_number <- as.integer(sub("Node ", "", shift_node_name))
      percent_complete <- round((i / length(sorted_candidates)) * 100, 2)
      .progress("Evaluating shift at node %d (%.2f%% complete)", shift_node_number, percent_complete)

      add_shift_result <- addShiftToModel(current_best_tree, shift_node_number, shift_id)
      shifted_tree <- add_shift_result$tree
      shift_id <- add_shift_result$shift_id

      if(plot == TRUE){
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
          model_fit_history <- list(
            model = model_with_shift,
            accepted = delta_ic >= shift_acceptance_threshold,
            delta_ic = delta_ic
          )
        }

        # Decision logic (unchanged)
        if (delta_ic >= shift_acceptance_threshold) {
          current_best_tree <- shifted_tree
          current_best_ic <- new_ic
          .progress(
            "Shift at node %d accepted. Updated %s: %.2f; Delta %s: %.2f",
            shift_node_number, IC, current_best_ic, IC, delta_ic
          )

          shift_vec[[length(shift_vec) + 1]] <- shift_node_number

          best_tree_no_uncertainty <- current_best_tree
          model_with_shift_no_uncertainty <- model_with_shift
        } else {
          .progress(
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
          model_fit_history<- list(
            model = NULL,
            accepted = FALSE,
            delta_ic = NA,
            error = e$message
          )
        }

      })
      if (isTRUE(store_model_fit_history) && !is.null(sub_dir)) {
        iteration_num <- i + 1L
        saveRDS(
          model_fit_history,
          file = file.path(sub_dir, paste0("iteration_", iteration_num, ".rds"))
        )
      }
      if(plot == TRUE){
        colorvec <- setNames(object = c('black', rainbow(length(unique(getStates(shifted_tree, type = 'both')))-1)),
                             nm = sort(as.numeric(unique(getStates(shifted_tree, type = 'both')))))
        plotSimmap(current_best_tree, colors=colorvec, ftype='off')
      }
    }

    #print(paste(shift_vec))
    shifts_no_uncertainty<-unlist(shift_vec)
    .progress("Shifts detected at nodes: %s", paste(shift_vec, collapse = ", "))

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
    .empty_ic_weights_df <- data.frame(
      node = integer(),
      ic_with_shift = numeric(),
      ic_without_shift = numeric(),
      delta_ic = numeric(),
      ic_weight_withshift = numeric(),
      ic_weight_withoutshift = numeric(),
      evidence_ratio = numeric()
    )

    ic_weights_df <- .empty_ic_weights_df  # default empty

    if (xor(uncertaintyweights, uncertaintyweights_par)) {

      # If no shifts, return empty df (consistent in both modes)
      if (length(unlist(shift_vec)) == 0) {
        .progress("%s", "No shifts were detected in the initial search; skipping IC weights calculation.")
        ic_weights_df <- .empty_ic_weights_df

      } else {

        # Retrieve the IC of the optimized model before uncertainty analysis
        original_ic <- if (IC == "GIC") {
          model_with_shift_no_uncertainty$GIC$GIC
        } else {
          model_with_shift_no_uncertainty$BIC$BIC
        }

        .progress("Considering %d shifts in the candidate set", length(shift_vec))
        .progress(
          "There are %d shifts in the mapped tree",
          length(unique(getStates(best_tree_no_uncertainty, type = "both"))) - 1
        )

        if (uncertaintyweights) {
          .progress("%s", "Calculating IC weights for initially identified shifts...")

          ic_weights_df <- .empty_ic_weights_df

          for (shift_node_number in unlist(shift_vec)) {
            .progress("Re-estimating model without shift at node %d", shift_node_number)

            tree_without_current_shift <- removeShiftFromTree(best_tree_no_uncertainty, shift_node_number)
            model_fun <- if (IC == "GIC") fitMvglsAndExtractGIC.formula else fitMvglsAndExtractBIC.formula
            model_without_current_shift <- model_fun(formula, tree_without_current_shift, trait_data, ...)

            ic_without_current_shift <- if (IC == "GIC") {
              model_without_current_shift$GIC$GIC
            } else {
              model_without_current_shift$BIC$BIC
            }

            delta_ic <- original_ic - ic_without_current_shift

            icw <- aicw(c(original_ic, ic_without_current_shift))$aicweights
            w_with <- icw[1]
            w_without <- icw[2]
            er <- w_with / w_without

            .progress("IC weight for the shift is %.2f", w_with)

            ic_weights_df <- rbind(
              ic_weights_df,
              data.frame(
                node = shift_node_number,
                ic_with_shift = original_ic,
                ic_without_shift = ic_without_current_shift,
                delta_ic = delta_ic,
                ic_weight_withshift = w_with,
                ic_weight_withoutshift = w_without,
                evidence_ratio = er
              )
            )
          }
        }

        if (uncertaintyweights_par) {
          .progress("%s", "Calculating IC weights for initially identified shifts in parallel...")

          ic_weights_df <- .empty_ic_weights_df

          shift_removed_trees <- lapply(unlist(shift_vec), function(shift_node_number) {
            removeShiftFromTree(best_tree_no_uncertainty, shift_node_number)
          })

          ic_results <- .run_future_lapply_safe(
            shift_removed_trees,
            function(tree) {
              model_fun <- if (IC == "GIC") fitMvglsAndExtractGIC.formula else fitMvglsAndExtractBIC.formula
              model_without_shift <- do.call(model_fun, c(list(formula, tree, trait_data), args_list))

              ic_without_shift <- if (IC == "GIC") model_without_shift$GIC$GIC else model_without_shift$BIC$BIC
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
            d_ic <- as.numeric(ic_res[["delta_ic"]])
            w_with <- as.numeric(ic_res[["ic_weight_withshift"]])
            w_without <- as.numeric(ic_res[["ic_weight_withoutshift"]])
            er <- w_with / w_without

            ic_weights_df <- rbind(
              ic_weights_df,
              data.frame(
                node = shift_node_number,
                ic_with_shift = original_ic,
                ic_without_shift = ic_without,
                delta_ic = d_ic,
                ic_weight_withshift = w_with,
                ic_weight_withoutshift = w_without,
                evidence_ratio = er
              )
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

    # Print statements for the optimal configuration and delta GIC/BIC
    if (IC == "GIC") {
      .progress("Optimal configuration found with GIC: %.2f", current_best_ic)
      .progress("Global Delta GIC: %.2f", baseline_ic - current_best_ic)
    } else if (IC == "BIC") {
      .progress("Optimal configuration found with BIC: %.2f", current_best_ic)
      .progress("Global Delta BIC: %.2f", baseline_ic - current_best_ic)
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
        model_no_uncertainty = if (length(shifts_no_uncertainty) > 0L) {
          model_with_shift_no_uncertainty$model
        } else {
          baseline_model$model
        },
        shift_nodes_no_uncertainty = shifts_no_uncertainty,
        optimal_ic = current_best_ic,
        baseline_ic = baseline_ic,
        IC_used = IC,
        num_candidates = length(sorted_candidates),
        model_fit_history = model_fit_history
      )

      # Create the IC and acceptance matrix from the model fit history
      if (isTRUE(store_model_fit_history) && !is.null(sub_dir)) {
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

