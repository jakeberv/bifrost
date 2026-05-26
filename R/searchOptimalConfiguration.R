#' Search for an Optimal Multi-Regime (Shift) Configuration on a Phylogeny
#'
#' @description
#' Greedy, stepwise search for evolutionary regime shifts on a phylogeny
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
#' are forwarded to \code{\link[mvMORPH]{mvgls}}. In practice, \code{method} and
#' \code{error} are often the most important of these: the package vignettes use
#' \code{method = "H&L"} for intercept-only, high-dimensional response matrices and
#' \code{method = "LL"} for formula-based searches with predictors, while
#' \code{error = TRUE} asks \code{mvgls()} to estimate a nuisance measurement-error
#' (intraspecific-variance) term from the data.
#'
#' @param baseline_tree A rooted \code{phylo} (or SIMMAP/\code{phylo}) object representing
#'   the starting tree. It does not need to already be painted: the function coerces the
#'   input to a \code{phylo} object and internally paints a single baseline state at the root
#'   before generating candidate shift configurations. Tip labels must match
#'   \code{trait_data}.
#' @param trait_data A \code{matrix} or \code{data.frame} of continuous trait values with row
#'   names matching \code{baseline_tree$tip.label} (same order). For the default
#'   \code{formula = "trait_data ~ 1"}, \code{trait_data} is typically supplied as a numeric
#'   matrix, but a numeric response-only \code{data.frame} is also accepted. When using
#'   more general formulas (e.g., pGLS-style models), a \code{data.frame} with named
#'   columns can be used instead.
#' @param formula Character string or formula object passed to \code{mvgls}. Defaults to
#'   \code{"trait_data ~ 1"}, which fits an intercept-only model treating the supplied
#'   multivariate trait matrix as the response. This is the appropriate choice for most
#'   morphometric data where there are no predictor variables. For more general models,
#'   \code{formula} can reference subsets of \code{trait_data} explicitly, for example
#'   \code{"trait_data[, 1:5] ~ 1"} to treat columns 1-5 as a multivariate response,
#'   \code{"trait_data[, 1:5] ~ trait_data[, 6]"} to fit a multivariate pGLS with an
#'   indexed predictor, or \code{cbind(y1, y2) ~ size + grp} to fit a named-column
#'   pGLS with numeric or factor predictors.
#' @param min_descendant_tips Integer (\eqn{\ge}1). Minimum number of tips required for an internal node
#'   to be considered as a candidate shift (forwarded to \code{generatePaintedTrees}). Larger values
#'   reduce the number of candidate shifts by excluding very small clades. For empirical datasets,
#'   values around \code{10} are a reasonable starting choice and can be tuned in sensitivity analyses.
#' @param num_cores Integer. Number of workers for candidate scoring. Uses plain
#'   serial evaluation when \code{num_cores = 1}. For \code{num_cores > 1}, uses
#'   \code{future::plan(multicore)} on Unix outside \code{RStudio}; otherwise uses
#'   \code{future::plan(multisession)}. During the parallel candidate-scoring blocks,
#'   BLAS/OpenMP threads are capped to 1 (per worker) to avoid CPU oversubscription.
#' @param ic_uncertainty_threshold Numeric (\eqn{\ge}0). Reserved for future development
#'   in post-search pruning and uncertainty analysis; currently not used by
#'   \code{searchOptimalConfiguration()}.
#' @param shift_acceptance_threshold Numeric (\eqn{\ge}0). Minimum IC improvement
#'   (baseline - new) required to accept a candidate shift during the forward search.
#'   Larger values yield more conservative models. For analyses based on the Generalized
#'   Information Criterion (\code{"GIC"}), a threshold on the order of \code{20} units is a
#'   conservative choice that tends to admit only strongly supported shifts. Simulation
#'   studies in Berv et al. (2026) suggest that this choice yields good balanced
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
#'   \code{penalty}, \code{target}, \code{error}, \code{REML}, etc.). In the workflows
#'   emphasized in the package vignettes, \code{method = "H&L"} is used for
#'   intercept-only searches on high-dimensional response matrices, whereas
#'   \code{method = "LL"} is used for formula-based searches with predictors and
#'   should also be used when \code{IC = "BIC"}. In \pkg{mvMORPH}, \code{method = "H&L"}
#'   is restricted to intercept-only models and the \code{"RidgeArch"} penalty.
#'   Setting \code{error = TRUE} asks \code{mvgls()} to estimate a nuisance
#'   measurement-error (intraspecific-variance) term from the data.
#'
#' @details
#' \strong{Input requirements.}
#' \itemize{
#'   \item \emph{Tree:} \code{baseline_tree} should be a rooted \code{phylo} tree
#'         with branch lengths interpreted in units of time. An ultrametric tree is not required.
#'         The starting tree does not need to already be painted; \code{searchOptimalConfiguration()}
#'         paints a single baseline regime internally before building shifted candidates.
#'   \item \emph{Trait data alignment:} \code{rownames(trait_data)} must match
#'         \code{baseline_tree$tip.label} in both names and order; any tips without data should be
#'         pruned beforehand.
#'   \item \emph{Data type:} \code{trait_data} is typically a numeric matrix of continuous traits;
#'         numeric response-only \code{data.frame}s are also supported for intercept-only
#'         searches, and named mixed-type \code{data.frame}s are supported for richer
#'         formulas. High-dimensional settings (p \eqn{\ge} n) are supported via
#'         penalized-likelihood \code{mvgls()} fits.
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
#' \strong{Parallelization.} When \code{num_cores = 1}, candidates are scored serially. For
#' larger values, candidate sub-model fits are distributed with \pkg{future} +
#' \pkg{future.apply}. On Unix outside \code{RStudio}, \code{multicore} is used; otherwise
#' \code{multisession} is used. A sequential plan is restored afterward.
#'
#' \strong{Plotting.} If \code{plot = TRUE}, trees are rendered with
#' \code{\link[phytools]{plotSimmap}()}; shift IDs are labeled with \code{\link[ape]{nodelabels}()}.
#'
#' \strong{Regime VCVs.} The returned \code{$VCVs} are extracted from the fitted multi-regime model via
#' \code{extractRegimeVCVs} and reflect regime-specific covariance
#' estimates (when \code{mvgls} is fitted under a PL/ML method).
#'
#' For high-dimensional trait datasets (p \eqn{\ge} n), penalized-likelihood settings in
#' \code{mvgls()} are often required for stable estimation. The package vignettes
#' distinguish two common workflows. For intercept-only searches on high-dimensional
#' response matrices (for example, GPA-aligned landmark data), the jaw-shape vignette
#' uses \code{method = "H&L"} with the default \code{"RidgeArch"} penalty; in
#' \pkg{mvMORPH}, this is a fast approximation to penalized LOOCV and is only available
#' for intercept-only models. For formula-based searches with predictors, the avian
#' skeleton vignette uses \code{method = "LL"} instead. When \code{IC = "BIC"},
#' \code{method = "LL"} should be used. Across empirical workflows, \code{error = TRUE}
#' is often a sensible default because it asks \code{mvgls()} to estimate a nuisance
#' measurement-error (intraspecific-variance) term from the data. Users should consult
#' the \pkg{mvMORPH} documentation for details on available methods and penalties and
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
#'
#' \dontrun{
#' # Intercept-only empirical-style search:
#' # high-dimensional response matrix with H&L + measurement error
#' res_hl <- searchOptimalConfiguration(
#'   baseline_tree              = as.phylo(simmap),
#'   trait_data                 = X,
#'   formula                    = "trait_data ~ 1",
#'   min_descendant_tips        = 10,
#'   num_cores                  = 2,
#'   shift_acceptance_threshold = 20,
#'   uncertaintyweights_par     = TRUE,
#'   IC                         = "GIC",
#'   plot                       = FALSE,
#'   method                     = "H&L",
#'   error                      = TRUE,
#'   store_model_fit_history    = TRUE,
#'   verbose                    = TRUE
#' )
#'
#' # Formula-based search with a predictor:
#' # use LL when the model includes predictors
#' dat <- data.frame(
#'   trait1    = X[, 1],
#'   trait2    = X[, 2],
#'   predictor = rnorm(nrow(X))
#' )
#' rownames(dat) <- simmap$tip.label
#'
#' res_ll <- searchOptimalConfiguration(
#'   baseline_tree              = as.phylo(simmap),
#'   trait_data                 = dat,
#'   formula                    = "trait_data[, 1:2] ~ trait_data[, 3]",
#'   min_descendant_tips        = 10,
#'   num_cores                  = 2,
#'   shift_acceptance_threshold = 20,
#'   IC                         = "GIC",
#'   plot                       = FALSE,
#'   method                     = "LL",
#'   error                      = TRUE,
#'   store_model_fit_history    = TRUE,
#'   verbose                    = TRUE
#' )
#' }
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
           IC = "GIC",
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

    if (!(inherits(formula, "formula") ||
          (is.character(formula) && length(formula) == 1L && !is.na(formula)))) {
      stop("formula must be a single character string or formula object.")
    }

    # Coerce to phylo and initialize a single baseline regime at the root
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
    .bifrost_search_validate_ic(IC)
    baseline_model <- .bifrost_search_fit_ic(IC, formula, candidate_trees[[1]], trait_data, ...)
    baseline_ic <- .bifrost_search_ic_value(baseline_model, IC)
    .progress("Baseline %s: %.2f", IC, baseline_ic)

    #evaluate all of the candidate trees under GIC or BIC
    # Capture additional arguments into a list
    args_list <- list(...)

    is_rstudio <- identical(Sys.getenv("RSTUDIO"), "1")

    .progress("%s", "Fitting sub-models in parallel...")

    candidate_scores <- .bifrost_search_score_candidates(
      candidate_trees_shifts = candidate_trees_shifts,
      baseline_ic = baseline_ic,
      IC = IC,
      formula = formula,
      trait_data = trait_data,
      args_list = args_list,
      num_cores = num_cores,
      is_rstudio = is_rstudio
    )

    .progress("%s", "Sorting and evaluating shifts...")
    sorted_candidates <- candidate_scores$sorted_candidates
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

    # Where to store on-disk history (CRAN-safe temp location)
    sub_dir <- .bifrost_search_history_dir(store_model_fit_history)

    forward_search <- .bifrost_search_forward(
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
      progress = .progress,
      ...
    )

    current_best_tree <- forward_search$current_best_tree
    current_best_ic <- forward_search$current_best_ic
    shift_vec <- forward_search$shift_vec
    model_with_shift_no_uncertainty <- forward_search$model_with_shift_no_uncertainty
    best_tree_no_uncertainty <- forward_search$best_tree_no_uncertainty
    warnings_list <- forward_search$warnings_list
    model_fit_history <- forward_search$model_fit_history

    #print(paste(shift_vec))
    shifts_no_uncertainty <- forward_search$shifts_no_uncertainty
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
    ic_weights_df <- .bifrost_search_calculate_ic_weights(
      uncertaintyweights = uncertaintyweights,
      uncertaintyweights_par = uncertaintyweights_par,
      shift_vec = shift_vec,
      best_tree_no_uncertainty = best_tree_no_uncertainty,
      model_with_shift_no_uncertainty = model_with_shift_no_uncertainty,
      IC = IC,
      formula = formula,
      trait_data = trait_data,
      args_list = args_list,
      num_cores = num_cores,
      is_rstudio = is_rstudio,
      progress = .progress
    )

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

      no_uncertainty <- .bifrost_search_no_uncertainty_components(
        shifts_no_uncertainty = shifts_no_uncertainty,
        model_with_shift_no_uncertainty = model_with_shift_no_uncertainty,
        best_tree_no_uncertainty = best_tree_no_uncertainty,
        baseline_model = baseline_model,
        baseline_candidate_tree = candidate_trees[[1]]
      )

      # Create the main list that will always be returned
      result_list <- list(
        user_input = user_input,
        tree_no_uncertainty_transformed = no_uncertainty$tree_transformed,
        tree_no_uncertainty_untransformed = no_uncertainty$tree_untransformed,
        model_no_uncertainty = no_uncertainty$model,
        shift_nodes_no_uncertainty = shifts_no_uncertainty,
        optimal_ic = current_best_ic,
        baseline_ic = baseline_ic,
        IC_used = IC,
        num_candidates = length(sorted_candidates),
        model_fit_history = model_fit_history
      )

      # Create the IC and acceptance matrix from the model fit history
      if (isTRUE(store_model_fit_history) && !is.null(sub_dir)) {
        result_list$model_fit_history <- .bifrost_search_load_history(sub_dir, IC)
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

    class(result_list) <- c("bifrost_search", class(result_list))
    return(result_list)

  }
#... is optional arguments to be passed to mvgls
