#' Generate Painted Subtrees for Eligible Nodes in a Phylogenetic Tree
#'
#' Identifies internal nodes in a (rooted) phylogenetic tree that have at least
#' a specified minimum number of descendant tips and returns a list of trees in
#' which each eligible node's downstream clade has been painted as a discrete
#' regime. Each painted tree represents a potential sub-model for further analysis.
#'
#' If `tree` is unrooted, it is first rooted using the first tip as an
#' outgroup (`ape::root(..., resolve.root = TRUE)`).
#'
#' @param tree An object of class \code{phylo}. If unrooted, it is rooted internally.
#' @param min_tips Integer (\eqn{\ge}1). Minimum number of descendant tips required for an
#'   internal node to be considered eligible.
#' @param state Character scalar. The regime label to paint on each eligible subtree.
#'   Defaults to \code{"shift"}.
#'
#' @return A named \code{list} of trees of class \code{simmap} (one per eligible node).
#'   Names are of the form \code{"Node <n>"} where \code{<n>} is the internal node
#'   number in the original \code{tree}. The function also prints counts of eligible
#'   nodes and generated sub-models via \code{cat()}.
#'
#' @details
#' Internally, the function:
#' \itemize{
#'   \item Ensures the tree is rooted.
#'   \item Scans internal nodes and keeps those with at least \code{min_tips} descendant tips
#'         (via \code{phytools::getDescendants} and \code{ape::Ntip} / \code{ape::Nnode}).
#'   \item Uses \code{phytools::paintSubTree} to paint the clade downstream of each eligible node
#'         with the specified \code{state}.
#' }
#'
#' @seealso \code{\link[phytools]{paintSubTree}}, \code{\link[phytools]{getDescendants}},
#'   \code{\link[ape]{Ntip}}, \code{\link[ape]{Nnode}}, \code{\link[ape]{root}}
#'
#' @examples
#' \dontrun{
#'   set.seed(123)
#'   tree <- phytools::pbtree(n = 30, scale = 1)
#'   painted <- generatePaintedTrees(tree, min_tips = 5, state = "shift")
#'   length(painted)
#'   names(painted)[1]
#' }
#'
#' @importFrom ape is.rooted root Ntip Nnode
#' @importFrom phytools getDescendants paintSubTree
#' @keywords internal
#' @noRd
generatePaintedTrees <- function(tree, min_tips, state = "shift") {
  if (!is.rooted(tree)) {
    tree <- root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
  }

  getEligibleNodes <- function(tree, min_tips) {
    eligible_nodes <- c()
    for (node in 1:(Nnode(tree))) {
      internal_node <- node + Ntip(tree)
      descendants <- getDescendants(tree, internal_node)
      tip_descendants <- tree$tip.label[descendants[descendants <= Ntip(tree)]]
      if (length(tip_descendants) >= min_tips) {
        eligible_nodes <- c(eligible_nodes, internal_node)
      }
    }
    return(eligible_nodes)
  }

  eligible_nodes <- getEligibleNodes(tree, min_tips)
  if (isTRUE(getOption("bifrost.verbose"))) {
    message(sprintf("%d eligible nodes are detected", length(eligible_nodes)))
  }

  painted_trees <- list()

  for (node in eligible_nodes) {
    tree_copy <- tree
    tree_copy <- paintSubTree(tree_copy, node, state)
    painted_trees[[paste("Node", node)]] <- tree_copy
  }

  if (isTRUE(getOption("bifrost.verbose"))) {
    message(sprintf("%d sub-models generated", length(painted_trees)))
  }
  return(painted_trees)
}

#' Fit mvgls Model to a Painted Tree and Extract GIC Score
#'
#' Fits a multivariate generalized least squares (mvgls) model using a
#' SIMMAP-formatted (painted) phylogenetic tree and a matrix of trait data,
#' returning both the fitted model object and its Generalized Information
#' Criterion (GIC) score.
#'
#' @param painted_tree An object of class \code{simmap} (see
#'   \code{\link[phytools]{paintSubTree}}) used as the phylogenetic tree.
#'   Must be fully dichotomous.
#' @param trait_data A numeric matrix of trait values (species × traits).
#'   Row names must match the tip labels of \code{painted_tree} exactly
#'   and be in the same order.
#'
#' @return A \code{list} with two components:
#' \describe{
#'   \item{\code{model}}{An object of class \code{mvgls}.}
#'   \item{\code{GIC}}{A single numeric value, the Generalized Information Criterion.}
#' }
#'
#' @details
#' This function always fits a multi-rates Brownian motion model
#' (\code{model = "BMM"}) under a standard intercept-only design
#' (\code{~ 1}). For trees with a single regime, the model is equivalent to
#' a standard Brownian motion fit.
#'
#' @seealso \code{\link[mvMORPH]{mvgls}}, \code{\link[mvMORPH]{GIC}},
#'   \code{\link[phytools]{paintSubTree}}
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   tree <- phytools::pbtree(n = 10, scale = 1)
#'   tree_painted <- phytools::paintSubTree(tree, node = 11, state = "A", anc.state = "A")
#'   trait_mat <- matrix(rnorm(10 * 2), ncol = 2)
#'   rownames(trait_mat) <- tree_painted$tip.label
#'   result <- fitMvglsAndExtractGIC(tree_painted, trait_mat)
#'   result$GIC
#' }
#'
#' @importFrom mvMORPH mvgls GIC
#' @keywords internal
#' @noRd
fitMvglsAndExtractGIC <- function(painted_tree, trait_data) {
  # Ensure trait_data is a matrix
  if (!is.matrix(trait_data)) {
    stop("trait_data must be a matrix.")
  }

  # Make sure the row names match the tip labels of the painted_tree
  if (!identical(rownames(trait_data), painted_tree$tip.label)) {
    stop("Row names of trait_data must exactly match the tip labels of the tree.")
  }

  # Fit the mvgls model directly using the matrix
  model <- mvgls(trait_data ~ 1, tree = painted_tree, model = "BMM", method='LL')
  gic_value <- GIC(model)

  # Return a list containing the model and the GIC
  return(list(model = model, GIC = gic_value))
}

#' Fit mvgls Model to a Painted Tree and Extract BIC Score
#'
#' Fits a multivariate generalized least squares (mvgls) model using a
#' SIMMAP-formatted (painted) phylogenetic tree and a matrix of trait data,
#' returning both the fitted model object and its Bayesian Information Criterion
#' (BIC) score.
#'
#' @param painted_tree An object of class \code{simmap} (see
#'   \code{\link[phytools]{paintSubTree}}) used as the phylogenetic tree.
#'   Must be fully dichotomous.
#' @param trait_data A numeric matrix of trait values (species × traits).
#'   Row names must exactly match the tip labels of \code{painted_tree} and
#'   be in the same order.
#'
#' @return A \code{list} with two components:
#' \describe{
#'   \item{\code{model}}{An object of class \code{mvgls}.}
#'   \item{\code{BIC}}{A single numeric value, the Bayesian Information Criterion.}
#' }
#'
#' @details
#' This function always fits a multi-rates Brownian motion model
#' (\code{model = "BMM"}) under a simple intercept-only design
#' (\code{~ 1}). For trees with a single painted regime, the model is
#' equivalent to a standard Brownian motion fit.
#'
#' @seealso \code{\link[mvMORPH]{mvgls}}, \code{\link[stats]{BIC}},
#'   \code{\link[phytools]{paintSubTree}}
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   tree <- phytools::pbtree(n = 10, scale = 1)
#'   tree_painted <- phytools::paintSubTree(tree, node = 11, state = "A", anc.state = "A")
#'   trait_mat <- matrix(rnorm(10 * 3), ncol = 3)
#'   rownames(trait_mat) <- tree_painted$tip.label
#'   result <- fitMvglsAndExtractBIC(tree_painted, trait_mat)
#'   result$BIC
#' }
#'
#' @importFrom mvMORPH mvgls
#' @importFrom stats BIC
#' @keywords internal
#' @noRd
fitMvglsAndExtractBIC <- function(painted_tree, trait_data) {
  # Ensure trait_data is a matrix
  if (!is.matrix(trait_data)) {
    stop("trait_data must be a matrix.")
  }

  # Make sure the row names match the tip labels of the painted_tree
  if (!identical(rownames(trait_data), painted_tree$tip.label)) {
    stop("Row names of trait_data must exactly match the tip labels of the tree.")
  }

  # Fit the mvgls model directly using the matrix
  model <- mvgls(trait_data ~ 1, tree = painted_tree, model = "BMM", method='LL')
  bic_value <- BIC(model)

  # Return a list containing the model and the GIC
  return(list(model = model, BIC = bic_value))
}

#' Fit mvgls Model Using a Formula and Extract GIC Score
#'
#' Fits a multivariate generalized least squares (mvgls) model to a
#' SIMMAP-formatted (painted) phylogenetic tree using a user-specified formula
#' and trait data. Automatically selects the evolutionary model:
#' \itemize{
#'   \item \code{model = "BM"} when the painted tree has a single regime.
#'   \item \code{model = "BMM"} when multiple regimes are present.
#' }
#' Returns the fitted model along with its Generalized Information Criterion (GIC) score.
#'
#' @param formula A character string specifying the model formula
#'   (e.g., \code{"cbind(trait1, trait2) ~ predictor"}).
#' @param painted_tree An object of class \code{simmap} (see
#'   \code{\link[phytools]{paintSubTree}}) used as the phylogenetic tree.
#' @param trait_data A \code{data.frame} or \code{matrix} of trait values with
#'   row names that exactly match the tip labels of \code{painted_tree}.
#' @param ... Additional arguments passed to \code{\link[mvMORPH]{mvgls}}.
#'
#' @return A \code{list} with two components:
#' \describe{
#'   \item{\code{model}}{An object of class \code{mvgls}.}
#'   \item{\code{GIC}}{A single numeric value with the Generalized Information Criterion.}
#' }
#'
#' @details
#' The function converts the character \code{formula} to an actual formula
#' object via \code{\link[stats]{as.formula}} before fitting. The number of regimes
#' in \code{painted_tree} is detected with \code{\link[phytools]{getStates}}.
#'
#' @seealso \code{\link[mvMORPH]{mvgls}}, \code{\link[mvMORPH]{GIC}},
#'   \code{\link[phytools]{getStates}}, \code{\link[stats]{as.formula}}
#'
#' @examples
#' \dontrun{
#'   set.seed(123)
#'   tree <- phytools::pbtree(n = 15, scale = 1)
#'   painted <- phytools::paintSubTree(tree, node = 16, state = "A", anc.state = "A")
#'   x <- rnorm(15)
#'   y1 <- 0.5 * x + rnorm(15)
#'   y2 <- -0.3 * x + rnorm(15)
#'   dat <- data.frame(x = x, y1 = y1, y2 = y2)
#'   rownames(dat) <- painted$tip.label
#'   result <- fitMvglsAndExtractGIC.formula("cbind(y1, y2) ~ x", painted, dat)
#'   result$GIC
#' }
#'
#' @importFrom mvMORPH mvgls GIC
#' @importFrom phytools getStates
#' @importFrom stats as.formula
#' @keywords internal
#' @noRd
fitMvglsAndExtractGIC.formula <- function(formula, painted_tree, trait_data, ...) {
  # Ensure trait_data is a matrix
  #if (!is.matrix(trait_data)) {
  #  stop("trait_data must be a matrix.")
  #}

  # Make sure the row names match the tip labels of the painted_tree
  if (!identical(rownames(trait_data), painted_tree$tip.label)) {
    stop("Row names of trait_data must exactly match the tip labels of the tree.")
  }

  # Validate that formula is provided and is a character
  if (missing(formula) || !is.character(formula)) {
    stop("A character formula must be provided.")
  }

  # Convert the string formula to an actual formula object
  formula_obj <- as.formula(formula)

  # Fit the mvgls model using the user-defined formula

  # On the baseline tree, we have to switch to model = BM bc there is only one painted regime
  # Then we switch back to BMM

  if(length(unique(getStates(tree=painted_tree))) == 1){
    model <- mvgls(formula_obj, tree = painted_tree, model = "BM", ...)
  } else {
    model <- mvgls(formula_obj, tree = painted_tree, model = "BMM", ...)
  }
  gic_value <- GIC(model)

  # Return a list containing the model and the GIC
  return(list(model = model, GIC = gic_value))
}

#' Fit mvgls Model Using a Formula and Extract BIC Score
#'
#' Fits a multivariate generalized least squares (mvgls) model to a
#' SIMMAP-formatted (painted) phylogenetic tree using a user-specified formula
#' and trait data. Automatically selects the evolutionary model:
#' \itemize{
#'   \item \code{model = "BM"} when the painted tree has a single regime.
#'   \item \code{model = "BMM"} when multiple regimes are present.
#' }
#' Returns the fitted model along with its Bayesian Information Criterion (BIC) score.
#'
#' @param formula A character string specifying the model formula
#'   (e.g., \code{"cbind(trait1, trait2) ~ predictor"}).
#' @param painted_tree An object of class \code{simmap} (see
#'   \code{\link[phytools]{paintSubTree}}) used as the phylogenetic tree.
#' @param trait_data A \code{data.frame} or \code{matrix} of trait values with
#'   row names that exactly match the tip labels of \code{painted_tree}.
#' @param ... Additional arguments passed to \code{\link[mvMORPH]{mvgls}}.
#'
#' @return A \code{list} with two components:
#' \describe{
#'   \item{\code{model}}{An object of class \code{mvgls}.}
#'   \item{\code{BIC}}{A single numeric value containing the Bayesian Information Criterion.}
#' }
#'
#' @details
#' The function converts the character \code{formula} to an actual formula
#' object via \code{\link[stats]{as.formula}} before fitting. The number of regimes
#' in \code{painted_tree} is detected with \code{\link[phytools]{getStates}}.
#'
#' @seealso \code{\link[mvMORPH]{mvgls}}, \code{\link[stats]{BIC}},
#'   \code{\link[phytools]{getStates}}, \code{\link[stats]{as.formula}}
#'
#' @examples
#' \dontrun{
#'   set.seed(321)
#'   tree <- phytools::pbtree(n = 12, scale = 1)
#'   painted <- phytools::paintSubTree(tree, node = 13, state = "A", anc.state = "A")
#'   x <- rnorm(12)
#'   y1 <- x + rnorm(12, sd = 0.5)
#'   y2 <- 0.8 * x + rnorm(12, sd = 0.3)
#'   dat <- data.frame(x = x, y1 = y1, y2 = y2)
#'   rownames(dat) <- painted$tip.label
#'   result <- fitMvglsAndExtractBIC.formula("cbind(y1, y2) ~ x", painted, dat)
#'   result$BIC
#' }
#'
#' @importFrom mvMORPH mvgls
#' @importFrom phytools getStates
#' @importFrom stats as.formula BIC
#' @keywords internal
#' @noRd
fitMvglsAndExtractBIC.formula <- function(formula, painted_tree, trait_data, ...) {
  # # Ensure trait_data is a matrix
  # if (!is.matrix(trait_data)) {
  #   stop("trait_data must be a matrix.")
  # }

  # Make sure the row names match the tip labels of the painted_tree
  if (!identical(rownames(trait_data), painted_tree$tip.label)) {
    stop("Row names of trait_data must exactly match the tip labels of the tree.")
  }

  # Validate that formula is provided and is a character
  if (missing(formula) || !is.character(formula)) {
    stop("A character formula must be provided.")
  }

  # Convert the string formula to an actual formula object
  formula_obj <- as.formula(formula)

  # Fit the mvgls model using the user-defined formula

  # On the baseline tree, we have to switch to model = BM bc there is only one painted regime
  # Then we switch back to BMM

  if(length(unique(getStates(tree=painted_tree))) == 1){
    model <- mvgls(formula_obj, tree = painted_tree, model = "BM", ...)
  } else {
    model <- mvgls(formula_obj, tree = painted_tree, model = "BMM", ...)
  }
  bic_value <- BIC(model)

  # Return a list containing the model and the GIC
  return(list(model = model, BIC = bic_value))
}

if (FALSE) {
#' Calculate \eqn{\Delta}GIC Scores Relative to a Baseline Model
#'
#' Computes the difference in Generalized Information Criterion (GIC) scores
#' between a baseline model (assumed to be the first element of `model_results`)
#' and one or more alternative models fitted to painted phylogenetic trees.
#' This allows rapid comparison of model fit across candidate regime-paintings.
#'
#' @param model_results A list of model fit results, typically the output of
#'   \code{fitMvglsAndExtractGIC.formula} or a similar internal helper.
#'   Each element must be a list containing at least a \code{GIC} component,
#'   where \code{GIC$GIC} is a numeric scalar.
#'   The first element is treated as the baseline model.
#'
#' @param painted_tree_list A named list of painted phylogenetic trees
#'   (objects of class \code{simmap}). The names are used as identifiers
#'   for the output vector of \eqn{\Delta}GIC values. Every element must have a
#'   non-\code{NULL} name.
#'
#' @return A named numeric vector of \eqn{\Delta}GIC scores, where each value is:
#'   \deqn{\Delta GIC_i = GIC_\mathrm{baseline} - GIC_i}
#'   Positive values indicate that the candidate model improves the fit
#'   relative to the baseline (lower GIC = better model).
#'
#' @details
#' This function is a convenience wrapper for comparing multiple models
#' on the same dataset. By subtracting each candidate model's GIC from the
#' baseline GIC, the resulting values are on a common, interpretable scale.
#'
#' @seealso
#' \code{\link[mvMORPH]{GIC}}, \code{\link[phytools]{paintSubTree}}
#'
#' @note
#' Related internal helper: \code{fitMvglsAndExtractGIC.formula} (not exported).
#'
#' @examples
#' \dontrun{
#'   # Assume we have a baseline and two alternative models
#'   baseline <- list(GIC = list(GIC = 100))
#'   alt1     <- list(GIC = list(GIC = 95))
#'   alt2     <- list(GIC = list(GIC = 105))
#'
#'   model_results <- list(baseline, alt1, alt2)
#'   painted_tree_list <- list(
#'     Baseline = "tree1",
#'     ShiftA   = "tree2",
#'     ShiftB   = "tree3"
#'   )
#'
#'   calculateAllDeltaGIC(model_results, painted_tree_list)
#'   # Named vector: Baseline=0, ShiftA=5, ShiftB=-5
#' }
#'
#' @keywords internal
#' @noRd
calculateAllDeltaGIC <- function(model_results, painted_tree_list) {
  # Check if the painted trees have names
  if (!all(sapply(painted_tree_list, function(x) !is.null(names(x))))) {
    stop("All trees in 'painted_tree_list' must have names.")
  }

  # Extract the names from the painted_tree_list
  tree_list_names <- names(painted_tree_list)

  # Retrieve the baseline GIC from the first model result
  # Ensure we are accessing the numeric GIC value correctly
  baseline_gic <- model_results[[1]]$GIC$GIC
  if (!is.numeric(baseline_gic)) {
    stop("The baseline GIC value must be numeric.")
  }

  # Use lapply to iterate over the list of model_results and calculate the delta GIC
  delta_gic <- lapply(seq_along(model_results), function(i) {
    # Ensure we are accessing the numeric GIC value correctly
    current_gic <- model_results[[i]]$GIC$GIC
    if (!is.numeric(current_gic)) {
      stop(paste("The GIC value for model", i, "must be numeric."))
    }
    # Calculate the difference in GIC between the baseline and the shift model
    return(baseline_gic - current_gic)
  })

  # Assign names to the delta GIC values and convert it to a named vector
  delta_gic <- setNames(unlist(delta_gic), tree_list_names)

  return(delta_gic)
}
}

#' Paint a Subtree in a Phylogenetic Tree with Optional Selective Overwriting
#'
#' Modifies a phylogenetic tree by painting the clade descending from a
#' specified node with a new discrete state and returns a SIMMAP-style tree.
#' You can either overwrite the entire subtree (\code{overwrite = TRUE}) or
#' selectively overwrite only edges currently in a \emph{target} state
#' (\code{overwrite = FALSE}). Optional stem painting splits the parent edge
#' into ancestral and derived states.
#'
#' @param tree An object of class \code{phylo}. If no mapped states are present,
#'   they are initialized (all edges set to \code{anc.state}).
#' @param node Integer node number (in \code{tree}) at which painting begins
#'   (i.e., the clade rooted at \code{node} will be painted).
#' @param state Character (or numeric) label of the new state to paint.
#' @param anc.state Character (or numeric) label of the ancestral (baseline)
#'   state used when initializing mappings or when splitting the stem.
#'   Default is \code{"1"}.
#' @param stem Logical or numeric. If \code{FALSE} (default), the incoming edge
#'   to \code{node} is left unchanged. If \code{TRUE}, the entire incoming edge
#'   is assigned to \code{state}. If a numeric value in \eqn{[0,1]}, the parent
#'   edge is split such that a fraction \code{stem} is assigned to \code{state}
#'   and \code{1 - stem} remains \code{anc.state}. For tip nodes, \code{stem}
#'   must be \code{TRUE} (cannot be \code{FALSE}).
#' @param overwrite Logical. If \code{TRUE} (default), overwrite the mappings
#'   on \emph{all} edges in the subtree with \code{state}. If \code{FALSE},
#'   only edges whose current mapping equals the \emph{target} state are
#'   overwritten (the target is the current state on the edge leading to
#'   \code{node}, or \code{anc.state} for tips).
#'
#' @return A modified tree of class \code{c("simmap", "phylo")} with updated
#'   \code{$maps} and \code{$mapped.edge}.
#'
#' @details
#' \itemize{
#'   \item If \code{tree$maps} is \code{NULL}, the function initializes a SIMMAP
#'         representation by assigning each edge length to \code{anc.state}.
#'   \item With \code{overwrite = TRUE}, every edge in the subtree is collapsed
#'         to a single segment labeled \code{state}.
#'   \item With \code{overwrite = FALSE}, only edges matching the target state
#'         (the state on the edge entering \code{node}, or \code{anc.state}
#'         for tips) are replaced; other existing states are preserved.
#'   \item If \code{stem} is numeric in \eqn{[0,1]}, the parent edge is split
#'         into two segments: \code{(1 - stem)} of \code{anc.state} and
#'         \code{stem} of \code{state}.
#' }
#'
#' @seealso \code{\link[phytools]{paintSubTree}}, \code{\link[phytools]{paintBranches}}
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   tr <- phytools::pbtree(n = 10, scale = 1)
#'   # Initialize mappings to "0" globally
#'   tr0 <- phytools::paintBranches(tr, edge = unique(tr$edge[,2]), state = "0", anc.state = "0")
#'
#'   # Paint an internal clade fully as "1"
#'   nd <- ape::Ntip(tr0) + 2L
#'   tr1 <- paintSubTree_mod(tr0, node = nd, state = "1", anc.state = "0",
#'     stem = TRUE, overwrite = TRUE)
#'
#'   # Selectively overwrite only edges already in state "0" within the subtree
#'   tr2 <- paintSubTree_mod(tr1, node = nd, state = "2", anc.state = "0",
#'     stem = 0.25, overwrite = FALSE)
#' }
#' @importFrom ape compute.brlen
#' @importFrom phytools getDescendants
#' @keywords internal
#' @noRd
paintSubTree_mod <- function(tree, node, state, anc.state="1", stem=FALSE, overwrite=TRUE) {
  if (!inherits(tree, "phylo")) stop("tree should be an object of class \"phylo\".")
  if (stem == 0 && node <= length(tree$tip)) stop("stem must be TRUE for node <= N")
  if (is.null(tree$edge.length)) tree <- compute.brlen(tree)

  if (is.null(tree$maps)) {
    maps <- as.list(tree$edge.length)
    for (i in 1:length(maps)) names(maps[[i]]) <- anc.state
  } else {
    maps <- tree$maps
  }

  if (overwrite) {
    # Original behavior: Overwrite entire subtree
    desc <- getDescendants(tree, node)
    z <- which(tree$edge[,2] %in% desc)
    for (i in z) {
      maps[[i]] <- sum(maps[[i]])
      names(maps[[i]]) <- state
    }
  } else {
    # Modified behavior: Selective overwriting
    target_state <- if (node > length(tree$tip)) names(maps[[which(tree$edge[,2] == node)]]) else anc.state
    desc <- getDescendants(tree, node)
    z <- which(tree$edge[,2] %in% desc)
    for (i in z) {
      if (names(maps[[i]]) == target_state) {
        maps[[i]] <- sum(maps[[i]])
        names(maps[[i]]) <- state
      }
    }
  }

  if (stem && node > length(tree$tip)) {
    stem_edge <- which(tree$edge[,2] == node)
    maps[[stem_edge]] <- sum(maps[[stem_edge]]) * c(1 - stem, stem)
    names(maps[[stem_edge]]) <- c(anc.state, state)
  }

  s <- vector()
  for (i in 1:nrow(tree$edge)) s <- c(s, names(maps[[i]]))
  s <- unique(s)
  mapped.edge <- matrix(0, length(tree$edge.length), length(s), dimnames=list(edge=apply(tree$edge, 1, function(x) paste(x, collapse=",")), state=s))
  for (i in 1:length(maps)) {
    for (j in 1:length(maps[[i]])) {
      mapped.edge[i, names(maps[[i]])[j]] <- mapped.edge[i, names(maps[[i]])[j]] + maps[[i]][j]
    }
  }

  tree$mapped.edge <- mapped.edge
  tree$maps <- maps
  class(tree) <- c("simmap", setdiff(class(tree), "simmap"))

  return(tree)
}

#' Remove a Painted Shift from a SIMMAP Tree
#'
#' Selectively removes a previously painted shift (regime) from a clade
#' descending from \code{shift_node} in a SIMMAP-style phylogenetic tree.
#' Edges in the clade that match the shift node's state are reassigned to
#' the ancestral state inherited from the parent node. Optionally, the
#' stem (incoming edge) can also be restored.
#'
#' @param tree A phylogenetic tree of class \code{"phylo"} or
#'   \code{c("simmap","phylo")}. The tree must have initialized
#'   \code{$edge.length} and \code{$maps}. If \code{$edge.length} is missing,
#'   it is computed using \code{\link[ape]{compute.brlen}}. If \code{$maps} is
#'   missing, all edges are assumed to be in state \code{"1"}.
#' @param shift_node Integer node number indicating where the shift
#'   begins (the clade rooted at this node will be repainted).
#' @param stem Logical or numeric. If \code{FALSE} (default), the
#'   incoming edge to \code{shift_node} is left unchanged. If \code{TRUE},
#'   the entire incoming edge is reassigned to the parental state.
#'   If a numeric value in \eqn{[0,1]}, the parent edge is split into
#'   two segments, both labeled with the parental state (effectively
#'   restoring it).
#'
#' @details
#' \itemize{
#'   \item The parental state is obtained from
#'   \code{\link[phytools]{getStates}(tree, type = "nodes")} for the parent of
#'   \code{shift_node}. If unavailable, defaults to \code{"1"}.
#'   \item Only edges whose single-segment state matches the
#'   state at \code{shift_node} are overwritten. Other states are preserved.
#'   \item If edges have multiple mapped segments, only exact single-state
#'   matches are considered for replacement.
#' }
#'
#' @return A modified tree of class \code{c("simmap","phylo")} with updated
#'   \code{$maps} and \code{$mapped.edge}, in which all affected edges are
#'   reassigned to the parent's state.
#'
#' @seealso
#' \code{paintSubTree_mod} (internal, not exported),
#' \code{\link[phytools]{paintSubTree}},
#' \code{\link[phytools]{getParent}},
#' \code{\link[phytools]{getDescendants}},
#' \code{\link[phytools]{getStates}},
#' \code{\link[ape]{compute.brlen}}
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   tr <- phytools::pbtree(n = 10, scale = 1)
#'   # Initialize to global state "0"
#'   tr0 <- phytools::paintBranches(tr, edge = unique(tr$edge[, 2]),
#'                                  state = "0", anc.state = "0")
#'   # Paint a shift (state "1") on an internal clade
#'   nd <- ape::Ntip(tr0) + 2L
#'   tr1 <- paintSubTree_mod(tr0, node = nd, state = "1",
#'                           anc.state = "0", stem = TRUE, overwrite = TRUE)
#'
#'   # Remove the shift: descendants revert to parental state "0"
#'   tr2 <- paintSubTree_removeShift(tr1, shift_node = nd, stem = FALSE)
#' }
#' @importFrom ape compute.brlen
#' @importFrom phytools getParent getDescendants getStates
#' @keywords internal
#' @noRd
paintSubTree_removeShift <- function(tree, shift_node, stem=FALSE) {
  if (!inherits(tree, "phylo")) stop("tree should be an object of class 'phylo'.")
  if (is.null(tree$edge.length)) tree <- compute.brlen(tree)

  if (is.null(tree$maps)) {
    maps <- as.list(tree$edge.length)
    for (i in 1:length(maps)) names(maps[[i]]) <- "1"  # Assuming '1' is the default ancestral state
  } else {
    maps <- tree$maps
  }

  # Get parent node and its state
  parent_node <- phytools::getParent(tree, shift_node)
  parent_state <- if (!is.na(parent_node)) getStates(tree, type = "nodes")[as.character(parent_node)] else "1"  # Default to '1' if parent is NA

  # Handle stem painting if applicable
  if (stem && shift_node > length(tree$tip)) {
    stem_edge <- which(tree$edge[,2] == shift_node)
    maps[[stem_edge]] <- sum(maps[[stem_edge]]) * c(1 - stem, stem)
    names(maps[[stem_edge]]) <- c(parent_state, parent_state)
  }

  # Get descendants and selectively overwrite branches
  desc <- getDescendants(tree, shift_node)
  shift_node_state <- getStates(tree, type = "nodes")[as.character(shift_node)]

  z <- which(tree$edge[,2] %in% desc)
  for (i in z) {
    if (names(maps[[i]]) == shift_node_state) {
      maps[[i]] <- sum(maps[[i]])
      names(maps[[i]]) <- parent_state
    }
  }

  # Update the tree with the new maps
  s <- vector()
  for (i in 1:nrow(tree$edge)) s <- c(s, names(maps[[i]]))
  s <- unique(s)
  mapped.edge <- matrix(0, length(tree$edge.length), length(s), dimnames = list(edge = apply(tree$edge, 1, function(x) paste(x, collapse = ",")), state = s))
  for (i in 1:length(maps)) {
    for (j in 1:length(maps[[i]])) {
      mapped.edge[i, names(maps[[i]])[j]] <- mapped.edge[i, names(maps[[i]])[j]] + maps[[i]][j]
    }
  }

  tree$mapped.edge <- mapped.edge
  tree$maps <- maps
  class(tree) <- c("simmap", setdiff(class(tree), "simmap"))

  return(tree)
}

#' Add a New Shift (Regime) to a SIMMAP Tree
#'
#' Increments a shift (regime) identifier and paints the subtree starting at
#' \code{shift_node} on a SIMMAP-style phylogenetic tree with the new state.
#' This enables stepwise construction of multi-regime models by applying
#' successive shifts with unique state IDs.
#'
#' @param tree A phylogenetic tree of class \code{phylo} or \code{simmap}.
#'   If not yet SIMMAP-initialized, ensure edge mappings exist upstream (e.g.,
#'   by painting a global baseline state) before calling.
#' @param shift_node Integer node number in \code{tree} indicating where the
#'   new regime should begin (the clade rooted at this node will be painted).
#' @param current_shift_id Integer ID of the last regime used. The function
#'   will assign \code{current_shift_id + 1} to the new shift.
#'
#' @return A \code{list} with:
#' \describe{
#'   \item{\code{tree}}{The updated SIMMAP tree with the new regime painted.}
#'   \item{\code{shift_id}}{The incremented shift ID (integer).}
#' }
#'
#' @details
#' The subtree is painted using \code{paintSubTree_mod()} with
#' \code{overwrite = FALSE} and \code{stem = FALSE}, meaning only edges
#' matching the target state (typically the current state along that path)
#' are selectively overwritten, and the parent (incoming) edge to
#' \code{shift_node} is not split or repainted. If you need to force a full
#' overwrite or paint the stem, call \code{paintSubTree_mod()} directly.
#'
#' @seealso
#' \code{paintSubTree_mod} (internal, not exported),
#' \code{\link[phytools]{paintSubTree}},
#' \code{\link[phytools]{paintBranches}}
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   tr <- phytools::pbtree(n = 12, scale = 1)
#'
#'   # Initialize global mapping to state "0" (baseline)
#'   tr0 <- phytools::paintBranches(tr, edge = unique(tr$edge[,2]),
#'                                  state = "0", anc.state = "0")
#'
#'   # Add first shift on an internal node
#'   nd <- ape::Ntip(tr0) + 2L
#'   res1 <- addShiftToModel(tr0, shift_node = nd, current_shift_id = 0L)
#'   res1$shift_id        # 1
#'   tr1 <- res1$tree     # SIMMAP with regime "1" painted on that clade
#'
#'   # Add a second shift somewhere else
#'   nd2 <- ape::Ntip(tr1) + 3L
#'   res2 <- addShiftToModel(tr1, shift_node = nd2, current_shift_id = res1$shift_id)
#'   res2$shift_id        # 2
#' }
#' @keywords internal
#' @noRd
addShiftToModel <- function(tree, shift_node, current_shift_id) {
  # Update the shift ID
  next_shift_id <- current_shift_id + 1

  # Paint the subtree with the new regime/shift id
  painted_tree <- paintSubTree_mod(tree, node = shift_node, state = as.character(next_shift_id), overwrite= FALSE, stem = FALSE)

  # Return a list with the updated tree and the new shift ID
  return(list(tree = painted_tree, shift_id = next_shift_id))
}

#' Remove a Painted Shift from a SIMMAP Tree
#'
#' Selectively removes a previously painted regime (shift) from the clade
#' descending from \code{shift_node} in a SIMMAP-style phylogenetic tree, and
#' restores the ancestral state inherited from the parent node. Optionally, the
#' stem (incoming edge) to \code{shift_node} can be adjusted as well.
#'
#' @param tree A phylogenetic tree of class \code{phylo} or \code{simmap}.
#'   If edge lengths are missing, they are initialized via
#'   \code{\link[ape]{compute.brlen}}. If \code{$maps} is \code{NULL}, a simple
#'   mapping is initialized with state \code{"1"} on all edges.
#' @param shift_node Integer node ID in \code{tree} at which the painted shift
#'   begins; the entire descendant clade is considered for removal.
#' @param stem Logical or numeric. If \code{FALSE} (default), the incoming edge
#'   to \code{shift_node} is left unchanged. If \code{TRUE}, the incoming edge
#'   is reassigned to the parental state. If a numeric value in \eqn{[0,1]},
#'   the parent edge is split into two segments: \code{1 - stem} and \code{stem}
#'   both labeled with the parental state (i.e., effectively restored).
#'
#' @return A modified tree of class \code{c("simmap","phylo")} with updated
#'   \code{$maps} and \code{$mapped.edge}, where edges in the target subtree
#'   that matched the shift state are overwritten by the parental state.
#'
#' @details
#' \itemize{
#'   \item The parental state is determined using \code{\link[phytools]{getParent}}
#'         and \code{\link[phytools]{getStates}} at the parent node of
#'         \code{shift_node}. If unavailable (e.g., missing parent), the default
#'         ancestral state \code{"1"} is used.
#'   \item Edges are selectively overwritten only when their current mapping
#'         equals the shift state found at \code{shift_node}. Existing mappings
#'         with other states are preserved.
#'   \item This implementation assumes edges to be represented as single-state
#'         segments when testing equality (i.e., \code{names(maps[[i]])} of
#'         length 1). If edges contain multiple segments, only exact single-state
#'         matches are overwritten.
#' }
#'
#' @seealso
#' \code{paintSubTree_mod},  # internal helper
#' \link[phytools]{paintSubTree},
#' \link[phytools]{paintBranches},
#' \link[phytools]{getParent},
#' \link[phytools]{getDescendants},
#' \link[phytools]{getStates},
#' \link[ape]{compute.brlen}
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   tr <- phytools::pbtree(n = 10, scale = 1)
#'   # Initialize a simple global mapping "0" on all edges
#'   tr0 <- phytools::paintBranches(tr, edge = unique(tr$edge[,2]), state = "0", anc.state = "0")
#'   # Paint a subtree as shift "1"
#'   nd <- ape::Ntip(tr0) + 2L
#'   tr1 <- paintSubTree_mod(tr0, node = nd, state = "1", anc.state = "0",
#'     stem = TRUE, overwrite = TRUE)
#'   # Remove that shift and restore parental state
#'   tr2 <- paintSubTree_removeShift(tr1, shift_node = nd, stem = FALSE)
#' }
#' @importFrom ape compute.brlen
#' @importFrom phytools getParent getDescendants getStates
#' @keywords internal
#' @noRd
removeShiftFromTree <- function(tree, shift_node, stem = FALSE) {
  #print(paste("Removing shift from node:", shift_node))

  # Retrieve node states; names of this vector are node indices
  node_states <- getStates(tree, type="nodes")
  #print("Current node states:")
  #print(node_states)

  # Using phytools' getParent to find the parent node index
  parent_node <- phytools::getParent(tree, shift_node)
  #print(paste("Parent node of", shift_node, "is:", parent_node))

  # Indicate the state of the shift node
  shift_state <- node_states[as.character(shift_node)]
  #print(paste("State of shift node", shift_node, "is:", shift_state))

  # Check if parent node index is valid
  if (!is.na(parent_node) && parent_node %in% names(node_states)) {
    # Get the state of the parent node
    parent_state <- node_states[as.character(parent_node)]
    #print(paste("State of parent node", parent_node, "is:", parent_state))

    # Check if parent state is NA
    if (!is.na(parent_state)) {
      # Paint the subtree at the shift node with the parent's state, without overwriting descendants
      #print(paste("Painting subtree at node", shift_node, "with state", parent_state, "to remove shift"))
      tree <- paintSubTree_removeShift(tree, shift_node, stem = stem)  # Using the specialized function for shift removal
    } else {
      #print(paste("State of parent node", parent_node, "is NA. Cannot remove shift."))
    }
  } else {
    #print("Invalid parent node. Cannot remove shift.")
  }

  return(tree)
}

#' Identify Internal Nodes Corresponding to Regime Shifts
#'
#' Scans a SIMMAP-style phylogenetic tree and returns the set of internal nodes
#' that correspond to regime shifts, as inferred from the most recent common
#' ancestor (MRCA) of tips sharing the same painted state.
#'
#' @param tree A SIMMAP-formatted phylogenetic tree (class \code{c("simmap","phylo")}).
#'   The tree must have discrete state mappings, accessible via
#'   \code{\link[phytools]{getStates}}.
#'
#' @return An integer vector of unique node numbers representing the MRCAs of
#'   clades in which a distinct state is expressed. If no shifts are found
#'   (e.g., all tips share the same state), returns an empty vector.
#'
#' @details
#' For each unique state painted on the tips, the function:
#' \enumerate{
#'   \item Identifies all tips carrying that state.
#'   \item Uses \code{\link[ape]{getMRCA}} to compute the MRCA node of those tips
#'         (if there are two or more tips).
#'   \item Collects all such MRCA nodes across states, ensuring uniqueness.
#' }
#'
#' @seealso
#' \code{\link[phytools]{getStates}},
#' \code{\link[ape]{getMRCA}}
#'
#' @note
#' Internally related helpers include \code{paintSubTree_mod}, \code{addShiftToModel},
#' and \code{paintSubTree_removeShift}; these are not exported in the current version.
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   tr <- phytools::pbtree(n = 12, scale = 1)
#'
#'   # Paint global state "0"
#'   tr0 <- phytools::paintBranches(tr, edge = unique(tr$edge[,2]), state = "0", anc.state = "0")
#'
#'   # Add a shift on a subtree as state "1"
#'   nd <- ape::Ntip(tr0) + 3L
#'   tr1 <- paintSubTree_mod(tr0, node = nd, state = "1", anc.state = "0", overwrite = TRUE)
#'
#'   # Identify shift nodes
#'   whichShifts(tr1)
#' }
#' @importFrom phytools getStates
#' @importFrom ape getMRCA
#' @keywords internal
#' @noRd
whichShifts <- function(tree) {
  tip_states <- getStates(tree, type = "tips")
  unique_states <- unique(tip_states)

  shift_nodes <- c()
  for (state in unique_states) {
    tips_with_state <- names(tip_states[tip_states == state])
    if (length(tips_with_state) > 1) {
      # Use ape::getMRCA to find the most recent common ancestor
      mrca_node <- ape::getMRCA(tree, tips_with_state)
      shift_nodes <- c(shift_nodes, mrca_node)
    }
  }

  return(unique(shift_nodes))
}

#' Extract Regime-Specific Variance–Covariance Matrices from a BMM mvgls Model
#'
#' Retrieves and scales the regime-specific variance–covariance (VCV) matrices
#' from a \code{\link[mvMORPH]{mvgls}} model fitted under a
#' Brownian Motion with multiple regimes (\code{"BMM"}) model.
#'
#' @param model_output An object returned by \code{\link[mvMORPH]{mvgls}}
#'   that was fitted with \code{model = "BMM"} and contains:
#'   \describe{
#'     \item{\code{param}}{A named numeric vector of regime-specific evolutionary rates.}
#'     \item{\code{sigma}}{A list that must contain \code{Pinv}, the precision matrix
#'       (inverse of the phylogenetic covariance matrix) used for model fitting.}
#'   }
#'
#' @return A named list of variance–covariance matrices (one per regime).
#' Each element is a numeric matrix with dimensions equal to the number of
#' traits in the model. The first regime is returned unscaled, and subsequent
#' regimes are scaled relative to the first regime's rate parameter.
#'
#' @details
#' The first regime's covariance matrix is taken directly from
#' \code{model_output$sigma$Pinv}. For each subsequent regime, the covariance
#' matrix is obtained by multiplying the base covariance matrix by the ratio
#' of the regime's rate parameter to the first regime's parameter.
#'
#' This function is particularly useful for extracting and comparing regime-specific
#' evolutionary rates from BMM models fitted with \code{\link[mvMORPH]{mvgls}}.
#'
#' @note
#' If \code{model_output} does not contain the required components
#' (\code{param}, \code{sigma}, and \code{sigma$Pinv}), the function
#' returns \code{NULL}.
#'
#' @seealso
#' \code{\link[mvMORPH]{mvgls}}
#'
#' @note
#' This function is typically used internally by model-fitting helpers such as
#' `fitMvglsAndExtractGIC.formula()` and `fitMvglsAndExtractBIC.formula()`.
#' These helpers are not exported in the current version.
#'
#' @examples
#' \dontrun{
#'   library(mvMORPH)
#'   set.seed(123)
#'   tree <- ape::rtree(5)
#'   dat <- data.frame(trait = rnorm(5))
#'   rownames(dat) <- tree$tip.label
#'
#'   # Fit a simple BMM model
#'   fit <- mvgls(trait ~ 1, tree = tree, model = "BMM")
#'
#'   # Extract regime-specific VCVs
#'   extractRegimeVCVs(fit)
#' }
#' @keywords internal
#' @noRd
extractRegimeVCVs <- function(model_output) {
  # Ensure the required components are in the model_output
  if (!"param" %in% names(model_output) || !"sigma" %in% names(model_output) || !"Pinv" %in% names(model_output$sigma)) {
    return(NULL)
    #stop("model_output does not contain the required components.")
  }

  # Extract the covariance matrix (Pinv) for the first regime (base VCV)
  base_Pinv <- model_output$sigma$Pinv

  # Get the parameter for the first regime (base rate)
  base_param <- model_output$param[1]

  # List to store VCVs for each regime
  vcv_list <- list()

  # Iterate through the parameters and calculate VCV for each regime
  param_names <- names(model_output$param)
  for (i in seq_along(param_names)) {
    regime_name <- param_names[i]
    regime_param <- model_output$param[regime_name]

    # Scale the precision matrix to get the regime's covariance matrix
    # For the first regime, no scaling is needed
    if (i == 1) {
      regime_vcv <- base_Pinv
    } else {
      regime_vcv <- base_Pinv * (regime_param / base_param)
    }

    # Add the VCV matrix to the list, using the regime's parameter name
    vcv_list[[regime_name]] <- regime_vcv
  }

  return(vcv_list)
}

#' Get All Descendants of a Node in a Phylogenetic Tree
#'
#' Recursively retrieves all descendant nodes (internal nodes and tips) from a specified
#' node in a phylogenetic tree of class \code{"phylo"}.
#'
#' @param tree An object of class \code{"phylo"} (see \code{\link[ape]{phylo}}).
#' @param node An integer specifying the node number from which to retrieve descendants.
#'   Internal nodes are numbered from \code{Ntip(tree) + 1} to
#'   \code{Ntip(tree) + Nnode(tree)} in the \code{phylo} object.
#' @param include.node Logical (default \code{FALSE}). If \code{TRUE}, the specified
#'   node itself is included in the returned vector.
#'
#' @return An integer vector of node numbers (tips and/or internal nodes)
#'   that descend from the specified node. The order of descendants is not guaranteed.
#'
#' @details
#' This function performs a depth-first traversal of the tree, recursively
#' identifying all nodes downstream of the specified starting node.
#'
#' @seealso
#' \code{\link[ape]{Ntip}}, \code{\link[ape]{Nnode}},
#' \code{\link[ape]{getMRCA}}, \code{\link[phytools]{getParent}}
#'
#' @examples
#' library(ape)
#' set.seed(123)
#' tree <- rtree(5)
#'
#' # Get descendants of the root node
#' root_node <- Ntip(tree) + 1
#' getDescendants(tree, root_node)
#'
#' # Include the starting node itself
#' getDescendants(tree, root_node, include.node = TRUE)
#' @keywords internal
#' @noRd
getDescendants <- function(tree, node, include.node = FALSE) {
  # Function to recursively find all descendants of a node
  descendants <- numeric(0)
  for (i in which(tree$edge[,1] == node)) {
    descendants <- c(descendants, tree$edge[i,2], getDescendants(tree, tree$edge[i,2]))
  }
  if (include.node) descendants <- c(node, descendants)
  return(unique(descendants))
}
