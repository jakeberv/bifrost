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
#' @param min_tips Integer (≥1). Minimum number of descendant tips required for an
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
#' \donttest{
#'   set.seed(123)
#'   tree <- phytools::pbtree(n = 30, scale = 1)
#'   painted <- generatePaintedTrees(tree, min_tips = 5, state = "shift")
#'   length(painted)
#'   names(painted)[1]
#' }
#'
#' @importFrom ape is.rooted root Ntip Nnode
#' @importFrom phytools getDescendants paintSubTree
#' @export
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
  cat(paste(length(eligible_nodes), "eligible nodes are detected", '\n'))

  painted_trees <- list()

  for (node in eligible_nodes) {
    tree_copy <- tree
    tree_copy <- paintSubTree(tree_copy, node, state)
    painted_trees[[paste("Node", node)]] <- tree_copy
  }

  cat(paste(length(painted_trees), "sub-models generated", '\n'))
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
#' \donttest{
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
#' @export
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
#' \donttest{
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
#' @export
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
#' \donttest{
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
#' @export
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
#' \donttest{
#'   set.seed(321)
#'   tree <- phytools::pbtree(n = 12, scale = 1)
#'   painted <- phytools::paintSubTree(tree, node = 13, state = "A", anc.state = "A")
#'   x <- rnorm(12)
#'   y <- x + rnorm(12)
#'   dat <- data.frame(x = x, y = y)
#'   rownames(dat) <- painted$tip.label
#'   result <- fitMvglsAndExtractBIC.formula("y ~ x", painted, dat)
#'   result$BIC
#' }
#'
#' @importFrom mvMORPH mvgls
#' @importFrom phytools getStates
#' @importFrom stats as.formula BIC
#' @export
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

#' Calculate ΔGIC Scores Relative to a Baseline Model
#'
#' Computes the difference in Generalized Information Criterion (GIC) scores
#' between a baseline model (assumed to be the first element of `model_results`)
#' and one or more alternative models fitted to painted phylogenetic trees.
#' This allows rapid comparison of model fit across candidate regime-paintings.
#'
#' @param model_results A list of model fit results, typically the output of
#'   \code{\link{fitMvglsAndExtractGIC.formula}} or a similar function.
#'   Each element must be a list containing at least a \code{GIC} component,
#'   where \code{GIC$GIC} is a numeric scalar.
#'   The first element is treated as the baseline model.
#'
#' @param painted_tree_list A named list of painted phylogenetic trees
#'   (objects of class \code{simmap}). The names are used as identifiers
#'   for the output vector of ΔGIC values. Every element must have a
#'   non-\code{NULL} name.
#'
#' @return A named numeric vector of ΔGIC scores, where each value is:
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
#' \code{\link[mvMORPH]{GIC}}, \code{\link{fitMvglsAndExtractGIC.formula}},
#' \code{\link[phytools]{paintSubTree}}
#'
#' @examples
#' \donttest{
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
#' @export
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

#' Paint a Subtree in a Phylogenetic Tree with Optional Selective Overwriting
#'
#' This function modifies a phylogenetic tree (of class "phylo") by painting a specified subtree
#' starting at a given node with a new state, optionally preserving existing state mappings unless overwritten.
#' It returns a SIMMAP-style tree with updated edge mappings and supports both full and selective painting,
#' as well as optional stem painting from the parent edge.
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

#' Remove a Painted Shift from a Subtree in a SIMMAP-Formatted Phylogenetic Tree
#'
#' This function removes a previously painted shift (regime change) from a specified node and its descendant branches
#' in a SIMMAP-style phylogenetic tree by selectively overwriting branches painted with the shift state.
#' The original ancestral state (inherited from the parent node) is restored, with optional stem painting of the edge leading to the node.
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

#' Add a New Shift to a Phylogenetic Tree Model
#'
#' This function adds a new evolutionary regime (shift) to a SIMMAP-style phylogenetic tree
#' by painting the subtree starting at a specified node with a new unique state identifier.
#' The shift ID is incremented and returned for tracking multiple shifts in downstream modeling.
addShiftToModel <- function(tree, shift_node, current_shift_id) {
  # Update the shift ID
  next_shift_id <- current_shift_id + 1

  # Paint the subtree with the new regime/shift id
  painted_tree <- paintSubTree_mod(tree, node = shift_node, state = as.character(next_shift_id), overwrite=F, stem = F)

  # Return a list with the updated tree and the new shift ID
  return(list(tree = painted_tree, shift_id = next_shift_id))
}

#' Remove a Shift from a Phylogenetic Tree by Reverting to Parent State
#'
#' This function removes a regime shift from a SIMMAP-style phylogenetic tree by repainting the subtree
#' starting at a specified node with the state of its parent node. It uses `paintSubTree_removeShift()`
#' to selectively revert the painted shift without affecting unrelated branches, with optional stem edge handling.
removeShiftFromTree <- function(tree, shift_node, stem=F) {
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
      tree <- paintSubTree_removeShift(tree, shift_node, stem=stem)  # Using the specialized function for shift removal
    } else {
      #print(paste("State of parent node", parent_node, "is NA. Cannot remove shift."))
    }
  } else {
    #print("Invalid parent node. Cannot remove shift.")
  }

  return(tree)
}

#' Identify Nodes Representing Shifts in a SIMMAP Tree
#'
#' This function identifies internal nodes in a SIMMAP-style phylogenetic tree
#' that represent regime shifts, based on the most recent common ancestors (MRCAs)
#' of tips sharing the same painted state.
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

#' Extract Regime-Specific Variance-Covariance Matrices from a BMM mvgls Model
#'
#' This function extracts and scales the regime-specific variance-covariance (VCV) matrices
#' from a mvgls model fitted under the BMM (Brownian Motion with multiple regimes) framework.
#' It returns a named list of VCV matrices, each corresponding to an evolutionary regime.
extractRegimeVCVs <- function(model_output) {
  # Ensure the required components are in the model_output
  if (!"param" %in% names(model_output) || !"sigma" %in% names(model_output) || !"Pinv" %in% names(model_output$sigma)) {
    return(NULL)
    stop("model_output does not contain the required components.")
  }

  # Extract the precision matrix (Pinv) for the first regime (base VCV)
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

# # Supporting function getDescendants (if not already defined in the environment)
#' Get All Descendants of a Node in a Phylogenetic Tree
#'
#' Recursively retrieves all descendant nodes (internal and/or tips) from a specified node
#' in a phylogenetic tree of class `"phylo"`. Optionally includes the starting node in the output.
getDescendants <- function(tree, node, include.node = FALSE) {
  # Function to recursively find all descendants of a node
  descendants <- numeric(0)
  for (i in which(tree$edge[,1] == node)) {
    descendants <- c(descendants, tree$edge[i,2], getDescendants(tree, tree$edge[i,2]))
  }
  if (include.node) descendants <- c(node, descendants)
  return(unique(descendants))
}
