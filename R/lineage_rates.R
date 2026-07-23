#' Compute Weighted Lineage Rates
#'
#' Compute the lineage-rate summary statistic described by Berv et al. (2026)
#' from mapped-state `phylo`/`simmap` tree objects (Paradis and Schliep 2019;
#' Revell 2012, 2024), including the mapped trees returned by
#' \emph{\code{bifrost}}. The function calculates present-biased weighted
#' summaries of the state values inherited at ancestral nodes along each
#' root-to-tip path. For \emph{\code{bifrost}} BMM fits, the state values are
#' the scalar regime-rate parameters stored in
#' `bifrost_search$model_no_uncertainty$param`; for generic input, set
#' `bifrost_search = NULL` and provide `tree` with a named `state_values`
#' vector. Results include shift counts, terminal-state values, and the summary
#' statistic, returned on the log scale by default.
#'
#' @param bifrost_search A `bifrost_search` object returned by
#'   \emph{\code{bifrost}}, a plain list with the required
#'   \emph{\code{bifrost}} components, or `NULL`. \emph{\code{bifrost}} inputs
#'   must include a fitted
#'   `model_no_uncertainty` component that identifies a multi-regime
#'   Brownian-motion model (`model = "BMM"`) and provides named numeric
#'   regime-rate parameters in `bifrost_search$model_no_uncertainty$param`.
#'   Use `NULL` for generic input mode, where `tree` and `state_values` are
#'   required. Supply either `bifrost_search` or both `tree` and `state_values`,
#'   but not both input modes at once. Plain-list \emph{\code{bifrost}} inputs
#'   must include
#'   `tree_no_uncertainty_untransformed` and `model_no_uncertainty`.
#' @param tree SIMMAP-style `phylo`/`simmap` tree (Paradis and Schliep 2019;
#'   Revell 2012, 2024) required when `bifrost_search = NULL` for generic input
#'   mode. Must be `NULL` when `bifrost_search` is supplied. The tree must
#'   include SIMMAP `maps` state data from which node and tip states can be
#'   recovered. For the age-decay interpretation, the tree should be
#'   time-scaled; non-ultrametric trees are accepted, with age calculations
#'   controlled by `age_reference`. Internal nodes must use the conventional
#'   `phylo` numbering in which the root is `Ntip(tree) + 1`.
#' @param state_values Named numeric vector required in generic input mode and
#'   invalid for \emph{\code{bifrost}} inputs. Names must include every SIMMAP
#'   state label in the tree; extra names are ignored. Values must be finite,
#'   and must be strictly positive when `log = TRUE`.
#' @param decay_base Numeric exponential-decay base `b` in the weighting kernel.
#'   Must be at least `1`. Larger values make weights decline more strongly
#'   as elements get farther from the present; `decay_base = 1` gives equal age
#'   weights regardless of `half_life`.
#' @param half_life Optional decay time scale `T`, in the same time units as
#'   the tree. With the default `decay_base = 2`, `T` is a true half-life: an
#'   element `T` time units farther from the present receives half the
#'   unnormalized weight. More generally, weights are multiplied by
#'   `1 / decay_base` over each interval of length `T`. Set `half_life = NULL`
#'   to use the unscaled kernel `b^(-age)`; this only changes the weights when
#'   `decay_base > 1`.
#' @param normalize_weights Logical; divide raw age-decay weights by their
#'   lineage-specific sum before averaging. Berv et al. (2026) used `TRUE`,
#'   which keeps the statistic as a weighted mean and makes lineages comparable
#'   when root-to-tip paths include different numbers of ancestral nodes. If
#'   every included ancestral node has the same state, the normalized value
#'   equals that state's supplied value or fitted rate.
#' @param age_reference Character; reference point for converting node heights
#'   into ages before the present. The default, `"tree"`, measures ages relative
#'   to the maximum tree height and preserves the original global-present
#'   behavior. Use `"tip"` for heterochronous or otherwise non-ultrametric trees
#'   when recency should be measured before each terminal tip's own sampling
#'   point. These choices are identical for ultrametric trees and mainly affect
#'   results when `normalize_weights = FALSE`.
#' @param log Logical; return log-scale lineage-rate columns by default.
#'   \emph{\code{bifrost}} BMM rates must always be finite and strictly positive.
#'   Generic `state_values` must be strictly positive when `log = TRUE`; when
#'   `log = FALSE`, generic values can remain on their original scale.
#' @param cores Integer number of future workers. Values greater than one use
#'   `future::multisession` through [future.apply::future_lapply()], which is
#'   platform agnostic.
#' @param progress Logical; show a local `progressr` text progress bar for
#'   the tip-wise computations. Defaults to `interactive()`, so progress is
#'   shown in interactive sessions and suppressed in non-interactive runs.
#' @return A `data.frame` with one row per tip. The stable top-level columns are
#'   `tip.label`,
#'   `shift_count`, the ancestral-node summary statistic, `tip_state`, and
#'   `tip_rate`. With `log = TRUE`, the summary statistic is named
#'   `log_lineage_rate`, and `log_tip_rate` is included. With `log = FALSE`,
#'   the summary statistic is named `lineage_rate`.
#'
#'   Additional diagnostic tables are stored as attributes. `branch_metrics`
#'   contains the per-tip intermediate table for the top-level ancestral-node
#'   statistic. Its columns are `Tip`, `Shift_Count`, `Total_Time`,
#'   `Shift_Rate_Per_Time`, `Shift_Rate_Per_Speciation`,
#'   `Weighted_Lineage_Value`, `Tip_State`, and `Tip_State_Value`.
#'   `shift_metrics` contains the same concise shift-count and time diagnostics
#'   without `Weighted_Lineage_Value`; its columns are `Tip`, `Shift_Count`,
#'   `Total_Time`, `Shift_Rate_Per_Time`, `Shift_Rate_Per_Speciation`,
#'   `Tip_State`, and `Tip_State_Value`. `Shift_Count` and both shift-rate
#'   columns are unweighted counts or ratios of detected parent-child node-state
#'   transitions. The `settings` attribute records the resolved input mode,
#'   model family, tree name, weighting parameters, age reference, log setting,
#'   and execution settings used for the calculation.
#' @details
#' \subsection{Overview}{
#' The summary statistic follows Berv et al. (2026) as a descriptive,
#' present-biased summary of each tip's inferred mapped-state history. In the
#' default BMM use case, it summarizes heterogeneous phenotypic tempo along a
#' tip-to-root path while giving more influence to recently inherited regimes.
#' The calculation uses the node and tip states recovered from the supplied
#' SIMMAP-style mapped-state tree. The parameters `decay_base`, `half_life`, and
#' `normalize_weights` define the weighting kernel used to explore how inferred
#' lineage summaries change as deeper history is weighted more or less
#' strongly.
#'
#' In the default \emph{\code{bifrost}} mode, `lineage_rates()` reads
#' regime-specific Brownian variance rates from
#' `bifrost_search$model_no_uncertainty$param` and summarizes how those rates
#' are inherited along mapped root-to-tip histories. That interpretation
#' requires a heterogeneous Brownian-motion model: a multi-regime BMM fit with
#' at least two mapped states.
#' }
#'
#' \subsection{Method and formula}{
#' Using the notation of Berv et al. (2026), generalized to an arbitrary decay
#' base `b` and state-associated value \eqn{r_i}, the normalized summary
#' statistic is
#'
#' \deqn{
#' \bar{r}_{\mathrm{WLR}} =
#' \frac{1}{Z} \sum_{i=1}^{L} b^{-a_i/T} \cdot r_i
#' }{
#' r_bar_WLR = (1 / Z) sum_{i = 1}^L b^(-a_i / T) * r_i
#' }
#'
#' where the normalizing constant is
#'
#' \deqn{
#' Z = \sum_{j=1}^{L} b^{-a_j/T}.
#' }{
#' Z = sum_{j = 1}^L b^(-a_j / T).
#' }
#'
#' Here, \eqn{\bar{r}_{\mathrm{WLR}}} is the weighted mean lineage value,
#' \eqn{r_i} is the state-associated value for node or path element \eqn{i},
#' \eqn{a_i} is that element's age before the selected reference point,
#' \eqn{b} is `decay_base`, \eqn{T} is the decay time scale, and \eqn{L} is the
#' number of nodes or path elements included in the statistic. For the
#' branch-based summary statistic aligned with Berv et al. (2026), \eqn{i}
#' indexes ancestral nodes along the lineage, including the root and excluding
#' the terminal tip state. This is the expected structure of \emph{\code{bifrost}}
#' mapped trees, where regime changes are represented by changes in node states
#' rather than integrated over within-branch SIMMAP segments.
#'
#' In \emph{\code{bifrost}} output, \eqn{r_i} is represented by the scalar BMM
#' regime-rate parameters in `bifrost_search$model_no_uncertainty$param`;
#' following Berv et al. (2026), these correspond to the arithmetic mean of the
#' diagonal variance terms in each regime's estimated evolutionary
#' variance-covariance matrix. This is the weighted phenotypic-rate
#' interpretation used in Berv et al. (2026). In generic input mode, \eqn{r_i}
#' is the corresponding value from `state_values`.
#' }
#'
#' \subsection{Weighting, normalization, and age reference}{
#' Berv et al. (2026) used the base-2 case, and the defaults
#' (`decay_base = 2`, `half_life = 5`, and `normalize_weights = TRUE`) preserve
#' that weighting scheme. With these defaults, a regime that is `half_life` time
#' units farther from the present receives half the unnormalized weight.
#' `decay_base = 1` gives every element raw weight 1 regardless of `half_life`,
#' so the summary statistic becomes an equal-weighted mean when
#' `normalize_weights = TRUE`. Values below `1` are rejected because they would
#' give greater weight to regimes farther from the present.
#'
#' The kernel terms \eqn{b^{-a_i/T}} are raw age-decay weights. When
#' `normalize_weights = TRUE`, the function divides those raw weights by
#' \eqn{Z} within each lineage before averaging. This keeps the statistic as a
#' weighted mean of state-associated values, rather than letting lineages with
#' more ancestral nodes accumulate larger values simply because more terms were
#' summed. If every included ancestral node is assigned to the same state, the
#' normalized value equals that state's supplied value or fitted rate regardless
#' of how many nodes occur along the path. When `normalize_weights = FALSE`, the
#' \eqn{1/Z} normalization is skipped.
#'
#' The mapped tree should be time-scaled for the age-decay weights to have a
#' direct temporal interpretation. The `age_reference` argument controls how the
#' ages \eqn{a_i} are measured. With the default `age_reference = "tree"`, ages
#' are measured relative to the maximum tree height, matching the original
#' global-present interpretation. For heterochronous or otherwise
#' non-ultrametric trees, `age_reference = "tip"` measures recency before each
#' terminal tip's own sampling point. These choices are identical for
#' ultrametric trees. When `normalize_weights = TRUE`, they give the same
#' top-level normalized values because all raw weights within a lineage are
#' rescaled by the same constant before normalization. The choice mainly
#' affects results when `normalize_weights = FALSE`.
#'
#' When `half_life = NULL`, the unscaled kernel \eqn{b^{-a_i}} is used instead
#' of \eqn{b^{-a_i/T}}. This removes the half-life divisor from the exponent and
#' only changes the weights when `decay_base > 1`. Parameter combinations that
#' make a mathematically positive decay weight unrepresentable in double
#' precision are rejected rather than returning a non-finite summary.
#' }
#'
#' \subsection{Interpretation}{
#' Interpret this as a flexible descriptive statistic for an inferred mapped
#' history, not as a model-derived estimator with a unique optimal weighting
#' scheme. It is not an instantaneous tip-rate estimate or a strict
#' recent-time-window average. The default statistic is node-based: it uses
#' ancestral node states along each root-to-tip path and does not integrate over
#' within-edge SIMMAP segment durations. For SIMMAP-style trees with mapped
#' changes inside branches, those segment durations are not part of the
#' top-level statistic, and within-edge changes are not counted as shifts unless
#' they are represented by parent-child node-state differences.
#' }
#'
#' \subsection{Lineage and shift diagnostics}{
#' The top-level `log_lineage_rate` or `lineage_rate` column is computed from
#' ancestral node states along each root-to-tip path. Because the statistic is
#' node-based and normalized within each lineage, its temporal resolution
#' depends on the ages and density of ancestral nodes. Lineages with many recent
#' nodes can downweight regimes farther from the present more rapidly, whereas
#' sparse lineages may retain more influence from regimes farther from the
#' present.
#'
#' The `branch_metrics` and `shift_metrics` attributes store concise unweighted
#' shift diagnostics: the number of detected parent-child node-state changes,
#' root-to-tip time, shifts per unit time, and shifts per speciation event.
#' `branch_metrics` additionally stores the weighted ancestral-node intermediate
#' used for the top-level lineage-rate column. `tip_rate` is the fitted rate or
#' supplied value for the terminal SIMMAP state. Compare the terminal value to
#' the lineage-rate columns when you want to assess how much the inherited state
#' history differs from the terminal regime state.
#' }
#'
#' \subsection{Generic input mode}{
#' `lineage_rates()` can also operate directly on a SIMMAP-style tree whose
#' mapped node states refer to externally estimated values supplied through
#' `state_values`. In this generic input mode, no fitted model is inspected:
#' the same node-based summary statistic is applied to the user-supplied state
#' values. This can be useful for summarizing arbitrary state-associated
#' quantities on a mapped tree when the node-state interpretation is
#' appropriate. The output column names retain the rate terminology used by the
#' \emph{\code{bifrost}} workflow for compatibility, but the values represent
#' the supplied `state_values` and their interpretation depends on what those
#' values measure. For generic SIMMAP trees, use this helper when ancestral
#' node states are the summary target; it is not a duration-weighted integrator
#' over mapped within-edge segments.
#' }
#' @references
#' Berv, J. S. et al. (2026). Rates of passerine body plan evolution in time
#' and space. \emph{Nature Ecology & Evolution}.
#' doi:10.1038/s41559-026-03110-5.
#'
#' Clavel, J., Escarguel, G., and Merceron, G. (2015). mvMORPH: an R package
#' for fitting multivariate evolutionary models to morphometric data.
#' \emph{Methods in Ecology and Evolution}, 6, 1311-1319.
#' doi:10.1111/2041-210X.12420.
#'
#' Clavel, J., Aristide, L., and Morlon, H. (2019). A penalized likelihood
#' framework for high-dimensional phylogenetic comparative methods and an
#' application to new-world monkeys brain evolution. \emph{Systematic Biology},
#' 68, 93-116. doi:10.1093/sysbio/syy045.
#'
#' Paradis, E., and Schliep, K. (2019). ape 5.0: an environment for modern
#' phylogenetics and evolutionary analyses in R. \emph{Bioinformatics}, 35,
#' 526-528. doi:10.1093/bioinformatics/bty633.
#'
#' Revell, L. J. (2012). phytools: an R package for phylogenetic comparative
#' biology (and other things). \emph{Methods in Ecology and Evolution}, 3,
#' 217-223. doi:10.1111/j.2041-210X.2011.00169.x.
#'
#' Revell, L. J. (2024). phytools 2.0: an updated R ecosystem for phylogenetic
#' comparative methods (and other things). \emph{PeerJ}, 12, e16505.
#' doi:10.7717/peerj.16505.
#' @examples
#' # Generic input mode: summarize arbitrary state-associated values on a mapped tree.
#' toy_tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
#' toy_tree <- phytools::paintSubTree(
#'   tree = toy_tree,
#'   node = ape::Ntip(toy_tree) + 1L,
#'   state = "0"
#' )
#' toy_tree <- phytools::paintSubTree(
#'   tree = toy_tree,
#'   node = ape::Ntip(toy_tree) + 2L,
#'   state = "1"
#' )
#' state_values <- c("0" = 1.2, "1" = 3.4)
#' lineage_rates(
#'   tree = toy_tree,
#'   state_values = state_values,
#'   log = FALSE,
#'   progress = FALSE
#' )
#'
#' \dontrun{
#' # Search-result mode: use the multi-regime BMM fit and SIMMAP tree.
#' lineage_summary <- lineage_rates(bifrost_search = search_result)
#' }
#' @export
lineage_rates <- function(
  bifrost_search = NULL,
  tree = NULL,
  state_values = NULL,
  decay_base = 2,
  half_life = 5,
  normalize_weights = TRUE,
  age_reference = c("tree", "tip"),
  log = TRUE,
  cores = 1L,
  progress = interactive()
) {
  .lineage_rates_check_packages()
  if (!is.numeric(decay_base) ||
        length(decay_base) != 1L ||
        !is.finite(decay_base) ||
        decay_base < 1) {
    stop("`decay_base` must be a finite numeric scalar greater than or equal to 1.", call. = FALSE)
  }
  if (!is.null(half_life) &&
        (!is.numeric(half_life) ||
           length(half_life) != 1L ||
           !is.finite(half_life) ||
           half_life <= 0)) {
    stop("`half_life` must be NULL or a positive finite numeric scalar.", call. = FALSE)
  }
  if (!is.logical(normalize_weights) ||
        length(normalize_weights) != 1L ||
        is.na(normalize_weights)) {
    stop("`normalize_weights` must be TRUE or FALSE.", call. = FALSE)
  }
  age_reference <- .lineage_rates_match_age_ref(age_reference)
  cores <- .lineage_rates_validate_cores(cores)
  if (!is.logical(log) || length(log) != 1L || is.na(log)) {
    stop("`log` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }

  resolved <- .lineage_rates_resolve_inputs(
    bifrost_search = bifrost_search,
    tree = tree,
    state_values = state_values,
    log = log
  )
  mapped_tree <- resolved$tree
  state_values <- resolved$state_values
  .lineage_rates_validate_root_numbering(mapped_tree)

  shift_metrics <- .lineage_rates_metric_table(
    mapped_tree = mapped_tree,
    mode = "shifts",
    decay_base = decay_base,
    state_values = state_values,
    normalize_weights = normalize_weights,
    age_reference = age_reference,
    half_life = half_life,
    cores = cores,
    progress = progress
  )
  branch_metrics <- .lineage_rates_metric_table(
    mapped_tree = mapped_tree,
    mode = "branches",
    decay_base = decay_base,
    state_values = state_values,
    normalize_weights = normalize_weights,
    age_reference = age_reference,
    half_life = half_life,
    cores = cores,
    progress = progress
  )

  if (isTRUE(log)) {
    out <- data.frame(
      tip.label = branch_metrics$Tip,
      shift_count = shift_metrics$Shift_Count,
      log_lineage_rate = base::log(branch_metrics$Weighted_Lineage_Value),
      tip_state = shift_metrics$Tip_State,
      tip_rate = shift_metrics$Tip_State_Value,
      log_tip_rate = base::log(shift_metrics$Tip_State_Value),
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(
      tip.label = branch_metrics$Tip,
      shift_count = shift_metrics$Shift_Count,
      lineage_rate = branch_metrics$Weighted_Lineage_Value,
      tip_state = shift_metrics$Tip_State,
      tip_rate = shift_metrics$Tip_State_Value,
      stringsAsFactors = FALSE
    )
  }

  attr(out, "branch_metrics") <- branch_metrics
  attr(out, "shift_metrics") <- shift_metrics
  attr(out, "settings") <- list(
    tree_name = resolved$tree_name,
    input_mode = resolved$input_mode,
    model_family = resolved$model_family,
    weighting = "exponential_age_decay",
    decay_base = decay_base,
    half_life = half_life,
    normalize_weights = normalize_weights,
    age_reference = age_reference,
    log = log,
    cores = cores,
    progress = progress
  )
  out
}

.lineage_rates_check_packages <- function() {
  pkgs <- c("ape", "phytools", "future", "future.apply", "progressr")
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0L) {
    stop("Missing required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.lineage_rates_validate_cores <- function(cores) {
  if (!is.numeric(cores) ||
        length(cores) != 1L ||
        !is.finite(cores) ||
        cores < 1 ||
        cores != floor(cores) ||
        cores > .Machine$integer.max) {
    stop(
      "`cores` must be a positive integer (a whole number) in R integer range.",
      call. = FALSE
    )
  }
  as.integer(cores)
}

.lineage_rates_match_age_ref <- function(age_reference) {
  tryCatch(
    match.arg(age_reference, choices = c("tree", "tip")),
    error = function(e) {
      stop("`age_reference` must be \"tree\" or \"tip\".", call. = FALSE)
    }
  )
}

.lineage_rates_as_character <- function(x) {
  tryCatch(
    as.character(x),
    error = function(e) character(0)
  )
}

.lineage_rates_model_family <- function(model) {
  if (is.null(model) || (!is.list(model) && !is.environment(model))) {
    return(NA_character_)
  }

  values <- character(0)

  if (!is.null(model$model)) {
    values <- c(values, .lineage_rates_as_character(model$model))
  }
  if (!is.null(model$call)) {
    call_parts <- tryCatch(as.list(model$call), error = function(e) list())
    if (!is.null(call_parts$model)) {
      values <- c(values, .lineage_rates_as_character(call_parts$model))
    }
  }

  values <- toupper(trimws(values))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0L) return(NA_character_)
  values[[1L]]
}

.lineage_rates_check_bmm <- function(model, caller = "lineage_rates()") {
  model_family <- .lineage_rates_model_family(model)
  if (is.na(model_family)) {
    stop(
      "`bifrost_search$model_no_uncertainty` must identify a multi-regime ",
      "BMM fitted model via `$model` or `$call$model`.",
      call. = FALSE
    )
  }

  if (!identical(model_family, "BMM")) {
    stop(
      "`", caller, "` requires a multi-regime BMM fit (`model = \"BMM\"`). ",
      "Single-regime BM/BM1 fits have no heterogeneous regime history to summarize. ",
      "Found model family `", model_family, "`.",
      call. = FALSE
    )
  }

  model_family
}

.lineage_rates_named_numeric <- function(values) {
  is.numeric(values) &&
    !is.null(names(values)) &&
    all(nzchar(names(values)))
}

.lineage_rates_validate_root_numbering <- function(tree) {
  root_nodes <- setdiff(unique(tree$edge[, 1]), unique(tree$edge[, 2]))
  expected_root <- ape::Ntip(tree) + 1L
  if (length(root_nodes) != 1L) {
    stop("Cannot identify a unique root node in `tree`.", call. = FALSE)
  }
  if (!identical(as.integer(root_nodes), expected_root)) {
    stop(
      "`lineage_rates()` requires the root node to be `Ntip(tree) + 1` ",
      "(expected ", expected_root, ", found ", root_nodes, "). Please ",
      "renumber the tree's internal nodes before calling this function.",
      call. = FALSE
    )
  }
  invisible(expected_root)
}

.lineage_rates_resolve_inputs <- function(
  bifrost_search,
  tree = NULL,
  state_values = NULL,
  log = TRUE,
  caller = "lineage_rates()"
) {
  if (is.null(bifrost_search)) {
    if (is.null(tree) && is.null(state_values)) {
      stop(
        "Provide either `bifrost_search` or both `tree` and `state_values`.",
        call. = FALSE
      )
    }
    return(.lineage_rates_generic_inputs(
      tree = tree,
      state_values = state_values,
      log = log
    ))
  }

  if (inherits(bifrost_search, "phylo")) {
    stop(
      "Generic tree inputs must be supplied with `tree`, not `bifrost_search`.",
      call. = FALSE
    )
  }

  if (!is.null(tree) || !is.null(state_values)) {
    stop(
      "`tree` and `state_values` are only used when `bifrost_search = NULL`; ",
      "`bifrost` inputs use `bifrost_search$tree_no_uncertainty_untransformed` ",
      "and `bifrost_search$model_no_uncertainty$param`.",
      call. = FALSE
    )
  }

  if (!inherits(bifrost_search, "bifrost_search") && !is.list(bifrost_search)) {
    stop(
      "`bifrost_search` must be a `bifrost_search` object or plain list ",
      "with required `bifrost` components.",
      call. = FALSE
    )
  }

  tree <- bifrost_search$tree_no_uncertainty_untransformed
  if (is.null(tree) || !inherits(tree, "phylo")) {
    stop(
      "`bifrost_search` must include `tree_no_uncertainty_untransformed`, ",
      "a SIMMAP/phylo tree.",
      call. = FALSE
    )
  }
  if (is.null(tree$maps)) {
    stop("The resolved tree must include SIMMAP `maps` state data.", call. = FALSE)
  }

  model <- bifrost_search$model_no_uncertainty
  model_family <- .lineage_rates_check_bmm(model, caller = caller)
  mapped_states <- unique(as.character(phytools::getStates(tree, type = "both")))
  if (length(mapped_states) < 2L) {
    stop(
      "`", caller, "` requires a multi-regime BMM fit with at least two mapped states; ",
      "single-regime models have no heterogeneous regime history to summarize.",
      call. = FALSE
    )
  }

  regime_rates <- model$param
  if (!.lineage_rates_named_numeric(regime_rates)) {
    stop(
      "`bifrost_search$model_no_uncertainty$param` must be a named numeric ",
      "vector of BMM regime rates.",
      call. = FALSE
    )
  }
  names(regime_rates) <- as.character(names(regime_rates))
  if (anyDuplicated(names(regime_rates))) {
    stop(
      "`bifrost_search$model_no_uncertainty$param` names must be unique.",
      call. = FALSE
    )
  }
  if (!all(is.finite(regime_rates)) || any(regime_rates <= 0)) {
    stop(
      "`bifrost_search$model_no_uncertainty$param` values must be finite ",
      "and strictly positive.",
      call. = FALSE
    )
  }

  missing_states <- setdiff(mapped_states, names(regime_rates))
  if (length(missing_states) > 0L) {
    stop(
      "`bifrost_search$model_no_uncertainty$param` is missing rate values ",
      "for mapped state(s): ",
      paste(missing_states, collapse = ", "),
      call. = FALSE
    )
  }

  list(
    tree = tree,
    state_values = regime_rates,
    tree_name = "tree_no_uncertainty_untransformed",
    input_mode = "bifrost",
    model_family = model_family
  )
}

.lineage_rates_generic_inputs <- function(
  tree,
  state_values,
  log = TRUE
) {
  if (is.null(tree) || !inherits(tree, "phylo")) {
    stop(
      "In generic input mode, `tree` must be a SIMMAP/phylo tree.",
      call. = FALSE
    )
  }
  if (is.null(tree$maps)) {
    stop("The resolved tree must include SIMMAP `maps` state data.", call. = FALSE)
  }

  if (!.lineage_rates_named_numeric(state_values)) {
    stop(
      "In generic input mode, `state_values` must be a named numeric vector ",
      "with names matching SIMMAP state labels.",
      call. = FALSE
    )
  }
  names(state_values) <- as.character(names(state_values))
  if (anyDuplicated(names(state_values))) {
    stop("`state_values` names must be unique.", call. = FALSE)
  }
  if (!all(is.finite(state_values))) {
    stop("`state_values` must be finite.", call. = FALSE)
  }
  if (isTRUE(log) && any(state_values <= 0)) {
    stop("`state_values` must be strictly positive when `log = TRUE`.", call. = FALSE)
  }

  mapped_states <- unique(as.character(phytools::getStates(tree, type = "both")))
  missing_states <- setdiff(mapped_states, names(state_values))
  if (length(missing_states) > 0L) {
    stop(
      "`state_values` is missing values for mapped state(s): ",
      paste(missing_states, collapse = ", "),
      call. = FALSE
    )
  }

  list(
    tree = tree,
    state_values = state_values,
    tree_name = "tree",
    input_mode = "generic_input",
    model_family = NA_character_
  )
}

.lineage_rates_tip_lapply <- function(n_tip, mode, cores, progress, FUN) {
  if (cores > 1L) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = cores)
  }

  .run_jobs <- function() {
    p <- NULL
    if (isTRUE(progress)) {
      p <- progressr::progressor(
        steps = n_tip,
        label = paste("lineage_rates", mode)
      )
    }

    .tip_fun <- function(k) {
      if (!is.null(p)) {
        p(sprintf("%s tips %d/%d", mode, k, n_tip))
      }
      FUN(k)
    }

    if (cores > 1L) {
      future.apply::future_lapply(seq_len(n_tip), .tip_fun, future.seed = FALSE)
    } else {
      lapply(seq_len(n_tip), .tip_fun)
    }
  }

  if (isTRUE(progress)) {
    progressr::with_progress(
      .run_jobs(),
      handlers = progressr::handler_txtprogressbar(style = 3L),
      enable = TRUE
    )
  } else {
    .run_jobs()
  }
}

.lineage_rates_metric_table <- function(
  mapped_tree,
  mode = c("shifts", "branches"),
  decay_base = 2,
  state_values,
  normalize_weights = TRUE,
  age_reference = c("tree", "tip"),
  half_life = 5,
  cores = 1L,
  progress = interactive()
) {
  mode <- match.arg(mode)

  n_tip <- length(mapped_tree$tip.label)
  tip_states <- phytools::getStates(mapped_tree, type = "tips")
  node_states <- phytools::getStates(mapped_tree, type = "nodes")
  names(node_states) <- as.character(names(node_states))
  root_node <- ape::Ntip(mapped_tree) + 1L

  .state_for_node <- function(node) {
    if (node <= n_tip) {
      return(unname(tip_states[[node]]))
    }
    unname(node_states[[as.character(node)]])
  }

  node_depths <- ape::node.depth.edgelength(mapped_tree)
  root_to_tip_distances <- node_depths[seq_len(n_tip)]
  tree_height <- max(phytools::nodeHeights(mapped_tree))
  parent_by_child <- stats::setNames(mapped_tree$edge[, 1], mapped_tree$edge[, 2])

  .ancestor_path <- function(node) {
    out <- integer(0)
    parent <- unname(parent_by_child[as.character(node)])
    while (!is.na(parent)) {
      out <- c(out, as.integer(parent))
      parent <- unname(parent_by_child[as.character(parent)])
    }
    out
  }

  .detect_shift <- function(parent_node, child_node) {
    parent_state <- .state_for_node(parent_node)
    child_state <- .state_for_node(child_node)
    !is.na(parent_state) && !is.na(child_state) && parent_state != child_state
  }

  .calculate_rates <- function(shift_count, total_time, speciation_events) {
    list(
      shift_rate_per_time = ifelse(
        total_time > 0,
        shift_count / total_time,
        NA_real_
      ),
      shift_rate_per_speciation = ifelse(
        speciation_events > 0,
        shift_count / speciation_events,
        NA_real_
      )
    )
  }

  .calculate_decay_weights <- function(ages) {
    if (!is.null(half_life)) {
      weights <- 1 / (decay_base ^ (ages / half_life))
    } else {
      weights <- 1 / (decay_base ^ ages)
    }
    if (
      length(weights) > 0L &&
        (any(!is.finite(weights)) || any(weights <= 0))
    ) {
      stop(
        "Age-decay weights experienced numerical underflow. ",
        "Increase `half_life` or reduce `decay_base`.",
        call. = FALSE
      )
    }
    if (isTRUE(normalize_weights) && length(weights) > 0L) {
      weight_sum <- sum(weights)
      if (!is.finite(weight_sum) || weight_sum <= 0) {
        stop(
          "Age-decay weights experienced numerical underflow. ",
          "Increase `half_life` or reduce `decay_base`.",
          call. = FALSE
        )
      }
      weights <- weights / weight_sum
    }
    weights
  }

  .calc_weighted_lineage_value <- function(node_states, node_ages) {
    values <- state_values[as.character(node_states)]
    weights <- .calculate_decay_weights(node_ages)
    sum(values * weights, na.rm = TRUE)
  }

  .process_tip <- function(k) {
    lineage_nodes <- c(k, .ancestor_path(k))
    shift_count <- 0L
    total_time <- unname(root_to_tip_distances[k])
    reference_height <- if (identical(age_reference, "tip")) total_time else tree_height
    speciation_events <- 0L

    for (i in seq_along(lineage_nodes)) {
      node <- lineage_nodes[i]
      if (i == 1L) next

      if (node != k && node != root_node) {
        speciation_events <- speciation_events + 1L
      }

      if (node != root_node) {
        parent_node <- unname(parent_by_child[as.character(node)])
        if (is.na(parent_node)) next # nocov

        shift_detected <- .detect_shift(
          parent_node = parent_node,
          child_node = node
        )
        if (isTRUE(shift_detected)) {
          shift_count <- shift_count + 1L
        }
      }
    }

    weighted_lineage_value <- NA_real_
    if (identical(mode, "branches")) {
      path_nodes <- .ancestor_path(k)
      rate_states <- vapply(path_nodes, .state_for_node, FUN.VALUE = character(1))
      rate_ages <- reference_height - unname(node_depths[path_nodes])
      weighted_lineage_value <- .calc_weighted_lineage_value(rate_states, rate_ages)
    }
    rate_summaries <- .calculate_rates(
      shift_count,
      total_time,
      speciation_events
    )

    c(
      shift_count = shift_count,
      total_time = total_time,
      shift_rate_per_time = rate_summaries$shift_rate_per_time,
      shift_rate_per_speciation = rate_summaries$shift_rate_per_speciation,
      weighted_lineage_value = weighted_lineage_value
    )
  }

  results <- .lineage_rates_tip_lapply(
    n_tip = n_tip,
    mode = mode,
    cores = cores,
    progress = progress,
    FUN = .process_tip
  )
  mat <- do.call(rbind, results)
  tip_state_value <- state_values[as.character(tip_states)]

  results_df <- data.frame(
    Tip = mapped_tree$tip.label,
    Shift_Count = as.numeric(mat[, "shift_count"]),
    Total_Time = as.numeric(mat[, "total_time"]),
    Shift_Rate_Per_Time = as.numeric(mat[, "shift_rate_per_time"]),
    Shift_Rate_Per_Speciation = as.numeric(mat[, "shift_rate_per_speciation"]),
    stringsAsFactors = FALSE
  )
  if (identical(mode, "branches")) {
    results_df$Weighted_Lineage_Value <- as.numeric(mat[, "weighted_lineage_value"])
  }
  results_df$Tip_State <- as.character(tip_states)
  results_df$Tip_State_Value <- as.numeric(tip_state_value)

  results_df
}
