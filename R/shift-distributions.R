#' Extract Regime Shift Transitions
#'
#' Extract node-state transitions from a mapped SIMMAP tree and annotate them
#' with state-associated numeric values. For `bifrost_search` objects, the
#' values are the fitted BMM regime rates in
#' `x$model_no_uncertainty$param`. `bifrost_search` inputs are accepted only
#' for multi-regime BMM fits, because other model families do not expose the
#' mapped heterogeneous-rate history summarized here. For generic inputs,
#' supply a SIMMAP tree and a named numeric `state_values` vector.
#'
#' @param x A `bifrost_search` object, a compatible list with
#'   `tree_no_uncertainty_untransformed` and `model_no_uncertainty`, a
#'   SIMMAP-style `phylo` tree when `state_values` is supplied, or `NULL` when
#'   using the `tree` argument.
#' @param tree Optional SIMMAP-style `phylo` tree for generic input mode.
#' @param state_values Named numeric vector mapping SIMMAP state labels to
#'   state-associated values for generic input mode.
#' @param include_root Logical; include a synthetic root row describing the
#'   root state and value.
#'
#' @return A `shift_transitions` data frame with one row per detected
#'   parent-child node-state transition, plus the optional root row. Columns
#'   include node identity, node height and age, parent and child states, parent
#'   and child values, rate delta, percentage change, log rate ratio, and the
#'   classified `rate_change`. The result carries `tree`, `state_values`, and
#'   `settings` attributes that record the mapped tree, resolved state-value
#'   vector, and input options used to create the table.
#' @details `shift_transitions()` implements the manuscript's node-based shift
#'   definition. Each positive-length edge map must therefore contain only one
#'   state. SIMMAP histories with state transitions inside an edge are rejected
#'   because assigning such a transition to the child node would give it the
#'   wrong time and branch location.
#' @examples
#' toy_tree <- ape::read.tree(text = "(((a:1,b:1):1,c:2):1,d:3);")
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
#' toy_tree <- phytools::paintSubTree(
#'   tree = toy_tree,
#'   node = ape::Ntip(toy_tree) + 3L,
#'   state = "2"
#' )
#' shift_transitions(
#'   tree = toy_tree,
#'   state_values = c("0" = 1, "1" = 4, "2" = 2)
#' )
#' @export
shift_transitions <- function(
  x = NULL,
  tree = NULL,
  state_values = NULL,
  include_root = TRUE
) {
  .shift_check_logical(include_root, "include_root")
  resolved <- .shift_resolve_inputs(
    x = x,
    tree = tree,
    state_values = state_values
  )

  mapped_tree <- resolved$tree
  values <- resolved$state_values
  .shift_validate_node_painted_tree(mapped_tree)
  states <- list(
    tips = phytools::getStates(mapped_tree, type = "tips"),
    nodes = phytools::getStates(mapped_tree, type = "nodes")
  )
  names(states$tips) <- as.character(names(states$tips))
  names(states$nodes) <- as.character(names(states$nodes))
  n_tip <- ape::Ntip(mapped_tree)
  root_node <- .shift_root_node(mapped_tree)
  node_heights <- .shift_node_heights(mapped_tree)
  tree_height <- max(node_heights, na.rm = TRUE)

  rows <- list()

  if (isTRUE(include_root)) {
    root_state <- .shift_state_at_node(mapped_tree, states, root_node)
    root_rate <- unname(values[as.character(root_state)])
    rows[[length(rows) + 1L]] <- data.frame(
      node = root_node,
      node_label = as.character(root_node),
      parent_node = NA_integer_,
      parent_state = NA_character_,
      parent_rate = NA_real_,
      child_node = root_node,
      child_state = as.character(root_state),
      child_rate = as.numeric(root_rate),
      height = as.numeric(node_heights[root_node]),
      age = as.numeric(tree_height - node_heights[root_node]),
      edge_length = NA_real_,
      rate_delta = NA_real_,
      percentage_change = NA_real_,
      log_ratio = NA_real_,
      rate_change = "root",
      stringsAsFactors = FALSE
    )
  }

  for (edge_index in seq_len(nrow(mapped_tree$edge))) {
    parent_node <- mapped_tree$edge[edge_index, 1]
    child_node <- mapped_tree$edge[edge_index, 2]
    parent_state <- .shift_state_at_node(mapped_tree, states, parent_node)
    child_state <- .shift_state_at_node(mapped_tree, states, child_node)

    if (is.na(parent_state) ||
        is.na(child_state) ||
        identical(as.character(parent_state), as.character(child_state))) {
      next
    }

    parent_rate <- as.numeric(values[as.character(parent_state)])
    child_rate <- as.numeric(values[as.character(child_state)])
    rate_delta <- child_rate - parent_rate
    percentage_change <- if (is.finite(parent_rate) && parent_rate != 0) {
      100 * rate_delta / parent_rate
    } else {
      NA_real_
    }
    log_ratio <- if (parent_rate > 0 && child_rate > 0) {
      base::log(child_rate / parent_rate)
    } else {
      NA_real_
    }
    node_label <- if (child_node <= n_tip) {
      mapped_tree$tip.label[child_node]
    } else {
      as.character(child_node)
    }

    rows[[length(rows) + 1L]] <- data.frame(
      node = child_node,
      node_label = node_label,
      parent_node = parent_node,
      parent_state = as.character(parent_state),
      parent_rate = parent_rate,
      child_node = child_node,
      child_state = as.character(child_state),
      child_rate = child_rate,
      height = as.numeric(node_heights[child_node]),
      age = as.numeric(tree_height - node_heights[child_node]),
      edge_length = as.numeric(mapped_tree$edge.length[edge_index]),
      rate_delta = rate_delta,
      percentage_change = percentage_change,
      log_ratio = log_ratio,
      rate_change = .shift_rate_change(rate_delta),
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0L) {
    out <- .shift_empty_transitions()
  } else {
    out <- do.call(rbind, rows)
    out <- out[order(out$height, out$node), , drop = FALSE]
    rownames(out) <- NULL
  }

  class(out) <- c("shift_transitions", class(out))
  attr(out, "tree") <- mapped_tree
  attr(out, "state_values") <- values
  attr(out, "settings") <- list(
    tree_name = resolved$tree_name,
    input_mode = resolved$input_mode,
    model_family = resolved$model_family,
    include_root = include_root
  )
  out
}

#' Prepare And Plot Shift Node Marks On A Tree
#'
#' `shift_node_marks()` prepares node markers for inferred regime shifts on an
#' already plotted tree. For `bifrost_search` inputs, transitions are computed
#' from the final multi-regime BMM mapped tree. Generic workflows can supply a
#' `shift_transitions()` table, or a SIMMAP-style tree plus a named
#' `state_values` vector.
#'
#' The workflow is intentionally plot-agnostic: first draw the tree with
#' [plot.rateMap()] or another tree plotting function, then call `plot()` on the
#' returned `shift_node_marks` object to add the node layer.
#'
#' @param x A `bifrost_search` object, a `shift_transitions()` table, a
#'   SIMMAP-style `phylo` tree when `state_values` is supplied, or an object
#'   previously returned by `shift_node_marks()`.
#' @param tree Optional SIMMAP-style `phylo` tree for generic input mode.
#' @param state_values Named numeric vector mapping SIMMAP state labels to
#'   state-associated values for generic input mode.
#' @param transitions Optional precomputed transition table. When supplied,
#'   `x` is used only as a possible source of `ic_weights`.
#' @param ic_weights Optional data frame with `node` and
#'   `ic_weight_withshift` columns. When omitted for a compatible
#'   `bifrost_search` object, `x$ic_weights` is used when available.
#' @param support_threshold Numeric threshold used to flag low-support shifts
#'   from `ic_weight_withshift`.
#' @param rate_changes Which directional shifts to draw.
#' @param letter_order Ordering used when assigning letters to increase nodes.
#'   The default, `"child_state"`, matches the Berv et al. (2026) manuscript
#'   convention by ordering increases by the destination regime state, with node
#'   number as a stable tie-breaker. Use `"node"` for node-number order or
#'   `"age"` for chronological order.
#' @param show_letters Logical; draw letter labels for increase nodes.
#' @param show_low_support Logical; draw an additional dot on shifts with
#'   `ic_weight_withshift < support_threshold`.
#' @param show_legend Logical; draw a compact legend for the node marks.
#' @param marker_base_cex,marker_scale_factor Controls for marker size. Size is
#'   `marker_base_cex + transformed(abs(percentage_change)) *
#'   marker_scale_factor`.
#' @param marker_transform Transformation applied to absolute percent change
#'   before scaling marker size.
#' @param marker_max_cex Maximum marker size.
#' @param marker_alpha Alpha used for filled node markers.
#' @param letter_cex Character expansion for letter labels.
#' @param increase_fill,decrease_fill Fill colors for increase and decrease
#'   markers.
#' @param marker_col Outline color for filled node markers.
#' @param low_support_cex,low_support_col Size and color for low-support dots.
#' @param legend_position Position passed to [graphics::legend()].
#' @param legend_cex Character expansion passed to [graphics::legend()].
#' @param row.names,optional Included for compatibility with
#'   [base::as.data.frame()]; ignored.
#' @param component Component to extract with [base::as.data.frame()]. Use
#'   `"marks"`, `"increase_key"`, or `"summary"`.
#' @param ... Reserved for future extensions. Supplying unused arguments to
#'   `plot()` or `as.data.frame()` methods is an error.
#'
#' @return `shift_node_marks()` returns an object of class `"shift_node_marks"`
#'   with `marks`, `increase_key`, `low_support_nodes`, `summary`, and
#'   `settings`. Plotting returns the same object invisibly. The
#'   `as.data.frame()` method returns the requested component as a plain
#'   `data.frame`.
#' @examples
#' toy_tree <- ape::read.tree(text = "(((a:1,b:1):1,c:2):1,d:3);")
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
#' toy_tree <- phytools::paintSubTree(
#'   tree = toy_tree,
#'   node = ape::Ntip(toy_tree) + 3L,
#'   state = "2"
#' )
#' marks <- shift_node_marks(
#'   tree = toy_tree,
#'   state_values = c("0" = 1, "1" = 4, "2" = 2)
#' )
#' marks$increase_key
#' @name shift_node_marks
NULL

#' @rdname shift_node_marks
#' @export
shift_node_marks <- function(
  x = NULL,
  tree = NULL,
  state_values = NULL,
  transitions = NULL,
  ic_weights = NULL,
  support_threshold = 0.9,
  letter_order = c("child_state", "node", "age"),
  marker_base_cex = 1.5,
  marker_scale_factor = 0.08,
  marker_transform = c("square_root", "log1p"),
  marker_max_cex = Inf
) {
  marker_transform <- match.arg(marker_transform)
  letter_order <- match.arg(letter_order)

  .shift_node_marks_resolve(
    x = x,
    tree = tree,
    state_values = state_values,
    transitions = transitions,
    ic_weights = ic_weights,
    support_threshold = support_threshold,
    letter_order = letter_order,
    marker_base_cex = marker_base_cex,
    marker_scale_factor = marker_scale_factor,
    marker_transform = marker_transform,
    marker_max_cex = marker_max_cex
  )
}

#' @rdname shift_node_marks
#' @method plot shift_node_marks
#' @export
plot.shift_node_marks <- function(
  x,
  rate_changes = c("increase", "decrease"),
  show_letters = TRUE,
  show_low_support = TRUE,
  show_legend = TRUE,
  marker_alpha = 0.75,
  letter_cex = 0.75,
  increase_fill = "white",
  decrease_fill = "#2b6cb0",
  marker_col = "#111827",
  low_support_cex = 0.70,
  low_support_col = grDevices::adjustcolor("black", alpha.f = 0.58),
  legend_position = "topleft",
  legend_cex = 0.72,
  ...
) {
  .distribution_check_unused_dots(list(...))
  rate_changes <- match.arg(
    rate_changes,
    c("increase", "decrease"),
    several.ok = TRUE
  )
  .shift_check_logical(show_letters, "show_letters")
  .shift_check_logical(show_low_support, "show_low_support")
  .shift_check_logical(show_legend, "show_legend")

  annotation <- .shift_node_marks_as_prepared(x)
  .shift_node_marks_draw(
    annotation = annotation,
    rate_changes = rate_changes,
    show_letters = show_letters,
    show_low_support = show_low_support,
    show_legend = show_legend,
    marker_alpha = marker_alpha,
    letter_cex = letter_cex,
    increase_fill = increase_fill,
    decrease_fill = decrease_fill,
    marker_col = marker_col,
    low_support_cex = low_support_cex,
    low_support_col = low_support_col,
    legend_position = legend_position,
    legend_cex = legend_cex
  )
}

#' @rdname shift_node_marks
#' @method as.data.frame shift_node_marks
#' @export
as.data.frame.shift_node_marks <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  component = c("marks", "increase_key", "summary"),
  ...
) {
  .distribution_check_unused_dots(list(...))
  component <- match.arg(component)
  x <- .shift_node_marks_as_prepared(x)
  out <- x[[component]]
  rownames(out) <- NULL
  out
}

#' Define Shift-Magnitude Analysis Groups
#'
#' Wrap named analysis groups so downstream shift-magnitude verbs know to run
#' once per group rather than pooling every input into one result. Each group
#' can be a transition table, BMM `bifrost_search` object, list of runs to pool,
#' or a precomputed shift-magnitude object accepted by the downstream verb. When
#' grouped inputs include precomputed `shift_magnitude_comparison` objects,
#' [compare_shift_magnitudes()] requires their stored settings to match the
#' settings requested for the grouped call.
#'
#' @param ... Named analysis groups. Use quoted names for labels that are not
#'   syntactic R names, such as `"GIC + BIC" = runs`.
#' @param .list Optional named list of analysis groups. Use this when groups
#'   are already stored in a list.
#' @param labels Optional replacement labels, one per analysis group. Defaults
#'   to names from the supplied groups, falling back to `"Analysis 1"`,
#'   `"Analysis 2"`, and so on when names are absent.
#'
#' @return A named list of class `shift_magnitude_groups`.
#' @examples
#' transitions <- data.frame(
#'   rate_change = c("increase", "increase", "decrease", "decrease"),
#'   rate_delta = c(2, 3, -1, -1.5)
#' )
#' groups <- shift_magnitude_groups(
#'   "GIC + BIC" = list(gic = transitions, bic = transitions),
#'   GIC = list(gic = transitions)
#' )
#' compare_shift_magnitudes(groups, ks_reps = 99)
#' @export
shift_magnitude_groups <- function(..., .list = NULL, labels = NULL) {
  dots <- list(...)
  if (!is.null(.list)) {
    if (length(dots) > 0L) {
      stop("Use either `...` or `.list`, not both.", call. = FALSE)
    }
    groups <- .list
  } else if (length(dots) == 1L &&
             is.null(names(dots)) &&
             .shift_magnitude_is_group_container(dots[[1L]])) {
    groups <- dots[[1L]]
  } else {
    groups <- dots
  }

  groups <- .shift_magnitude_group_inputs(groups, labels)
  class(groups) <- c("shift_magnitude_groups", "list")
  groups
}

#' Compare Rate-Shift Magnitudes
#'
#' Compare the magnitude of fitted rate increases and decreases from a
#' `bifrost_search` BMM result or a transition table returned by
#' `shift_transitions()`. This is the package-level analogue of the Berv et al.
#' (2026) workflow that summarized rate changes and tested whether increase
#' and decrease magnitudes differed. When `x` is a `bifrost_search` object, it
#' must be a multi-regime BMM fit; use a transition table or `tree` plus
#' `state_values` for generic mapped-tree inputs.
#'
#' @param x A `bifrost_search` object, compatible list, SIMMAP tree, data
#'   frame with `rate_change` and the selected `measure`, a list of transition
#'   tables or BMM `bifrost_search` objects to pool before testing, a
#'   `shift_magnitude_groups` object to compare each named group separately,
#'   or `NULL` when using the `tree` argument. Fitted `bifrost_search` inputs
#'   are accepted only for multi-regime BMM fits.
#' @param tree Optional SIMMAP-style `phylo` tree for generic input mode.
#' @param state_values Named numeric state-value vector for generic input mode.
#' @param measure Which transition column to compare. `rate_delta` compares
#'   fitted child-minus-parent rate differences, `percentage_change` compares
#'   percent changes relative to the parent rate, and `log_ratio` compares log
#'   child/parent rate ratios.
#' @param transform How to transform the selected measure before testing.
#'   `"absolute"` compares magnitudes, `"log_absolute"` compares log
#'   magnitudes, and `"signed"` compares signed values.
#' @param tests Character vector of tests to run. Supported values are `"t"`
#'   for a Welch t-test, `"wilcox"` for an unpaired Wilcoxon rank-sum test, and
#'   `"ks"` for a two-sample Kolmogorov-Smirnov test. Use `"all"` for all
#'   three. The default reproduces the Berv et al. (2026) rate-delta magnitude
#'   comparison: a KS test on `log(abs(rate_delta))`.
#' @param alternative Alternative hypothesis passed to the tests.
#' @param ks_simulate_p_value Logical; passed to [stats::ks.test()] when
#'   `"ks"` is requested. The default matches the Berv et al. (2026) utility.
#' @param ks_B Integer number of replicates for the simulated KS p-value.
#' @param ks_reps Optional clearer alias for `ks_B`.
#' @param bootstrap_p_value Logical; estimate the mean p-value from stratified
#'   bootstrap resampling within the increase and decrease groups. This
#'   reproduces the bootstrap p-value annotation used in the Berv et al. (2026)
#'   magnitude-density panels when `"ks"` is requested.
#' @param bootstrap_R Integer number of bootstrap replicates when
#'   `bootstrap_p_value = TRUE`.
#' @param bootstrap_reps Optional clearer alias for `bootstrap_R`.
#' @param seed Optional random seed for stochastic KS simulation and
#'   bootstrap p-value resampling. The Berv et al. (2026) defaults are unchanged;
#'   when a seed is supplied, the previous global RNG state is restored before
#'   returning.
#'
#' @return A list of class `shift_magnitude_comparison` with the transformed
#'   values used for testing, group summaries, modal summaries, test results,
#'   and settings. When `x` is a `shift_magnitude_groups` object, a named
#'   `shift_magnitude_comparison_set` is returned instead.
#'
#' @details
#' Plain lists of transition tables or BMM search objects are pooled into one
#' comparison before testing. Use [shift_magnitude_groups()] when list elements
#' are analysis groups that should be compared separately. If a group is already
#' a `shift_magnitude_comparison`, its `measure`, `transform`, `tests`,
#' `alternative`, KS settings, and bootstrap settings must exactly match the
#' grouped call; otherwise, recompute that group from raw transitions or call
#' `compare_shift_magnitudes()` with matching settings.
#' @examples
#' transitions <- data.frame(
#'   rate_change = c("increase", "increase", "decrease", "decrease"),
#'   rate_delta = c(2, 3, -1, -1.5)
#' )
#' compare_shift_magnitudes(transitions, ks_reps = 99)
#' groups <- shift_magnitude_groups(A = transitions, B = transitions)
#' compare_shift_magnitudes(groups, ks_reps = 99)
#' @export
compare_shift_magnitudes <- function(
  x = NULL,
  tree = NULL,
  state_values = NULL,
  measure = c("rate_delta", "percentage_change", "log_ratio"),
  transform = c("log_absolute", "absolute", "signed"),
  tests = "ks",
  alternative = c("two.sided", "less", "greater"),
  ks_simulate_p_value = TRUE,
  ks_B = 10000L,
  bootstrap_p_value = FALSE,
  bootstrap_R = 100L,
  ks_reps = NULL,
  bootstrap_reps = NULL,
  seed = NULL
) {
  measure <- match.arg(measure)
  transform <- match.arg(transform)
  alternative <- match.arg(alternative)
  ks_B <- .shift_magnitude_resolve_rep_alias(ks_reps, ks_B, "ks_reps")
  bootstrap_R <- .shift_magnitude_resolve_rep_alias(
    bootstrap_reps,
    bootstrap_R,
    "bootstrap_reps"
  )
  .shift_magnitude_validate_ks(ks_simulate_p_value, ks_B)
  .shift_magnitude_validate_bootstrap(bootstrap_p_value, bootstrap_R)
  .distribution_check_seed(seed)
  ks_B <- as.integer(ks_B)
  bootstrap_R <- as.integer(bootstrap_R)
  seed <- if (is.null(seed)) NULL else as.integer(seed)
  tests <- .shift_magnitude_normalize_tests(tests)
  settings <- list(
    measure = measure,
    transform = transform,
    tests = tests,
    alternative = alternative,
    ks_simulate_p_value = ks_simulate_p_value,
    ks_B = ks_B,
    bootstrap_p_value = bootstrap_p_value,
    bootstrap_R = bootstrap_R,
    seed = seed
  )
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  if (inherits(x, "shift_magnitude_groups")) {
    if (!is.null(tree) || !is.null(state_values)) {
      stop(
        "`tree` and `state_values` are not used when `x` is a shift_magnitude_groups object.",
        call. = FALSE
      )
    }
    groups <- .shift_magnitude_group_inputs(x)
    comparisons <- Map(function(group_name, group) {
      if (inherits(group, "shift_magnitude_comparison")) {
        mismatches <- .shift_magnitude_comparison_setting_mismatches(group, settings)
        if (length(mismatches) > 0L) {
          stop(
            "Precomputed shift_magnitude_comparison in analysis group `",
            group_name,
            "` was created with settings that do not match this ",
            "`compare_shift_magnitudes()` call. Mismatched setting(s): ",
            paste(mismatches, collapse = ", "),
            ". Recompute the comparison from raw transitions or call ",
            "`compare_shift_magnitudes()` with matching settings.",
            call. = FALSE
          )
        }
        return(group)
      }
      tryCatch(
        compare_shift_magnitudes(
          x = group,
          measure = measure,
          transform = transform,
          tests = tests,
          alternative = alternative,
          ks_simulate_p_value = ks_simulate_p_value,
          ks_B = ks_B,
          bootstrap_p_value = bootstrap_p_value,
          bootstrap_R = bootstrap_R,
          seed = seed
        ),
        error = function(err) {
          stop(
            "Error in analysis group `",
            group_name,
            "`: ",
            conditionMessage(err),
            call. = FALSE
          )
        }
      )
    }, names(groups), groups)

    attr(comparisons, "settings") <- settings
    class(comparisons) <- c("shift_magnitude_comparison_set", "list")
    return(comparisons)
  }

  transitions <- .shift_magnitude_resolve_transitions(
    x = x,
    tree = tree,
    state_values = state_values
  )
  values <- .shift_magnitude_values(
    transitions = transitions,
    measure = measure,
    transform = transform
  )
  summary <- .shift_magnitude_summary(values)
  modal <- .shift_magnitude_modal(values, transform)
  test_results <- .shift_magnitude_tests(
    values = values,
    tests = tests,
    alternative = alternative,
    ks_simulate_p_value = ks_simulate_p_value,
    ks_B = ks_B,
    bootstrap_p_value = bootstrap_p_value,
    bootstrap_R = bootstrap_R
  )

  out <- list(
    values = values,
    summary = summary,
    modal = modal,
    tests = test_results,
    settings = settings
  )
  class(out) <- c("shift_magnitude_comparison", "list")
  out
}

#' Compute Waiting Times Between Regime Shifts
#'
#' Compute chronological whole-tree waiting times, between-shift intervals
#' along root-to-tip lineages, and acceleration-focused increase-to-increase
#' waiting times from a mapped SIMMAP tree or from the output of
#' `shift_transitions()`. Lineage intervals exclude the censored
#' root-to-first-shift and last-shift-to-tip boundaries. Fitted
#' `bifrost_search` inputs are accepted only for multi-regime BMM fits; generic
#' SIMMAP trees remain supported with `tree` and `state_values`.
#' Tied or simultaneous transitions are preserved as zero-length waiting times
#' in the descriptive tables.
#'
#' @param x A BMM `bifrost_search` object, compatible list, SIMMAP tree,
#'   `shift_transitions` data frame, or `NULL` when using `tree`.
#' @param tree Optional SIMMAP-style `phylo` tree for generic input mode, or
#'   the tree corresponding to a user-supplied transition table.
#' @param state_values Named numeric state-value vector for generic input mode.
#' @param scope Which components to compute. `"all"` computes global, lineage,
#'   and acceleration summaries.
#'
#' @return A list of class `shift_waiting_times` with `transitions`, `global`,
#'   `lineage`, `lineage_by_tip`, `acceleration`, and `summaries` components.
#'   Components outside the requested `scope` are returned as `NULL`.
#' @examples
#' toy_tree <- ape::read.tree(text = "(((a:1,b:1):1,c:2):1,d:3);")
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
#' toy_tree <- phytools::paintSubTree(
#'   tree = toy_tree,
#'   node = ape::Ntip(toy_tree) + 3L,
#'   state = "2"
#' )
#' shift_waiting_times(
#'   tree = toy_tree,
#'   state_values = c("0" = 1, "1" = 4, "2" = 2)
#' )
#' @export
shift_waiting_times <- function(
  x = NULL,
  tree = NULL,
  state_values = NULL,
  scope = c("all", "global", "lineage", "acceleration")
) {
  scope <- match.arg(scope)
  resolved <- .shift_waiting_resolve_transitions(
    x = x,
    tree = tree,
    state_values = state_values
  )
  transitions <- resolved$transitions
  mapped_tree <- resolved$tree

  needs_global <- scope %in% c("all", "global")
  needs_lineage <- scope %in% c("all", "lineage")
  needs_acceleration <- scope %in% c("all", "acceleration")
  if ((needs_lineage || needs_acceleration) && is.null(mapped_tree)) {
    stop(
      "`tree` is required for lineage and acceleration waiting times.",
      call. = FALSE
    )
  }

  global <- if (needs_global) .shift_global_waiting_times(transitions) else NULL
  lineage_result <- if (needs_lineage) {
    .shift_lineage_waiting_times(transitions, mapped_tree)
  } else {
    list(combined = NULL, by_tip = NULL)
  }
  acceleration <- if (needs_acceleration) {
    .shift_acceleration_waiting_times(transitions, mapped_tree)
  } else {
    NULL
  }

  summaries <- list(
    global = if (needs_global) .shift_waiting_summary(global) else NULL,
    lineage = if (needs_lineage) .shift_waiting_summary(lineage_result$combined) else NULL,
    acceleration = if (needs_acceleration) acceleration$summary else NULL
  )

  out <- list(
    transitions = transitions,
    global = global,
    lineage = lineage_result$combined,
    lineage_by_tip = lineage_result$by_tip,
    acceleration = acceleration,
    summaries = summaries,
    scope = scope
  )
  class(out) <- c("shift_waiting_times", "list")
  out
}

#' Fit Candidate Distributions To Rates
#'
#' Fit and rank candidate continuous distributions for a numeric vector of
#' rates, a lineage-rate table, or the fitted regime rates in a `bifrost_search`
#' object. `bifrost_search` inputs must be multi-regime BMM fits; other model
#' families are rejected before fitting. The function is generic with respect
#' to the candidate distribution: no distribution, including Gumbel, receives
#' special handling unless requested through `models`.
#'
#' @param x Numeric vector, lineage-rate table, or `bifrost_search` object.
#' @param log Logical; log-transform positive rate values before fitting.
#'   Lineage-rate tables with a `log_lineage_rate` column are treated as already
#'   log transformed when `log = TRUE`. When `log = FALSE`, lineage-rate
#'   tables must include a raw rate column rather than only
#'   `log_lineage_rate`.
#' @param models Character vector of `univariateML` model identifiers. When
#'   `NULL`, a broad real-valued candidate set is used.
#' @param select_by Criterion used to select `selected_model`.
#'
#' @return A list of class `rate_distribution_fit` with the transformed data,
#'   model rankings, a long-form `parameters` table of fitted parameter
#'   estimates, individual fitted objects, selected model, and failed model
#'   attempts. The raw `univariateML` objects remain available in `fits` for
#'   workflows that need package-specific methods.
#' @examples
#' if (requireNamespace("univariateML", quietly = TRUE)) {
#'   rates <- c(0.8, 1.1, 1.4, 1.9, 2.5, 3.1)
#'   fit_rate_distribution(rates, models = c("norm", "gumbel"))
#' }
#' @export
fit_rate_distribution <- function(
  x,
  log = TRUE,
  models = NULL,
  select_by = c("AIC", "BIC")
) {
  .shift_check_logical(log, "log")
  select_by <- match.arg(select_by)
  if (is.null(models)) {
    models <- .rate_distribution_default_models()
  }

  data_info <- .rate_distribution_data(x, log = log)
  fit <- .distribution_fit_models(
    values = data_info$values,
    models = models,
    select_by = select_by
  )

  out <- list(
    data = data_info$values,
    original_data = data_info$original_values,
    source = data_info$source,
    source_column = data_info$source_column,
    transformation = data_info$transformation,
    log = log,
    models = fit$models,
    rankings = fit$rankings,
    parameters = fit$parameters,
    fits = fit$fits,
    bootstrap = list(),
    selected = fit$selected,
    selected_model = fit$selected_model,
    select_by = select_by,
    failed = fit$failed
  )
  class(out) <- c("rate_distribution_fit", "list")
  out
}

#' Bootstrap A Fitted Rate Distribution
#'
#' Precompute parametric bootstrap parameter draws for a fitted rate
#' distribution. This keeps stochastic uncertainty estimation out of the plotting
#' method while preserving the bootstrap draws on the returned
#' `rate_distribution_fit` object.
#'
#' @details
#' For Gumbel fits, the default `kl = NULL` also computes the
#' Kullback-Leibler divergence between the empirical kernel density of the
#' fitted data and the fitted Gumbel density, then repeats that calculation for
#' each bootstrap parameter draw. Extract the compact diagnostic table with
#' `as.data.frame(x, component = "goodness")`.
#'
#' @param x A `rate_distribution_fit` object from [fit_rate_distribution()].
#' @param model Distribution model to bootstrap. If `NULL`, the selected model
#'   stored in `x$selected_model` is used.
#' @param reps Number of bootstrap parameter draws.
#' @param seed Optional random seed for the bootstrap draw. The previous global
#'   RNG state is restored after the bootstrap when a seed is supplied.
#' @param kl `NULL` or logical. When `NULL`, compute the bootstrap
#'   Kullback-Leibler divergence summary for Gumbel fits and skip it for other
#'   models. Set to `FALSE` to store only the bootstrap parameter draws.
#'
#' @return A `rate_distribution_fit` object with a
#'   `rate_distribution_bootstrap` entry stored in `x$bootstrap[[model]]`.
#'   The entry contains the model name, number of draws, seed, and parameter
#'   matrix returned by [univariateML::bootstrapml()]. For Gumbel fits, the
#'   entry also stores a bootstrapped Kullback-Leibler divergence summary unless
#'   `kl = FALSE`.
#' @examples
#' if (requireNamespace("univariateML", quietly = TRUE) &&
#'     requireNamespace("evd", quietly = TRUE)) {
#'   fit <- fit_rate_distribution(c(1, 1.4, 1.9, 2.8, 4.1, 6.5), models = "gumbel")
#'   fit <- bootstrap_rate_distribution(fit, model = "gumbel", reps = 25, seed = 1)
#' }
#' @export
bootstrap_rate_distribution <- function(
  x,
  model = NULL,
  reps = 1000L,
  seed = NULL,
  kl = NULL
) {
  model <- .distribution_rate_fit_model(
    x = x,
    model = model,
    caller = "bootstrap_rate_distribution()"
  )
  .distribution_check_count(reps, "reps", minimum = 1L)
  .distribution_check_seed(seed)
  compute_kl <- .distribution_bootstrap_compute_kl(model, kl)
  if (!requireNamespace("univariateML", quietly = TRUE)) {
    stop("Package `univariateML` is required to bootstrap rate distributions.", call. = FALSE)
  }
  if (compute_kl && !requireNamespace("evd", quietly = TRUE)) {
    stop(
      "Package `evd` is required to compute the Gumbel Kullback-Leibler summary. ",
      "Install it or set `kl = FALSE`.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  boot_parameters <- univariateML::bootstrapml(
    x$fits[[model]],
    reps = as.integer(reps),
    map = identity,
    reducer = identity
  )
  boot_parameters <- as.matrix(boot_parameters)

  if (is.null(rownames(boot_parameters)) || any(!nzchar(rownames(boot_parameters)))) {
    parameter_names <- x$parameters$parameter[x$parameters$model == model]
    if (length(parameter_names) == nrow(boot_parameters)) {
      rownames(boot_parameters) <- parameter_names
    }
  }

  if (is.null(x$bootstrap) || !is.list(x$bootstrap)) {
    x$bootstrap <- list()
  }

  kl_summary <- if (compute_kl) {
    .distribution_gumbel_kl(
      values = x$data,
      fit = x$fits[[model]],
      boot_parameters = boot_parameters
    )
  } else {
    NULL
  }

  x$bootstrap[[model]] <- structure(
    list(
      model = model,
      reps = ncol(boot_parameters),
      seed = if (is.null(seed)) NA_integer_ else as.integer(seed),
      parameters = boot_parameters,
      kl = kl_summary,
      call = match.call()
    ),
    class = c("rate_distribution_bootstrap", "list")
  )
  x
}

#' Fit Candidate Distributions To Shift Waiting Times
#'
#' Fit and rank candidate positive continuous distributions for shift waiting
#' times. Inputs can be numeric waiting times, data frames with a `TimeToNext`
#' column, BMM `bifrost_search` output, or the output of
#' `shift_waiting_times()`. Fitted `bifrost_search` inputs are accepted only
#' for multi-regime BMM fits; numeric vectors and data frames remain generic.
#' Because the candidate distributions have strictly positive support, inputs
#' containing zero waits from tied or simultaneous events are rejected with an
#' explanatory error; the zero waits remain available from
#' [shift_waiting_times()] for descriptive analysis.
#'
#' @param waiting_times Numeric vector, waiting-time data frame,
#'   `shift_waiting_times` object, or BMM `bifrost_search` output.
#' @param by_lineage Logical; when `TRUE`, use pooled between-shift lineage
#'   intervals for the top-level fit instead of chronological whole-tree
#'   intervals, and also fit candidate distributions separately for each
#'   lineage with at least two waiting times.
#' @param global_include_root_wait Logical; when fitting whole-tree global waits
#'   from a `shift_waiting_times` object, `bifrost_search` object, or global
#'   waiting-time data frame, include the initial root-to-first-shift interval.
#'   The default `TRUE` matches the Berv et al. (2026) reported whole-tree
#'   Weibull summary. Set to `FALSE` to exclude the row labeled
#'   `PreviousShift == "root"`. Multiple such rows are ambiguous and produce an
#'   error; data without that label are unchanged.
#' @param models Character vector of positive-support `univariateML` model
#'   identifiers.
#' @param select_by Criterion used to select best models.
#'
#' @return A list of class `waiting_time_distribution_fit` with pooled model
#'   rankings and a long-form `parameters` table of fitted parameter estimates.
#'   With `by_lineage = FALSE`, the top-level pool contains chronological
#'   whole-tree intervals; with `by_lineage = TRUE`, it contains pooled
#'   between-shift lineage intervals and the result additionally includes
#'   per-lineage best-model, parameter, and exponential-rate summaries. The raw
#'   `univariateML` objects remain available in `fits` and `lineage_fits`.
#' @examples
#' if (requireNamespace("univariateML", quietly = TRUE)) {
#'   waits <- c(0.5, 0.8, 1.2, 1.9, 2.4, 3.0)
#'   fit_waiting_time_distribution(waits, models = c("exp", "gamma", "weibull"))
#' }
#' @export
fit_waiting_time_distribution <- function(
  waiting_times,
  by_lineage = FALSE,
  global_include_root_wait = TRUE,
  models = c("exp", "gamma", "weibull", "lnorm", "invgauss"),
  select_by = c("AIC", "BIC")
) {
  .shift_check_logical(by_lineage, "by_lineage")
  .shift_check_logical(global_include_root_wait, "global_include_root_wait")
  select_by <- match.arg(select_by)
  data_info <- .waiting_time_data(
    waiting_times,
    prefer_lineage = by_lineage,
    global_include_root_wait = global_include_root_wait
  )
  fit <- .distribution_fit_models(
    values = data_info$values,
    models = models,
    select_by = select_by
  )

  lineage <- NULL
  lineage_fits <- NULL
  lambda_summary <- NULL
  if (isTRUE(by_lineage)) {
    lineage <- .waiting_time_fit_by_lineage(
      data = data_info$data_frame,
      models = models,
      select_by = select_by
    )
    lineage_fits <- lineage$fits
    lambda_summary <- .waiting_time_lambda_summary(lineage$summary$exp_lambda)
  }

  out <- list(
    data = data_info$values,
    source = data_info$source,
    models = fit$models,
    rankings = fit$rankings,
    parameters = fit$parameters,
    fits = fit$fits,
    selected = fit$selected,
    selected_model = fit$selected_model,
    select_by = select_by,
    failed = fit$failed,
    by_lineage = by_lineage,
    global_include_root_wait = global_include_root_wait,
    lineage = if (isTRUE(by_lineage)) lineage$summary else NULL,
    lineage_parameters = if (isTRUE(by_lineage)) lineage$parameters else NULL,
    lineage_fits = lineage_fits,
    lambda_summary = lambda_summary
  )
  class(out) <- c("waiting_time_distribution_fit", "list")
  out
}

#' Coerce Shift And Distribution Summaries To Data Frames
#'
#' Extract concise, vignette-friendly tables from fitted shift-magnitude and
#' distribution objects. These methods do not change the underlying objects;
#' they only select and lightly round commonly inspected columns.
#'
#' @param x A `shift_magnitude_comparison`,
#'   `shift_magnitude_comparison_set`, `shift_magnitude_counts`,
#'   `shift_magnitude_count_set`, `rate_distribution_fit`, or
#'   `waiting_time_distribution_fit` object.
#' @param row.names,optional Included for compatibility with
#'   [base::as.data.frame()]; ignored.
#' @param component Summary table to extract. Ignored for shift-magnitude count
#'   objects. For shift-magnitude comparisons and
#'   comparison sets, use `"summary"`, `"tests"`, `"modal"`, or `"ks"`. For
#'   distribution fits, use `"rankings"` or `"parameters"`;
#'   rate-distribution fits with bootstrapped Gumbel diagnostics also support
#'   `"goodness"`, which currently reports the Gumbel Kullback-Leibler
#'   divergence summary. For waiting-time fits with `by_lineage = TRUE`, use
#'   `"lineage"` or `"lineage_parameters"`.
#' @param analysis Optional label used when `component = "ks"` for
#'   `shift_magnitude_comparison` objects.
#' @param models Optional model identifiers used to filter parameter tables.
#' @param n Maximum number of rows to return.
#' @param digits Number of significant digits used for numeric columns.
#' @param ... Reserved for future extensions. Supplying unused arguments is an
#'   error.
#'
#' @return A `data.frame` whose columns depend on `component` and the input class.
#'   Shift-magnitude summary tables report directional sample sizes and moments;
#'   test tables report the requested test statistics and p-values; modal tables
#'   report KDE modes and ratios; and `component = "ks"` returns one KS-focused row
#'   per comparison. Distribution-fit ranking tables report model information
#'   criteria, parameter tables report model parameters and estimates, Gumbel
#'   goodness tables report bootstrapped KL-divergence summaries, and
#'   waiting-time lineage tables report per-lineage fit summaries.
#' @name distribution-fit-data-frames
NULL

#' @rdname distribution-fit-data-frames
#' @method as.data.frame shift_magnitude_comparison
#' @export
as.data.frame.shift_magnitude_comparison <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  component = c("summary", "tests", "modal", "ks"),
  analysis = NULL,
  digits = 4L,
  ...
) {
  .distribution_check_unused_dots(list(...))
  component <- match.arg(component)
  if (identical(component, "summary")) {
    out <- x$summary[, c(
      "rate_change",
      "n",
      "mean",
      "median",
      "mean_abs_raw",
      "median_abs_raw"
    ), drop = FALSE]
  } else if (identical(component, "tests")) {
    out <- x$tests[, c(
      "test",
      "statistic",
      "p_value",
      "mean_difference",
      "median_difference",
      "method"
    ), drop = FALSE]
  } else if (identical(component, "modal")) {
    out <- x$modal
  } else {
    ks_row <- x$tests[x$tests$test == "ks", , drop = FALSE]
    if (nrow(ks_row) == 0L) {
      stop("`x` does not contain a KS test result.", call. = FALSE)
    }
    out <- data.frame(
      analysis = if (is.null(analysis)) NA_character_ else as.character(analysis),
      increases = x$summary$n[x$summary$rate_change == "increase"],
      decreases = x$summary$n[x$summary$rate_change == "decrease"],
      ks_D = ks_row$statistic,
      ks_p = ks_row$p_value,
      bootstrap_mean_p = ks_row$bootstrap_mean_p,
      mode_increase = x$modal$mode_increase,
      mode_decrease = x$modal$mode_decrease,
      linear_ratio = x$modal$linear_ratio,
      stringsAsFactors = FALSE
    )
  }

  .distribution_round_numeric(out, digits = digits)
}

#' @rdname distribution-fit-data-frames
#' @method as.data.frame shift_magnitude_comparison_set
#' @export
as.data.frame.shift_magnitude_comparison_set <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  component = c("ks", "summary", "tests", "modal"),
  digits = 4L,
  ...
) {
  .distribution_check_unused_dots(list(...))
  component <- match.arg(component)
  rows <- Map(function(analysis, comparison) {
    if (identical(component, "ks")) {
      return(as.data.frame(
        comparison,
        component = "ks",
        analysis = analysis,
        digits = digits
      ))
    }
    data <- as.data.frame(comparison, component = component, digits = digits)
    data.frame(
      analysis = analysis,
      data,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }, names(x), unclass(x))

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @rdname distribution-fit-data-frames
#' @method as.data.frame shift_magnitude_counts
#' @export
as.data.frame.shift_magnitude_counts <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  digits = 4L,
  ...
) {
  .distribution_check_unused_dots(list(...))
  out <- x
  class(out) <- "data.frame"
  attr(out, "analysis") <- NULL
  rownames(out) <- NULL
  .distribution_round_numeric(out, digits = digits)
}

#' @rdname distribution-fit-data-frames
#' @method as.data.frame shift_magnitude_count_set
#' @export
as.data.frame.shift_magnitude_count_set <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  digits = 4L,
  ...
) {
  .distribution_check_unused_dots(list(...))
  rows <- Map(function(analysis, counts) {
    data <- as.data.frame(counts, digits = digits)
    data.frame(
      analysis = analysis,
      data,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }, names(x), unclass(x))

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @rdname distribution-fit-data-frames
#' @method as.data.frame rate_distribution_fit
#' @export
as.data.frame.rate_distribution_fit <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  component = c("rankings", "parameters", "goodness"),
  models = NULL,
  n = Inf,
  digits = 4L,
  ...
) {
  .distribution_check_unused_dots(list(...))
  component <- match.arg(component)
  out <- .distribution_fit_table(
    x = x,
    component = component,
    models = models,
    n = n
  )
  .distribution_round_numeric(out, digits = digits)
}

#' @rdname distribution-fit-data-frames
#' @method as.data.frame waiting_time_distribution_fit
#' @export
as.data.frame.waiting_time_distribution_fit <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  component = c("rankings", "parameters", "lineage", "lineage_parameters"),
  models = NULL,
  n = Inf,
  digits = 4L,
  ...
) {
  .distribution_check_unused_dots(list(...))
  component <- match.arg(component)
  out <- .distribution_fit_table(
    x = x,
    component = component,
    models = models,
    n = n
  )
  .distribution_round_numeric(out, digits = digits)
}

#' Plot A Fitted Rate Distribution
#'
#' Plot a fitted candidate distribution from [fit_rate_distribution()]. The
#' current method supports the Gumbel candidate used by Berv et al. (2026). It
#' can overlay a pointwise bootstrap band and bootstrap density curves when
#' [bootstrap_rate_distribution()] has already stored bootstrap draws on `x`.
#'
#' @param x A `rate_distribution_fit` object.
#' @param model Distribution model to plot. Currently only `"gumbel"` is
#'   supported.
#' @param bootstrap `NULL` or logical. If `NULL`, draw the bootstrap overlay
#'   when precomputed draws for `model` are available. If `TRUE`, require
#'   precomputed draws from [bootstrap_rate_distribution()]. If `FALSE`, draw no
#'   bootstrap overlay.
#' @param bootstrap_curves Number of precomputed bootstrap density curves to
#'   draw over the band. Set to `0` to draw only the band.
#' @param col Color for the fitted density, bootstrap band, and bootstrap
#'   curves.
#' @param hist_col Histogram fill color. A vector of colors is passed through
#'   to [graphics::hist()] and can color bins individually.
#' @param band_alpha,curve_alpha Alpha values for the bootstrap band and
#'   individual bootstrap curves.
#' @param data_density Logical; overlay the empirical kernel density.
#' @param rug Logical; add a rug for the fitted data values.
#' @param breaks Histogram breaks passed to [graphics::hist()].
#' @param main,xlab,ylab Plot labels.
#' @param ... Additional arguments passed to [graphics::hist()].
#'
#' @return Invisibly returns a list with the x grid, fitted density,
#'   bootstrap densities, bootstrap band, and plotted parameter table.
#' @examples
#' if (requireNamespace("univariateML", quietly = TRUE) &&
#'     requireNamespace("evd", quietly = TRUE)) {
#'   fit <- fit_rate_distribution(c(1, 1.4, 1.9, 2.8, 4.1, 6.5), models = "gumbel")
#'   fit <- bootstrap_rate_distribution(fit, model = "gumbel", reps = 25, seed = 1)
#'   plot(fit, bootstrap_curves = 5)
#' }
#' @method plot rate_distribution_fit
#' @export
plot.rate_distribution_fit <- function(
  x,
  model = c("gumbel"),
  bootstrap = NULL,
  bootstrap_curves = 120L,
  col = "#b01f2e",
  hist_col = "grey90",
  band_alpha = 0.14,
  curve_alpha = 0.08,
  data_density = TRUE,
  rug = TRUE,
  breaks = "FD",
  main = NULL,
  xlab = NULL,
  ylab = "Density",
  ...
) {
  model <- match.arg(model)
  .distribution_plot_check_rate_fit(x, model)
  if (!is.null(bootstrap) &&
      (!is.logical(bootstrap) || length(bootstrap) != 1L || is.na(bootstrap))) {
    stop("`bootstrap` must be NULL, TRUE, or FALSE.", call. = FALSE)
  }
  if (!requireNamespace("evd", quietly = TRUE)) {
    stop("Package `evd` is required to plot the Gumbel density.", call. = FALSE)
  }

  .distribution_check_count(bootstrap_curves, "bootstrap_curves", minimum = 0L)
  bootstrap_result <- .distribution_bootstrap_result(x, model)
  draw_bootstrap <- if (is.null(bootstrap)) {
    !is.null(bootstrap_result)
  } else {
    isTRUE(bootstrap)
  }
  if (isTRUE(draw_bootstrap) && is.null(bootstrap_result)) {
    stop(
      "No precomputed bootstrap draws found for model `", model, "`. ",
      "Call `bootstrap_rate_distribution(x, model = \"", model, "\")` before plotting, ",
      "or set `bootstrap = FALSE`.",
      call. = FALSE
    )
  }

  values <- x$data
  fitted_model <- x$fits[[model]]
  mu_hat <- as.numeric(fitted_model["mu"])
  sigma_hat <- as.numeric(fitted_model["sigma"])
  empirical_density <- stats::density(values)
  x_values <- seq(
    min(c(empirical_density$x, mu_hat - 4 * sigma_hat)),
    max(c(empirical_density$x, mu_hat + 4 * sigma_hat)),
    length.out = 300L
  )
  fitted_density <- evd::dgumbel(x_values, loc = mu_hat, scale = sigma_hat)

  boot_density <- NULL
  band <- NULL
  if (isTRUE(draw_bootstrap)) {
    boot_density <- .distribution_bootstrap_density(
      bootstrap_result = bootstrap_result,
      model = model,
      x_values = x_values
    )
    band <- apply(
      boot_density,
      1L,
      stats::quantile,
      probs = c(0.025, 0.975),
      na.rm = TRUE
    )
  }

  y_values <- c(empirical_density$y, fitted_density, if (!is.null(band)) band)
  if (is.null(main)) {
    main <- "Gumbel fit to fitted regime rates"
  }
  if (is.null(xlab)) {
    xlab <- .distribution_plot_xlab(x)
  }

  graphics::hist(
    values,
    breaks = breaks,
    freq = FALSE,
    border = "white",
    col = hist_col,
    xlab = xlab,
    ylab = ylab,
    main = main,
    ylim = range(y_values, finite = TRUE),
    las = 1,
    ...
  )

  if (!is.null(band)) {
    graphics::polygon(
      c(x_values, rev(x_values)),
      c(band[1L, ], rev(band[2L, ])),
      col = grDevices::adjustcolor(col, alpha.f = band_alpha),
      border = NA
    )
  }

  if (!is.null(boot_density) && bootstrap_curves > 0L) {
    curve_count <- min(as.integer(bootstrap_curves), ncol(boot_density))
    for (i in seq_len(curve_count)) {
      graphics::lines(
        x_values,
        boot_density[, i],
        col = grDevices::adjustcolor(col, alpha.f = curve_alpha),
        lwd = 0.6
      )
    }
  }

  if (isTRUE(data_density)) {
    graphics::lines(empirical_density, col = "grey35", lwd = 1.4, lty = 3)
  }
  graphics::lines(x_values, fitted_density, col = col, lwd = 2.2)
  if (isTRUE(rug)) {
    graphics::rug(values, col = grDevices::adjustcolor("grey20", alpha.f = 0.45))
  }

  legend_labels <- c(
    sprintf("Gumbel: mu = %.2f, sigma = %.2f", mu_hat, sigma_hat),
    if (isTRUE(data_density)) "Empirical density" else NULL,
    if (!is.null(band)) "95% bootstrap band" else NULL
  )
  graphics::legend(
    "topright",
    bty = "n",
    lwd = c(2.2, if (isTRUE(data_density)) 1.4 else NULL, if (!is.null(band)) NA else NULL),
    lty = c(1, if (isTRUE(data_density)) 3 else NULL, if (!is.null(band)) NA else NULL),
    pch = c(NA, if (isTRUE(data_density)) NA else NULL, if (!is.null(band)) 15 else NULL),
    pt.cex = c(NA, if (isTRUE(data_density)) NA else NULL, if (!is.null(band)) 1.4 else NULL),
    col = c(
      col,
      if (isTRUE(data_density)) "grey35" else NULL,
      if (!is.null(band)) grDevices::adjustcolor(col, alpha.f = 0.2) else NULL
    ),
    legend = legend_labels
  )

  invisible(list(
    x = x_values,
    fitted_density = fitted_density,
    bootstrap_density = boot_density,
    band = band,
    parameters = x$parameters[x$parameters$model == model, , drop = FALSE]
  ))
}

#' Plot One Shift-Magnitude Comparison
#'
#' Draw one density panel from a [compare_shift_magnitudes()] result. By
#' default, the method reproduces the Berv et al. (2026) Supplementary Figure 5
#' scaling: each density integrates to that direction's observed frequency, so
#' panel height reflects both shift magnitude and relative abundance. Set
#' `scale_by_frequency = FALSE` for conventional unit-area densities.
#'
#' @param x A single `shift_magnitude_comparison` object.
#' @param scale_by_frequency Logical; multiply each density by its observed
#'   increase/decrease frequency.
#' @param density_args Optional named list of additional arguments passed to
#'   [stats::density()].
#' @param colors Named colors for `"increase"` and `"decrease"`.
#' @param fill_alpha Alpha value for density fills.
#' @param lwd Line width for density curves.
#' @param xlim,ylim Optional axis limits.
#' @param xlab,ylab,main Plot labels.
#' @param annotate Logical; draw panel annotations. This is a master switch for
#'   `show_sample_size`, `show_test_label`, and `show_modal_ratio`.
#' @param show_modes Logical; draw dashed vertical lines at the KDE modes for
#'   increase and decrease magnitudes.
#' @param show_sample_size Logical; include increase/decrease sample sizes in
#'   panel annotations when `annotate = TRUE`.
#' @param show_test_label Logical; include KS statistics and p-values in panel
#'   annotations when `annotate = TRUE` and a KS test is available.
#' @param show_modal_ratio Logical; include the back-transformed modal ratio in
#'   panel annotations when `annotate = TRUE` and it is finite.
#' @param legend Logical; draw an increase/decrease legend.
#' @param rug Logical; add rugs for the transformed values.
#' @param mar,oma Graphical margins passed to [graphics::par()].
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns a list with the comparison object, summary tables,
#'   density estimates, and plot settings.
#' @examples
#' transitions <- data.frame(
#'   rate_change = c("increase", "increase", "increase", "decrease", "decrease", "decrease"),
#'   rate_delta = c(3, 4, 5, -1, -1.5, -2)
#' )
#' comparison <- compare_shift_magnitudes(transitions, ks_reps = 99)
#' plot(comparison)
#' @method plot shift_magnitude_comparison
#' @export
plot.shift_magnitude_comparison <- function(
  x,
  scale_by_frequency = TRUE,
  density_args = list(),
  colors = c(increase = "#377eb8", decrease = "#e41a1c"),
  fill_alpha = 0.22,
  lwd = 2,
  xlim = NULL,
  ylim = NULL,
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  annotate = TRUE,
  show_modes = TRUE,
  show_sample_size = TRUE,
  show_test_label = TRUE,
  show_modal_ratio = TRUE,
  legend = TRUE,
  rug = FALSE,
  mar = c(4.5, 4.5, 2.4, 0.8),
  oma = c(0, 0, 0, 0),
  ...
) {
  if (!inherits(x, "shift_magnitude_comparison")) {
    stop("`x` must be a shift_magnitude_comparison object.", call. = FALSE)
  }
  .shift_check_logical(scale_by_frequency, "scale_by_frequency")
  .shift_check_logical(annotate, "annotate")
  .shift_check_logical(show_modes, "show_modes")
  .shift_check_logical(show_sample_size, "show_sample_size")
  .shift_check_logical(show_test_label, "show_test_label")
  .shift_check_logical(show_modal_ratio, "show_modal_ratio")
  .shift_check_logical(legend, "legend")
  .shift_check_logical(rug, "rug")
  .shift_magnitude_plot_check_alpha(fill_alpha, "fill_alpha")
  .shift_magnitude_plot_check_density_args(density_args)
  colors <- .shift_magnitude_plot_colors(colors)

  prepared <- .shift_magnitude_distribution_prepare(
    comparison = x,
    scale_by_frequency = scale_by_frequency,
    density_args = density_args,
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    main = main,
    annotate = annotate,
    show_modes = show_modes,
    show_sample_size = show_sample_size,
    show_test_label = show_test_label,
    show_modal_ratio = show_modal_ratio
  )

  .shift_magnitude_distribution_draw(
    prepared,
    colors = colors,
    fill_alpha = fill_alpha,
    lwd = lwd,
    legend = legend,
    rug = rug,
    mar = mar,
    oma = oma,
    ...
  )
}

#' Plot A Shift-Magnitude Comparison Set
#'
#' Comparison sets contain multiple named comparisons. Select one named
#' comparison, for example `x[["GIC"]]`, before plotting.
#'
#' @param x A `shift_magnitude_comparison_set` object.
#' @param ... Ignored.
#'
#' @return Stops with an informative message.
#' @method plot shift_magnitude_comparison_set
#' @export
plot.shift_magnitude_comparison_set <- function(x, ...) {
  stop(
    "`x` contains multiple magnitude comparisons; select one with `x[[\"name\"]]` before plotting.",
    call. = FALSE
  )
}

#' Count Inferred Rate Increases And Decreases
#'
#' Count inferred rate increases and decreases for a transition table,
#' `bifrost_search` object, `shift_magnitude_comparison` object, or list of
#' runs. This is the tabular input to `plot()` for shift-magnitude counts and
#' reproduces the count/frequency summaries used by Berv et al. (2026) before
#' plotting. Fitted `bifrost_search` inputs are accepted only for multi-regime
#' BMM fits; generic mapped-tree inputs should use `tree` plus `state_values`.
#'
#' @param x A transition table, `bifrost_search` object, list of
#'   transition/search objects, `shift_magnitude_comparison` object,
#'   `shift_magnitude_groups` object, or `NULL` when using the `tree` argument.
#' @param tree Optional SIMMAP-style `phylo` tree for generic input mode.
#' @param state_values Named numeric state-value vector for generic input mode.
#'
#' @return A `shift_magnitude_counts` data frame with one row per input run.
#'   Columns include `source`, increase/decrease counts, total directional
#'   shifts, and increase/decrease frequencies. When `x` is a
#'   `shift_magnitude_groups` object, a named `shift_magnitude_count_set` is
#'   returned instead.
#' @examples
#' transitions <- data.frame(
#'   rate_change = c("root", "increase", "increase", "decrease", "decrease"),
#'   rate_delta = c(NA, 3, 4, -1, -2)
#' )
#' shift_magnitude_counts(list(run1 = transitions, run2 = transitions))
#' shift_magnitude_counts(shift_magnitude_groups(A = transitions, B = transitions))
#' @export
shift_magnitude_counts <- function(
  x = NULL,
  tree = NULL,
  state_values = NULL
) {
  if (inherits(x, "shift_magnitude_groups")) {
    if (!is.null(tree) || !is.null(state_values)) {
      stop(
        "`tree` and `state_values` are not used when `x` is a shift_magnitude_groups object.",
        call. = FALSE
      )
    }
    groups <- .shift_magnitude_group_inputs(x)
    counts <- lapply(names(groups), function(group_name) {
      out <- tryCatch(
        shift_magnitude_counts(groups[[group_name]]),
        error = function(err) {
          stop(
            "Error in analysis group `",
            group_name,
            "`: ",
            conditionMessage(err),
            call. = FALSE
          )
        }
      )
      attr(out, "analysis") <- group_name
      out
    })
    names(counts) <- names(groups)
    class(counts) <- c("shift_magnitude_count_set", "list")
    return(counts)
  }

  out <- .shift_magnitude_count_rows(
    x = x,
    tree = tree,
    state_values = state_values
  )
  rownames(out) <- NULL
  class(out) <- c("shift_magnitude_counts", "data.frame")
  out
}

#' Plot One Shift Count Or Frequency Comparison
#'
#' Summarize inferred rate increases and decreases for one set of search
#' outputs. The default paired frequency plot mirrors the count/frequency
#' comparison used for Berv et al. (2026) Supplementary Figure 5.
#'
#' @param x A `shift_magnitude_counts` object.
#' @param statistic Plot `"frequency"` or raw `"count"` values.
#' @param colors Named colors for `"increase"` and `"decrease"`.
#' @param point_cex Point size for run-level observations.
#' @param line_col Color for paired run-level segments.
#' @param boxplot Logical; overlay compact boxplots.
#' @param paired_lines Logical; draw paired decrease-to-increase segments for
#'   each run.
#' @param paired_test Logical; annotate the plot with a paired Wilcoxon test
#'   when at least two runs are available.
#' @param reference_line Optional horizontal reference line. Defaults to `0.5`
#'   for frequencies and no line for counts.
#' @param ylim Optional common y-axis limits.
#' @param ylab,main Plot labels.
#' @param mar,oma Graphical margins passed to [graphics::par()].
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns a list with the count data, paired tests, and plot
#'   settings.
#' @examples
#' transitions <- data.frame(
#'   rate_change = c("root", "increase", "increase", "decrease", "decrease"),
#'   rate_delta = c(NA, 3, 4, -1, -2)
#' )
#' counts <- shift_magnitude_counts(list(run1 = transitions, run2 = transitions))
#' plot(counts)
#' @method plot shift_magnitude_counts
#' @export
plot.shift_magnitude_counts <- function(
  x,
  statistic = c("frequency", "count"),
  colors = c(increase = "#377eb8", decrease = "#e41a1c"),
  point_cex = 1.15,
  line_col = grDevices::adjustcolor("grey35", alpha.f = 0.35),
  boxplot = TRUE,
  paired_lines = TRUE,
  paired_test = TRUE,
  reference_line = NULL,
  ylim = NULL,
  ylab = NULL,
  main = NULL,
  mar = c(4.5, 4.2, 2.4, 0.8),
  oma = c(0, 0, 0, 0),
  ...
) {
  if (!inherits(x, "shift_magnitude_counts")) {
    stop("`x` must be a shift_magnitude_counts object.", call. = FALSE)
  }
  statistic <- match.arg(statistic)
  .shift_check_logical(boxplot, "boxplot")
  .shift_check_logical(paired_lines, "paired_lines")
  .shift_check_logical(paired_test, "paired_test")
  colors <- .shift_magnitude_plot_colors(colors)

  count_data <- x
  count_tests <- .shift_magnitude_plot_count_tests(count_data, statistic)
  if (is.null(ylim)) {
    ylim <- .shift_magnitude_plot_count_ylim(count_data, statistic)
  }
  if (is.null(ylab)) {
    ylab <- if (identical(statistic, "frequency")) "Shift frequency" else "Shift count"
  }
  if (is.null(reference_line) && identical(statistic, "frequency")) {
    reference_line <- 0.5
  }
  if (is.null(main)) {
    main <- attr(count_data, "analysis", exact = TRUE)
    if (is.null(main)) {
      main <- ""
    }
  } else {
    main <- .shift_magnitude_plot_single_title(main)
  }
  old_mar <- graphics::par("mar")
  old_oma <- graphics::par("oma")
  on.exit({
    graphics::par(mar = old_mar, oma = old_oma)
  }, add = TRUE)
  graphics::par(mar = mar, oma = oma)

  values <- .shift_magnitude_plot_count_values(count_data, statistic)
  graphics::plot(
    NA,
    NA,
    xlim = c(0.55, 2.45),
    ylim = ylim,
    xaxt = "n",
    xlab = "",
    ylab = ylab,
    main = main,
    las = 1,
    bty = "n",
    ...
  )
  if (!is.null(reference_line)) {
    graphics::abline(
      h = reference_line,
      lty = 3,
      col = grDevices::adjustcolor("grey20", alpha.f = 0.7)
    )
  }
  if (isTRUE(boxplot)) {
    graphics::boxplot(
      list(decrease = values$decrease, increase = values$increase),
      at = c(1, 2),
      add = TRUE,
      axes = FALSE,
      outline = FALSE,
      boxwex = 0.26,
      col = grDevices::adjustcolor(
        c(colors[["decrease"]], colors[["increase"]]),
        alpha.f = 0.16
      ),
      border = c(colors[["decrease"]], colors[["increase"]])
    )
  }
  offsets <- .shift_magnitude_plot_offsets(nrow(count_data))
  if (isTRUE(paired_lines)) {
    graphics::segments(
      x0 = 1 + offsets,
      y0 = values$decrease,
      x1 = 2 + offsets,
      y1 = values$increase,
      col = line_col
    )
  }
  graphics::points(
    1 + offsets,
    values$decrease,
    pch = 21,
    bg = colors[["decrease"]],
    col = "white",
    cex = point_cex
  )
  graphics::points(
    2 + offsets,
    values$increase,
    pch = 21,
    bg = colors[["increase"]],
    col = "white",
    cex = point_cex
  )
  graphics::axis(1, at = c(1, 2), labels = c("decrease", "increase"))
  graphics::box(bty = "L")
  if (isTRUE(paired_test)) {
    .shift_magnitude_plot_count_annotation(count_tests)
  }

  invisible(list(
    data = count_data,
    tests = count_tests,
    settings = list(
      statistic = statistic,
      analysis = attr(count_data, "analysis", exact = TRUE)
    )
  ))
}

#' Plot A Shift-Magnitude Count Set
#'
#' Count sets contain multiple named count tables. Select one named count
#' table, for example `x[["GIC"]]`, before plotting.
#'
#' @param x A `shift_magnitude_count_set` object.
#' @param ... Ignored.
#'
#' @return Stops with an informative message.
#' @method plot shift_magnitude_count_set
#' @export
plot.shift_magnitude_count_set <- function(x, ...) {
  stop(
    "`x` contains multiple count groups; select one with `x[[\"name\"]]` before plotting.",
    call. = FALSE
  )
}

.shift_magnitude_distribution_prepare <- function(
  comparison,
  scale_by_frequency,
  density_args,
  xlim,
  ylim,
  xlab,
  ylab,
  main,
  annotate,
  show_modes,
  show_sample_size,
  show_test_label,
  show_modal_ratio
) {
  measure <- comparison$settings$measure
  transform <- comparison$settings$transform
  densities <- .shift_magnitude_plot_density_pair(
    comparison,
    scale_by_frequency = scale_by_frequency,
    density_args = density_args
  )

  if (is.null(xlim)) {
    xlim <- .shift_magnitude_plot_density_xlim(list(densities))
  }
  if (is.null(ylim)) {
    ylim <- .shift_magnitude_plot_density_ylim(list(densities))
  }
  if (is.null(xlab)) {
    xlab <- .shift_magnitude_plot_xlab(measure, transform, comparison)
  }
  if (is.null(ylab)) {
    ylab <- if (isTRUE(scale_by_frequency)) "Relative frequency density" else "Density"
  }
  main <- .shift_magnitude_plot_single_title(main)

  out <- list(
    comparison = comparison,
    summary = comparison$summary,
    tests = comparison$tests,
    modal = comparison$modal,
    values = comparison$values,
    densities = densities,
    settings = list(
      scale_by_frequency = scale_by_frequency,
      measure = measure,
      transform = transform,
      density_args = density_args,
      xlim = xlim,
      ylim = ylim,
      xlab = xlab,
      ylab = ylab,
      main = main,
      annotate = annotate,
      show_modes = show_modes,
      show_sample_size = show_sample_size,
      show_test_label = show_test_label,
      show_modal_ratio = show_modal_ratio
    )
  )
  class(out) <- c("shift_magnitude_density_plot", "list")
  out
}

.shift_magnitude_distribution_draw <- function(
  prepared,
  colors,
  fill_alpha,
  lwd,
  legend,
  rug,
  mar,
  oma,
  ...
) {
  comparison <- prepared$comparison
  densities <- prepared$densities
  settings <- prepared$settings
  old_mar <- graphics::par("mar")
  old_oma <- graphics::par("oma")
  on.exit({
    graphics::par(mar = old_mar, oma = old_oma)
  }, add = TRUE)
  graphics::par(mar = mar, oma = oma)
  graphics::plot(
    NA,
    NA,
    xlim = settings$xlim,
    ylim = settings$ylim,
    xlab = settings$xlab,
    ylab = settings$ylab,
    main = settings$main,
    las = 1,
    bty = "n",
    ...
  )
  .shift_magnitude_plot_draw_density(
    densities$decrease,
    col = colors[["decrease"]],
    fill_alpha = fill_alpha,
    lwd = lwd
  )
  .shift_magnitude_plot_draw_density(
    densities$increase,
    col = colors[["increase"]],
    fill_alpha = fill_alpha,
    lwd = lwd
  )
  if (isTRUE(settings$show_modes)) {
    .shift_magnitude_plot_mode_lines(comparison, colors)
  }
  graphics::box(bty = "L")

  if (isTRUE(rug)) {
    values <- comparison$values
    graphics::rug(
      values$value[values$rate_change == "decrease"],
      col = grDevices::adjustcolor(colors[["decrease"]], alpha.f = 0.45)
    )
    graphics::rug(
      values$value[values$rate_change == "increase"],
      col = grDevices::adjustcolor(colors[["increase"]], alpha.f = 0.45)
    )
  }
  if (isTRUE(legend)) {
    .shift_magnitude_plot_density_legend(colors, fill_alpha)
  }
  if (isTRUE(settings$annotate)) {
    .shift_magnitude_plot_density_annotation(
      comparison,
      show_sample_size = settings$show_sample_size,
      show_test_label = settings$show_test_label,
      show_modal_ratio = settings$show_modal_ratio
    )
  }

  invisible(prepared)
}

.shift_magnitude_plot_density_pair <- function(
  comparison,
  scale_by_frequency,
  density_args
) {
  values <- comparison$values
  increases <- values$value[values$rate_change == "increase"]
  decreases <- values$value[values$rate_change == "decrease"]
  density_increase <- .shift_magnitude_plot_density(increases, density_args)
  density_decrease <- .shift_magnitude_plot_density(decreases, density_args)

  if (isTRUE(scale_by_frequency)) {
    total <- length(increases) + length(decreases)
    density_increase$y <- density_increase$y * length(increases) / total
    density_decrease$y <- density_decrease$y * length(decreases) / total
  }

  list(increase = density_increase, decrease = density_decrease)
}

.shift_magnitude_plot_density <- function(x, density_args) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    stop("At least two finite values are required to draw a density.", call. = FALSE)
  }
  density_args$x <- NULL
  density_args$na.rm <- NULL
  do.call(stats::density, c(list(x = x, na.rm = TRUE), density_args))
}

.shift_magnitude_plot_density_xlim <- function(densities) {
  x_values <- unlist(lapply(densities, function(panel) {
    c(panel$increase$x, panel$decrease$x)
  }), use.names = FALSE)
  range(x_values, finite = TRUE)
}

.shift_magnitude_plot_density_ylim <- function(densities) {
  y_values <- unlist(lapply(densities, function(panel) {
    c(panel$increase$y, panel$decrease$y)
  }), use.names = FALSE)
  c(0, max(y_values, na.rm = TRUE) * 1.18)
}

.shift_magnitude_plot_xlab <- function(measure, transform, comparison) {
  if (inherits(comparison, "shift_magnitude_comparison") &&
      !is.null(comparison$settings$measure)) {
    measure <- comparison$settings$measure
  }
  if (inherits(comparison, "shift_magnitude_comparison") &&
      !is.null(comparison$settings$transform)) {
    transform <- comparison$settings$transform
  }
  label <- switch(
    measure,
    rate_delta = "rate_delta",
    percentage_change = "percentage_change",
    log_ratio = "log_ratio"
  )
  switch(
    transform,
    log_absolute = paste0("log(abs(", label, "))"),
    absolute = paste0("abs(", label, ")"),
    signed = label
  )
}

.shift_magnitude_plot_single_title <- function(main) {
  if (is.null(main)) {
    return("")
  }
  if (!is.character(main) || length(main) != 1L || is.na(main)) {
    stop("`main` must be NULL or a single string.", call. = FALSE)
  }
  main
}

.shift_magnitude_plot_draw_density <- function(density, col, fill_alpha, lwd) {
  graphics::polygon(
    c(density$x, rev(density$x)),
    c(density$y, rep(0, length(density$y))),
    col = grDevices::adjustcolor(col, alpha.f = fill_alpha),
    border = NA
  )
  graphics::lines(density$x, density$y, col = col, lwd = lwd)
}

.shift_magnitude_plot_mode_lines <- function(comparison, colors) {
  modal <- comparison$modal
  if (is.finite(modal$mode_decrease)) {
    graphics::abline(v = modal$mode_decrease, col = colors[["decrease"]], lty = 2)
  }
  if (is.finite(modal$mode_increase)) {
    graphics::abline(v = modal$mode_increase, col = colors[["increase"]], lty = 2)
  }
}

.shift_magnitude_plot_density_legend <- function(colors, fill_alpha) {
  graphics::legend(
    "topright",
    bty = "n",
    legend = c("Increases", "Decreases"),
    col = c(colors[["increase"]], colors[["decrease"]]),
    fill = grDevices::adjustcolor(
      c(colors[["increase"]], colors[["decrease"]]),
      alpha.f = fill_alpha
    ),
    lwd = 2
  )
}

.shift_magnitude_plot_density_annotation <- function(
  comparison,
  show_sample_size,
  show_test_label,
  show_modal_ratio
) {
  label <- .shift_magnitude_plot_density_label(
    comparison,
    show_sample_size = show_sample_size,
    show_test_label = show_test_label,
    show_modal_ratio = show_modal_ratio
  )
  if (!nzchar(label)) {
    return(invisible(NULL))
  }
  usr <- graphics::par("usr")
  x <- usr[2L] - 0.04 * diff(usr[1:2])
  y <- usr[3L] + 0.72 * diff(usr[3:4])
  graphics::text(
    x,
    y,
    labels = label,
    adj = c(1, 0.5),
    cex = 0.95
  )
}

.shift_magnitude_plot_density_label <- function(
  comparison,
  show_sample_size = TRUE,
  show_test_label = TRUE,
  show_modal_ratio = TRUE
) {
  summary <- comparison$summary
  inc_n <- summary$n[summary$rate_change == "increase"]
  dec_n <- summary$n[summary$rate_change == "decrease"]
  label <- character()
  if (isTRUE(show_sample_size)) {
    label <- c(label, sprintf("n inc/dec = %s/%s", inc_n, dec_n))
  }
  ks <- comparison$tests[comparison$tests$test == "ks", , drop = FALSE]
  if (isTRUE(show_test_label) && nrow(ks) > 0L) {
    label <- c(
      label,
      sprintf(
        "KS D = %.3f, P = %s",
        ks$statistic,
        .shift_magnitude_plot_p_value(ks$p_value)
      )
    )
    if (is.finite(ks$bootstrap_mean_p)) {
      label <- c(
        label,
        sprintf("boot mean P = %s", .shift_magnitude_plot_p_value(ks$bootstrap_mean_p))
      )
    }
  }
  if (isTRUE(show_modal_ratio) && is.finite(comparison$modal$linear_ratio)) {
    label <- c(label, sprintf("mode ratio = %.2f", comparison$modal$linear_ratio))
  }
  paste(label, collapse = "\n")
}

.shift_magnitude_count_rows <- function(x, tree, state_values) {
  if (inherits(x, "shift_magnitude_counts")) {
    if (!is.null(tree) || !is.null(state_values)) {
      stop("`tree` and `state_values` are not used when `x` is a count table.", call. = FALSE)
    }
    return(as.data.frame.shift_magnitude_counts(x))
  }

  if (inherits(x, "shift_magnitude_count_set")) {
    stop(
      "`x` contains multiple count groups; select one with `x[[\"name\"]]` or use `as.data.frame(x)`.",
      call. = FALSE
    )
  }

  if (inherits(x, "shift_magnitude_comparison")) {
    if (!is.null(tree) || !is.null(state_values)) {
      stop(
        "`tree` and `state_values` are not used when `x` is a shift-magnitude comparison.",
        call. = FALSE
      )
    }
    return(.shift_magnitude_plot_count_row_from_comparison(x, "comparison"))
  }

  if (.shift_magnitude_is_input_list(x)) {
    if (!is.null(tree) || !is.null(state_values)) {
      stop("`tree` and `state_values` are not used when `x` is a list of inputs.", call. = FALSE)
    }
    if (length(x) == 0L) {
      stop("`x` must contain at least one transition table or BMM search object.", call. = FALSE)
    }

    input_names <- names(x)
    rows <- lapply(seq_along(x), function(i) {
      source <- if (!is.null(input_names) && nzchar(input_names[[i]])) {
        input_names[[i]]
      } else {
        as.character(i)
      }
      input <- x[[i]]
      if (inherits(input, "shift_magnitude_comparison")) {
        return(.shift_magnitude_plot_count_row_from_comparison(input, source))
      }
      if (.shift_magnitude_is_input_list(input)) {
        stop(
          "`x` appears to contain analysis groups. Wrap grouped inputs with `shift_magnitude_groups()`.",
          call. = FALSE
        )
      }
      transitions <- .shift_magnitude_resolve_transitions(input)
      .shift_magnitude_plot_count_row(transitions, source)
    })
    return(do.call(rbind, rows))
  }

  transitions <- .shift_magnitude_resolve_transitions(
    x = x,
    tree = tree,
    state_values = state_values
  )
  .shift_magnitude_plot_count_row(transitions, "input")
}

.shift_magnitude_plot_count_row_from_comparison <- function(comparison, source) {
  summary <- comparison$summary
  increase <- summary$n[summary$rate_change == "increase"]
  decrease <- summary$n[summary$rate_change == "decrease"]
  .shift_magnitude_plot_count_row_from_counts(increase, decrease, source)
}

.shift_magnitude_plot_count_row <- function(transitions, source) {
  counts <- table(factor(transitions$rate_change, levels = c("increase", "decrease")))
  .shift_magnitude_plot_count_row_from_counts(
    increase = unname(counts[["increase"]]),
    decrease = unname(counts[["decrease"]]),
    source = source
  )
}

.shift_magnitude_plot_count_row_from_counts <- function(
  increase,
  decrease,
  source
) {
  total <- increase + decrease
  data.frame(
    source = source,
    decrease = as.integer(decrease),
    increase = as.integer(increase),
    total = as.integer(total),
    decrease_frequency = if (total > 0L) decrease / total else NA_real_,
    increase_frequency = if (total > 0L) increase / total else NA_real_,
    stringsAsFactors = FALSE
  )
}

.shift_magnitude_plot_count_tests <- function(count_data, statistic) {
  values <- .shift_magnitude_plot_count_values(count_data, statistic)
  finite <- is.finite(values$decrease) & is.finite(values$increase)
  decrease <- values$decrease[finite]
  increase <- values$increase[finite]
  analysis <- attr(count_data, "analysis", exact = TRUE)
  if (is.null(analysis)) {
    analysis <- NA_character_
  }
  if (length(increase) < 2L) {
    return(data.frame(
      analysis = analysis,
      n = length(increase),
      wilcox_statistic = NA_real_,
      wilcox_p_value = NA_real_,
      t_statistic = NA_real_,
      t_p_value = NA_real_,
      mean_difference = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  wilcox <- suppressWarnings(stats::wilcox.test(
    decrease,
    increase,
    paired = TRUE,
    exact = FALSE
  ))
  t_test <- tryCatch(
    stats::t.test(decrease, increase, paired = TRUE),
    error = function(err) NULL
  )
  data.frame(
    analysis = analysis,
    n = length(increase),
    wilcox_statistic = unname(wilcox$statistic[[1L]]),
    wilcox_p_value = wilcox$p.value,
    t_statistic = if (is.null(t_test)) NA_real_ else unname(t_test$statistic[[1L]]),
    t_p_value = if (is.null(t_test)) NA_real_ else t_test$p.value,
    mean_difference = mean(increase - decrease),
    stringsAsFactors = FALSE
  )
}

.shift_magnitude_plot_count_values <- function(panel, statistic) {
  if (identical(statistic, "frequency")) {
    return(list(
      decrease = panel$decrease_frequency,
      increase = panel$increase_frequency
    ))
  }
  list(decrease = panel$decrease, increase = panel$increase)
}

.shift_magnitude_plot_count_ylim <- function(count_data, statistic) {
  values <- .shift_magnitude_plot_count_values(count_data, statistic)
  y <- c(values$decrease, values$increase)
  if (identical(statistic, "frequency")) {
    return(c(0, 1))
  }
  c(0, max(y, na.rm = TRUE) * 1.18)
}

.shift_magnitude_plot_offsets <- function(n) {
  if (n <= 1L) {
    return(0)
  }
  seq(-0.07, 0.07, length.out = n)
}

.shift_magnitude_plot_count_annotation <- function(test_row) {
  if (nrow(test_row) == 0L || !is.finite(test_row$wilcox_p_value)) {
    return(invisible(NULL))
  }
  usr <- graphics::par("usr")
  graphics::text(
    mean(usr[1:2]),
    usr[4L] - 0.08 * diff(usr[3:4]),
    labels = sprintf(
      "paired Wilcoxon V = %s, P = %s",
      .shift_magnitude_plot_statistic(test_row$wilcox_statistic),
      .shift_magnitude_plot_p_value(test_row$wilcox_p_value)
    ),
    cex = 0.95
  )
  invisible(NULL)
}

.shift_magnitude_plot_statistic <- function(x) {
  if (!is.finite(x)) {
    return("NA")
  }
  format(signif(x, 4), scientific = FALSE, trim = TRUE)
}

.shift_magnitude_plot_p_value <- function(p) {
  if (!is.finite(p)) {
    return("NA")
  }
  if (p < 0.001) {
    return("<0.001")
  }
  format(signif(p, 3), scientific = FALSE, trim = TRUE)
}

.shift_magnitude_plot_colors <- function(colors) {
  if (!is.character(colors) || length(colors) < 2L) {
    stop("`colors` must provide colors for increases and decreases.", call. = FALSE)
  }
  if (is.null(names(colors)) ||
      !all(c("increase", "decrease") %in% names(colors))) {
    colors <- c(increase = colors[[1L]], decrease = colors[[2L]])
  }
  colors[c("increase", "decrease")]
}

.shift_check_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

.shift_magnitude_plot_check_alpha <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
    stop("`", name, "` must be a number between 0 and 1.", call. = FALSE)
  }
  invisible(TRUE)
}

.shift_magnitude_plot_check_density_args <- function(density_args) {
  if (!is.list(density_args)) {
    stop("`density_args` must be a named list.", call. = FALSE)
  }
  if (length(density_args) > 0L &&
      (is.null(names(density_args)) || any(!nzchar(names(density_args))))) {
    stop("`density_args` must be a named list.", call. = FALSE)
  }
  invisible(TRUE)
}

.shift_validate_node_painted_tree <- function(tree) {
  within_edge <- which(vapply(tree$maps, function(edge_map) {
    positive <- is.finite(edge_map) & edge_map > 0
    states <- names(edge_map)[positive]
    length(unique(states[!is.na(states)])) > 1L
  }, logical(1)))

  if (length(within_edge) > 0L) {
    stop(
      "`shift_transitions()` is node-based and cannot represent within-edge ",
      "SIMMAP transitions. Found within-edge transitions on edge row(s): ",
      paste(within_edge, collapse = ", "),
      ". Supply a node-painted history whose positive-length edge maps each ",
      "contain one state.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.shift_resolve_inputs <- function(x, tree = NULL, state_values = NULL) {
  .lineage_rates_check_packages()
  if (inherits(x, "phylo")) {
    if (!is.null(tree)) {
      stop("Supply a generic tree through either `x` or `tree`, not both.", call. = FALSE)
    }
    tree <- x
    x <- NULL
  }

  .lineage_rates_resolve_inputs(
    bifrost_search = x,
    tree = tree,
    state_values = state_values,
    log = FALSE,
    caller = "shift_transitions()"
  )
}

.shift_root_node <- function(tree) {
  root_node <- setdiff(unique(tree$edge[, 1]), unique(tree$edge[, 2]))
  if (length(root_node) != 1L) {
    stop("Cannot identify a unique root node in `tree`.", call. = FALSE)
  }
  as.integer(root_node)
}

.shift_node_heights <- function(tree) {
  node_count <- ape::Ntip(tree) + tree$Nnode
  heights <- rep(NA_real_, node_count)
  heights[.shift_root_node(tree)] <- 0
  node_heights <- phytools::nodeHeights(tree)
  for (i in seq_len(nrow(node_heights))) {
    heights[tree$edge[i, 2]] <- node_heights[i, 2]
  }
  heights
}

.shift_state_at_node <- function(tree, states, node) {
  if (node <= ape::Ntip(tree)) {
    state <- states$tips[tree$tip.label[node]]
    if (!is.na(state)) {
      return(unname(state))
    }
  }
  unname(states$nodes[as.character(node)])
}

.shift_rate_change <- function(rate_delta) {
  if (!is.finite(rate_delta) || rate_delta == 0) {
    return("no_change")
  }
  if (rate_delta > 0) "increase" else "decrease"
}

.shift_empty_transitions <- function() {
  data.frame(
    node = integer(),
    node_label = character(),
    parent_node = integer(),
    parent_state = character(),
    parent_rate = numeric(),
    child_node = integer(),
    child_state = character(),
    child_rate = numeric(),
    height = numeric(),
    age = numeric(),
    edge_length = numeric(),
    rate_delta = numeric(),
    percentage_change = numeric(),
    log_ratio = numeric(),
    rate_change = character(),
    stringsAsFactors = FALSE
  )
}

.shift_node_marks_resolve <- function(x,
                                      tree = NULL,
                                      state_values = NULL,
                                      transitions = NULL,
                                      ic_weights = NULL,
                                      support_threshold = 0.9,
                                      letter_order = "child_state",
                                      marker_base_cex = 1.5,
                                      marker_scale_factor = 0.08,
                                      marker_transform = "square_root",
                                      marker_max_cex = Inf) {
  if (.shift_node_marks_is_prepared(x)) {
    if (!is.null(tree) ||
        !is.null(state_values) ||
        !is.null(transitions) ||
        !is.null(ic_weights)) {
      stop(
        "`tree`, `state_values`, `transitions`, and `ic_weights` are not used ",
        "when `x` is a prepared shift-node mark object.",
        call. = FALSE
      )
    }
    return(.shift_node_marks_as_prepared(x))
  }

  if (is.null(transitions)) {
    if (inherits(x, "shift_transitions") ||
        (is.data.frame(x) && "rate_change" %in% names(x))) {
      if (!is.null(tree) || !is.null(state_values)) {
        stop(
          "`tree` and `state_values` are not used when `x` is a transition table.",
          call. = FALSE
        )
      }
      transitions <- x
    } else {
      transitions <- shift_transitions(
        x = x,
        tree = tree,
        state_values = state_values,
        include_root = FALSE
      )
    }
  } else if (!is.data.frame(transitions)) {
    stop("`transitions` must be a data frame.", call. = FALSE)
  } else if (!is.null(tree) || !is.null(state_values)) {
    stop(
      "`tree` and `state_values` are not used when `transitions` is supplied.",
      call. = FALSE
    )
  }

  if (is.null(ic_weights) &&
      .shift_is_bifrost_output(x) &&
      is.data.frame(x$ic_weights)) {
    ic_weights <- x$ic_weights
  }

  .shift_node_marks_prepare(
    transitions = transitions,
    ic_weights = ic_weights,
    support_threshold = support_threshold,
    letter_order = letter_order,
    marker_base_cex = marker_base_cex,
    marker_scale_factor = marker_scale_factor,
    marker_transform = marker_transform,
    marker_max_cex = marker_max_cex
  )
}

.shift_node_marks_is_prepared <- function(x) {
  inherits(x, "shift_node_marks") ||
    (
      is.list(x) &&
        !is.data.frame(x) &&
        is.data.frame(x$marks) &&
        !is.null(x$summary)
    )
}

.shift_node_marks_as_prepared <- function(x) {
  if (!.shift_node_marks_is_prepared(x)) {
    stop(
      "`x` must be a shift_node_marks object or a compatible prepared ",
      "shift-node mark list.",
      call. = FALSE
    )
  }
  if (!inherits(x, "shift_node_marks")) {
    class(x) <- c("shift_node_marks", class(x))
  }
  x
}

.shift_node_marks_prepare <- function(transitions,
                                      ic_weights = NULL,
                                      support_threshold = 0.9,
                                      letter_order = "child_state",
                                      marker_base_cex = 1.5,
                                      marker_scale_factor = 0.08,
                                      marker_transform = "square_root",
                                      marker_max_cex = Inf) {
  .shift_node_marks_check_number(
    support_threshold,
    "support_threshold",
    minimum = 0,
    maximum = 1
  )
  .shift_node_marks_check_number(marker_base_cex, "marker_base_cex", minimum = 0)
  .shift_node_marks_check_number(marker_scale_factor, "marker_scale_factor", minimum = 0)
  .shift_node_marks_check_number(marker_max_cex, "marker_max_cex", minimum = 0, allow_inf = TRUE)

  marks <- .shift_node_marks_transition_table(transitions)
  marks <- marks[marks$rate_change %in% c("increase", "decrease"), , drop = FALSE]

  if (!is.null(ic_weights)) {
    marks <- .shift_node_marks_merge_weights(marks, ic_weights)
  } else if (!"ic_weight_withshift" %in% names(marks)) {
    marks$ic_weight_withshift <- rep(NA_real_, nrow(marks))
  }

  transformed_change <- switch(
    marker_transform,
    square_root = sqrt(abs(marks$percentage_change)),
    log1p = log1p(abs(marks$percentage_change))
  )
  marks$node_cex <- marker_base_cex + transformed_change * marker_scale_factor
  marks$node_cex[!is.finite(marks$node_cex)] <- marker_base_cex
  marks$node_cex <- pmin(marks$node_cex, marker_max_cex)

  marks$letter <- rep(NA_character_, nrow(marks))
  increase_rows <- which(marks$rate_change == "increase")
  letter_rows <- .shift_node_marks_letter_rows(
    marks = marks,
    rows = increase_rows,
    letter_order = letter_order
  )
  marks$letter[letter_rows] <- .shift_node_marks_letters(length(letter_rows))

  low_support_nodes <- marks$node[
    !is.na(marks$ic_weight_withshift) &
      marks$ic_weight_withshift < support_threshold
  ]

  increase_key <- marks[
    marks$rate_change == "increase",
    c(
      "letter",
      "child_state",
      "node",
      "age",
      "percentage_change",
      "ic_weight_withshift"
    ),
    drop = FALSE
  ]
  increase_key <- increase_key[order(increase_key$letter), , drop = FALSE]
  rownames(increase_key) <- NULL

  summary <- data.frame(
    directional_shifts = nrow(marks),
    increases = sum(marks$rate_change == "increase"),
    decreases = sum(marks$rate_change == "decrease"),
    low_support_shifts = length(low_support_nodes),
    stringsAsFactors = FALSE
  )

  out <- list(
    marks = marks,
    increase_key = increase_key,
    low_support_nodes = low_support_nodes,
    summary = summary,
    settings = list(
      support_threshold = support_threshold,
      letter_order = letter_order,
      marker_base_cex = marker_base_cex,
      marker_scale_factor = marker_scale_factor,
      marker_transform = marker_transform,
      marker_max_cex = marker_max_cex
    )
  )
  class(out) <- c("shift_node_marks", "list")
  out
}

.shift_node_marks_letter_rows <- function(marks, rows, letter_order) {
  if (length(rows) == 0L) {
    return(integer())
  }

  if (identical(letter_order, "child_state")) {
    rows[order(as.character(marks$child_state[rows]), marks$node[rows])]
  } else if (identical(letter_order, "age")) {
    rows[order(marks$age[rows], marks$node[rows])]
  } else {
    rows[order(marks$node[rows])]
  }
}

.shift_node_marks_transition_table <- function(transitions) {
  transitions <- as.data.frame(transitions)
  required <- c("node", "rate_change")
  missing <- setdiff(required, names(transitions))
  if (length(missing) > 0L) {
    stop(
      "`transitions` must include columns: ",
      paste(required, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!"percentage_change" %in% names(transitions)) {
    if (all(c("parent_rate", "child_rate") %in% names(transitions))) {
      transitions$percentage_change <- ifelse(
        is.finite(transitions$parent_rate) & transitions$parent_rate != 0,
        100 * (transitions$child_rate - transitions$parent_rate) /
          transitions$parent_rate,
        NA_real_
      )
    } else {
      transitions$percentage_change <- rep(NA_real_, nrow(transitions))
    }
  }

  if (!"child_state" %in% names(transitions)) {
    transitions$child_state <- if ("node_label" %in% names(transitions)) {
      as.character(transitions$node_label)
    } else {
      as.character(transitions$node)
    }
  }
  if (!"age" %in% names(transitions)) {
    transitions$age <- rep(NA_real_, nrow(transitions))
  }
  transitions
}

.shift_node_marks_merge_weights <- function(marks, ic_weights) {
  if (!is.data.frame(ic_weights) ||
      !all(c("node", "ic_weight_withshift") %in% names(ic_weights))) {
    stop(
      "`ic_weights` must be a data frame with `node` and ",
      "`ic_weight_withshift` columns.",
      call. = FALSE
    )
  }

  if ("ic_weight_withshift" %in% names(marks)) {
    marks$ic_weight_withshift <- NULL
  }
  marks$.shift_node_mark_order <- seq_len(nrow(marks))
  weights <- ic_weights[, c("node", "ic_weight_withshift"), drop = FALSE]
  weights <- weights[!duplicated(weights$node), , drop = FALSE]
  marks <- merge(marks, weights, by = "node", all.x = TRUE, sort = FALSE)
  marks <- marks[order(marks$.shift_node_mark_order), , drop = FALSE]
  marks$.shift_node_mark_order <- NULL
  rownames(marks) <- NULL
  marks
}

.shift_node_marks_letters <- function(n) {
  if (n == 0L) {
    return(character())
  }
  if (n <= length(LETTERS)) {
    return(LETTERS[seq_len(n)])
  }
  c(LETTERS, paste0("N", seq_len(n - length(LETTERS))))
}

.shift_node_marks_draw <- function(annotation,
                                   rate_changes = c("increase", "decrease"),
                                   show_letters = TRUE,
                                   show_low_support = TRUE,
                                   show_legend = TRUE,
                                   marker_alpha = 0.75,
                                   letter_cex = 0.75,
                                   increase_fill = "white",
                                   decrease_fill = "#2b6cb0",
                                   marker_col = "#111827",
                                   low_support_cex = 0.70,
                                   low_support_col = grDevices::adjustcolor("black", alpha.f = 0.58),
                                   legend_position = "topleft",
                                   legend_cex = 0.72) {
  .shift_node_marks_check_number(
    marker_alpha,
    "marker_alpha",
    minimum = 0,
    maximum = 1
  )
  .shift_node_marks_check_number(letter_cex, "letter_cex", minimum = 0)
  .shift_node_marks_check_number(low_support_cex, "low_support_cex", minimum = 0)
  .shift_node_marks_check_number(legend_cex, "legend_cex", minimum = 0)

  marks <- annotation$marks
  if (nrow(marks) == 0L) {
    return(invisible(annotation))
  }

  for (change in intersect(c("decrease", "increase"), rate_changes)) {
    rows <- marks$rate_change == change
    if (!any(rows)) {
      next
    }

    node_fill <- if (change == "increase") increase_fill else decrease_fill
    ape::nodelabels(
      node = marks$node[rows],
      pch = 21,
      cex = marks$node_cex[rows],
      bg = grDevices::adjustcolor(node_fill, alpha.f = marker_alpha),
      col = grDevices::adjustcolor(marker_col, alpha.f = marker_alpha)
    )
  }

  letter_rows <- marks$rate_change %in% rate_changes & !is.na(marks$letter)
  if (show_letters && any(letter_rows)) {
    ape::nodelabels(
      node = marks$node[letter_rows],
      text = marks$letter[letter_rows],
      frame = "none",
      cex = letter_cex,
      col = marker_col
    )
  }

  low_support_nodes <- intersect(
    annotation$low_support_nodes,
    marks$node[marks$rate_change %in% rate_changes]
  )
  if (show_low_support && length(low_support_nodes) > 0L) {
    ape::nodelabels(
      node = low_support_nodes,
      pch = 16,
      cex = low_support_cex,
      col = low_support_col
    )
  }

  if (show_legend) {
    legend_items <- rate_changes
    legend_pch <- rep(21, length(legend_items))
    legend_bg <- ifelse(legend_items == "increase", increase_fill, decrease_fill)
    legend_col <- rep(marker_col, length(legend_items))
    legend_pt_cex <- rep(1.15, length(legend_items))

    if (show_low_support && length(low_support_nodes) > 0L) {
      legend_items <- c(legend_items, paste0("IC weight < ", annotation$settings$support_threshold))
      legend_pch <- c(legend_pch, 16)
      legend_bg <- c(legend_bg, NA)
      legend_col <- c(legend_col, low_support_col)
      legend_pt_cex <- c(legend_pt_cex, low_support_cex)
    }

    graphics::legend(
      legend_position,
      legend = legend_items,
      pch = legend_pch,
      pt.bg = legend_bg,
      col = legend_col,
      pt.cex = legend_pt_cex,
      cex = legend_cex,
      bty = "n"
    )
  }

  invisible(annotation)
}

.shift_node_marks_check_number <- function(x,
                                           name,
                                           minimum = -Inf,
                                           maximum = Inf,
                                           allow_inf = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x)
  if (!allow_inf) {
    valid <- valid && is.finite(x)
  }
  valid <- valid && x >= minimum && x <= maximum
  if (!valid) {
    stop("`", name, "` must be a single numeric value.", call. = FALSE)
  }
  invisible(TRUE)
}

.shift_magnitude_resolve_transitions <- function(x, tree = NULL, state_values = NULL) {
  if (inherits(x, "shift_transitions") ||
      (is.data.frame(x) && "rate_change" %in% names(x))) {
    if (!is.null(tree) || !is.null(state_values)) {
      stop(
        "`tree` and `state_values` are not used when `x` is a transition table.",
        call. = FALSE
      )
    }
    return(x)
  }

  if (.shift_magnitude_is_input_list(x)) {
    if (!is.null(tree) || !is.null(state_values)) {
      stop(
        "`tree` and `state_values` are not used when `x` is a list of inputs.",
        call. = FALSE
      )
    }
    return(.shift_magnitude_bind_transition_list(x))
  }

  shift_transitions(
    x = x,
    tree = tree,
    state_values = state_values,
    include_root = TRUE
  )
}

.shift_magnitude_is_input_list <- function(x) {
  is.list(x) &&
    !is.data.frame(x) &&
    !inherits(x, c(
      "shift_transitions",
      "shift_magnitude_groups",
      "shift_magnitude_comparison",
      "shift_magnitude_comparison_set",
      "shift_magnitude_count_set",
      "bifrost_search",
      "phylo",
      "multiPhylo"
    )) &&
    !.shift_is_bifrost_output(x)
}

.shift_magnitude_is_group_container <- function(x) {
  is.list(x) &&
    !is.data.frame(x) &&
    !inherits(x, c(
      "shift_transitions",
      "shift_magnitude_groups",
      "shift_magnitude_comparison",
      "shift_magnitude_counts",
      "shift_magnitude_comparison_set",
      "shift_magnitude_count_set",
      "bifrost_search",
      "phylo",
      "multiPhylo"
    )) &&
    !.shift_is_bifrost_output(x)
}

.shift_magnitude_group_inputs <- function(x, labels = NULL) {
  if (inherits(x, "shift_magnitude_groups")) {
    x <- unclass(x)
  }
  if (!is.list(x) ||
      is.data.frame(x) ||
      inherits(x, c(
        "shift_magnitude_groups",
        "shift_magnitude_comparison",
        "shift_magnitude_comparison_set",
        "shift_magnitude_counts",
        "shift_magnitude_count_set",
        "shift_transitions",
        "bifrost_search",
        "phylo",
        "multiPhylo"
      )) ||
      .shift_is_bifrost_output(x)) {
    stop("`x` must be a non-empty list of analysis groups.", call. = FALSE)
  }
  if (length(x) == 0L) {
    stop("`x` must contain at least one analysis group.", call. = FALSE)
  }

  if (is.null(labels)) {
    group_names <- names(x)
    if (is.null(group_names)) {
      group_names <- rep("", length(x))
    }
    group_names <- as.character(group_names)
    group_names[is.na(group_names)] <- ""
    empty <- !nzchar(group_names)
    group_names[empty] <- paste("Analysis", which(empty))
  } else {
    group_names <- as.character(labels)
    if (length(group_names) != length(x) ||
        anyNA(group_names) ||
        any(!nzchar(group_names))) {
      stop("`labels` must contain one non-empty label per analysis group.", call. = FALSE)
    }
  }
  if (anyDuplicated(group_names)) {
    stop("Analysis group labels must be unique.", call. = FALSE)
  }
  names(x) <- group_names
  x
}

.shift_is_bifrost_output <- function(x) {
  inherits(x, "bifrost_search") ||
    (
      is.list(x) &&
        !is.null(x$tree_no_uncertainty_untransformed) &&
        !is.null(x$model_no_uncertainty)
    )
}

.shift_magnitude_bind_transition_list <- function(x) {
  if (length(x) == 0L) {
    stop("`x` must contain at least one transition table or BMM search object.", call. = FALSE)
  }

  input_names <- names(x)
  rows <- lapply(seq_along(x), function(i) {
    source <- if (!is.null(input_names) && nzchar(input_names[[i]])) {
      input_names[[i]]
    } else {
      as.character(i)
    }

    item <- x[[i]]
    if (.shift_magnitude_is_input_list(item)) {
      stop(
        "`x` appears to contain analysis groups. Wrap grouped inputs with `shift_magnitude_groups()`.",
        call. = FALSE
      )
    }

    data <- tryCatch(
      {
        if (inherits(item, "shift_transitions") ||
            (is.data.frame(item) && "rate_change" %in% names(item))) {
          as.data.frame(item)
        } else {
          as.data.frame(shift_transitions(item, include_root = TRUE))
        }
      },
      error = function(err) {
        stop(
          "Could not resolve `x[[", i, "]]` as a transition table or BMM search object: ",
          conditionMessage(err),
          call. = FALSE
        )
      }
    )
    data$source <- source
    data
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.shift_magnitude_normalize_tests <- function(tests) {
  if (!is.character(tests) || length(tests) == 0L || anyNA(tests)) {
    stop("`tests` must be a non-empty character vector.", call. = FALSE)
  }
  tests <- tolower(trimws(tests))
  if (all(!nzchar(tests))) {
    stop("`tests` must be a non-empty character vector.", call. = FALSE)
  }
  if (any(!nzchar(tests))) {
    stop("`tests` must not contain blank values.", call. = FALSE)
  }
  tests <- unique(tests)
  if ("all" %in% tests && length(tests) > 1L) {
    stop(
      "The special `all` value cannot be combined with other test names.",
      call. = FALSE
    )
  }
  supported <- c("t", "wilcox", "ks", "all")
  unknown <- setdiff(tests, supported)
  if (length(unknown) > 0L) {
    stop(
      "`tests` contains unsupported value(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  if (identical(tests, "all")) {
    tests <- c("t", "wilcox", "ks")
  }
  tests
}

.shift_magnitude_comparison_setting_mismatches <- function(comparison, settings) {
  existing <- comparison$settings
  if (!is.list(existing)) {
    return(names(settings))
  }

  names(settings)[!vapply(names(settings), function(name) {
    identical(existing[[name]], settings[[name]])
  }, logical(1))]
}

.shift_magnitude_validate_ks <- function(ks_simulate_p_value, ks_B) {
  if (!is.logical(ks_simulate_p_value) ||
      length(ks_simulate_p_value) != 1L ||
      is.na(ks_simulate_p_value)) {
    stop("`ks_simulate_p_value` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(ks_B) ||
      length(ks_B) != 1L ||
      !is.finite(ks_B) ||
      ks_B < 1 ||
      ks_B != floor(ks_B) ||
      ks_B > .Machine$integer.max) {
    stop("`ks_B` must be a positive integer in R integer range.", call. = FALSE)
  }
  invisible(TRUE)
}

.shift_magnitude_validate_bootstrap <- function(bootstrap_p_value, bootstrap_R) {
  if (!is.logical(bootstrap_p_value) ||
      length(bootstrap_p_value) != 1L ||
      is.na(bootstrap_p_value)) {
    stop("`bootstrap_p_value` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(bootstrap_R) ||
      length(bootstrap_R) != 1L ||
      !is.finite(bootstrap_R) ||
      bootstrap_R < 1 ||
      bootstrap_R != floor(bootstrap_R) ||
      bootstrap_R > .Machine$integer.max) {
    stop("`bootstrap_R` must be a positive integer in R integer range.", call. = FALSE)
  }
  invisible(TRUE)
}

.shift_magnitude_resolve_rep_alias <- function(alias_value, current, alias) {
  if (is.null(alias_value)) {
    return(current)
  }
  if (!is.numeric(alias_value) ||
      length(alias_value) != 1L ||
      !is.finite(alias_value) ||
      alias_value < 1 ||
      alias_value != floor(alias_value) ||
      alias_value > .Machine$integer.max) {
    stop("`", alias, "` must be a positive integer in R integer range.", call. = FALSE)
  }
  as.integer(alias_value)
}

.shift_magnitude_values <- function(transitions, measure, transform) {
  required <- c("rate_change", measure)
  missing <- setdiff(required, names(transitions))
  if (length(missing) > 0L) {
    stop(
      "`x` is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  keep <- transitions$rate_change %in% c("increase", "decrease")
  data <- data.frame(
    rate_change = as.character(transitions$rate_change[keep]),
    raw_value = as.numeric(transitions[[measure]][keep]),
    stringsAsFactors = FALSE
  )
  if ("source" %in% names(transitions)) {
    data$source <- as.character(transitions$source[keep])
  }
  data <- data[is.finite(data$raw_value), , drop = FALSE]

  if (nrow(data) == 0L) {
    stop("No finite increase or decrease values are available to compare.", call. = FALSE)
  }

  if (identical(transform, "signed")) {
    data$value <- data$raw_value
  } else {
    data$value <- abs(data$raw_value)
    if (identical(transform, "log_absolute")) {
      data <- data[data$value > 0, , drop = FALSE]
      data$value <- base::log(data$value)
    }
  }

  counts <- table(factor(data$rate_change, levels = c("increase", "decrease")))
  if (any(counts < 2L)) {
    stop(
      "At least two finite increase and decrease values are required. ",
      "Found increase = ", counts[["increase"]],
      ", decrease = ", counts[["decrease"]], ".",
      call. = FALSE
    )
  }

  data
}

.shift_magnitude_summary <- function(values) {
  groups <- c("increase", "decrease")
  out <- lapply(groups, function(group) {
    x <- values$value[values$rate_change == group]
    raw <- values$raw_value[values$rate_change == group]
    data.frame(
      rate_change = group,
      n = length(x),
      mean = mean(x),
      median = stats::median(x),
      sd = if (length(x) > 1L) stats::sd(x) else NA_real_,
      min = min(x),
      max = max(x),
      mean_signed_raw = mean(raw),
      mean_abs_raw = mean(abs(raw)),
      median_abs_raw = stats::median(abs(raw)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

.shift_magnitude_modal <- function(values, transform) {
  increases <- values$value[values$rate_change == "increase"]
  decreases <- values$value[values$rate_change == "decrease"]
  mode_increase <- .shift_density_mode(increases)
  mode_decrease <- .shift_density_mode(decreases)

  linear_ratio <- NA_real_
  if (identical(transform, "log_absolute") &&
      is.finite(mode_increase) &&
      is.finite(mode_decrease)) {
    linear_ratio <- exp(mode_increase) / exp(mode_decrease)
  } else if (identical(transform, "absolute") &&
             is.finite(mode_increase) &&
             is.finite(mode_decrease) &&
             mode_decrease != 0) {
    linear_ratio <- mode_increase / mode_decrease
  }

  data.frame(
    mode_increase = mode_increase,
    mode_decrease = mode_decrease,
    mode_difference = mode_increase - mode_decrease,
    linear_ratio = linear_ratio,
    stringsAsFactors = FALSE
  )
}

.shift_density_mode <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  if (length(unique(x)) == 1L) {
    return(x[1L])
  }
  dens <- stats::density(x, na.rm = TRUE)
  dens$x[which.max(dens$y)]
}

.shift_magnitude_tests <- function(
  values,
  tests,
  alternative,
  ks_simulate_p_value = TRUE,
  ks_B = 10000L,
  bootstrap_p_value = FALSE,
  bootstrap_R = 100L
) {
  increases <- values$value[values$rate_change == "increase"]
  decreases <- values$value[values$rate_change == "decrease"]

  rows <- lapply(tests, function(test) {
    result <- .shift_magnitude_run_test(
      test = test,
      increases = increases,
      decreases = decreases,
      alternative = alternative,
      ks_simulate_p_value = ks_simulate_p_value,
      ks_B = ks_B,
      conf_int = TRUE
    )
    row <- .shift_magnitude_test_row(test, result, increases, decreases)
    if (isTRUE(bootstrap_p_value)) {
      boot_p <- .shift_magnitude_bootstrap_p_values(
        test = test,
        increases = increases,
        decreases = decreases,
        alternative = alternative,
        ks_simulate_p_value = ks_simulate_p_value,
        ks_B = ks_B,
        bootstrap_R = bootstrap_R
      )
      row$bootstrap_mean_p <- mean(boot_p, na.rm = TRUE)
      row$bootstrap_R <- as.integer(bootstrap_R)
    } else {
      row$bootstrap_mean_p <- NA_real_
      row$bootstrap_R <- NA_integer_
    }
    row
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.shift_magnitude_run_test <- function(
  test,
  increases,
  decreases,
  alternative,
  ks_simulate_p_value,
  ks_B,
  conf_int = FALSE
) {
  switch(
    test,
    t = stats::t.test(
      increases,
      decreases,
      var.equal = FALSE,
      alternative = alternative
    ),
    wilcox = stats::wilcox.test(
      increases,
      decreases,
      paired = FALSE,
      exact = FALSE,
      conf.int = conf_int,
      alternative = alternative
    ),
    ks = suppressWarnings(stats::ks.test(
      increases,
      decreases,
      alternative = alternative,
      simulate.p.value = ks_simulate_p_value,
      B = as.integer(ks_B)
    ))
  )
}

.shift_magnitude_bootstrap_p_values <- function(
  test,
  increases,
  decreases,
  alternative,
  ks_simulate_p_value,
  ks_B,
  bootstrap_R
) {
  n_increase <- length(increases)
  n_decrease <- length(decreases)

  replicate(as.integer(bootstrap_R), {
    sampled_increases <- sample(increases, n_increase, replace = TRUE)
    sampled_decreases <- sample(decreases, n_decrease, replace = TRUE)
    result <- tryCatch(
      .shift_magnitude_run_test(
        test = test,
        increases = sampled_increases,
        decreases = sampled_decreases,
        alternative = alternative,
        ks_simulate_p_value = ks_simulate_p_value,
        ks_B = ks_B,
        conf_int = FALSE
      ),
      error = function(err) NULL
    )
    if (is.null(result)) {
      NA_real_
    } else {
      result$p.value
    }
  })
}

.shift_magnitude_test_row <- function(test, result, increases, decreases) {
  conf <- result$conf.int
  if (is.null(conf)) {
    conf <- c(NA_real_, NA_real_)
  }
  data.frame(
    test = test,
    statistic = unname(result$statistic[[1L]]),
    parameter = if (!is.null(result$parameter)) unname(result$parameter[[1L]]) else NA_real_,
    p_value = result$p.value,
    mean_increase = mean(increases),
    mean_decrease = mean(decreases),
    mean_difference = mean(increases) - mean(decreases),
    median_increase = stats::median(increases),
    median_decrease = stats::median(decreases),
    median_difference = stats::median(increases) - stats::median(decreases),
    conf_low = unname(conf[1L]),
    conf_high = unname(conf[2L]),
    method = result$method,
    stringsAsFactors = FALSE
  )
}

.shift_waiting_resolve_transitions <- function(x, tree = NULL, state_values = NULL) {
  if (inherits(x, "shift_transitions") ||
      (is.data.frame(x) && all(c("node", "height", "rate_change") %in% names(x)))) {
    if (!is.null(state_values)) {
      stop("`state_values` is not used when `x` is a transition table.", call. = FALSE)
    }
    mapped_tree <- if (!is.null(tree)) tree else attr(x, "tree", exact = TRUE)
    return(list(transitions = x, tree = mapped_tree))
  }

  transitions <- shift_transitions(
    x = x,
    tree = tree,
    state_values = state_values,
    include_root = TRUE
  )
  list(
    transitions = transitions,
    tree = attr(transitions, "tree", exact = TRUE)
  )
}

.shift_waiting_category <- function(previous_shift, next_shift) {
  if (previous_shift %in% c("increase", "decrease") &&
      next_shift %in% c("increase", "decrease")) {
    return(paste(previous_shift, next_shift, sep = "_to_"))
  }
  "other"
}

.shift_empty_global_waiting <- function() {
  data.frame(
    PreviousNode = integer(),
    NextNode = integer(),
    PreviousState = character(),
    NextState = character(),
    PreviousShift = character(),
    NextShift = character(),
    TimeToNext = numeric(),
    Category = character(),
    stringsAsFactors = FALSE
  )
}

.shift_empty_lineage_waiting <- function() {
  data.frame(
    Lineage = integer(),
    tip_label = character(),
    PreviousNode = integer(),
    NextNode = integer(),
    PreviousState = character(),
    NextState = character(),
    PreviousShift = character(),
    NextShift = character(),
    TimeToNext = numeric(),
    Category = character(),
    stringsAsFactors = FALSE
  )
}

.shift_global_waiting_times <- function(transitions) {
  transitions <- transitions[order(transitions$height, transitions$node), , drop = FALSE]
  if (nrow(transitions) < 2L) {
    return(.shift_empty_global_waiting())
  }

  previous <- transitions[-nrow(transitions), , drop = FALSE]
  next_rows <- transitions[-1L, , drop = FALSE]
  categories <- mapply(
    .shift_waiting_category,
    previous$rate_change,
    next_rows$rate_change,
    USE.NAMES = FALSE
  )

  data.frame(
    PreviousNode = previous$node,
    NextNode = next_rows$node,
    PreviousState = previous$child_state,
    NextState = next_rows$child_state,
    PreviousShift = previous$rate_change,
    NextShift = next_rows$rate_change,
    TimeToNext = next_rows$height - previous$height,
    Category = categories,
    stringsAsFactors = FALSE
  )
}

.shift_lineage_waiting_times <- function(transitions, tree) {
  shifts <- transitions[transitions$rate_change != "root", , drop = FALSE]
  if (nrow(shifts) < 2L) {
    return(list(combined = .shift_empty_lineage_waiting(), by_tip = list()))
  }

  root_node <- .shift_root_node(tree)
  by_tip <- vector("list", ape::Ntip(tree))
  names(by_tip) <- tree$tip.label

  for (tip in seq_len(ape::Ntip(tree))) {
    path <- ape::nodepath(tree, from = root_node, to = tip)
    path_shifts <- shifts[shifts$node %in% path, , drop = FALSE]
    path_shifts <- path_shifts[order(path_shifts$height, path_shifts$node), , drop = FALSE]
    if (nrow(path_shifts) < 2L) {
      by_tip[[tip]] <- .shift_empty_lineage_waiting()
      next
    }

    previous <- path_shifts[-nrow(path_shifts), , drop = FALSE]
    next_rows <- path_shifts[-1L, , drop = FALSE]
    categories <- mapply(
      .shift_waiting_category,
      previous$rate_change,
      next_rows$rate_change,
      USE.NAMES = FALSE
    )

    by_tip[[tip]] <- data.frame(
      Lineage = tip,
      tip_label = tree$tip.label[tip],
      PreviousNode = previous$node,
      NextNode = next_rows$node,
      PreviousState = previous$child_state,
      NextState = next_rows$child_state,
      PreviousShift = previous$rate_change,
      NextShift = next_rows$rate_change,
      TimeToNext = next_rows$height - previous$height,
      Category = categories,
      stringsAsFactors = FALSE
    )
  }

  non_empty <- by_tip[vapply(by_tip, nrow, integer(1)) > 0L]
  combined <- if (length(non_empty) == 0L) {
    .shift_empty_lineage_waiting()
  } else {
    out <- do.call(rbind, non_empty)
    rownames(out) <- NULL
    out
  }

  list(combined = combined, by_tip = by_tip)
}

.shift_empty_acceleration_pairs <- function() {
  data.frame(
    Lineage = integer(),
    tip_label = character(),
    older_node = integer(),
    younger_node = integer(),
    gap = numeric(),
    stringsAsFactors = FALSE
  )
}

.shift_acceleration_waiting_times <- function(transitions, tree) {
  shifts <- transitions[
    transitions$rate_change != "root" & transitions$rate_change == "increase",
    ,
    drop = FALSE
  ]
  if (nrow(shifts) < 2L) {
    pairs <- .shift_empty_acceleration_pairs()
  } else {
    root_node <- .shift_root_node(tree)
    pair_list <- vector("list", ape::Ntip(tree))
    for (tip in seq_len(ape::Ntip(tree))) {
      path <- ape::nodepath(tree, from = root_node, to = tip)
      inc <- shifts[shifts$node %in% path, , drop = FALSE]
      inc <- inc[order(inc$height, inc$node), , drop = FALSE]
      if (nrow(inc) < 2L) {
        pair_list[[tip]] <- NULL
        next
      }
      pair_list[[tip]] <- data.frame(
        Lineage = tip,
        tip_label = tree$tip.label[tip],
        older_node = inc$node[-nrow(inc)],
        younger_node = inc$node[-1L],
        gap = diff(inc$height),
        stringsAsFactors = FALSE
      )
    }
    pair_list <- Filter(Negate(is.null), pair_list)
    pairs <- if (length(pair_list) == 0L) {
      .shift_empty_acceleration_pairs()
    } else {
      out <- do.call(rbind, pair_list)
      rownames(out) <- NULL
      out
    }
  }

  collapsed <- unique(pairs[, c("older_node", "younger_node", "gap"), drop = FALSE])
  summary <- rbind(
    lineage_weighted = .shift_gap_summary(pairs$gap),
    collapsed = .shift_gap_summary(collapsed$gap)
  )
  summary <- data.frame(
    summary = rownames(summary),
    summary,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  list(
    lineage_weighted = pairs,
    collapsed = collapsed,
    summary = summary
  )
}

.shift_waiting_summary <- function(waiting_times) {
  if (is.null(waiting_times) || nrow(waiting_times) == 0L) {
    return(data.frame(
      Category = character(),
      Mean = numeric(),
      Median = numeric(),
      Variance = numeric(),
      Count = integer(),
      stringsAsFactors = FALSE
    ))
  }

  categories <- unique(waiting_times$Category)
  out <- lapply(categories, function(category) {
    values <- waiting_times$TimeToNext[waiting_times$Category == category]
    data.frame(
      Category = category,
      Mean = mean(values),
      Median = stats::median(values),
      Variance = if (length(values) > 1L) stats::var(values) else NA_real_,
      Count = length(values),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

.shift_gap_summary <- function(values) {
  if (length(values) == 0L) {
    return(c(n = 0, mean = NA_real_, median = NA_real_, sd = NA_real_))
  }
  c(
    n = length(values),
    mean = mean(values),
    median = stats::median(values),
    sd = if (length(values) > 1L) stats::sd(values) else NA_real_
  )
}

.distribution_check_univariateML <- function() {
  if (!requireNamespace("univariateML", quietly = TRUE)) {
    stop(
      "`univariateML` is required for distribution fitting. ",
      "Install it or skip fitting candidate distributions.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.distribution_normalize_models <- function(models) {
  if (!is.character(models) || length(models) == 0L || anyNA(models)) {
    stop("`models` must be a non-empty character vector.", call. = FALSE)
  }
  models <- tolower(trimws(models))
  models <- sub("^ml", "", models)
  models <- models[nzchar(models)]
  if (length(models) == 0L) {
    stop("`models` must include at least one non-empty model name.", call. = FALSE)
  }
  unique(models)
}

.rate_distribution_default_models <- function() {
  c(
    "norm",
    "gumbel",
    "logis",
    "laplace",
    "cauchy",
    "std",
    "sstd",
    "snorm",
    "ged",
    "sged"
  )
}

.rate_distribution_data <- function(x, log = TRUE) {
  original_values <- NULL
  source <- NULL
  source_column <- NA_character_
  already_logged <- FALSE

  if (is.numeric(x) && !is.data.frame(x)) {
    original_values <- as.numeric(x)
    source <- "numeric"
  } else if (is.data.frame(x)) {
    extracted <- .rate_distribution_from_data_frame(x, log = log)
    original_values <- extracted$values
    source <- "data.frame"
    source_column <- extracted$column
    already_logged <- extracted$already_logged
  } else if (inherits(x, "bifrost_search") || is.list(x)) {
    resolved <- .lineage_rates_resolve_inputs(
      bifrost_search = x,
      tree = NULL,
      state_values = NULL,
      log = FALSE,
      caller = "fit_rate_distribution()"
    )
    original_values <- as.numeric(resolved$state_values)
    source <- "bifrost_search$model_no_uncertainty$param"
  } else {
    stop(
      "`x` must be a numeric vector, data frame, `bifrost_search` object, ",
      "or compatible list.",
      call. = FALSE
    )
  }

  if (!is.numeric(original_values) ||
      length(original_values) == 0L ||
      any(!is.finite(original_values))) {
    stop("Rate values must be a non-empty finite numeric vector.", call. = FALSE)
  }

  if (isTRUE(log) && !isTRUE(already_logged)) {
    if (any(original_values <= 0)) {
      stop("Rate values must be strictly positive when `log = TRUE`.", call. = FALSE)
    }
    values <- base::log(original_values)
    transformation <- "log"
  } else {
    values <- original_values
    transformation <- if (isTRUE(already_logged)) "already_log" else "none"
  }

  list(
    values = values,
    original_values = original_values,
    source = source,
    source_column = source_column,
    transformation = transformation
  )
}

.rate_distribution_from_data_frame <- function(x, log = TRUE) {
  raw_candidates <- c("lineage_rate", "Weighted_Lineage_Value", "tip_rate", "rate")
  candidates <- if (isTRUE(log) && "log_lineage_rate" %in% names(x)) {
    "log_lineage_rate"
  } else {
    raw_candidates
  }
  column <- candidates[candidates %in% names(x)][1L]
  if (!isTRUE(log) && is.na(column) && "log_lineage_rate" %in% names(x)) {
    stop(
      "`log = FALSE` requires a raw rate column; `log_lineage_rate` is already logged.",
      call. = FALSE
    )
  }
  if (is.na(column)) {
    stop(
      "Rate data frames must include one of: `log_lineage_rate`, ",
      "`lineage_rate`, `Weighted_Lineage_Value`, `tip_rate`, or `rate`.",
      call. = FALSE
    )
  }
  list(
    values = as.numeric(x[[column]]),
    column = column,
    already_logged = identical(column, "log_lineage_rate")
  )
}

.distribution_fit_models <- function(values, models, select_by = "AIC") {
  .distribution_check_univariateML()
  models <- .distribution_normalize_models(models)
  if (!is.numeric(values) || length(values) == 0L || any(!is.finite(values))) {
    stop("`values` must be a non-empty finite numeric vector.", call. = FALSE)
  }

  rows <- list()
  failed <- .distribution_empty_failures()
  had_non_finite_statistics <- FALSE
  for (model in models) {
    fit <- tryCatch(
      univariateML::model_select(
        x = values,
        models = model,
        criterion = select_by,
        type = "continuous",
        return = "all"
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      failed <- rbind(
        failed,
        data.frame(model = model, reason = conditionMessage(fit), stringsAsFactors = FALSE)
      )
      next
    }
    fit <- as.data.frame(fit)
    if (nrow(fit) == 0L) {
      failed <- rbind(
        failed,
        data.frame(model = model, reason = "No fitted model returned.", stringsAsFactors = FALSE)
      )
      next
    }
    candidate <- fit[1L, , drop = FALSE]
    statistic_names <- c("logLik", "AIC", "BIC")
    statistics <- if (all(statistic_names %in% names(candidate))) {
      suppressWarnings(as.numeric(unlist(
        candidate[1L, statistic_names, drop = FALSE],
        use.names = FALSE
      )))
    } else {
      rep(NA_real_, length(statistic_names))
    }
    if (length(statistics) != length(statistic_names) || any(!is.finite(statistics))) {
      had_non_finite_statistics <- TRUE
      statistic_text <- paste0(
        statistic_names,
        "=",
        format(statistics, trim = TRUE),
        collapse = ", "
      )
      failed <- rbind(
        failed,
        data.frame(
          model = model,
          reason = paste0("non-finite fit statistics (", statistic_text, ")."),
          stringsAsFactors = FALSE
        )
      )
      next
    }
    rows[[model]] <- candidate
  }

  if (length(rows) == 0L) {
    if (isTRUE(had_non_finite_statistics)) {
      stop(
        "No candidate distributions produced finite fit statistics. First failure: ",
        failed$reason[1L],
        call. = FALSE
      )
    }
    stop(
      "No candidate distributions could be fit. First failure: ",
      failed$reason[1L],
      call. = FALSE
    )
  }

  rankings <- do.call(rbind, rows)
  rownames(rankings) <- NULL
  rankings$d_AIC <- rankings$AIC - min(rankings$AIC, na.rm = TRUE)
  rankings$d_BIC <- rankings$BIC - min(rankings$BIC, na.rm = TRUE)
  rankings$d_logLik <- max(rankings$logLik, na.rm = TRUE) - rankings$logLik
  rankings <- rankings[order(rankings[[select_by]], rankings$AIC, rankings$BIC, rankings$ml), , drop = FALSE]
  rownames(rankings) <- NULL

  fits <- stats::setNames(rankings$univariateML, rankings$ml)
  parameters <- .distribution_parameter_table(rankings)
  list(
    models = models,
    rankings = rankings,
    parameters = parameters,
    fits = fits,
    selected = rankings[1L, , drop = FALSE],
    selected_model = rankings$ml[1L],
    failed = failed
  )
}

.distribution_empty_failures <- function() {
  data.frame(
    model = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

.distribution_fit_table <- function(
  x,
  component,
  models = NULL,
  n = Inf
) {
  if (identical(component, "rankings")) {
    out <- x$rankings[, c("model", "ml", "AIC", "BIC", "d_AIC", "d_BIC"), drop = FALSE]
  } else if (identical(component, "parameters")) {
    out <- x$parameters
    if (!is.null(models)) {
      out <- out[out$model %in% models, , drop = FALSE]
    }
    out <- out[, c(
      "model",
      "distribution",
      "parameter",
      "estimate",
      "selected",
      "selected_AIC",
      "selected_BIC"
    ), drop = FALSE]
  } else if (identical(component, "goodness")) {
    out <- .distribution_goodness_table(x, models = models)
  } else if (identical(component, "lineage")) {
    out <- x$lineage
    if (is.null(out)) {
      out <- data.frame()
    }
  } else if (identical(component, "lineage_parameters")) {
    out <- x$lineage_parameters
    if (is.null(out)) {
      out <- data.frame()
    } else {
      if (!is.null(models)) {
        out <- out[out$model %in% models, , drop = FALSE]
      }
      out <- out[, c(
        "Lineage",
        "model",
        "distribution",
        "parameter",
        "estimate",
        "selected_AIC",
        "selected_BIC"
      ), drop = FALSE]
    }
  } else {
    stop("Unsupported table component: ", component, call. = FALSE)
  }

  if (is.finite(n)) {
    out <- utils::head(out, n)
  }
  rownames(out) <- NULL
  out
}

.distribution_goodness_table <- function(x, models = NULL) {
  boot <- x$bootstrap
  if (is.null(boot) || !is.list(boot) || length(boot) == 0L) {
    return(data.frame(
      model = character(),
      statistic = character(),
      estimate = numeric(),
      conf_low = numeric(),
      conf_high = numeric(),
      conf_level = numeric(),
      bootstrap_reps = integer(),
      stringsAsFactors = FALSE
    ))
  }

  if (!is.null(models)) {
    boot <- boot[names(boot) %in% models]
  }

  tables <- lapply(boot, function(entry) {
    kl <- entry$kl
    if (is.null(kl) || is.null(kl$summary)) {
      return(NULL)
    }
    kl$summary
  })
  tables <- Filter(Negate(is.null), tables)

  if (length(tables) == 0L) {
    return(data.frame(
      model = character(),
      statistic = character(),
      estimate = numeric(),
      conf_low = numeric(),
      conf_high = numeric(),
      conf_level = numeric(),
      bootstrap_reps = integer(),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, tables)
  rownames(out) <- NULL
  out
}

.distribution_round_numeric <- function(x, digits = 4L) {
  if (!is.data.frame(x) || nrow(x) == 0L) {
    return(x)
  }
  digits <- as.integer(digits)
  if (!is.finite(digits) || digits < 1L) {
    stop("`digits` must be a positive integer.", call. = FALSE)
  }
  numeric_columns <- vapply(x, is.numeric, logical(1))
  x[numeric_columns] <- lapply(x[numeric_columns], signif, digits = digits)
  rownames(x) <- NULL
  x
}

.distribution_check_unused_dots <- function(dots) {
  if (length(dots) > 0L) {
    names <- names(dots)
    missing_names <- is.na(names) | !nzchar(names)
    names[missing_names] <- paste0("..", which(missing_names))
    stop("Unused argument(s): ", paste(names, collapse = ", "), call. = FALSE)
  }
}

.distribution_plot_check_rate_fit <- function(x, model) {
  if (!inherits(x, "rate_distribution_fit")) {
    stop("`x` must be a `rate_distribution_fit` object.", call. = FALSE)
  }
  if (!model %in% names(x$fits)) {
    stop(
      "`x` does not include a fitted `", model, "` distribution. ",
      "Include `models = \"", model, "\"` in `fit_rate_distribution()`.",
      call. = FALSE
    )
  }
  parameters <- x$parameters[x$parameters$model == model, , drop = FALSE]
  if (!all(c("mu", "sigma") %in% parameters$parameter)) {
    stop("The fitted Gumbel model must include `mu` and `sigma` parameters.", call. = FALSE)
  }
  invisible(TRUE)
}

.distribution_rate_fit_model <- function(x, model = NULL, caller = "rate distribution function") {
  if (!inherits(x, "rate_distribution_fit")) {
    stop("`x` must be a `rate_distribution_fit` object.", call. = FALSE)
  }

  if (is.null(model)) {
    model <- x$selected_model
    if (is.null(model) || length(model) != 1L || is.na(model) || !nzchar(model)) {
      model <- names(x$fits)[[1L]]
    }
  }

  if (!is.character(model) || length(model) != 1L || is.na(model) || !nzchar(model)) {
    stop("`model` must be a single fitted model name.", call. = FALSE)
  }
  if (!model %in% names(x$fits)) {
    stop(
      "`x` does not include a fitted `", model, "` distribution. ",
      "Include `models = \"", model, "\"` in `fit_rate_distribution()` before calling ",
      caller,
      ".",
      call. = FALSE
    )
  }
  model
}

.distribution_check_count <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      x < minimum ||
      x != floor(x) ||
      x > .Machine$integer.max) {
    qualifier <- if (minimum <= 0L) "a non-negative integer" else "a positive integer"
    stop("`", name, "` must be ", qualifier, " in R integer range.", call. = FALSE)
  }
  invisible(TRUE)
}

.distribution_check_seed <- function(seed) {
  if (is.null(seed)) {
    return(invisible(TRUE))
  }
  if (!is.numeric(seed) ||
      length(seed) != 1L ||
      !is.finite(seed) ||
      seed != floor(seed) ||
      seed < -.Machine$integer.max ||
      seed > .Machine$integer.max) {
    stop(
      "`seed` must be NULL or a single finite integer in the R integer range.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.distribution_bootstrap_compute_kl <- function(model, kl) {
  if (is.null(kl)) {
    return(identical(model, "gumbel"))
  }

  if (!is.logical(kl) || length(kl) != 1L || is.na(kl)) {
    stop("`kl` must be NULL or a single logical value.", call. = FALSE)
  }

  if (isTRUE(kl) && !identical(model, "gumbel")) {
    stop("Kullback-Leibler summaries are currently supported only for `model = \"gumbel\"`.", call. = FALSE)
  }

  isTRUE(kl)
}

.distribution_bootstrap_result <- function(x, model) {
  boot <- x$bootstrap
  if (is.null(boot) || !is.list(boot) || is.null(boot[[model]])) {
    return(NULL)
  }
  boot[[model]]
}

.distribution_gumbel_kl <- function(values, fit, boot_parameters, conf_level = 0.95) {
  empirical <- stats::density(values)
  fit_parameters <- as.numeric(fit[c("mu", "sigma")])
  names(fit_parameters) <- c("mu", "sigma")

  estimate <- .distribution_gumbel_kl_value(
    empirical = empirical,
    mu = fit_parameters[["mu"]],
    sigma = fit_parameters[["sigma"]]
  )

  if (!is.matrix(boot_parameters) ||
      nrow(boot_parameters) < 2L ||
      !all(c("mu", "sigma") %in% rownames(boot_parameters))) {
    stop(
      "Precomputed bootstrap draws for `gumbel` must contain `mu` and `sigma` rows.",
      call. = FALSE
    )
  }

  bootstrap_values <- vapply(seq_len(ncol(boot_parameters)), function(i) {
    .distribution_gumbel_kl_value(
      empirical = empirical,
      mu = unname(boot_parameters["mu", i]),
      sigma = unname(boot_parameters["sigma", i])
    )
  }, numeric(1))

  finite_bootstrap <- bootstrap_values[is.finite(bootstrap_values)]
  alpha <- (1 - conf_level) / 2
  interval <- if (length(finite_bootstrap) > 0L) {
    stats::quantile(
      finite_bootstrap,
      probs = c(alpha, 1 - alpha),
      na.rm = TRUE,
      names = FALSE
    )
  } else {
    c(NA_real_, NA_real_)
  }

  out <- list(
    summary = data.frame(
      model = "gumbel",
      statistic = "DKL",
      estimate = estimate,
      conf_low = interval[[1L]],
      conf_high = interval[[2L]],
      conf_level = conf_level,
      bootstrap_reps = length(finite_bootstrap),
      stringsAsFactors = FALSE
    ),
    values = bootstrap_values
  )
  class(out) <- c("rate_distribution_kl", "list")
  out
}

.distribution_gumbel_kl_value <- function(empirical, mu, sigma) {
  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
    return(NA_real_)
  }

  dx <- diff(empirical$x)
  if (length(dx) == 0L || !all(is.finite(dx)) || any(dx <= 0)) {
    return(NA_real_)
  }
  dx <- mean(dx)

  empirical_density <- empirical$y
  theoretical_density <- evd::dgumbel(empirical$x, loc = mu, scale = sigma)
  empirical_area <- sum(empirical_density * dx)
  theoretical_area <- sum(theoretical_density * dx)

  if (!is.finite(empirical_area) || empirical_area <= 0 ||
      !is.finite(theoretical_area) || theoretical_area <= 0) {
    return(NA_real_)
  }

  empirical_density <- empirical_density / empirical_area
  theoretical_density <- theoretical_density / theoretical_area

  keep <- is.finite(empirical_density) &
    is.finite(theoretical_density) &
    empirical_density > 0 &
    theoretical_density > 0

  # nocov start
  if (!any(keep)) {
    return(NA_real_)
  }
  # nocov end

  sum(
    empirical_density[keep] *
      log(empirical_density[keep] / theoretical_density[keep]) *
      dx
  )
}

.distribution_bootstrap_density <- function(bootstrap_result, model, x_values) {
  parameters <- bootstrap_result$parameters
  if (!is.matrix(parameters) ||
      nrow(parameters) < 2L ||
      !all(c("mu", "sigma") %in% rownames(parameters))) {
    stop(
      "Precomputed bootstrap draws for `gumbel` must contain `mu` and `sigma` rows.",
      call. = FALSE
    )
  }

  vapply(seq_len(ncol(parameters)), function(i) {
    evd::dgumbel(
      x_values,
      loc = unname(parameters["mu", i]),
      scale = unname(parameters["sigma", i])
    )
  }, numeric(length(x_values)))
}

.distribution_plot_xlab <- function(x) {
  prefix <- if (identical(x$transformation, "log") ||
      identical(x$transformation, "already_log")) {
    "Log "
  } else {
    ""
  }

  column <- x$source_column
  if (is.null(column) || length(column) != 1L || is.na(column)) {
    column <- ""
  }

  noun <- if (identical(x$source, "bifrost_search$model_no_uncertainty$param")) {
    "fitted BMM regime rate"
  } else if (identical(x$source, "data.frame")) {
    switch(
      column,
      log_lineage_rate = "weighted lineage rate",
      lineage_rate = "weighted lineage rate",
      Weighted_Lineage_Value = "weighted lineage rate",
      tip_rate = "tip rate",
      rate = "rate",
      "rate"
    )
  } else {
    "rate"
  }

  paste0(prefix, noun)
}

.distribution_empty_parameters <- function() {
  data.frame(
    model = character(),
    distribution = character(),
    parameter = character(),
    estimate = numeric(),
    n = integer(),
    logLik = numeric(),
    AIC = numeric(),
    BIC = numeric(),
    d_AIC = numeric(),
    d_BIC = numeric(),
    d_logLik = numeric(),
    selected = logical(),
    selected_AIC = logical(),
    selected_BIC = logical(),
    density = character(),
    support_lower = numeric(),
    support_upper = numeric(),
    stringsAsFactors = FALSE
  )
}

.distribution_parameter_table <- function(rankings) {
  if (is.null(rankings) || nrow(rankings) == 0L || !"univariateML" %in% names(rankings)) {
    return(.distribution_empty_parameters())
  }

  min_aic <- min(rankings$AIC, na.rm = TRUE)
  min_bic <- min(rankings$BIC, na.rm = TRUE)
  selected_model <- rankings$ml[[1L]]

  rows <- lapply(seq_len(nrow(rankings)), function(i) {
    .distribution_fit_parameter_rows(
      fit = rankings$univariateML[[i]],
      ranking = rankings[i, , drop = FALSE],
      selected_model = selected_model,
      min_aic = min_aic,
      min_bic = min_bic
    )
  })

  rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
  if (length(rows) == 0L) {
    return(.distribution_empty_parameters())
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.distribution_fit_parameter_rows <- function(
  fit,
  ranking,
  selected_model,
  min_aic,
  min_bic
) {
  estimates <- .distribution_fit_parameter_estimates(fit)
  if (length(estimates) == 0L) {
    return(.distribution_empty_parameters())
  }

  support <- attr(fit, "support", exact = TRUE)
  if (is.null(support) || length(support) < 2L) {
    support <- c(NA_real_, NA_real_)
  }
  model <- as.character(ranking$ml[[1L]])
  distribution <- attr(fit, "model", exact = TRUE)
  if (is.null(distribution) || length(distribution) == 0L) {
    distribution <- as.character(ranking$model[[1L]])
  }
  density <- attr(fit, "density", exact = TRUE)
  if (is.null(density) || length(density) == 0L) {
    density <- NA_character_
  }
  n <- attr(fit, "n", exact = TRUE)
  if (is.null(n) || length(n) == 0L) {
    n <- NA_integer_
  }
  fit_logLik <- attr(fit, "logLik", exact = TRUE)
  if (is.null(fit_logLik) || length(fit_logLik) == 0L) {
    fit_logLik <- ranking$logLik[[1L]]
  }

  data.frame(
    model = model,
    distribution = as.character(distribution[[1L]]),
    parameter = names(estimates),
    estimate = unname(estimates),
    n = as.integer(n[[1L]]),
    logLik = as.numeric(fit_logLik[[1L]]),
    AIC = as.numeric(ranking$AIC[[1L]]),
    BIC = as.numeric(ranking$BIC[[1L]]),
    d_AIC = as.numeric(ranking$d_AIC[[1L]]),
    d_BIC = as.numeric(ranking$d_BIC[[1L]]),
    d_logLik = as.numeric(ranking$d_logLik[[1L]]),
    selected = identical(model, as.character(selected_model)),
    selected_AIC = isTRUE(all.equal(as.numeric(ranking$AIC[[1L]]), min_aic)),
    selected_BIC = isTRUE(all.equal(as.numeric(ranking$BIC[[1L]]), min_bic)),
    density = as.character(density[[1L]]),
    support_lower = as.numeric(support[[1L]]),
    support_upper = as.numeric(support[[2L]]),
    stringsAsFactors = FALSE
  )
}

.distribution_fit_parameter_estimates <- function(fit) {
  raw <- tryCatch(unclass(fit), error = function(e) numeric())
  values <- suppressWarnings(as.numeric(raw))
  if (length(values) == 0L) {
    return(stats::setNames(numeric(), character()))
  }

  parameter_names <- names(raw)
  if (is.null(parameter_names) || length(parameter_names) != length(values)) {
    parameter_names <- paste0("parameter", seq_along(values))
  }
  parameter_names[is.na(parameter_names) | !nzchar(parameter_names)] <-
    paste0("parameter", which(is.na(parameter_names) | !nzchar(parameter_names)))

  stats::setNames(values, parameter_names)
}

.waiting_time_data <- function(waiting_times,
                               prefer_lineage = FALSE,
                               global_include_root_wait = TRUE) {
  if (inherits(waiting_times, "shift_waiting_times")) {
    if (isTRUE(prefer_lineage) &&
        !is.null(waiting_times$lineage) &&
        nrow(waiting_times$lineage) > 0L) {
      data_frame <- waiting_times$lineage
      source <- "shift_waiting_times$lineage"
    } else if (!is.null(waiting_times$global) && nrow(waiting_times$global) > 0L) {
      data_frame <- waiting_times$global
      source <- "shift_waiting_times$global"
    } else if (!is.null(waiting_times$lineage) && nrow(waiting_times$lineage) > 0L) {
      data_frame <- waiting_times$lineage
      source <- "shift_waiting_times$lineage"
    } else {
      stop("`waiting_times` does not contain any waiting-time rows.", call. = FALSE)
    }
  } else if (.shift_is_bifrost_output(waiting_times)) {
    computed <- shift_waiting_times(waiting_times)
    data_info <- .waiting_time_data(
      computed,
      prefer_lineage = prefer_lineage,
      global_include_root_wait = global_include_root_wait
    )
    data_info$source <- paste0(
      "bifrost_search via ",
      data_info$source
    )
    return(data_info)
  } else if (is.data.frame(waiting_times)) {
    data_frame <- waiting_times
    source <- "data.frame"
  } else if (is.numeric(waiting_times)) {
    values <- as.numeric(waiting_times)
    .waiting_time_validate_values(values)
    return(list(
      values = values,
      data_frame = data.frame(TimeToNext = values, stringsAsFactors = FALSE),
      source = "numeric"
    ))
  } else {
    stop(
      "`waiting_times` must be numeric, a data frame, or a `shift_waiting_times` object.",
      call. = FALSE
    )
  }

  value_column <- .waiting_time_value_column(data_frame)
  if (!isTRUE(prefer_lineage) && !isTRUE(global_include_root_wait)) {
    data_frame <- .waiting_time_drop_root_wait(data_frame)
  }
  values <- as.numeric(data_frame[[value_column]])
  .waiting_time_validate_values(values)
  data_frame$TimeToNext <- values
  list(
    values = values,
    data_frame = data_frame,
    source = source
  )
}

.waiting_time_drop_root_wait <- function(data_frame) {
  if (!"PreviousShift" %in% names(data_frame) || nrow(data_frame) == 0L) {
    return(data_frame)
  }
  previous_shift <- as.character(data_frame$PreviousShift)
  root_rows <- which(!is.na(previous_shift) & previous_shift == "root")
  if (length(root_rows) > 1L) {
    stop("`PreviousShift` identifies multiple root waiting intervals.", call. = FALSE)
  }
  if (length(root_rows) == 0L) {
    return(data_frame)
  }
  data_frame[-root_rows, , drop = FALSE]
}

.waiting_time_value_column <- function(data) {
  candidates <- c("TimeToNext", "waiting_time", "wait", "gap")
  column <- candidates[candidates %in% names(data)][1L]
  if (is.na(column)) {
    stop(
      "Waiting-time data frames must include `TimeToNext`, `waiting_time`, ",
      "`wait`, or `gap`.",
      call. = FALSE
    )
  }
  column
}

.waiting_time_validate_values <- function(values) {
  if (!is.numeric(values) ||
      length(values) == 0L ||
      any(!is.finite(values))) {
    stop(
      "Waiting times must be a non-empty finite numeric vector.",
      call. = FALSE
    )
  }
  if (any(values < 0)) {
    stop("Waiting times cannot be negative.", call. = FALSE)
  }
  if (any(values == 0)) {
    stop(
      "Waiting times include zero waiting times from tied or simultaneous events; ",
      "the available continuous distributions require strictly positive values.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.waiting_time_lineage_column <- function(data) {
  candidates <- c("Lineage", "lineage", "tip_label", "tip.label")
  column <- candidates[candidates %in% names(data)][1L]
  if (is.na(column)) {
    stop(
      "`by_lineage = TRUE` requires a lineage column named `Lineage`, ",
      "`lineage`, `tip_label`, or `tip.label`.",
      call. = FALSE
    )
  }
  column
}

.waiting_time_fit_by_lineage <- function(data, models, select_by = "AIC") {
  lineage_column <- .waiting_time_lineage_column(data)
  split_data <- split(data$TimeToNext, data[[lineage_column]])
  summary_rows <- vector("list", length(split_data))
  fits <- vector("list", length(split_data))
  names(summary_rows) <- names(split_data)
  names(fits) <- names(split_data)

  for (lineage in names(split_data)) {
    values <- split_data[[lineage]]
    values <- values[is.finite(values) & values > 0]
    if (length(values) < 2L) {
      summary_rows[[lineage]] <- data.frame(
        Lineage = lineage,
        n = length(values),
        fit_success = FALSE,
        best_AIC = NA_character_,
        best_BIC = NA_character_,
        exp_lambda = NA_real_,
        error = "Fewer than two positive waiting times.",
        stringsAsFactors = FALSE
      )
      next
    }

    fit <- tryCatch(
      .distribution_fit_models(values, models = models, select_by = select_by),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      summary_rows[[lineage]] <- data.frame(
        Lineage = lineage,
        n = length(values),
        fit_success = FALSE,
        best_AIC = NA_character_,
        best_BIC = NA_character_,
        exp_lambda = NA_real_,
        error = conditionMessage(fit),
        stringsAsFactors = FALSE
      )
      next
    }

    fits[[lineage]] <- fit
    summary_rows[[lineage]] <- data.frame(
      Lineage = lineage,
      n = length(values),
      fit_success = TRUE,
      best_AIC = fit$rankings$ml[which.min(fit$rankings$AIC)],
      best_BIC = fit$rankings$ml[which.min(fit$rankings$BIC)],
      exp_lambda = .waiting_time_exp_lambda(fit),
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL
  fits <- fits[!vapply(fits, is.null, logical(1))]
  parameters <- .waiting_time_lineage_parameters(fits)
  list(summary = summary, fits = fits, parameters = parameters)
}

.waiting_time_lineage_parameters <- function(fits) {
  if (length(fits) == 0L) {
    out <- .distribution_empty_parameters()
    return(data.frame(Lineage = character(), out, check.names = FALSE))
  }

  rows <- lapply(names(fits), function(lineage) {
    parameters <- fits[[lineage]]$parameters
    if (is.null(parameters) || nrow(parameters) == 0L) {
      return(NULL)
    }
    cbind(Lineage = lineage, parameters)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    out <- .distribution_empty_parameters()
    return(data.frame(Lineage = character(), out, check.names = FALSE))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.waiting_time_exp_lambda <- function(fit) {
  if (!"exp" %in% names(fit$fits)) {
    return(NA_real_)
  }
  exp_fit <- fit$fits[["exp"]]
  if (!"rate" %in% names(exp_fit)) {
    return(NA_real_)
  }
  as.numeric(exp_fit[["rate"]])
}

.waiting_time_lambda_summary <- function(lambda) {
  lambda <- lambda[is.finite(lambda)]
  if (length(lambda) == 0L) {
    return(data.frame(
      n = 0L,
      mean = NA_real_,
      median = NA_real_,
      sd = NA_real_,
      hdi_95_low = NA_real_,
      hdi_95_high = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  hdi <- .waiting_time_hdi(lambda, prob = 0.95)
  data.frame(
    n = length(lambda),
    mean = mean(lambda),
    median = stats::median(lambda),
    sd = if (length(lambda) > 1L) stats::sd(lambda) else NA_real_,
    hdi_95_low = hdi[["low"]],
    hdi_95_high = hdi[["high"]],
    stringsAsFactors = FALSE
  )
}

.waiting_time_hdi <- function(x, prob = 0.95) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (n == 0L) {
    return(c(low = NA_real_, high = NA_real_))
  }
  if (n == 1L) {
    return(c(low = x[[1L]], high = x[[1L]]))
  }

  window_size <- max(1L, ceiling(prob * n))
  if (window_size >= n) {
    return(c(low = x[[1L]], high = x[[n]]))
  }

  starts <- seq_len(n - window_size + 1L)
  ends <- starts + window_size - 1L
  widths <- x[ends] - x[starts]
  best <- starts[which.min(widths)]
  c(low = x[[best]], high = x[[best + window_size - 1L]])
}
