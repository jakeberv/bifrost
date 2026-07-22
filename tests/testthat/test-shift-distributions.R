# tests/testthat/test-shift-distributions.R

.shift_distributions_skip_tree_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  testthat::skip_if_not_installed("progressr")
}

.shift_distributions_skip_fit_deps <- function() {
  testthat::skip_if_not_installed("univariateML")
}

.shift_distributions_fixture <- function() {
  tr <- ape::read.tree(text = "(((a:1,b:1):1,c:2):1,d:3);")
  tr <- phytools::paintSubTree(
    tree = tr,
    node = ape::Ntip(tr) + 1L,
    state = "0"
  )
  tr <- phytools::paintSubTree(
    tree = tr,
    node = ape::Ntip(tr) + 2L,
    state = "1"
  )
  phytools::paintSubTree(
    tree = tr,
    node = ape::Ntip(tr) + 3L,
    state = "2"
  )
}

.shift_distributions_fit <- function(tree, rates = c("0" = 1, "1" = 4, "2" = 2), model = "BMM") {
  structure(
    list(
      tree_no_uncertainty_untransformed = tree,
      model_no_uncertainty = list(
        model = model,
        call = list(model = model),
        param = rates
      )
    ),
    class = c("bifrost_search", "list")
  )
}

.shift_distributions_values <- function() {
  c("0" = 1, "1" = 4, "2" = 2)
}

.shift_distributions_transitions <- function() {
  data.frame(
    rate_change = c(
      "root",
      "increase",
      "increase",
      "increase",
      "decrease",
      "decrease",
      "decrease"
    ),
    rate_delta = c(NA, 3, 4, 5, -1, -1.5, -2),
    percentage_change = c(NA, 300, 250, 200, -50, -60, -70),
    log_ratio = c(NA, log(4), log(3.5), log(3), log(0.5), log(0.4), log(0.3)),
    stringsAsFactors = FALSE
  )
}

.shift_distributions_root_only_tree <- function() {
  phytools::paintSubTree(
    ape::read.tree(text = "((a:1,b:1):1,c:2);"),
    node = 4L,
    state = "0"
  )
}

test_that("shift_transitions extracts and labels BMM transitions", {
  .shift_distributions_skip_tree_deps()
  tree <- .shift_distributions_fixture()
  fit <- .shift_distributions_fit(tree)

  out <- shift_transitions(fit)

  testthat::expect_s3_class(out, "shift_transitions")
  testthat::expect_equal(out$rate_change, c("root", "increase", "decrease"))
  testthat::expect_equal(out$node, c(5, 6, 7))
  testthat::expect_equal(out$height, c(0, 1, 2))
  testthat::expect_equal(out$child_state, c("0", "1", "2"))
  testthat::expect_equal(out$child_rate, c(1, 4, 2))
  testthat::expect_equal(out$rate_delta, c(NA, 3, -2))
  testthat::expect_equal(out$percentage_change, c(NA, 300, -50))
  testthat::expect_equal(out$log_ratio, c(NA, log(4), log(0.5)))
  testthat::expect_s3_class(attr(out, "tree"), "phylo")
  testthat::expect_equal(attr(out, "settings")$input_mode, "bifrost")

  generic <- shift_transitions(
    tree = tree,
    state_values = c("0" = 1, "1" = 4, "2" = 2)
  )
  testthat::expect_equal(
    unclass(generic[, names(out)]),
    unclass(out),
    ignore_attr = TRUE
  )

  no_root <- shift_transitions(
    tree = tree,
    state_values = c("0" = 1, "1" = 4, "2" = 2),
    include_root = FALSE
  )
  testthat::expect_equal(no_root$rate_change, c("increase", "decrease"))

  generic_counts <- shift_magnitude_counts(
    tree = tree,
    state_values = c("0" = 1, "1" = 4, "2" = 2)
  )
  testthat::expect_s3_class(generic_counts, "shift_magnitude_counts")
  testthat::expect_equal(generic_counts$source, "input")
  testthat::expect_equal(generic_counts$increase, 1L)
  testthat::expect_equal(generic_counts$decrease, 1L)
})

test_that("shift_transitions rejects within-edge SIMMAP transitions", {
  .shift_distributions_skip_tree_deps()
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  tree <- phytools::paintSubTree(
    tree,
    node = ape::Ntip(tree) + 1L,
    state = "0"
  )
  tree <- phytools::paintSubTree(
    tree,
    node = ape::Ntip(tree) + 2L,
    state = "1",
    anc.state = "0",
    stem = 0.5
  )

  testthat::expect_error(
    shift_transitions(
      tree = tree,
      state_values = c("0" = 1, "1" = 2)
    ),
    "node-based.*within-edge"
  )
})

test_that("shift_transitions keeps numeric tip labels separate from node numbers", {
  .shift_distributions_skip_tree_deps()
  tree <- ape::read.tree(text = "((4:1,2:1):1,3:2);")
  tree <- phytools::paintSubTree(tree, node = 4L, state = "0")
  tree <- phytools::paintSubTree(tree, node = 5L, state = "1")

  transitions <- shift_transitions(
    tree = tree,
    state_values = c("0" = 1, "1" = 2)
  )

  testthat::expect_equal(transitions$node, c(4L, 5L))
  testthat::expect_equal(transitions$child_state, c("0", "1"))
  testthat::expect_equal(transitions$rate_change, c("root", "increase"))
})

test_that("shift_node_marks prepares and draws generic shift annotations", {
  .shift_distributions_skip_tree_deps()
  tree <- .shift_distributions_fixture()
  fit <- .shift_distributions_fit(tree)
  fit$ic_weights <- data.frame(
    node = c(6L, 7L),
    ic_weight_withshift = c(0.95, 0.40)
  )

  marks <- shift_node_marks(fit)

  testthat::expect_s3_class(marks, "shift_node_marks")
  testthat::expect_named(
    marks,
    c("marks", "increase_key", "low_support_nodes", "summary", "settings")
  )
  testthat::expect_equal(marks$summary$directional_shifts, 2L)
  testthat::expect_equal(marks$summary$increases, 1L)
  testthat::expect_equal(marks$summary$decreases, 1L)
  testthat::expect_equal(marks$summary$low_support_shifts, 1L)
  testthat::expect_equal(marks$increase_key$letter, "A")
  testthat::expect_equal(marks$low_support_nodes, 7L)
  malformed <- structure(
    list(tree_no_uncertainty_untransformed = tree),
    class = "list"
  )
  testthat::expect_error(
    shift_node_marks(malformed),
    "^`bifrost_search\\$model_no_uncertainty` must identify a multi-regime BMM fitted model via `\\$model` or `\\$call\\$model`\\.$"
  )
  testthat::expect_equal(marks$settings$letter_order, "child_state")
  testthat::expect_equal(as.data.frame(marks), marks$marks)
  testthat::expect_equal(
    as.data.frame(marks, component = "increase_key"),
    marks$increase_key
  )
  testthat::expect_equal(as.data.frame(marks, component = "summary"), marks$summary)

  generic_marks <- shift_node_marks(
    tree = tree,
    state_values = c("0" = 1, "1" = 4, "2" = 2)
  )
  testthat::expect_s3_class(generic_marks, "shift_node_marks")
  testthat::expect_true(all(is.na(generic_marks$marks$ic_weight_withshift)))

  prepared_again <- shift_node_marks(marks)
  testthat::expect_s3_class(prepared_again, "shift_node_marks")
  testthat::expect_equal(prepared_again$marks, marks$marks)

  ordering_transitions <- data.frame(
    node = c(12L, 10L, 11L),
    rate_change = rep("increase", 3L),
    child_state = c("2", "3", "1"),
    age = c(3, 2, 1),
    percentage_change = c(100, 400, 900),
    stringsAsFactors = FALSE
  )
  child_order <- shift_node_marks(ordering_transitions)
  node_order <- shift_node_marks(
    ordering_transitions,
    letter_order = "node",
    marker_max_cex = 2
  )
  age_order <- shift_node_marks(
    ordering_transitions,
    letter_order = "age"
  )
  testthat::expect_equal(child_order$increase_key$node, c(11L, 12L, 10L))
  testthat::expect_equal(node_order$increase_key$node, c(10L, 11L, 12L))
  testthat::expect_equal(age_order$increase_key$node, c(11L, 10L, 12L))
  testthat::expect_true(all(node_order$marks$node_cex <= 2))

  png_path <- tempfile(fileext = ".png")
  grDevices::png(png_path)
  on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
  utils::capture.output(plot(tree, fsize = 0))
  drawn <- plot(
    marks,
    rate_changes = "increase",
    show_low_support = FALSE,
    show_legend = FALSE
  )
  grDevices::dev.off()

  testthat::expect_s3_class(drawn, "shift_node_marks")
  testthat::expect_true(file.exists(png_path))
  testthat::expect_true(file.info(png_path)$size > 0)
})

test_that("shift_node_marks is the only public shift-node mark entry point", {
  testthat::expect_true("shift_node_marks" %in% getNamespaceExports("bifrost"))
  testthat::expect_false("plot_shift_node_marks" %in% getNamespaceExports("bifrost"))
  testthat::expect_false(exists(
    "plot_shift_node_marks",
    envir = asNamespace("bifrost"),
    inherits = FALSE
  ))
})

test_that("compare_shift_magnitudes duplicates Berv et al. (2026)-style magnitude test", {
  transitions <- data.frame(
    rate_change = c(
      "root",
      "increase",
      "increase",
      "increase",
      "decrease",
      "decrease",
      "decrease"
    ),
    rate_delta = c(NA, 3, 4, 5, -1, -1.5, -2),
    percentage_change = c(NA, 300, 250, 200, -50, -60, -70),
    log_ratio = c(NA, log(4), log(3.5), log(3), log(0.5), log(0.4), log(0.3))
  )

  berv_style <- compare_shift_magnitudes(
    transitions,
    ks_B = 99
  )

  testthat::expect_s3_class(berv_style, "shift_magnitude_comparison")
  testthat::expect_equal(berv_style$settings$measure, "rate_delta")
  testthat::expect_equal(berv_style$settings$transform, "log_absolute")
  testthat::expect_equal(berv_style$settings$tests, "ks")
  testthat::expect_true(berv_style$settings$ks_simulate_p_value)
  testthat::expect_equal(berv_style$values$value, log(abs(transitions$rate_delta[-1L])))
  testthat::expect_equal(berv_style$summary$rate_change, c("increase", "decrease"))
  testthat::expect_equal(berv_style$summary$n, c(3, 3))
  testthat::expect_named(
    berv_style$modal,
    c("mode_increase", "mode_decrease", "mode_difference", "linear_ratio")
  )
  testthat::expect_true(is.finite(berv_style$modal$linear_ratio))
  testthat::expect_equal(berv_style$tests$test, "ks")
  testthat::expect_match(berv_style$tests$method, "Kolmogorov")
  testthat::expect_true(is.finite(berv_style$tests$p_value))
  testthat::expect_equal(
    as.data.frame(berv_style, component = "ks", analysis = "focal")$analysis,
    "focal"
  )
  testthat::expect_equal(
    as.data.frame(berv_style, component = "summary")$rate_change,
    c("increase", "decrease")
  )

  set.seed(123)
  seed_before <- .Random.seed
  seeded_a <- compare_shift_magnitudes(
    transitions,
    ks_reps = 99,
    bootstrap_p_value = TRUE,
    bootstrap_reps = 5,
    seed = 2026
  )
  testthat::expect_equal(.Random.seed, seed_before)
  seeded_b <- compare_shift_magnitudes(
    transitions,
    ks_reps = 99,
    bootstrap_p_value = TRUE,
    bootstrap_reps = 5,
    seed = 2026
  )
  testthat::expect_equal(seeded_a$settings$seed, 2026L)
  testthat::expect_equal(seeded_a$tests$p_value, seeded_b$tests$p_value)
  testthat::expect_equal(
    seeded_a$tests$bootstrap_mean_p,
    seeded_b$tests$bootstrap_mean_p
  )
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  compare_shift_magnitudes(transitions, ks_reps = 99, seed = 11)
  testthat::expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

  optional_tests <- compare_shift_magnitudes(
    transitions,
    transform = "absolute",
    tests = c("t", "wilcox", "ks"),
    ks_simulate_p_value = FALSE
  )
  testthat::expect_equal(optional_tests$summary$mean_abs_raw, c(4, 1.5))
  testthat::expect_equal(optional_tests$tests$test, c("t", "wilcox", "ks"))

  bootstrapped <- compare_shift_magnitudes(
    transitions,
    ks_reps = 9,
    bootstrap_p_value = TRUE,
    bootstrap_reps = 3
  )
  testthat::expect_equal(bootstrapped$settings$ks_B, 9L)
  testthat::expect_true(bootstrapped$settings$bootstrap_p_value)
  testthat::expect_equal(bootstrapped$settings$bootstrap_R, 3L)
  testthat::expect_equal(bootstrapped$tests$bootstrap_R, 3L)
  testthat::expect_true(is.finite(bootstrapped$tests$bootstrap_mean_p))

  pooled_runs <- compare_shift_magnitudes(
    list(gic = transitions, bic = transitions),
    ks_B = 99
  )
  testthat::expect_equal(pooled_runs$summary$n, c(6, 6))
  testthat::expect_setequal(pooled_runs$values$source, c("gic", "bic"))

  grouped_inputs <- shift_magnitude_groups(
    combined = list(gic = transitions, bic = transitions),
    focal = transitions,
    precomputed = berv_style
  )
  testthat::expect_s3_class(grouped_inputs, "shift_magnitude_groups")

  list_wrapped_groups <- shift_magnitude_groups(list(
    A = transitions,
    B = transitions
  ))
  testthat::expect_s3_class(list_wrapped_groups, "shift_magnitude_groups")
  testthat::expect_equal(names(list_wrapped_groups), c("A", "B"))

  grouped_runs <- compare_shift_magnitudes(
    grouped_inputs,
    ks_B = 99
  )
  testthat::expect_s3_class(grouped_runs, "shift_magnitude_comparison_set")
  testthat::expect_equal(names(grouped_runs), c("combined", "focal", "precomputed"))
  testthat::expect_s3_class(grouped_runs$combined, "shift_magnitude_comparison")
  testthat::expect_equal(grouped_runs$combined$summary$n, c(6, 6))
  testthat::expect_equal(grouped_runs$focal$summary$n, c(3, 3))
  testthat::expect_identical(grouped_runs$precomputed, berv_style)
  testthat::expect_equal(attr(grouped_runs, "settings")$ks_B, 99L)

  alias_precomputed <- compare_shift_magnitudes(transitions, ks_reps = 99)
  alias_grouped <- compare_shift_magnitudes(
    shift_magnitude_groups(alias = alias_precomputed),
    ks_B = 99
  )
  testthat::expect_identical(alias_grouped$alias, alias_precomputed)

  bootstrap_alias_precomputed <- compare_shift_magnitudes(
    transitions,
    ks_reps = 9,
    bootstrap_p_value = TRUE,
    bootstrap_reps = 2
  )
  bootstrap_alias_grouped <- compare_shift_magnitudes(
    shift_magnitude_groups(alias = bootstrap_alias_precomputed),
    ks_B = 9,
    bootstrap_p_value = TRUE,
    bootstrap_R = 2
  )
  testthat::expect_identical(bootstrap_alias_grouped$alias, bootstrap_alias_precomputed)

  mismatched_precomputed <- compare_shift_magnitudes(
    transitions,
    measure = "percentage_change",
    transform = "absolute",
    tests = "t"
  )
  testthat::expect_error(
    compare_shift_magnitudes(
      shift_magnitude_groups(mismatch = mismatched_precomputed),
      ks_B = 99
    ),
    "mismatch.*Mismatched setting"
  )
  testthat::expect_error(
    compare_shift_magnitudes(grouped_inputs, tree = list(), ks_B = 99),
    "shift_magnitude_groups"
  )

  grouped_ks <- as.data.frame(grouped_runs, component = "ks")
  testthat::expect_equal(grouped_ks$analysis, c("combined", "focal", "precomputed"))
  testthat::expect_true(all(c("ks_D", "linear_ratio") %in% names(grouped_ks)))
  testthat::expect_equal(
    as.data.frame(grouped_runs, component = "summary")$analysis,
    rep(names(grouped_runs), each = 2L)
  )
  testthat::expect_true("analysis" %in% names(as.data.frame(grouped_runs, component = "tests")))
  testthat::expect_true("analysis" %in% names(as.data.frame(grouped_runs, component = "modal")))

  labelled_groups <- compare_shift_magnitudes(
    shift_magnitude_groups(
      .list = list(transitions, transitions),
      labels = c("A", "B")
    ),
    ks_B = 99
  )
  testthat::expect_equal(names(labelled_groups), c("A", "B"))
  default_label_group <- compare_shift_magnitudes(
    shift_magnitude_groups(transitions),
    ks_B = 99
  )
  testthat::expect_equal(names(default_label_group), "Analysis 1")

  testthat::expect_error(
    compare_shift_magnitudes(transitions[1:3, ]),
    "At least two finite increase and decrease values"
  )
  testthat::expect_error(
    compare_shift_magnitudes(transitions, tests = "anova"),
    "unsupported"
  )
  testthat::expect_error(
    compare_shift_magnitudes(transitions, tests = ""),
    "non-empty"
  )
  testthat::expect_error(
    compare_shift_magnitudes(transitions, tests = c("ks", "")),
    "blank"
  )
  testthat::expect_error(
    compare_shift_magnitudes(transitions, tests = c(" ", "t")),
    "blank"
  )
  testthat::expect_error(
    compare_shift_magnitudes(transitions, tests = c("all", "anova")),
    "`all`"
  )
  testthat::expect_error(
    compare_shift_magnitudes(
      shift_magnitude_groups(Good = transitions, Bad = transitions[1:3, ]),
      ks_B = 99
    ),
    "analysis group `Bad`"
  )
  testthat::expect_error(
    shift_magnitude_counts(
      shift_magnitude_groups(
        Good = transitions,
        Bad = list(run = list(nested = transitions))
      )
    ),
    "analysis group `Bad`"
  )
  testthat::expect_error(
    shift_magnitude_groups(.list = transitions),
    "non-empty list of analysis groups"
  )
  testthat::expect_error(
    shift_magnitude_groups(),
    "at least one analysis group"
  )
  testthat::expect_error(
    shift_magnitude_groups(transitions, labels = c("A", "B")),
    "one non-empty"
  )
  testthat::expect_error(
    shift_magnitude_groups(
      .list = list(transitions, transitions),
      labels = c("A", "A")
    ),
    "unique"
  )
  testthat::expect_error(
    shift_magnitude_groups(A = transitions, .list = list(B = transitions)),
    "either `...` or `.list`"
  )
})

test_that("shift magnitude count and comparison plot methods use selected objects", {
  transitions_a <- data.frame(
    rate_change = c(
      "root",
      "increase",
      "increase",
      "increase",
      "decrease",
      "decrease",
      "decrease"
    ),
    rate_delta = c(NA, 3, 4, 5, -1, -1.5, -2),
    stringsAsFactors = FALSE
  )
  transitions_b <- data.frame(
    rate_change = c(
      "root",
      "increase",
      "increase",
      "increase",
      "increase",
      "decrease",
      "decrease",
      "decrease"
    ),
    rate_delta = c(NA, 4, 5, 6, 7, -0.8, -1.1, -1.4),
    stringsAsFactors = FALSE
  )
  grouped <- shift_magnitude_groups(
    "GIC + BIC" = list(gic = transitions_a, bic = transitions_b),
    GIC = list(gic = transitions_a)
  )

  counts <- shift_magnitude_counts(grouped[["GIC + BIC"]])
  testthat::expect_s3_class(counts, "shift_magnitude_counts")
  testthat::expect_equal(nrow(counts), 2)
  testthat::expect_equal(names(counts)[seq_len(6)], c(
    "source",
    "decrease",
    "increase",
    "total",
    "decrease_frequency",
    "increase_frequency"
  ))
  testthat::expect_equal(counts$source, c("gic", "bic"))
  testthat::expect_equal(counts$increase, c(3L, 4L))
  testthat::expect_equal(counts$decrease, c(3L, 3L))

  grouped_counts <- shift_magnitude_counts(grouped)
  testthat::expect_s3_class(grouped_counts, "shift_magnitude_count_set")
  testthat::expect_equal(names(grouped_counts), c("GIC + BIC", "GIC"))
  testthat::expect_s3_class(grouped_counts[["GIC + BIC"]], "shift_magnitude_counts")
  testthat::expect_equal(attr(grouped_counts[["GIC + BIC"]], "analysis"), "GIC + BIC")
  testthat::expect_error(
    shift_magnitude_counts(grouped, tree = list()),
    "shift_magnitude_groups"
  )
  grouped_count_table <- as.data.frame(grouped_counts)
  testthat::expect_equal(grouped_count_table$analysis, c("GIC + BIC", "GIC + BIC", "GIC"))
  testthat::expect_equal(grouped_count_table$source, c("gic", "bic", "gic"))

  labeled_counts <- shift_magnitude_counts(
    shift_magnitude_groups(
      .list = unclass(grouped),
      labels = c("Combined", "GIC only")
    )
  )
  testthat::expect_equal(names(labeled_counts), c("Combined", "GIC only"))

  default_labels <- shift_magnitude_counts(shift_magnitude_groups(transitions_a))
  testthat::expect_equal(names(default_labels), "Analysis 1")
  testthat::expect_error(
    shift_magnitude_groups(.list = transitions_a),
    "list of analysis groups"
  )
  testthat::expect_error(
    shift_magnitude_counts(list(A = list(run = transitions_a))),
    "shift_magnitude_groups"
  )

  count_path <- tempfile(fileext = ".png")
  grDevices::png(count_path)
  on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
  count_plot <- plot(
    grouped_counts[["GIC + BIC"]],
    statistic = "frequency",
    main = "GIC + BIC"
  )
  grDevices::dev.off()

  testthat::expect_true(file.exists(count_path))
  testthat::expect_true(file.info(count_path)$size > 0)
  testthat::expect_named(count_plot, c("data", "tests", "settings"))
  testthat::expect_equal(nrow(count_plot$data), 2)
  testthat::expect_equal(count_plot$settings$analysis, "GIC + BIC")
  testthat::expect_equal(count_plot$settings$statistic, "frequency")
  testthat::expect_true(is.finite(count_plot$tests$wilcox_statistic[1]))
  testthat::expect_error(plot(grouped_counts), "select one")
  testthat::expect_error(plot.shift_magnitude_counts(list()), "shift_magnitude_counts")

  comparisons <- compare_shift_magnitudes(grouped, ks_reps = 9)
  testthat::expect_error(plot(comparisons), "select one")

  distribution_path <- tempfile(fileext = ".png")
  grDevices::png(distribution_path)
  distribution_plot <- plot(
    comparisons[["GIC + BIC"]],
    annotate = FALSE,
    legend = FALSE
  )
  grDevices::dev.off()

  testthat::expect_true(file.exists(distribution_path))
  testthat::expect_true(file.info(distribution_path)$size > 0)
  testthat::expect_s3_class(distribution_plot, "shift_magnitude_density_plot")
  testthat::expect_named(
    distribution_plot,
    c(
      "comparison",
      "summary",
      "tests",
      "modal",
      "values",
      "densities",
      "settings"
    )
  )
  testthat::expect_s3_class(distribution_plot$comparison, "shift_magnitude_comparison")
  testthat::expect_equal(distribution_plot$summary$rate_change, c("increase", "decrease"))
  testthat::expect_equal(distribution_plot$comparison$settings$ks_B, 9L)
  testthat::expect_true(distribution_plot$settings$scale_by_frequency)

  distribution_method_path <- tempfile(fileext = ".png")
  grDevices::png(distribution_method_path)
  method_plot <- plot(
    comparisons[["GIC + BIC"]],
    scale_by_frequency = FALSE,
    annotate = FALSE,
    legend = FALSE,
    rug = TRUE,
    show_modes = FALSE,
    show_sample_size = FALSE,
    show_test_label = FALSE,
    show_modal_ratio = FALSE
  )
  grDevices::dev.off()
  testthat::expect_true(file.exists(distribution_method_path))
  testthat::expect_false(method_plot$settings$scale_by_frequency)
  testthat::expect_false(method_plot$settings$show_modes)
  testthat::expect_false(method_plot$settings$show_sample_size)
  testthat::expect_false(method_plot$settings$show_test_label)
  testthat::expect_false(method_plot$settings$show_modal_ratio)

  count_label <- NULL
  testthat::local_mocked_bindings(
    text = function(..., labels) {
      count_label <<- labels
      invisible(NULL)
    },
    .package = "graphics"
  )
  default_plot_path <- file.path(getwd(), "Rplots.pdf")
  unlink(default_plot_path)
  annotation_count_path <- tempfile(fileext = ".png")
  grDevices::png(annotation_count_path)
  on.exit(
    if (grDevices::dev.cur() > 1L) {
      grDevices::dev.off()
    },
    add = TRUE
  )
  graphics::plot(1, 1, type = "n")
  .shift_magnitude_plot_count_annotation(count_plot$tests[1, , drop = FALSE])
  grDevices::dev.off()
  testthat::expect_true(file.exists(annotation_count_path))
  testthat::expect_false(file.exists(default_plot_path))
  testthat::expect_match(count_label, "Wilcoxon V =")
  testthat::expect_match(count_label, "P =")
})

test_that("shift_transitions validates BMM inputs and empty-transition cases", {
  .shift_distributions_skip_tree_deps()
  tree <- .shift_distributions_fixture()
  values <- .shift_distributions_values()

  tip_tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  tip_tree <- phytools::paintSubTree(tip_tree, node = ape::Ntip(tip_tree) + 1L, state = "0")
  tip_tree <- phytools::paintSubTree(tip_tree, node = 1L, state = "1", stem = TRUE)
  tip_transition <- shift_transitions(
    tree = tip_tree,
    state_values = c("0" = 1, "1" = 2)
  )
  testthat::expect_true("a" %in% tip_transition$node_label)

  x_tree_transition <- shift_transitions(tree, state_values = values)
  testthat::expect_equal(x_tree_transition$rate_change, c("root", "increase", "decrease"))

  empty_transitions <- shift_transitions(
    tree = .shift_distributions_root_only_tree(),
    state_values = c("0" = 1),
    include_root = FALSE
  )
  testthat::expect_equal(nrow(empty_transitions), 0L)
  testthat::expect_error(
    shift_transitions(tree, tree = tree, state_values = values),
    "either `x` or `tree`"
  )
  testthat::expect_error(shift_transitions(tree, include_root = NA), "include_root")
  bad_root <- structure(
    list(edge = matrix(c(3L, 1L, 4L, 2L), ncol = 2L, byrow = TRUE)),
    class = "phylo"
  )
  testthat::expect_error(.shift_root_node(bad_root), "unique root")
  testthat::expect_equal(.shift_rate_change(0), "no_change")
})

test_that("shift_waiting_times handles empty histories and scoped summaries", {
  .shift_distributions_skip_tree_deps()
  root_only_tree <- .shift_distributions_root_only_tree()
  values <- .shift_distributions_values()

  empty_waits <- shift_waiting_times(
    tree = root_only_tree,
    state_values = c("0" = 1)
  )
  testthat::expect_equal(nrow(empty_waits$global), 0L)
  testthat::expect_equal(nrow(empty_waits$lineage), 0L)
  testthat::expect_equal(nrow(empty_waits$summaries$global), 0L)
  testthat::expect_null(shift_waiting_times(empty_waits$transitions, scope = "global")$lineage)
  testthat::expect_null(
    shift_waiting_times(empty_waits$transitions, tree = root_only_tree, scope = "lineage")$global
  )
  testthat::expect_null(
    shift_waiting_times(empty_waits$transitions, tree = root_only_tree, scope = "lineage")$acceleration
  )
  transitions_without_tree <- empty_waits$transitions
  attr(transitions_without_tree, "tree") <- NULL
  testthat::expect_error(
    shift_waiting_times(transitions_without_tree, scope = "lineage"),
    "`tree` is required"
  )
  testthat::expect_error(
    shift_waiting_times(empty_waits$transitions, state_values = values),
    "`state_values` is not used"
  )
  testthat::expect_equal(.shift_waiting_category("root", "increase"), "other")

  one_shift_table <- data.frame(
    node = c(1L, 2L),
    height = c(0, 1),
    child_state = c("0", "1"),
    rate_change = c("root", "increase"),
    stringsAsFactors = FALSE
  )
  one_shift <- shift_waiting_times(
    one_shift_table,
    tree = root_only_tree,
    scope = "global"
  )
  testthat::expect_equal(nrow(one_shift$global), 1L)

  separate_tip_shifts <- data.frame(
    node = c(1L, 2L),
    height = c(2, 2),
    child_state = c("1", "2"),
    rate_change = c("increase", "increase"),
    stringsAsFactors = FALSE
  )
  no_lineage_pairs <- .shift_lineage_waiting_times(separate_tip_shifts, root_only_tree)
  no_acceleration_pairs <- .shift_acceleration_waiting_times(separate_tip_shifts, root_only_tree)
  testthat::expect_equal(nrow(no_lineage_pairs$combined), 0L)
  testthat::expect_equal(nrow(no_acceleration_pairs$collapsed), 0L)
})

test_that("shift_node_marks validates minimal and prepared annotations", {
  .shift_distributions_skip_tree_deps()
  tree <- .shift_distributions_fixture()

  minimal_marks <- shift_node_marks(
    data.frame(node = c(10L, 11L), rate_change = c("increase", "decrease"))
  )
  testthat::expect_true(all(c("percentage_change", "child_state", "age") %in% names(minimal_marks$marks)))
  node_label_marks <- shift_node_marks(
    data.frame(node = 12L, node_label = "state-a", rate_change = "increase")
  )
  testthat::expect_equal(node_label_marks$marks$child_state, "state-a")
  rate_marks <- shift_node_marks(
    data.frame(
      node = 12L,
      rate_change = "increase",
      parent_rate = 2,
      child_rate = 4,
      stringsAsFactors = FALSE
    )
  )
  testthat::expect_equal(rate_marks$marks$percentage_change, 100)
  weighted_marks <- shift_node_marks(
    data.frame(
      node = 10L,
      rate_change = "increase",
      percentage_change = 1,
      ic_weight_withshift = 0.2
    ),
    ic_weights = data.frame(node = 10L, ic_weight_withshift = 0.8)
  )
  testthat::expect_equal(weighted_marks$marks$ic_weight_withshift, 0.8)
  prepared_list <- list(marks = minimal_marks$marks, summary = minimal_marks$summary)
  testthat::expect_s3_class(shift_node_marks(prepared_list), "shift_node_marks")
  testthat::expect_error(
    shift_node_marks(minimal_marks$marks, marker_base_cex = NA),
    "marker_base_cex"
  )
  testthat::expect_error(shift_node_marks(minimal_marks, tree = tree), "prepared")
  testthat::expect_error(
    shift_node_marks(minimal_marks$marks, tree = tree),
    "transition table"
  )
  testthat::expect_error(shift_node_marks(tree, transitions = list()), "data frame")
  testthat::expect_error(
    shift_node_marks(tree, transitions = minimal_marks$marks, tree = tree),
    "not used when `transitions`"
  )
  testthat::expect_error(
    shift_node_marks(NULL, transitions = data.frame(node = 1L)),
    "must include columns"
  )
  testthat::expect_error(
    shift_node_marks(minimal_marks$marks, ic_weights = data.frame(node = 10L)),
    "`ic_weights`"
  )
  testthat::expect_error(
    plot(minimal_marks, show_legend = NA),
    "`show_legend` must be TRUE or FALSE"
  )
  testthat::expect_error(
    plot(minimal_marks, marker_alpha = 1.01),
    "`marker_alpha` must be a single numeric value"
  )
  testthat::expect_error(plot.shift_node_marks(list()), "shift_node_marks")
  testthat::expect_error(as.data.frame.shift_node_marks(list()), "shift_node_marks")
  testthat::expect_equal(.shift_node_marks_letters(0), character())
  testthat::expect_equal(tail(.shift_node_marks_letters(28), 2), c("N1", "N2"))
  testthat::expect_s3_class(
    plot(shift_node_marks(data.frame(node = 1L, rate_change = "root"))),
    "shift_node_marks"
  )

  fit <- .shift_distributions_fit(tree)
  fit$ic_weights <- data.frame(
    node = c(6L, 7L),
    ic_weight_withshift = c(0.95, 0.40)
  )
  marks <- shift_node_marks(fit)
  mark_path <- tempfile(fileext = ".png")
  grDevices::png(mark_path)
  on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
  utils::capture.output(plot(tree, fsize = 0))
  plot(marks, show_legend = TRUE, show_low_support = TRUE)
  grDevices::dev.off()
  testthat::expect_true(file.exists(mark_path))

  increase_only <- shift_node_marks(
    data.frame(node = 6L, rate_change = "increase", percentage_change = 100)
  )
  increase_path <- tempfile(fileext = ".png")
  grDevices::png(increase_path)
  utils::capture.output(plot(tree, fsize = 0))
  plot(increase_only, rate_changes = c("decrease", "increase"), show_legend = FALSE)
  grDevices::dev.off()
  testthat::expect_true(file.exists(increase_path))
})

test_that("compare_shift_magnitudes validates extraction and comparison inputs", {
  .shift_distributions_skip_tree_deps()
  tree <- .shift_distributions_fixture()
  transitions <- .shift_distributions_transitions()

  all_tests <- compare_shift_magnitudes(
    transitions,
    transform = "signed",
    tests = "all",
    ks_reps = 9
  )
  testthat::expect_equal(all_tests$tests$test, c("t", "wilcox", "ks"))
  testthat::expect_true(is.na(all_tests$modal$linear_ratio))
  testthat::expect_equal(as.data.frame(all_tests, component = "tests")$test, c("t", "wilcox", "ks"))
  testthat::expect_named(as.data.frame(all_tests, component = "modal"), names(all_tests$modal))
  testthat::expect_error(as.data.frame(all_tests, type = "tests"), "Unused argument.*type")
  t_only <- compare_shift_magnitudes(transitions, tests = "t")
  testthat::expect_error(as.data.frame(t_only, component = "ks"), "KS test")
  unnamed <- compare_shift_magnitudes(list(transitions, transitions), ks_reps = 9)
  testthat::expect_setequal(unnamed$values$source, c("1", "2"))
  testthat::expect_error(
    .shift_magnitude_resolve_transitions(list(run = transitions), tree = tree),
    "list of inputs"
  )
  testthat::expect_error(compare_shift_magnitudes(list(), ks_reps = 9), "at least one")
  testthat::expect_error(
    compare_shift_magnitudes(list(A = list(run = transitions)), ks_reps = 9),
    "shift_magnitude_groups"
  )
  testthat::expect_error(
    compare_shift_magnitudes(list(list(foo = 1)), ks_reps = 9),
    "shift_magnitude_groups"
  )
  testthat::expect_error(compare_shift_magnitudes(transitions, tests = character()), "non-empty")
  testthat::expect_error(compare_shift_magnitudes(transitions, ks_simulate_p_value = NA), "ks_simulate")
  testthat::expect_error(compare_shift_magnitudes(transitions, ks_B = 0), "`ks_B`")
  testthat::expect_error(compare_shift_magnitudes(transitions, bootstrap_p_value = NA), "bootstrap_p_value")
  testthat::expect_error(compare_shift_magnitudes(transitions, bootstrap_R = 0), "`bootstrap_R`")
  testthat::expect_error(compare_shift_magnitudes(transitions, ks_reps = 1.5), "`ks_reps`")
  outside_int <- as.double(.Machine$integer.max) + 1
  testthat::expect_error(
    compare_shift_magnitudes(transitions, ks_B = outside_int),
    "ks_B.*R integer range"
  )
  testthat::expect_error(
    compare_shift_magnitudes(transitions, bootstrap_R = outside_int),
    "bootstrap_R.*R integer range"
  )
  testthat::expect_error(
    compare_shift_magnitudes(transitions, ks_reps = outside_int),
    "ks_reps.*R integer range"
  )
  testthat::expect_error(compare_shift_magnitudes(transitions, seed = 1.5), "`seed`")
  testthat::expect_error(compare_shift_magnitudes(transitions, seed = 3e9), "R integer range")
  testthat::expect_error(
    compare_shift_magnitudes(transitions[, "rate_change", drop = FALSE]),
    "missing required"
  )
  testthat::expect_error(
    compare_shift_magnitudes(data.frame(rate_change = c("increase", "decrease"), rate_delta = c(NA, NA))),
    "No finite"
  )
  testthat::expect_error(compare_shift_magnitudes(transitions, tree = tree), "transition table")
  testthat::expect_equal(.shift_density_mode(numeric()), NA_real_)
  testthat::expect_equal(.shift_density_mode(c(2, 2)), 2)
})

test_that("shift magnitude comparison plots expose stable settings and summaries", {
  transitions <- .shift_distributions_transitions()
  comparison <- compare_shift_magnitudes(
    transitions,
    ks_reps = 9,
    bootstrap_p_value = TRUE,
    bootstrap_reps = 2
  )
  density_path <- tempfile(fileext = ".png")
  grDevices::png(density_path)
  density_result <- plot(
    comparison,
    rug = TRUE,
    legend = TRUE,
    annotate = TRUE,
    scale_by_frequency = FALSE,
    main = "Magnitude comparison"
  )
  grDevices::dev.off()
  testthat::expect_true(file.exists(density_path))
  testthat::expect_s3_class(density_result, "shift_magnitude_density_plot")
  custom_label_path <- tempfile(fileext = ".png")
  grDevices::png(custom_label_path)
  custom_label_result <- plot(comparison, xlab = "Custom magnitude")
  grDevices::dev.off()
  testthat::expect_equal(custom_label_result$settings$xlab, "Custom magnitude")
  testthat::expect_equal(custom_label_result$settings$measure, "rate_delta")
  testthat::expect_equal(custom_label_result$settings$transform, "log_absolute")
  testthat::expect_error(plot.shift_magnitude_comparison(list()), "shift_magnitude_comparison")
  testthat::expect_match(.shift_magnitude_plot_density_label(comparison), "boot mean P")
  testthat::expect_equal(
    .shift_magnitude_plot_density_label(
      comparison,
      show_sample_size = FALSE,
      show_test_label = FALSE,
      show_modal_ratio = FALSE
    ),
    ""
  )
  annotation_path <- tempfile(fileext = ".png")
  grDevices::png(annotation_path)
  graphics::plot(1, 1, type = "n")
  .shift_magnitude_plot_density_annotation(
    comparison,
    show_sample_size = FALSE,
    show_test_label = FALSE,
    show_modal_ratio = FALSE
  )
  grDevices::dev.off()
  testthat::expect_true(file.exists(annotation_path))
  testthat::expect_equal(.shift_magnitude_plot_p_value(NA_real_), "NA")
  testthat::expect_equal(.shift_magnitude_plot_p_value(1e-6), "<0.001")
  testthat::expect_error(.shift_magnitude_plot_density(c(1), list()), "At least two")
  testthat::expect_equal(.shift_magnitude_plot_xlab("percentage_change", "absolute", list()), "abs(percentage_change)")
  testthat::expect_equal(.shift_magnitude_plot_xlab("log_ratio", "signed", list()), "log_ratio")
  testthat::expect_equal(.shift_magnitude_plot_single_title(NULL), "")
  testthat::expect_equal(.shift_magnitude_plot_single_title("GIC"), "GIC")
  testthat::expect_error(.shift_magnitude_plot_single_title(c("A", "B")), "`main`")
  testthat::expect_equal(.shift_magnitude_plot_colors(c("red", "blue")), c(increase = "red", decrease = "blue"))
  testthat::expect_error(.shift_magnitude_plot_colors("red"), "`colors`")
  testthat::expect_error(.shift_check_logical(NA, "flag"), "`flag`")
  testthat::expect_error(.shift_magnitude_plot_check_alpha(2, "alpha"), "`alpha`")
  testthat::expect_error(.shift_magnitude_plot_check_density_args(1), "`density_args`")
  testthat::expect_error(.shift_magnitude_plot_check_density_args(list(1)), "named list")
})

test_that("shift_magnitude_counts summarizes count and grouped inputs", {
  .shift_distributions_skip_tree_deps()
  tree <- .shift_distributions_fixture()
  transitions <- .shift_distributions_transitions()
  comparison <- compare_shift_magnitudes(
    transitions,
    ks_reps = 9,
    bootstrap_p_value = TRUE,
    bootstrap_reps = 2
  )

  counts <- shift_magnitude_counts(list(run1 = transitions, run2 = transitions))
  testthat::expect_s3_class(counts, "shift_magnitude_counts")
  testthat::expect_equal(shift_magnitude_counts(counts)$source, counts$source)
  testthat::expect_false("analysis" %in% names(as.data.frame(counts)))
  testthat::expect_error(shift_magnitude_counts(counts, tree = tree), "count table")
  unnamed_counts <- shift_magnitude_counts(list(comparison, transitions))
  testthat::expect_equal(unnamed_counts$source, c("1", "2"))
  direct_comparison_counts <- shift_magnitude_counts(comparison)
  testthat::expect_equal(direct_comparison_counts$source, "comparison")
  testthat::expect_error(
    shift_magnitude_counts(comparison, tree = tree),
    "shift-magnitude comparison"
  )
  comparison_counts <- shift_magnitude_counts(list(comp = comparison))
  testthat::expect_equal(comparison_counts$increase, 3L)
  testthat::expect_error(shift_magnitude_counts(list(run = transitions), tree = tree), "list of inputs")
  testthat::expect_error(shift_magnitude_counts(list()), "at least one")
  testthat::expect_error(
    shift_magnitude_counts(list(A = list(run = transitions))),
    "shift_magnitude_groups"
  )
  count_groups <- shift_magnitude_counts(
    shift_magnitude_groups(A = list(run = transitions), B = transitions)
  )
  testthat::expect_error(shift_magnitude_counts(count_groups), "multiple count groups")
  testthat::expect_equal(as.data.frame(count_groups)$analysis, c("A", "B"))
  one_count_test <- .shift_magnitude_plot_count_tests(comparison_counts, "frequency")
  testthat::expect_equal(one_count_test$n, 1L)
  testthat::expect_true(is.na(one_count_test$wilcox_p_value))
  testthat::expect_equal(.shift_magnitude_plot_offsets(1L), 0)
  testthat::expect_null(.shift_magnitude_plot_count_annotation(one_count_test))
  count_path <- tempfile(fileext = ".png")
  grDevices::png(count_path)
  count_result <- plot(
    counts,
    statistic = "count",
    boxplot = FALSE,
    paired_lines = FALSE,
    paired_test = FALSE,
    colors = c("black", "orange")
  )
  grDevices::dev.off()
  testthat::expect_equal(count_result$settings$statistic, "count")
  testthat::expect_equal(.shift_magnitude_plot_statistic(NA_real_), "NA")
  testthat::expect_equal(.shift_magnitude_plot_count_values(counts, "count")$increase, counts$increase)
})

test_that("rate-distribution helpers validate rate, bootstrap, and plotting contracts", {
  .shift_distributions_skip_fit_deps()
  testthat::skip_if_not_installed("evd")

  testthat::expect_error(fit_rate_distribution(TRUE, models = "norm"), "`x` must be")
  testthat::expect_error(fit_rate_distribution(c(1, Inf), models = "norm"), "finite")
  testthat::expect_error(fit_rate_distribution(c(0, 1), models = "norm"), "strictly positive")
  testthat::expect_error(fit_rate_distribution(c(1, 2), log = NA, models = "norm"), "`log`")
  testthat::expect_error(fit_rate_distribution(c(1, 2), models = character()), "`models`")
  testthat::expect_error(fit_rate_distribution(c(1, 2), models = c(" ", "")), "at least one")
  default_fit <- fit_rate_distribution(seq(1, 3, length.out = 12), models = NULL)
  testthat::expect_true("gumbel" %in% default_fit$models)
  data_rate <- fit_rate_distribution(data.frame(rate = c(1, 2, 3, 4)), models = "norm")
  testthat::expect_equal(data_rate$source_column, "rate")
  log_only <- data.frame(log_lineage_rate = log(c(0.8, 1.1, 1.4, 1.9, 2.5, 3.1)))
  testthat::expect_error(
    fit_rate_distribution(log_only, log = FALSE, models = "norm"),
    "log = FALSE.*raw rate column.*log_lineage_rate"
  )
  raw <- data.frame(lineage_rate = c(0.8, 1.1, 1.4, 1.9, 2.5, 3.1))
  raw_fit <- fit_rate_distribution(raw, log = FALSE, models = "norm")
  testthat::expect_identical(raw_fit$source_column, "lineage_rate")
  testthat::expect_identical(raw_fit$transformation, "none")
  testthat::expect_equal(raw_fit$data, raw$lineage_rate)
  testthat::expect_error(.rate_distribution_from_data_frame(data.frame(x = 1)), "Rate data frames")
  testthat::expect_error(.distribution_fit_models(c(1, NA), models = "norm"), "`values`")
  testthat::expect_error(.distribution_fit_models(c(1, 2, 3), models = "notamodel"), "No candidate")
  finite_constant_fit <- fit_rate_distribution(
    rep(2, 8),
    log = FALSE,
    models = c("norm", "cauchy")
  )
  testthat::expect_equal(finite_constant_fit$selected_model, "cauchy")
  testthat::expect_equal(finite_constant_fit$rankings$ml, "cauchy")
  testthat::expect_equal(finite_constant_fit$failed$model, "norm")
  testthat::expect_match(finite_constant_fit$failed$reason, "non-finite")
  testthat::expect_error(
    fit_rate_distribution(rep(2, 8), log = FALSE, models = "norm"),
    "No candidate distributions produced finite fit statistics"
  )

  gumbel <- fit_rate_distribution(c(1, 1.4, 1.9, 2.8, 4.1, 6.5), models = "gumbel")
  testthat::expect_equal(nrow(as.data.frame(gumbel, component = "goodness")), 0L)
  boot_no_kl <- bootstrap_rate_distribution(gumbel, model = "gumbel", reps = 3, seed = 11, kl = FALSE)
  testthat::expect_null(boot_no_kl$bootstrap$gumbel$kl)
  local({
    testthat::local_mocked_bindings(
      bootstrapml = function(...) matrix(seq_len(4), nrow = 2),
      .package = "univariateML"
    )
    repaired_names <- bootstrap_rate_distribution(gumbel, model = "gumbel", reps = 2, kl = FALSE)
    testthat::expect_equal(rownames(repaired_names$bootstrap$gumbel$parameters), c("mu", "sigma"))
  })
  testthat::expect_equal(nrow(as.data.frame(boot_no_kl, component = "goodness")), 0L)
  testthat::expect_equal(nrow(as.data.frame(boot_no_kl, component = "goodness", models = "absent")), 0L)
  gumbel_no_boot_list <- gumbel
  gumbel_no_boot_list$bootstrap <- NULL
  boot_repaired <- bootstrap_rate_distribution(gumbel_no_boot_list, model = "gumbel", reps = 2, seed = 12, kl = FALSE)
  testthat::expect_type(boot_repaired$bootstrap, "list")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  if (had_seed) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  bootstrap_rate_distribution(gumbel, model = "gumbel", reps = 2, seed = 13, kl = FALSE)
  testthat::expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  testthat::expect_error(.distribution_bootstrap_compute_kl("norm", TRUE), "Kullback-Leibler")
  testthat::expect_error(.distribution_plot_check_rate_fit(list(), "gumbel"), "`x`")
  testthat::expect_error(.distribution_plot_check_rate_fit(gumbel, "norm"), "does not include")
  norm_fit <- fit_rate_distribution(c(1, 2, 3, 4, 5), models = "norm")
  bad_gumbel <- gumbel
  bad_gumbel$parameters <- bad_gumbel$parameters[bad_gumbel$parameters$parameter != "sigma", , drop = FALSE]
  testthat::expect_error(.distribution_plot_check_rate_fit(bad_gumbel, "gumbel"), "mu.*sigma")
  testthat::expect_error(plot(gumbel, bootstrap = NA), "`bootstrap`")
  plot_path <- tempfile(fileext = ".png")
  grDevices::png(plot_path)
  plot(boot_no_kl, bootstrap = NULL, bootstrap_curves = 0, data_density = FALSE, rug = FALSE)
  grDevices::dev.off()
  testthat::expect_true(file.exists(plot_path))
  broken_boot <- boot_no_kl$bootstrap$gumbel
  rownames(broken_boot$parameters) <- NULL
  testthat::expect_error(.distribution_bootstrap_density(broken_boot, "gumbel", 1:3), "mu")
  testthat::expect_error(
    .distribution_gumbel_kl(c(1, 2, 3), gumbel$fits$gumbel, matrix(1:4, nrow = 2)),
    "mu"
  )
  bad_kl_draws <- matrix(c(0, -1), nrow = 2, dimnames = list(c("mu", "sigma"), NULL))
  bad_kl <- .distribution_gumbel_kl(c(1, 2, 3), gumbel$fits$gumbel, bad_kl_draws)
  testthat::expect_equal(bad_kl$summary$bootstrap_reps, 0L)
  empirical <- stats::density(c(1, 2, 3))
  testthat::expect_true(is.na(.distribution_gumbel_kl_value(empirical, 0, -1)))
  testthat::expect_true(is.na(.distribution_gumbel_kl_value(list(x = 1, y = 1), 0, 1)))
  testthat::expect_true(is.na(.distribution_gumbel_kl_value(list(x = c(1, 2), y = c(0, 0)), 0, 1)))
  no_selected <- gumbel
  no_selected$selected_model <- NULL
  testthat::expect_equal(.distribution_rate_fit_model(no_selected), "gumbel")
  testthat::expect_error(.distribution_rate_fit_model(list()), "`x`")
  testthat::expect_error(.distribution_rate_fit_model(gumbel, NA), "`model`")
  testthat::expect_error(.distribution_fit_table(gumbel, component = "unsupported"), "Unsupported")
  testthat::expect_equal(.distribution_plot_xlab(list(transformation = "none", source_column = NA, source = "numeric")), "rate")
  testthat::expect_equal(
    .distribution_plot_xlab(list(transformation = "none", source_column = "tip_rate", source = "data.frame")),
    "tip rate"
  )
  testthat::expect_equal(
    .distribution_plot_xlab(list(transformation = "none", source_column = "Weighted_Lineage_Value", source = "data.frame")),
    "weighted lineage rate"
  )
  testthat::expect_equal(
    .distribution_plot_xlab(list(transformation = "none", source_column = "lineage_rate", source = "data.frame")),
    "weighted lineage rate"
  )
  testthat::expect_equal(
    .distribution_plot_xlab(list(transformation = "none", source_column = "rate", source = "data.frame")),
    "rate"
  )
  testthat::expect_equal(
    .distribution_plot_xlab(list(transformation = "none", source_column = "unknown", source = "data.frame")),
    "rate"
  )
  testthat::expect_equal(nrow(.distribution_empty_parameters()), 0L)
  testthat::expect_equal(nrow(.distribution_parameter_table(NULL)), 0L)
  ranking <- data.frame(
    ml = "fake",
    model = "fake",
    logLik = -1,
    AIC = 2,
    BIC = 3,
    d_AIC = 0,
    d_BIC = 0,
    d_logLik = 0,
    stringsAsFactors = FALSE
  )
  fake_fit <- structure(c(1, 2), support = c(0, Inf))
  fake_rows <- .distribution_fit_parameter_rows(fake_fit, ranking, "fake", 2, 3)
  testthat::expect_equal(fake_rows$parameter, c("parameter1", "parameter2"))
  no_support_rows <- .distribution_fit_parameter_rows(structure(c(alpha = 1)), ranking, "fake", 2, 3)
  testthat::expect_true(is.na(no_support_rows$support_lower))
  empty_ranking <- ranking
  empty_ranking$univariateML <- I(list(structure(character())))
  testthat::expect_equal(nrow(.distribution_parameter_table(empty_ranking)), 0L)
  testthat::expect_equal(.distribution_fit_parameter_estimates(structure(character())) |> length(), 0L)
  testthat::expect_error(.distribution_round_numeric(data.frame(x = 1), digits = 0), "`digits`")
})

test_that("waiting-time distribution helpers validate pooled and lineage inputs", {
  .shift_distributions_skip_fit_deps()

  numeric_waits <- fit_waiting_time_distribution(c(0.5, 1, 2), models = "exp")
  testthat::expect_equal(numeric_waits$source, "numeric")
  testthat::expect_equal(nrow(as.data.frame(numeric_waits, component = "rankings", n = 1)), 1L)
  testthat::expect_equal(nrow(as.data.frame(numeric_waits, component = "lineage")), 0L)
  testthat::expect_equal(nrow(as.data.frame(numeric_waits, component = "lineage_parameters")), 0L)
  testthat::expect_error(
    fit_waiting_time_distribution(c(0.5, 1, 2), by_lineage = NA, models = "exp"),
    "`by_lineage`"
  )
  lineage_only <- structure(
    list(
      global = .shift_empty_global_waiting(),
      lineage = data.frame(
        Lineage = c("a", "a"),
        TimeToNext = c(1, 2),
        PreviousShift = c("increase", "decrease"),
        stringsAsFactors = FALSE
      )
    ),
    class = c("shift_waiting_times", "list")
  )
  lineage_only_fit <- fit_waiting_time_distribution(lineage_only, models = "exp")
  testthat::expect_equal(lineage_only_fit$source, "shift_waiting_times$lineage")
  testthat::expect_error(
    fit_waiting_time_distribution(structure(list(global = data.frame(), lineage = data.frame()), class = "shift_waiting_times")),
    "does not contain"
  )
  testthat::expect_error(fit_waiting_time_distribution(list(), models = "exp"), "`waiting_times`")
  testthat::expect_equal(.waiting_time_drop_root_wait(data.frame(TimeToNext = 1)), data.frame(TimeToNext = 1))
  testthat::expect_equal(.waiting_time_value_column(data.frame(wait = 1)), "wait")
  testthat::expect_equal(.waiting_time_value_column(data.frame(gap = 1)), "gap")
  testthat::expect_error(.waiting_time_value_column(data.frame(x = 1)), "Waiting-time data frames")
  testthat::expect_equal(.waiting_time_lineage_column(data.frame(tip.label = "a")), "tip.label")
  testthat::expect_error(.waiting_time_lineage_column(data.frame(x = 1)), "lineage column")
  failed_lineage <- .waiting_time_fit_by_lineage(
    data.frame(Lineage = "a", TimeToNext = c(1, 2)),
    models = "notamodel"
  )
  testthat::expect_false(failed_lineage$summary$fit_success)
  testthat::expect_equal(nrow(.waiting_time_lineage_parameters(list())), 0L)
  testthat::expect_equal(nrow(.waiting_time_lineage_parameters(list(a = list(parameters = data.frame())))), 0L)
  testthat::expect_equal(.waiting_time_exp_lambda(list(fits = list(norm = 1))), NA_real_)
  testthat::expect_equal(.waiting_time_exp_lambda(list(fits = list(exp = c(scale = 1)))), NA_real_)
  testthat::expect_equal(.waiting_time_lambda_summary(NA_real_)$n, 0L)
  testthat::expect_equal(.waiting_time_hdi(numeric()), c(low = NA_real_, high = NA_real_))
  testthat::expect_equal(.waiting_time_hdi(2), c(low = 2, high = 2))
  testthat::expect_equal(.waiting_time_hdi(1:100, prob = 0.95), c(low = 1, high = 95))
})

test_that("global root waits are removed by label rather than position", {
  ordered <- data.frame(
    PreviousShift = c("root", "a", "b"),
    TimeToNext = c(9, 2, 4),
    stringsAsFactors = FALSE
  )
  shuffled <- ordered[c(2L, 1L, 3L), , drop = FALSE]

  ordered_out <- .waiting_time_data(
    ordered,
    global_include_root_wait = FALSE
  )
  shuffled_out <- .waiting_time_data(
    shuffled,
    global_include_root_wait = FALSE
  )
  testthat::expect_equal(ordered_out$values, c(2, 4))
  testthat::expect_equal(shuffled_out$values, c(2, 4))

  no_root <- ordered[-1L, , drop = FALSE]
  testthat::expect_equal(
    .waiting_time_data(no_root, global_include_root_wait = FALSE)$values,
    c(2, 4)
  )

  duplicated_root <- rbind(ordered, ordered[1L, , drop = FALSE])
  testthat::expect_error(
    .waiting_time_data(duplicated_root, global_include_root_wait = FALSE),
    "multiple root waiting intervals"
  )
})

test_that("compact avian skeleton sensitivity bundle reproduces pooled shift counts", {
  runs <- .avian_skeleton_compact_sensitivity()
  testthat::expect_s3_class(runs, "bifrost_search_compact_bundle")
  testthat::expect_length(runs, 12)
  testthat::expect_type(attr(runs, "provenance"), "list")
  testthat::expect_match(attr(runs[[1L]], "provenance")$source_file, ".RDS", fixed = TRUE)

  transitions <- lapply(runs, shift_transitions)
  gic <- transitions[grepl("gic$", names(transitions))]
  bic <- transitions[grepl("bic$", names(transitions))]

  pooled_count <- function(x, direction) {
    sum(vapply(
      x,
      function(data) sum(data$rate_change == direction),
      integer(1)
    ))
  }

  testthat::expect_equal(
    c(
      pooled_count(transitions, "increase"),
      pooled_count(transitions, "decrease")
    ),
    c(110L, 278L)
  )
  testthat::expect_equal(
    c(pooled_count(gic, "increase"), pooled_count(gic, "decrease")),
    c(59L, 160L)
  )
  testthat::expect_equal(
    c(pooled_count(bic, "increase"), pooled_count(bic, "decrease")),
    c(51L, 118L)
  )

  compare_no_sim <- function(x) {
    compare_shift_magnitudes(
      x,
      ks_simulate_p_value = FALSE,
      bootstrap_p_value = FALSE
    )
  }

  all_fit <- compare_no_sim(transitions)
  gic_fit <- compare_no_sim(gic)
  bic_fit <- compare_no_sim(bic)

  testthat::expect_equal(
    round(c(
      all_fit$tests$statistic,
      gic_fit$tests$statistic,
      bic_fit$tests$statistic
    ), 3),
    c(0.240, 0.240, 0.234)
  )
  testthat::expect_equal(
    round(c(
      all_fit$modal$linear_ratio,
      gic_fit$modal$linear_ratio,
      bic_fit$modal$linear_ratio
    ), 2),
    c(2.76, 2.53, 2.92)
  )

  grouped_fit <- compare_no_sim(shift_magnitude_groups(
    "GIC + BIC" = transitions,
    GIC = gic,
    BIC = bic
  ))
  grouped_ks <- as.data.frame(grouped_fit, component = "ks", digits = 8)
  testthat::expect_equal(
    round(grouped_ks$mode_increase, 3),
    c(-7.821, -7.866, -7.793)
  )
  testthat::expect_equal(
    round(grouped_ks$mode_decrease, 3),
    c(-8.836, -8.793, -8.865)
  )
})

test_that("shift magnitude input resolution reports unresolved list entries", {
  testthat::expect_error(
    compare_shift_magnitudes(list(bad = 1)),
    "Could not resolve"
  )
  testthat::expect_equal(
    .shift_magnitude_comparison_setting_mismatches(
      list(settings = "legacy"),
      list(log = TRUE, alternative = "two.sided")
    ),
    c("log", "alternative")
  )
})

test_that("shift_transitions and fit_rate_distribution require BMM for bifrost output", {
  .shift_distributions_skip_tree_deps()
  tree <- .shift_distributions_fixture()
  non_bmm <- .shift_distributions_fit(tree, model = "BM1")

  testthat::expect_error(
    shift_transitions(non_bmm),
    "`shift_transitions\\(\\)` requires a multi-regime BMM fit"
  )

  .shift_distributions_skip_fit_deps()
  testthat::expect_error(
    fit_rate_distribution(non_bmm, models = c("norm", "gumbel")),
    "`fit_rate_distribution\\(\\)` requires a multi-regime BMM fit"
  )
  testthat::expect_error(
    fit_waiting_time_distribution(non_bmm, models = "exp"),
    "`shift_transitions\\(\\)` requires a multi-regime BMM fit"
  )
  testthat::expect_error(
    compare_shift_magnitudes(non_bmm),
    "`shift_transitions\\(\\)` requires a multi-regime BMM fit"
  )
})

test_that("shift_waiting_times returns global, lineage, and acceleration summaries", {
  .shift_distributions_skip_tree_deps()
  tree <- .shift_distributions_fixture()

  waits <- shift_waiting_times(
    tree = tree,
    state_values = c("0" = 1, "1" = 4, "2" = 2)
  )

  testthat::expect_s3_class(waits, "shift_waiting_times")
  testthat::expect_equal(waits$global$TimeToNext, c(1, 1))
  testthat::expect_equal(waits$global$Category, c("other", "increase_to_decrease"))
  testthat::expect_equal(nrow(waits$lineage), 2)
  testthat::expect_equal(waits$lineage$tip_label, c("a", "b"))
  testthat::expect_equal(waits$lineage$TimeToNext, c(1, 1))
  testthat::expect_equal(
    waits$summaries$lineage$Category,
    "increase_to_decrease"
  )
  testthat::expect_equal(waits$acceleration$summary$n, c(0, 0))

  increasing_waits <- shift_waiting_times(
    tree = tree,
    state_values = c("0" = 1, "1" = 2, "2" = 4)
  )
  testthat::expect_equal(increasing_waits$acceleration$summary$n, c(2, 1))
  testthat::expect_equal(increasing_waits$acceleration$summary$mean, c(1, 1))

  tied_transitions <- data.frame(
    node = c(5L, 6L),
    height = c(1, 1),
    child_state = c("1", "2"),
    rate_change = c("increase", "decrease"),
    stringsAsFactors = FALSE
  )
  tied_waits <- shift_waiting_times(tied_transitions, scope = "global")
  testthat::expect_equal(tied_waits$global$TimeToNext, 0)
})

test_that("fit_rate_distribution is generic and can fit explicit Gumbel requests", {
  .shift_distributions_skip_fit_deps()
  .shift_distributions_skip_tree_deps()

  values <- c(-1.3, -0.8, -0.2, 0.1, 0.4, 0.9, 1.2)
  out <- fit_rate_distribution(
    values,
    log = FALSE,
    models = c("norm", "gumbel")
  )

  testthat::expect_s3_class(out, "rate_distribution_fit")
  testthat::expect_setequal(out$models, c("norm", "gumbel"))
  testthat::expect_type(out$bootstrap, "list")
  testthat::expect_equal(.distribution_plot_xlab(out), "rate")
  testthat::expect_equal(out$selected_model, out$rankings$ml[which.min(out$rankings$AIC)])
  testthat::expect_true(all(c("AIC", "BIC", "ml", "univariateML") %in% names(out$rankings)))
  testthat::expect_true(all(c("model", "parameter", "estimate", "selected") %in% names(out$parameters)))
  testthat::expect_setequal(
    out$parameters$parameter[out$parameters$model == "gumbel"],
    c("mu", "sigma")
  )
  testthat::expect_equal(
    out$parameters$estimate[
      out$parameters$model == "gumbel" &
        out$parameters$parameter == "mu"
    ],
    as.numeric(out$fits$gumbel["mu"]),
    tolerance = 1e-12
  )

  gumbel <- fit_rate_distribution(
    c(1, 1.4, 1.9, 2.8, 4.1, 6.5),
    models = "gumbel"
  )
  testthat::expect_equal(.distribution_plot_xlab(gumbel), "Log rate")
  testthat::expect_equal(gumbel$models, "gumbel")
  testthat::expect_equal(gumbel$selected_model, "gumbel")
  testthat::expect_equal(
    gumbel$parameters[, c("model", "parameter")],
    data.frame(
      model = c("gumbel", "gumbel"),
      parameter = c("mu", "sigma"),
      stringsAsFactors = FALSE
    )
  )
  testthat::expect_equal(
    as.data.frame(gumbel, component = "parameters")$parameter,
    c("mu", "sigma")
  )
  testthat::expect_equal(
    as.data.frame(gumbel, component = "rankings")$ml,
    "gumbel"
  )

  fit <- .shift_distributions_fit(.shift_distributions_fixture())
  from_bifrost <- fit_rate_distribution(fit, models = c("norm", "gumbel"))
  testthat::expect_equal(from_bifrost$source, "bifrost_search$model_no_uncertainty$param")
  testthat::expect_equal(from_bifrost$original_data, c(1, 4, 2))
  testthat::expect_equal(
    .distribution_plot_xlab(from_bifrost),
    "Log fitted BMM regime rate"
  )

  lineage_fit <- fit_rate_distribution(
    data.frame(log_lineage_rate = values),
    models = c("norm", "gumbel")
  )
  testthat::expect_equal(lineage_fit$source_column, "log_lineage_rate")
  testthat::expect_equal(lineage_fit$transformation, "already_log")
  testthat::expect_equal(
    .distribution_plot_xlab(lineage_fit),
    "Log weighted lineage rate"
  )
})

test_that("fit_rate_distribution excludes candidate fits without finite statistics", {
  .shift_distributions_skip_fit_deps()

  testthat::local_mocked_bindings(
    model_select = function(...) data.frame(ml = "norm", stringsAsFactors = FALSE),
    .package = "univariateML"
  )

  testthat::expect_error(
    fit_rate_distribution(c(1, 2, 3, 4), log = FALSE, models = "norm"),
    "No candidate distributions produced finite fit statistics"
  )
})

test_that("fit_rate_distribution reports candidate fits that return no rows", {
  .shift_distributions_skip_fit_deps()

  testthat::local_mocked_bindings(
    model_select = function(...) data.frame(),
    .package = "univariateML"
  )

  testthat::expect_error(
    fit_rate_distribution(c(1, 2, 3, 4), log = FALSE, models = "norm"),
    "No candidate distributions could be fit.*No fitted model returned"
  )
})

test_that("distribution optional-package guards identify the missing dependency", {
  .shift_distributions_skip_fit_deps()

  check_fitting_dependency <- .distribution_check_univariateML
  environment(check_fitting_dependency) <- list2env(
    list(requireNamespace = function(...) FALSE),
    parent = environment(.distribution_check_univariateML)
  )
  testthat::expect_error(
    check_fitting_dependency(),
    "univariateML.*required for distribution fitting"
  )

  fit <- fit_rate_distribution(
    c(1, 1.4, 1.9, 2.8, 4.1, 6.5),
    models = "gumbel"
  )
  bootstrap_without_fitter <- bootstrap_rate_distribution
  environment(bootstrap_without_fitter) <- list2env(
    list(requireNamespace = function(package, ...) !identical(package, "univariateML")),
    parent = environment(bootstrap_rate_distribution)
  )
  testthat::expect_error(
    bootstrap_without_fitter(fit, model = "gumbel", reps = 2L, kl = FALSE),
    "univariateML.*required to bootstrap"
  )

  bootstrap_without_kl <- bootstrap_rate_distribution
  environment(bootstrap_without_kl) <- list2env(
    list(requireNamespace = function(package, ...) !identical(package, "evd")),
    parent = environment(bootstrap_rate_distribution)
  )
  testthat::expect_error(
    bootstrap_without_kl(fit, model = "gumbel", reps = 2L, kl = TRUE),
    "evd.*required to compute.*Kullback-Leibler"
  )

  plot_without_density <- plot.rate_distribution_fit
  environment(plot_without_density) <- list2env(
    list(requireNamespace = function(...) FALSE),
    parent = environment(plot.rate_distribution_fit)
  )
  testthat::expect_error(
    plot_without_density(fit, model = "gumbel"),
    "evd.*required to plot"
  )
})

test_that("bootstrap_rate_distribution precomputes fitted rate uncertainty", {
  .shift_distributions_skip_fit_deps()
  testthat::skip_if_not_installed("evd")

  gumbel <- fit_rate_distribution(
    c(1, 1.4, 1.9, 2.8, 4.1, 6.5),
    models = "gumbel"
  )
  booted <- bootstrap_rate_distribution(
    gumbel,
    model = "gumbel",
    reps = 5,
    seed = 2026
  )

  testthat::expect_s3_class(booted, "rate_distribution_fit")
  testthat::expect_s3_class(booted$bootstrap$gumbel, "rate_distribution_bootstrap")
  testthat::expect_equal(booted$bootstrap$gumbel$model, "gumbel")
  testthat::expect_equal(booted$bootstrap$gumbel$reps, 5)
  testthat::expect_equal(booted$bootstrap$gumbel$seed, 2026)
  testthat::expect_equal(dim(booted$bootstrap$gumbel$parameters), c(2L, 5L))
  testthat::expect_setequal(
    rownames(booted$bootstrap$gumbel$parameters),
    c("mu", "sigma")
  )
  testthat::expect_s3_class(booted$bootstrap$gumbel$kl, "rate_distribution_kl")
  testthat::expect_equal(booted$bootstrap$gumbel$kl$summary$statistic, "DKL")
  testthat::expect_true(is.finite(booted$bootstrap$gumbel$kl$summary$estimate))
  testthat::expect_equal(
    as.data.frame(booted, component = "goodness")$bootstrap_reps,
    5L
  )

  testthat::expect_error(
    bootstrap_rate_distribution(gumbel, model = "norm", reps = 5),
    "does not include a fitted `norm` distribution"
  )
  testthat::expect_error(
    bootstrap_rate_distribution(gumbel, model = "gumbel", reps = 0),
    "`reps` must be a positive integer"
  )
  outside_int <- as.double(.Machine$integer.max) + 1
  testthat::expect_error(
    bootstrap_rate_distribution(gumbel, model = "gumbel", reps = outside_int),
    "reps.*R integer range"
  )
  testthat::expect_error(
    bootstrap_rate_distribution(gumbel, model = "gumbel", seed = 1.5),
    "`seed` must be NULL or a single finite integer"
  )
  testthat::expect_error(
    bootstrap_rate_distribution(gumbel, model = "gumbel", seed = 3e9),
    "R integer range"
  )
  testthat::expect_error(
    bootstrap_rate_distribution(gumbel, model = "gumbel", kl = NA),
    "`kl` must be NULL or a single logical value"
  )
})

test_that("plot.rate_distribution_fit draws the Gumbel overlay", {
  .shift_distributions_skip_fit_deps()
  testthat::skip_if_not_installed("evd")

  gumbel <- fit_rate_distribution(
    c(1, 1.4, 1.9, 2.8, 4.1, 6.5),
    models = "gumbel"
  )

  testthat::expect_error(
    plot(gumbel, bootstrap = TRUE),
    "No precomputed bootstrap draws found"
  )

  gumbel <- bootstrap_rate_distribution(
    gumbel,
    model = "gumbel",
    reps = 5,
    seed = 2026
  )

  outside_int <- as.double(.Machine$integer.max) + 1
  testthat::expect_error(
    plot(gumbel, bootstrap = TRUE, bootstrap_curves = outside_int),
    "bootstrap_curves.*R integer range"
  )

  png_path <- tempfile(fileext = ".png")
  grDevices::png(png_path)
  on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
  plot_result <- plot(
    gumbel,
    bootstrap = TRUE,
    bootstrap_curves = 2
  )
  grDevices::dev.off()

  testthat::expect_true(file.exists(png_path))
  testthat::expect_true(file.info(png_path)$size > 0)
  testthat::expect_named(
    plot_result,
    c("x", "fitted_density", "bootstrap_density", "band", "parameters")
  )
  testthat::expect_equal(ncol(plot_result$bootstrap_density), 5L)
  testthat::expect_equal(nrow(plot_result$band), 2L)

  colored_path <- tempfile(fileext = ".png")
  grDevices::png(colored_path)
  plot(
    gumbel,
    bootstrap = FALSE,
    breaks = seq(min(gumbel$data), max(gumbel$data), length.out = 6L),
    hist_col = grDevices::rainbow(5L)
  )
  grDevices::dev.off()
  testthat::expect_true(file.exists(colored_path))
  testthat::expect_true(file.info(colored_path)$size > 0)
})

test_that("fit_waiting_time_distribution handles pooled and lineage inputs", {
  .shift_distributions_skip_fit_deps()

  waits <- data.frame(
    Lineage = c("a", "a", "a", "b", "b", "b", "c"),
    TimeToNext = c(0.8, 1.1, 1.6, 1.2, 1.9, 2.4, 0.9),
    stringsAsFactors = FALSE
  )

  out <- fit_waiting_time_distribution(
    waits,
    by_lineage = TRUE,
    models = c("exp", "gamma", "weibull")
  )

  testthat::expect_s3_class(out, "waiting_time_distribution_fit")
  testthat::expect_equal(out$selected_model, out$rankings$ml[which.min(out$rankings$AIC)])
  testthat::expect_setequal(
    out$parameters$parameter[out$parameters$model == "exp"],
    "rate"
  )
  testthat::expect_equal(nrow(out$lineage), 3)
  testthat::expect_equal(out$lineage$fit_success, c(TRUE, TRUE, FALSE))
  testthat::expect_true(all(is.finite(out$lineage$exp_lambda[1:2])))
  testthat::expect_true(all(c("Lineage", "model", "parameter", "estimate") %in% names(out$lineage_parameters)))
  testthat::expect_setequal(out$lineage_parameters$Lineage, c("a", "b"))
  testthat::expect_equal(
    out$lineage_parameters$estimate[
      out$lineage_parameters$Lineage == "a" &
        out$lineage_parameters$model == "exp" &
        out$lineage_parameters$parameter == "rate"
    ],
    out$lineage$exp_lambda[out$lineage$Lineage == "a"],
    tolerance = 1e-12
  )
  testthat::expect_equal(out$lambda_summary$n, 2)
  testthat::expect_equal(
    as.data.frame(out, component = "parameters", models = "exp")$parameter,
    "rate"
  )
  testthat::expect_equal(
    as.data.frame(out, component = "lineage_parameters", models = "exp")$model,
    c("exp", "exp")
  )

  testthat::expect_error(
    fit_waiting_time_distribution(c(1, 0, 2), models = "exp"),
    paste0(
      "zero waiting times from tied or simultaneous events.*",
      "continuous distributions require strictly positive values"
    )
  )
  testthat::expect_error(
    fit_waiting_time_distribution(c(1, -1, 2), models = "exp"),
    "Waiting times cannot be negative"
  )
  testthat::expect_error(
    fit_waiting_time_distribution(c(1, Inf, 2), models = "exp"),
    "non-empty finite numeric vector"
  )

  .shift_distributions_skip_tree_deps()
  fit <- .shift_distributions_fit(.shift_distributions_fixture())
  from_bifrost <- fit_waiting_time_distribution(fit, models = "exp")
  testthat::expect_equal(
    from_bifrost$source,
    "bifrost_search via shift_waiting_times$global"
  )
  testthat::expect_equal(from_bifrost$selected_model, "exp")
  malformed <- structure(
    list(tree_no_uncertainty_untransformed = fit$tree_no_uncertainty_untransformed),
    class = "list"
  )
  testthat::expect_error(
    fit_waiting_time_distribution(malformed, models = "exp"),
    "^`waiting_times` must be numeric, a data frame, or a `shift_waiting_times` object\\.$"
  )

  generated_waits <- shift_waiting_times(fit)
  with_root <- fit_waiting_time_distribution(generated_waits, models = "exp")
  without_root <- fit_waiting_time_distribution(
    generated_waits,
    global_include_root_wait = FALSE,
    models = "exp"
  )
  without_root_frame <- fit_waiting_time_distribution(
    generated_waits$global,
    global_include_root_wait = FALSE,
    models = "exp"
  )
  testthat::expect_equal(length(with_root$data), 2L)
  testthat::expect_equal(length(without_root$data), 1L)
  testthat::expect_equal(length(without_root_frame$data), 1L)
  testthat::expect_true(with_root$global_include_root_wait)
  testthat::expect_false(without_root$global_include_root_wait)
})

test_that("compact avian skeleton focal regressions reproduce manuscript summaries", {
  .shift_distributions_skip_fit_deps()
  .shift_distributions_skip_tree_deps()
  testthat::skip_if_not_installed("evd")

  search <- .avian_skeleton_compact_search()
  testthat::expect_type(attr(search, "provenance"), "list")

  default_rate_fit <- fit_rate_distribution(search)
  testthat::expect_equal(default_rate_fit$selected_model, "gumbel")

  rate_fit <- default_rate_fit
  rate_fit <- bootstrap_rate_distribution(
    rate_fit,
    model = "gumbel",
    reps = 1000,
    seed = 2026
  )
  gumbel_parameters <- as.data.frame(
    rate_fit,
    component = "parameters", models = "gumbel"
  )
  gumbel_goodness <- as.data.frame(rate_fit, component = "goodness", digits = 6)
  testthat::expect_equal(
    round(gumbel_parameters$estimate[gumbel_parameters$parameter == "mu"], 3),
    -8.838
  )
  testthat::expect_equal(
    round(gumbel_parameters$estimate[gumbel_parameters$parameter == "sigma"], 3),
    0.830
  )
  testthat::expect_equal(round(gumbel_goodness$estimate, 3), 0.031)
  testthat::expect_equal(round(gumbel_goodness$conf_low, 3), 0.013)
  testthat::expect_equal(round(gumbel_goodness$conf_high, 3), 0.142)

  waits <- shift_waiting_times(search)
  default_waiting_fit <- fit_waiting_time_distribution(waits)
  testthat::expect_equal(default_waiting_fit$selected_model, "weibull")

  with_root <- default_waiting_fit
  without_root <- fit_waiting_time_distribution(
    waits,
    models = "weibull",
    global_include_root_wait = FALSE
  )
  with_root_parameters <- as.data.frame(
    with_root,
    component = "parameters", models = "weibull"
  )
  without_root_parameters <- as.data.frame(without_root, component = "parameters")
  testthat::expect_equal(
    round(with_root_parameters$estimate[with_root_parameters$parameter == "shape"], 3),
    0.770
  )
  testthat::expect_equal(
    round(with_root_parameters$estimate[with_root_parameters$parameter == "scale"], 3),
    0.450
  )
  testthat::expect_equal(
    round(without_root_parameters$estimate[without_root_parameters$parameter == "shape"], 3),
    0.803
  )
  testthat::expect_equal(
    round(without_root_parameters$estimate[without_root_parameters$parameter == "scale"], 3),
    0.423
  )

  default_lineage <- fit_waiting_time_distribution(waits, by_lineage = TRUE)
  default_lineage_summary <- as.data.frame(default_lineage, component = "lineage")
  default_best_aic <- table(default_lineage_summary$best_AIC, useNA = "ifany")
  testthat::expect_equal(unname(default_best_aic[["exp"]]), 1603L)
  testthat::expect_equal(sum(default_lineage_summary$fit_success), 1708L)

  lineage_exp <- fit_waiting_time_distribution(
    waits,
    by_lineage = TRUE,
    models = "exp"
  )
  testthat::expect_equal(lineage_exp$lambda_summary$n, 1708L)
  testthat::expect_equal(round(lineage_exp$lambda_summary$mean, 3), 0.201)
  testthat::expect_equal(round(lineage_exp$lambda_summary$hdi_95_low, 3), 0.079)
  testthat::expect_equal(round(lineage_exp$lambda_summary$hdi_95_high, 3), 0.403)
})
