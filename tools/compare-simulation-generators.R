.comparison_generators <- c("original", "empirical")
.comparison_required_packages <- c(
  "ape", "future", "future.apply", "mvMORPH", "phytools", "pkgload",
  "progressr", "RRphylo", "withr"
)

.comparison_integer <- function(value, name, minimum = 1L) {
  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) != 1L || is.na(parsed) || !is.finite(parsed) ||
      parsed != floor(parsed) || parsed < minimum ||
      parsed > .Machine$integer.max) {
    stop("`", name, "` must be an integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(parsed)
}

.comparison_parse_args <- function(args) {
  defaults <- list(
    replicates = 5L,
    output_dir = NULL,
    workers = as.integer(max(
      1L,
      min(4L, parallel::detectCores(logical = FALSE))
    )),
    diagnostic_draws = 100L,
    tree_tip_count = 200L,
    num_shifts = 5L
  )
  value_options <- c(
    "--replicates",
    "--output-dir",
    "--workers",
    "--diagnostic-draws",
    "--tree-tip-count",
    "--num-shifts"
  )
  i <- 1L
  while (i <= length(args)) {
    option <- args[[i]]
    if (!option %in% value_options) {
      stop("Unknown option: ", option, call. = FALSE)
    }
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      stop("Option ", option, " requires a value.", call. = FALSE)
    }
    value <- args[[i + 1L]]
    if (option == "--replicates") {
      defaults$replicates <- .comparison_integer(value, "replicates")
    } else if (option == "--output-dir") {
      defaults$output_dir <- value
    } else if (option == "--workers") {
      defaults$workers <- .comparison_integer(value, "workers")
    } else if (option == "--diagnostic-draws") {
      defaults$diagnostic_draws <- .comparison_integer(
        value,
        "diagnostic_draws"
      )
    } else if (option == "--tree-tip-count") {
      defaults$tree_tip_count <- .comparison_integer(value, "tree_tip_count", 20L)
    } else if (option == "--num-shifts") {
      defaults$num_shifts <- .comparison_integer(value, "num_shifts")
    }
    i <- i + 2L
  }

  defaults
}

.comparison_covariance_row <- function(sigma) {
  correlation <- stats::cov2cor(sigma)
  values <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
  data.frame(
    mean_absolute_correlation = if (nrow(sigma) == 1L) {
      0
    } else {
      mean(abs(correlation[lower.tri(correlation)]))
    },
    effective_dimensionality = sum(values)^2 / sum(values^2),
    condition_number = max(values) / min(values),
    mean_marginal_variance = mean(diag(sigma)),
    stringsAsFactors = FALSE
  )
}

.comparison_build_template <- function() {
  tree_path <- file.path(
    "inst", "extdata", "avian-skeleton", "passerine_bodyplan_tree.tre"
  )
  data_path <- file.path(
    "inst", "extdata", "avian-skeleton", "passerine_bodyplan_data.RDS"
  )
  if (!file.exists(tree_path) || !file.exists(data_path)) {
    stop("Run this script from the bifrost package root.", call. = FALSE)
  }

  tree <- ape::read.tree(tree_path)
  trait_data <- readRDS(data_path)
  trait_data <- as.matrix(trait_data[tree$tip.label, , drop = FALSE])
  skeletal_columns <- setdiff(colnames(trait_data), "vertnet_mass")
  trait_data <- trait_data[, c(skeletal_columns, "vertnet_mass"), drop = FALSE]
  formula_string <- sprintf(
    "trait_data[, 1:%d] ~ trait_data[, %d]",
    length(skeletal_columns),
    ncol(trait_data)
  )

  createSimulationTemplate(
    baseline_tree = tree,
    trait_data = trait_data,
    formula = formula_string,
    response_columns = seq_along(skeletal_columns),
    predictor_columns = ncol(trait_data),
    method = "LL",
    error = TRUE
  )
}

.comparison_covariance_draws <- function(template, generator, n_draws, seed) {
  draw_covariance <- getFromNamespace(
    ".simulation_draw_covariance",
    "bifrost"
  )
  set.seed(seed)
  rows <- lapply(seq_len(n_draws), function(i) {
    draw <- draw_covariance(template, generator, NULL)
    cbind(
      data.frame(
        simulation_generator = generator,
        draw = i,
        covariance_df = draw$covariance_df,
        stringsAsFactors = FALSE
      ),
      .comparison_covariance_row(draw$sigma)
    )
  })
  do.call(rbind, rows)
}

.comparison_search_options <- function() {
  list(
    formula = "trait_data ~ 1",
    min_descendant_tips = 10L,
    shift_acceptance_threshold = 10,
    IC = "GIC",
    num_cores = 1L,
    method = "LL",
    error = TRUE,
    uncertaintyweights_par = FALSE,
    plot = FALSE,
    store_model_fit_history = FALSE,
    verbose = FALSE,
    progress = FALSE
  )
}

.comparison_run_study <- function(template,
                                  generator,
                                  scenario,
                                  options,
                                  seed,
                                  warning_log) {
  warning_messages <- character(0)
  value <- withCallingHandlers(
    tryCatch({
      if (scenario == "null") {
        runFalsePositiveSimulationStudy(
          template = template,
          n_replicates = options$replicates,
          tree_tip_count = options$tree_tip_count,
          simulation_options = list(simulation_generator = generator),
          search_options = .comparison_search_options(),
          num_cores = options$workers,
          seed = seed
        )
      } else {
        runShiftRecoverySimulationStudy(
          template = template,
          n_replicates = options$replicates,
          tree_tip_count = options$tree_tip_count,
          simulation_options = list(
            num_shifts = options$num_shifts,
            min_shift_tips = 10L,
            max_shift_tips = 20L,
            scale_mode = scenario,
            scale_factor_range = c(0.1, 2),
            exclude_range = c(0.5, 1.5),
            simulation_generator = generator,
            integration_power_range = c(0.5, 1.25),
            integration_exclude_range = c(0.8, 1.1),
            buffer = 3L
          ),
          search_options = .comparison_search_options(),
          fuzzy_distance = 2L,
          weighted = FALSE,
          num_cores = options$workers,
          seed = seed
        )
      }
    }, error = function(e) {
      structure(list(error = conditionMessage(e)), class = "comparison_error")
    }),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(warning_messages)) {
    warning_log[[paste(generator, scenario, sep = "/")]] <- warning_messages
  }
  list(value = value, warning_log = warning_log)
}

.comparison_scalar <- function(x, name) {
  value <- x[[name]]
  if (is.null(value) || length(value) == 0L) NA_real_ else as.numeric(value[[1L]])
}

.comparison_per_replicate <- function(study, generator, scenario) {
  if (inherits(study, "comparison_error")) {
    return(data.frame(
      simulation_generator = generator,
      scenario = scenario,
      replicate = NA_integer_,
      status = "study_error",
      error = study$error,
      n_true_shifts = if (scenario == "null") 0L else NA_integer_,
      n_inferred_shifts = NA_integer_,
      strict_precision = NA_real_,
      strict_recall = NA_real_,
      strict_f1 = NA_real_,
      fuzzy_precision = NA_real_,
      fuzzy_recall = NA_real_,
      fuzzy_f1 = NA_real_,
      ancestral_mean_absolute_correlation = NA_real_,
      ancestral_effective_dimensionality = NA_real_,
      ancestral_condition_number = NA_real_,
      mean_integration_power = NA_real_,
      mean_variance_scale = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(seq_along(study$simdata), function(i) {
    sim <- study$simdata[[i]]
    search <- study$results[[i]]
    status <- if (is.null(search$error)) "ok" else "error"
    ancestral <- if (scenario == "null") {
      sim$covariance_matrix
    } else {
      sim$VCVs[["ancestral"]]
    }
    covariance_metrics <- .comparison_covariance_row(ancestral)
    if (scenario == "null") {
      strict <- fuzzy <- list(precision = NA_real_, recall = NA_real_, f1 = NA_real_)
      n_true <- 0L
    } else {
      evaluation <- evaluateShiftRecovery(
        simdata = list(sim),
        simresults = list(search),
        fuzzy_distance = 2L,
        weighted = FALSE,
        verbose = FALSE
      )
      strict <- evaluation$strict
      fuzzy <- evaluation$fuzzy
      n_true <- length(sim$shiftNodes)
    }
    integration_powers <- sim$sampledIntegrationPowers
    integration_powers <- integration_powers[is.finite(integration_powers)]
    variance_scales <- sim$sampledVarianceScaleFactors
    variance_scales <- variance_scales[is.finite(variance_scales)]

    data.frame(
      simulation_generator = generator,
      scenario = scenario,
      replicate = i,
      status = status,
      error = if (is.null(search$error)) NA_character_ else as.character(search$error),
      n_true_shifts = n_true,
      n_inferred_shifts = if (identical(status, "ok")) {
        length(search$shift_nodes_no_uncertainty)
      } else {
        NA_integer_
      },
      strict_precision = .comparison_scalar(strict, "precision"),
      strict_recall = .comparison_scalar(strict, "recall"),
      strict_f1 = .comparison_scalar(strict, "f1"),
      fuzzy_precision = .comparison_scalar(fuzzy, "precision"),
      fuzzy_recall = .comparison_scalar(fuzzy, "recall"),
      fuzzy_f1 = .comparison_scalar(fuzzy, "f1"),
      ancestral_mean_absolute_correlation =
        covariance_metrics$mean_absolute_correlation,
      ancestral_effective_dimensionality =
        covariance_metrics$effective_dimensionality,
      ancestral_condition_number = covariance_metrics$condition_number,
      mean_integration_power = if (length(integration_powers)) {
        mean(integration_powers)
      } else {
        NA_real_
      },
      mean_variance_scale = if (length(variance_scales)) {
        mean(variance_scales)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.comparison_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

.comparison_summarize <- function(per_replicate) {
  groups <- split(
    per_replicate,
    interaction(
      per_replicate$simulation_generator,
      per_replicate$scenario,
      drop = TRUE
    )
  )
  rows <- lapply(groups, function(x) {
    data.frame(
      simulation_generator = x$simulation_generator[[1L]],
      scenario = x$scenario[[1L]],
      n_rows = nrow(x),
      n_ok = sum(x$status == "ok"),
      n_error = sum(x$status != "ok"),
      mean_inferred_shifts = .comparison_mean(x$n_inferred_shifts),
      mean_strict_precision = .comparison_mean(x$strict_precision),
      mean_strict_recall = .comparison_mean(x$strict_recall),
      mean_strict_f1 = .comparison_mean(x$strict_f1),
      mean_fuzzy_precision = .comparison_mean(x$fuzzy_precision),
      mean_fuzzy_recall = .comparison_mean(x$fuzzy_recall),
      mean_fuzzy_f1 = .comparison_mean(x$fuzzy_f1),
      mean_ancestral_absolute_correlation =
        .comparison_mean(x$ancestral_mean_absolute_correlation),
      mean_ancestral_effective_dimensionality =
        .comparison_mean(x$ancestral_effective_dimensionality),
      mean_ancestral_condition_number =
        .comparison_mean(x$ancestral_condition_number),
      mean_integration_power = .comparison_mean(x$mean_integration_power),
      mean_variance_scale = .comparison_mean(x$mean_variance_scale),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result[order(result$scenario, result$simulation_generator), , drop = FALSE]
}

.comparison_warning_lines <- function(warning_log) {
  if (length(warning_log) == 0L) {
    return(character(0))
  }

  unlist(lapply(names(warning_log), function(name) {
    paste0("[", name, "] ", warning_log[[name]])
  }), use.names = FALSE)
}

.comparison_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- .comparison_parse_args(args)
  if (is.null(options$output_dir) || !nzchar(options$output_dir)) {
    stop("`--output-dir` is required.", call. = FALSE)
  }
  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(options$output_dir)) {
    stop("Could not create output directory: ", options$output_dir, call. = FALSE)
  }

  required <- .comparison_required_packages
  missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  pkgload::load_all(".", quiet = TRUE)
  withr::local_envvar(c(
    RSTUDIO = "1",
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  ))
  progressr::handlers(global = FALSE)

  message("Building the 2,057-tip passerine calibration template...")
  template <- .comparison_build_template()
  message(
    "Template: ", template$n_tips, " tips, ", template$n_response_traits,
    " response traits, residual df = ", template$residual_df
  )

  covariance_draws <- rbind(
    .comparison_covariance_draws(
      template, "original", options$diagnostic_draws, 9001L
    ),
    .comparison_covariance_draws(
      template, "empirical", options$diagnostic_draws, 9002L
    )
  )

  generators <- .comparison_generators
  scenarios <- c("null", "proportional", "correlation")
  scenario_seeds <- c(null = 1101L, proportional = 2101L, correlation = 3101L)
  studies <- setNames(vector("list", length(generators)), generators)
  warning_log <- list()
  per_replicate <- list()
  row_index <- 0L
  for (generator in generators) {
    studies[[generator]] <- setNames(vector("list", length(scenarios)), scenarios)
    for (scenario in scenarios) {
      message(
        "Running ", options$replicates, " replicate(s): ", generator, "/",
        scenario, " with multisession workers = ", options$workers
      )
      run <- .comparison_run_study(
        template = template,
        generator = generator,
        scenario = scenario,
        options = options,
        seed = scenario_seeds[[scenario]],
        warning_log = warning_log
      )
      studies[[generator]][[scenario]] <- run$value
      warning_log <- run$warning_log
      row_index <- row_index + 1L
      per_replicate[[row_index]] <- .comparison_per_replicate(
        run$value,
        generator,
        scenario
      )
    }
  }
  per_replicate <- do.call(rbind, per_replicate)
  rownames(per_replicate) <- NULL
  summary <- .comparison_summarize(per_replicate)

  utils::write.csv(
    covariance_draws,
    file.path(options$output_dir, "covariance_draws.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    per_replicate,
    file.path(options$output_dir, "per_replicate.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summary,
    file.path(options$output_dir, "comparison_summary.csv"),
    row.names = FALSE
  )
  warning_lines <- .comparison_warning_lines(warning_log)
  writeLines(warning_lines, file.path(options$output_dir, "warnings.txt"))
  saveRDS(
    list(
      options = options,
      template_summary = list(
        n_tips = template$n_tips,
        n_response_traits = template$n_response_traits,
        residual_df = template$residual_df,
        empirical_covariance = .comparison_covariance_row(
          template$residual_covariance
        )
      ),
      covariance_draws = covariance_draws,
      per_replicate = per_replicate,
      summary = summary,
      warnings = warning_log,
      studies = studies
    ),
    file.path(options$output_dir, "comparison_results.rds"),
    compress = "xz"
  )

  print(summary, row.names = FALSE)
  message("Wrote comparison outputs to ", normalizePath(options$output_dir))
  invisible(summary)
}

if (identical(environment(), globalenv()) && !interactive()) {
  .comparison_main()
}
