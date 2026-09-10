#!/usr/bin/env Rscript

# Local-only companion to the manuscript runners. All package work uses the
# installed development bifrost; sourcing the companion only reuses run helpers.
.tuning_path <- local({
  explicit <- Sys.getenv("BIFROST_TUNING_RUNNER_PATH", "")
  if (nzchar(explicit)) normalizePath(explicit, mustWork = TRUE) else {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(arg) != 1L) stop("Set BIFROST_TUNING_RUNNER_PATH when sourcing this file.")
    normalizePath(sub("^--file=", "", arg), mustWork = TRUE)
  }
})
.tuning <- new.env(parent = globalenv())
.tuning$runner_path <- .tuning_path
.tuning$companion_path <- file.path(dirname(.tuning_path), "run-manuscript-simulation-suite.R")
sys.source(.tuning$companion_path, envir = .tuning)

evalq({
  original <- mget(c("parse_cli", "load_development_package", "build_template",
                    "collect_provenance", "set_single_threaded_math",
                    "summarize_config", "aggregate_complete_groups", "run_suite",
                    "simulation_options", "search_job_dataset"))
  suite_name <- "bifrost-vignette-tuning-grid-12d-clades10-40"
  selection_policy <- list(
    max_false_positive_rate = 0.10, max_any_false_positive = 0.05,
    min_evaluable_fraction = 0.50, primary_metric = "fuzzy_balanced_accuracy",
    scenario_weights = c(proportional = 0.5, correlation = 0.5),
    tie_break = "conservative", allow_infeasible = FALSE
  )

  usage <- function() cat(paste0(
    "Portable 12-dimensional vignette tuning grid (local-only runner).\n",
    "Usage: Rscript run-vignette-tuning-grid.R [options]\n",
    "  --execute          Run/resume searches; default only prints the design.\n",
    "  --dry-run          Print the design without loading bifrost or downloading data.\n",
    "  --mode=full        500 replicates per scenario; 18,000 searches (default).\n",
    "  --mode=smoke       1 replicate per scenario; 36 searches, separate outputs.\n",
    "  --cores=N          Single-core replicate workers; defaults to min(63, physical cores).\n",
    "  --only=GROUP_ID    null-n250, proportional-n250-s5, integration-rate-n250-s5.\n",
    "  --job-index=N      Run one numbered replicate (12 sequential searches).\n",
    "  --output-dir=PATH  Separate output root; default vignette-tuning-grid-12d-clades10-40.\n",
    "  --repo=PATH        Optional checkout for provenance/output placement, not package loading.\n",
    "  --data-dir=PATH    Optional offline input override; default bifrost_example_file().\n",
    "  --summarize-only   Rebuild available summaries without fitting or simulating.\n",
    "  --list / --list-jobs / --help\n",
    "Existing results are validated and resumed; --overwrite is deliberately unsupported.\n",
    "Keep run-manuscript-simulation-suite.R next to this file. No scheduler is required.\n"
  ))

  parse_cli <- function(args) {
    summarize_only <- "--summarize-only" %in% args
    if ("--overwrite" %in% args) stop_cli("Use a new output directory instead of --overwrite.")
    if (summarize_only && any(c("--execute", "--dry-run") %in% args)) {
      stop_cli("--summarize-only cannot be combined with execution or dry-run options.")
    }
    cli <- original$parse_cli(setdiff(args, "--summarize-only"))
    if (summarize_only && (!is.null(cli$only) || !is.null(cli$job_index) || cli$list || cli$list_jobs)) {
      stop_cli("--summarize-only summarizes all available complete groups.")
    }
    cli$summarize_only <- summarize_only
    if (is.null(cli$cores)) cli$cores <- min(63L, detect_workers())
    cli
  }

  build_groups <- function(mode = "full") {
    n <- if (mode == "full") 500L else 1L
    groups <- rbind(
      new_group("null-n250", "null", 250L, 0L, 5L, 12L, 500L),
      new_group("proportional-n250-s5", "proportional", 250L, 5L, 5L, 12L, 500L),
      new_group("integration-rate-n250-s5", "integration-rate", 250L, 5L, 5L, 12L, 500L)
    )
    groups$index <- seq_len(nrow(groups))
    groups$n_replicates <- n
    groups$total_searches <- n * groups$searches_per_replicate
    groups
  }

  search_configs <- function(group) {
    configs <- expand.grid(shift_acceptance_threshold = c(10, 20, 30),
                           min_descendant_tips = c(10L, 20L), IC = c("GIC", "BIC"),
                           stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    configs$setting_id <- rep(1:6, 2)
    configs$config_id <- paste0(group$id, "-", tolower(configs$IC), "-d",
                                configs$shift_acceptance_threshold, "-m", configs$min_descendant_tips)
    configs$weighted <- group$scenario != "null"
    configs
  }

  resolve_output_dir <- function(repo, supplied = NULL) {
    path <- if (!is.null(supplied)) path.expand(supplied) else {
      if (is.null(repo)) "vignette-tuning-grid-12d-clades10-40" else file.path("local-cache", "vignette-tuning-grid-12d-clades10-40")
    }
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
      path <- file.path(if (is.null(repo)) getwd() else repo, path)
    }
    normalizePath(path, mustWork = FALSE)
  }

  print_design <- function(repo, data_dir, mode, workers, groups, jobs, selected,
                           output_dir, executing = FALSE) {
    cat("Bifrost vignette tuning grid\n",
        "Mode: ", mode, "\nOutput root: ", output_dir,
        "\nReplicate jobs: ", nrow(jobs), "\nTotal searches: ", sum(jobs$n_searches),
        "\nSelected replicate jobs: ", nrow(selected), "\nSelected searches: ", sum(selected$n_searches),
        "\nWorkers: ", workers, "\nResponse traits: 12; tree tips: 250",
        "\nGIC and BIC; delta IC 10/20/30; minimum clade size 10/20",
        "\nEmpirical calibration error=TRUE; simulated searches error=FALSE",
        "\nFive planted shifts in 10-40-tip clades; true shifts below a search cutoff remain in recovery accounting",
        "\nEach simulated dataset is shared across all 12 settings",
        "\nParallelism: dynamic single-core replicate jobs; searches serial within each worker",
        "\nNull any-false-positive tolerance: 0.05; equal-weight scenario balanced accuracy ranking",
        "\nData inputs: ", if (is.null(data_dir)) "bifrost_example_file() at execution" else data_dir,
        "\n", sep = "")
    print_group_list(groups[groups$id %in% selected$group_id, ])
    if (!executing) cat("No simulations were run. Supply --execute to start.\n")
  }

  set_single_threaded_math <- function() {
    original$set_single_threaded_math()
    current <- getOption("future.globals.maxSize", 0)
    if (!is.numeric(current) || length(current) != 1 || is.na(current)) current <- 0
    options(future.globals.maxSize = max(current, 16 * 1024^3))
  }

  check_evaluator <- function(namespace) {
    evaluator <- get("evaluateShiftRecovery", namespace)
    tree <- ape::read.tree(text = "(((a,b),(c,d)),((e,f),(g,h)));")
    result <- evaluator(
      list(list(paintedTree = tree, shiftNodes = 10L)),
      list(list(shift_nodes_no_uncertainty = NULL, num_candidates = 1L, candidate_nodes = 11L)),
      weighted = FALSE, verbose = FALSE)
    if (!identical(result$n_evaluable_replicates, 1L) ||
        !isTRUE(all.equal(unname(result$counts$fuzzy), c(0, 0, 1, 1)))) {
      stop_cli("Install development bifrost with both zero-shift and candidate-aware recovery fixes.")
    }
    if (!"allow_infeasible" %in% names(formals(get("selectTunedSearchParameters", namespace)))) {
      stop_cli("Installed bifrost must support allow_infeasible=FALSE in the tuning selector.")
    }
    if (!exists(".simulation_tuning_study_seeds", namespace, inherits = FALSE)) {
      stop_cli("Install development bifrost with paired runSearchTuningGrid() support.")
    }
  }

  load_development_package <- function() {
    namespace <- original$load_development_package()
    check_evaluator(namespace)
    namespace
  }

  runtime_identity <- function(inputs, namespace) {
    path <- getNamespaceInfo(namespace, "path")
    files <- file.path(path, c("DESCRIPTION", "NAMESPACE", "R/bifrost.rdb", "R/bifrost.rdx"))
    list(package_md5 = unname(tools::md5sum(files)), data_md5 = unname(tools::md5sum(inputs)),
         R = R.version.string, platform = R.version$platform)
  }

  build_template <- function(inputs, namespace) {
    path <- file.path(active_output, active_mode, "calibration-template.rds")
    identity <- runtime_identity(inputs, namespace)
    if (file.exists(path)) {
      saved <- readRDS(path)
      if (!identical(saved$identity, identity)) stop_cli("Calibration package/data changed; use a new output directory.")
      return(saved$template)
    }
    template <- with_local_seed(5L, function() original$build_template(inputs, namespace))
    atomic_save_rds(list(schema_version = 1L, identity = identity, template = template), path)
    template
  }

  collect_provenance <- function(repo, data_dir, inputs, template, namespace) {
    provenance <- original$collect_provenance(repo, data_dir, inputs, template, namespace)
    provenance$runner_files_sha256 <- vapply(c(runner_path, companion_path),
      function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE), character(1))
    names(provenance$runner_files_sha256) <- basename(c(runner_path, companion_path))
    provenance$source_fingerprint <- digest::digest(list(
      package_data = provenance$source_fingerprint, runner = provenance$runner_files_sha256,
      template = provenance$template_covariance_sha256, R = R.version.string), algo = "sha256")
    path <- file.path(active_output, active_mode, "provenance.rds")
    if (file.exists(path)) {
      saved <- readRDS(path)$provenance
      if (!identical(saved$source_fingerprint, provenance$source_fingerprint)) {
        stop_cli("Package/data/runner changed; use a new output directory.")
      }
      return(saved)
    }
    atomic_save_rds(list(schema_version = 1L, provenance = provenance), path)
    provenance
  }

  simulation_options <- function(group) {
    options <- original$simulation_options(group)
    if (group$scenario != "null") {
      options$min_shift_tips <- 10L
      options$max_shift_tips <- 40L
    }
    options
  }

  suite_contract <- function(groups, mode) list(
    schema_version = 1L, suite = suite_name, mode = mode, groups = groups,
    configs = lapply(seq_len(nrow(groups)), function(i) search_configs(groups[i, ])),
    simulations = lapply(seq_len(nrow(groups)), function(i) simulation_options(groups[i, ])),
    formula = "trait_data ~ 1", method = "LL", error = FALSE, template_error = TRUE,
    response_traits = 12L, fuzzy_distance = 2L, selection_policy = selection_policy,
    seed_scheme = "runSearchTuningGrid-paired-v1", paired_settings = TRUE
  )

  base_search_options <- function(config) list(
    formula = "trait_data ~ 1", min_descendant_tips = config$min_descendant_tips,
    num_cores = 1L, IC = config$IC, shift_acceptance_threshold = config$shift_acceptance_threshold,
    ic_uncertainty_threshold = 10, uncertaintyweights = FALSE,
    uncertaintyweights_par = isTRUE(config$weighted), plot = FALSE,
    store_model_fit_history = FALSE, verbose = FALSE, progress = FALSE, method = "LL", error = FALSE
  )

  job_seeds <- function(group, replicate) {
    namespace <- asNamespace("bifrost")
    scenario <- if (group$scenario == "integration-rate") "correlation" else group$scenario
    study_seed <- get(".simulation_tuning_study_seeds", namespace)(group$seed)[[scenario]]
    with_seed <- get(".simulation_seed_runner", namespace)()
    with_seed(study_seed, {
      simulation <- get(".simulation_replicate_seeds", namespace)(group$n_replicates, TRUE)
      search <- get(".simulation_replicate_seeds", namespace)(group$n_replicates, TRUE)
      list(study = study_seed, simulation = simulation[[replicate]], search = search[[replicate]])
    })
  }

  simulate_job_dataset <- function(group, seed, namespace, template) {
    name <- if (group$scenario == "null") "simulateNullDataset" else "simulateShiftedDataset"
    args <- c(list(template = template, tree_tip_count = group$tree_tip_count), simulation_options(group))
    with_seed <- get(".simulation_seed_runner", namespace)()
    sim <- with_seed(seed, do.call(get(name, namespace), args))
    sim$user_input <- get(".simulation_study_compact_call", namespace)(name, args)
    sim
  }

  search_job_dataset <- function(simulation, group, config, seed, namespace) {
    with_seed <- get(".simulation_seed_runner", namespace)()
    with_seed(seed, original$search_job_dataset(simulation, group, config, seed, namespace))
  }

  validate_search_result <- function(result, true_nodes = integer(0)) {
    candidates <- result$candidate_nodes
    inferred <- result$shift_nodes_no_uncertainty
    if (!"shift_nodes_no_uncertainty" %in% names(result) || !is.numeric(candidates) ||
        anyNA(candidates) || any(!is.finite(candidates)) || any(candidates != floor(candidates)) ||
        anyDuplicated(candidates) || !identical(as.integer(length(candidates)), as.integer(result$num_candidates)) ||
        any(!inferred %in% candidates)) stop_cli("Invalid saved candidate or inferred-node vector.")
    # A true shift below the search cutoff is intentionally allowed. The
    # candidate-aware evaluator retains it in the recovery denominator.
    invisible(TRUE)
  }

  validate_job_wrapper <- function(wrapper, job, group, configs, provenance, design_fingerprint) {
    if (!identical(wrapper$suite, suite_name) || !identical(wrapper$schema_version, 1L) ||
        !identical(wrapper$job, as_scalar_list(job)) || !identical(wrapper$group, as_scalar_list(group)) ||
        !identical(wrapper$search_configs, configs) ||
        !identical(wrapper$design_fingerprint, design_fingerprint) ||
        !identical(wrapper$source_fingerprint, provenance$source_fingerprint) ||
        !identical(names(wrapper$results), configs$config_id)) stop_cli("Saved job provenance/design mismatch: ", job$job_id)
    if (ncol(as.matrix(wrapper$simulation$trait_data)) != 12L ||
        !identical(wrapper$simulation_hash, simulation_content_hash(wrapper$simulation))) stop_cli("Saved simulation mismatch.")
    for (result in wrapper$results) validate_search_result(result, wrapper$simulation$shiftNodes)
    invisible(TRUE)
  }

  run_replicate_job <- function(job, group, configs, output_dir, mode, namespace,
                               template, provenance, design_fingerprint, overwrite = FALSE) {
    tryCatch({
      set_single_threaded_math()
      path <- job_output_path(output_dir, mode, job)
      if (file.exists(path)) {
        validate_job_wrapper(readRDS(path), job, group, configs, provenance, design_fingerprint)
        return(list(status = "skipped", job_id = job$job_id))
      }
      checkpoint <- file.path(output_dir, mode, "checkpoints", group$id, job$job_id)
      key <- list(job = as_scalar_list(job), design = design_fingerprint, source = provenance$source_fingerprint)
      sim_path <- file.path(checkpoint, "simulation.rds")
      seeds <- job_seeds(group, job$replicate)
      if (file.exists(sim_path)) {
        saved <- readRDS(sim_path)
        if (!identical(saved$key, key) || !identical(saved$hash, simulation_content_hash(saved$simulation))) {
          stop_cli("Simulation checkpoint mismatch: ", job$job_id)
        }
        simulation <- saved$simulation
      } else {
        simulation <- simulate_job_dataset(group, seeds$simulation, namespace, template)
        atomic_save_rds(list(schema_version = 1L, key = key, simulation = simulation,
                            hash = simulation_content_hash(simulation)), sim_path)
      }
      started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
      results <- setNames(vector("list", nrow(configs)), configs$config_id)
      for (i in seq_len(nrow(configs))) {
        config <- configs[i, ]
        config_path <- file.path(checkpoint, paste0(config$config_id, ".rds"))
        if (file.exists(config_path)) {
          saved <- readRDS(config_path)
          if (!identical(saved$key, key) || !identical(saved$config, config) ||
              !identical(saved$simulation_hash, simulation_content_hash(simulation))) stop_cli("Search checkpoint mismatch.")
          result <- saved$result
        } else {
          result <- search_job_dataset(simulation, group, config, seeds$search, namespace)
          validate_search_result(result, simulation$shiftNodes)
          atomic_save_rds(list(schema_version = 1L, key = key, config = config,
                              simulation_hash = simulation_content_hash(simulation), result = result), config_path)
        }
        validate_search_result(result, simulation$shiftNodes)
        results[[i]] <- result
        invisible(gc(verbose = FALSE))
      }
      wrapper <- list(schema_version = 1L, suite = suite_name, mode = mode,
        job = as_scalar_list(job), group = as_scalar_list(group), search_configs = configs,
        seeds = seeds, design_fingerprint = design_fingerprint, source_fingerprint = provenance$source_fingerprint,
        provenance = provenance, started_at = started, finished_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
        simulation_hash = simulation_content_hash(simulation), simulation = simulation, results = results)
      validate_job_wrapper(wrapper, job, group, configs, provenance, design_fingerprint)
      atomic_save_rds(wrapper, path)
      message("Completed ", job$job_id)
      list(status = "written", job_id = job$job_id, path = path)
    }, error = function(e) list(status = "error", job_id = job$job_id, error = conditionMessage(e)))
  }

  compare_candidate_and_legacy <- function(evaluator, simdata, results) {
    # Always use the current, candidate-node-aware evaluator.
    evaluator(simdata, results, fuzzy_distance = 2L, weighted = TRUE, verbose = FALSE)
  }

  summarize_config <- function(group, config, wrappers, namespace, provenance, design_fingerprint, mode) {
    summary <- original$summarize_config(group, config, wrappers, namespace, provenance, design_fingerprint, mode)
    summary$suite <- suite_name
    summary$schema_version <- 1L
    summary$summary_row$min_descendant_tips <- config$min_descendant_tips
    summary$summary_row$setting_id <- config$setting_id
    summary$summary_row$study_seed <- job_seeds(group, 1L)$study
    p <- summary$per_replicate
    size_counts <- t(vapply(wrappers, function(wrapper) {
      simulation <- wrapper$simulation
      sizes <- vapply(simulation$shiftNodes, function(node) {
        sum(phytools::getDescendants(simulation$paintedTree, node) <=
              ape::Ntip(simulation$paintedTree))
      }, integer(1))
      c(n_true_shifts_10_19 = sum(sizes >= 10L & sizes <= 19L),
        n_true_shifts_20_40 = sum(sizes >= 20L & sizes <= 40L))
    }, c(n_true_shifts_10_19 = 0L, n_true_shifts_20_40 = 0L)))
    for (field in colnames(size_counts)) {
      summary$per_replicate[[field]] <- size_counts[, field]
      summary$summary_row[[field]] <- sum(size_counts[, field])
    }
    summary$summary_row$evaluable_fraction <- safe_mean(p$n_candidates[p$status == "ok"] > 0)
    summary$summary_row$n_evaluable_replicates <- if (is.null(summary$evaluation)) {
      sum(p$status == "ok" & p$n_candidates > 0)
    } else summary$evaluation$n_evaluable_replicates
    summary
  }

  build_tuning_grid <- function(table, ic) {
    configs <- search_configs(build_groups()[1, ])
    settings <- configs[configs$IC == ic, c("setting_id", "IC", "shift_acceptance_threshold", "min_descendant_tips")]
    for (i in seq_len(nrow(settings))) {
      for (scenario in c("null", "proportional", "integration-rate")) {
        row <- table[table$IC == ic & table$setting_id == settings$setting_id[i] & table$scenario == scenario, ]
        if (nrow(row) != 1L) stop_cli("Incomplete tuning grid; no recommendation yet.")
        prefix <- if (scenario == "integration-rate") "correlation" else scenario
        settings[i, paste0(prefix, "_seed")] <- row$study_seed
        common <- c(completion_rate = "completion_rate", failure_rate = "failure_rate",
                    mean_inferred_shifts = "mean_inferred_shifts", evaluable_fraction = "evaluable_fraction")
        row$failure_rate <- 1 - row$completion_rate
        for (name in names(common)) settings[i, paste0(prefix, "_", name)] <- row[[common[[name]]]]
        if (scenario == "null") {
          settings[i, "null_mean_false_positive_rate"] <- row$mean_false_positive_rate
          settings[i, "null_fraction_any_false_positive"] <- row$fraction_any_false_positive
        } else {
          for (metric in grep("^(strict_|fuzzy_|weighted_fuzzy_f1$)", names(row), value = TRUE)) {
            settings[i, paste0(prefix, "_", metric)] <- row[[metric]]
          }
        }
      }
    }
    rownames(settings) <- NULL
    structure(list(IC = ic, grid = settings[, 1:4], summary_table = settings,
      base_search_options = list(formula = "trait_data ~ 1", method = "LL", error = FALSE),
      null_replicates = unique(table$n_replicates[table$scenario == "null"]),
      recovery_replicates = unique(table$n_replicates[table$scenario != "null"]),
      simulation_generators = c(null = "empirical", proportional = "empirical", correlation = "empirical"),
      fuzzy_distance = 2L, weighted = TRUE, studies = NULL,
      paired_settings = TRUE, store_studies = FALSE,
      study_seeds = stats::setNames(as.integer(settings[1L, c("null_seed", "proportional_seed", "correlation_seed")]),
                                    c("null", "proportional", "correlation")),
      replicate_outputs = "Saved separately under replicates and checkpoints"),
      class = c("bifrost_search_tuning_grid", "list"))
  }

  select_grid <- function(grid, namespace) {
    tryCatch(list(status = "selected", selection = do.call(get("selectTunedSearchParameters", namespace),
      c(list(tuning_grid = grid), selection_policy))), error = function(e) {
        if (!startsWith(conditionMessage(e), "No rankable tuning settings") &&
            !startsWith(conditionMessage(e), "No tuning settings have finite")) stop(e)
        list(status = "no_recommendation", reason = conditionMessage(e), selection = NULL)
      })
  }

  aggregate_complete_groups <- function(groups, jobs, output_dir, mode, namespace, provenance, design_fingerprint) {
    result <- original$aggregate_complete_groups(groups, jobs, output_dir, mode, namespace, provenance, design_fingerprint)
    root <- file.path(output_dir, mode, "summaries")
    path <- file.path(root, "suite-summary.rds")
    if (!file.exists(path)) return(result)
    saved <- readRDS(path)
    saved$suite <- suite_name
    saved$schema_version <- 1L
    saved$selection_policy <- selection_policy
    atomic_save_rds(saved, path, overwrite = TRUE)
    if (!isTRUE(saved$all_groups_complete)) return(result)
    for (ic in c("GIC", "BIC")) {
      grid <- build_tuning_grid(saved$summary_table, ic)
      choice <- if (mode == "smoke") list(status = "smoke_only_no_recommendation", selection = NULL) else select_grid(grid, namespace)
      atomic_save_rds(list(schema_version = 1L, grid = grid, policy = selection_policy,
                          recommendation = choice, provenance = provenance),
                      file.path(root, paste0(tolower(ic), "-tuning.rds")), overwrite = TRUE)
      atomic_write_csv(grid$summary_table, file.path(root, paste0(tolower(ic), "-tuning.csv")))
      message(ic, ": ", choice$status)
    }
    result
  }

  run_suite <- function(cli) {
    active_output <<- resolve_output_dir(find_repo(cli$repo), cli$output_dir)
    active_mode <<- cli$mode
    if (!cli$summarize_only) return(original$run_suite(cli))
    namespace <- load_development_package()
    path <- file.path(active_output, active_mode, "provenance.rds")
    if (!file.exists(path)) stop_cli("No run provenance exists in this output directory.")
    saved <- readRDS(path)$provenance
    if (!file.exists(file.path(active_output, active_mode, "calibration-template.rds"))) {
      stop_cli("No saved calibration template; --summarize-only will not fit one.")
    }
    inputs <- resolve_data_inputs(find_data_dir(cli$data_dir), namespace)
    template <- build_template(inputs, namespace)
    current <- collect_provenance(find_repo(cli$repo), find_data_dir(cli$data_dir), inputs, template, namespace)
    if (!identical(current$source_fingerprint, saved$source_fingerprint)) stop_cli("Package/data/runner changed; cannot silently reassess this run.")
    groups <- build_groups(cli$mode)
    aggregate_complete_groups(groups, build_jobs(groups), active_output, cli$mode, namespace, saved,
                              digest::digest(suite_contract(groups, cli$mode), algo = "sha256"))
  }
}, envir = .tuning)

if (sys.nframe() == 0L) {
  tryCatch(.tuning$main(), error = function(e) {
    message("Error: ", conditionMessage(e))
    quit(save = "no", status = 1L, runLast = FALSE)
  })
}
