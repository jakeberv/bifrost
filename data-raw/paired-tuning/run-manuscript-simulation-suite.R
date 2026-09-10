#!/usr/bin/env Rscript

# Local, portable runner for the 12-dimensional manuscript-aligned simulation
# suite. This file is intentionally excluded from Git tracking.

usage <- function() {
  cat(paste0(
    "Run the 12-dimensional bifrost manuscript simulation suite.\n\n",
    "Usage:\n",
    "  Rscript run-manuscript-simulation-suite.R [options]\n\n",
    "Safety:\n",
    "  The runner only prints its design unless --execute is supplied.\n\n",
    "Options:\n",
    "  --execute                  Run the selected replicate jobs.\n",
    "  --dry-run                  Print the selected design without running it.\n",
    "  --mode=full                Run 100 replicates per group (default).\n",
    "  --mode=smoke               Run one replicate per group.\n",
    "  --repo=PATH                Optional checkout for Git provenance and output placement.\n",
    "  --data-dir=PATH            Optional local override for the empirical input files.\n",
    "  --output-dir=PATH          Output root; relative paths resolve from repo or cwd.\n",
    "  --cores=N                  Single-core replicate workers on this machine.\n",
    "  --only=GROUP_ID            Run one group or comma-separated group IDs.\n",
    "  --job-index=N              Run one numbered replicate job.\n",
    "  --overwrite                Replace existing selected replicate outputs.\n",
    "  --list                     List the seven replicate groups and exit.\n",
    "  --list-jobs                List every numbered replicate job and exit.\n",
    "  --help                     Show this message and exit.\n\n",
    "Single-machine example:\n",
    "  Rscript run-manuscript-simulation-suite.R --list\n",
    "  Rscript run-manuscript-simulation-suite.R --execute --mode=full\n\n",
    "The development version of bifrost must already be available in the active\n",
    "R library. By default, empirical inputs are resolved with\n",
    "bifrost_example_file(); --data-dir provides an offline override. --repo\n",
    "never loads package code from a source checkout.\n\n",
    "The local physical core count is detected unless --cores is supplied. A\n",
    "dynamically scheduled fork queue remains within one machine. Each worker\n",
    "simulates one dataset and performs all paired searches sequentially on one\n",
    "core.\n"
  ))
}

stop_cli <- function(...) stop(paste0(...), call. = FALSE)

parse_positive_integer <- function(value, option) {
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed < 1L ||
      !identical(as.character(parsed), sub("^\\+", "", value))) {
    stop_cli(option, " must be a positive integer.")
  }
  parsed
}

parse_cli <- function(args) {
  out <- list(
    execute = FALSE, dry_run = FALSE, mode = "full", repo = NULL,
    data_dir = NULL,
    output_dir = NULL, cores = NULL, only = NULL, job_index = NULL,
    overwrite = FALSE, list = FALSE, list_jobs = FALSE, help = FALSE
  )
  for (arg in args) {
    if (identical(arg, "--execute")) out$execute <- TRUE
    else if (identical(arg, "--dry-run")) out$dry_run <- TRUE
    else if (identical(arg, "--overwrite")) out$overwrite <- TRUE
    else if (identical(arg, "--list")) out$list <- TRUE
    else if (identical(arg, "--list-jobs")) out$list_jobs <- TRUE
    else if (arg %in% c("--help", "-h")) out$help <- TRUE
    else if (startsWith(arg, "--mode=")) out$mode <- sub("^--mode=", "", arg)
    else if (startsWith(arg, "--repo=")) out$repo <- sub("^--repo=", "", arg)
    else if (startsWith(arg, "--data-dir=")) {
      out$data_dir <- sub("^--data-dir=", "", arg)
    }
    else if (startsWith(arg, "--output-dir=")) {
      out$output_dir <- sub("^--output-dir=", "", arg)
    } else if (startsWith(arg, "--cores=")) {
      out$cores <- parse_positive_integer(sub("^--cores=", "", arg), "--cores")
    } else if (startsWith(arg, "--only=")) {
      requested <- strsplit(sub("^--only=", "", arg), ",", fixed = TRUE)[[1L]]
      requested <- trimws(requested)
      if (length(requested) < 1L || any(!nzchar(requested))) {
        stop_cli("--only must contain at least one group ID.")
      }
      out$only <- requested
    } else if (startsWith(arg, "--job-index=")) {
      out$job_index <- parse_positive_integer(
        sub("^--job-index=", "", arg), "--job-index"
      )
    } else stop_cli("Unknown option: ", arg)
  }
  if (!out$mode %in% c("full", "smoke")) {
    stop_cli("--mode must be one of: full, smoke.")
  }
  if (out$execute && out$dry_run) stop_cli("Use only one of --execute and --dry-run.")
  if (!is.null(out$only) && !is.null(out$job_index)) {
    stop_cli("Use only one of --only and --job-index.")
  }
  if (out$list && out$list_jobs) stop_cli("Use only one of --list and --list-jobs.")
  out
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) != 1L) return(NULL)
  normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE)
}

is_bifrost_repo <- function(path) {
  description <- file.path(path, "DESCRIPTION")
  if (!file.exists(description)) return(FALSE)
  package_name <- tryCatch(
    read.dcf(description, fields = "Package")[[1L]],
    error = function(e) NA_character_
  )
  identical(package_name, "bifrost")
}

ancestor_paths <- function(path) {
  current <- normalizePath(path, mustWork = FALSE)
  paths <- character(0)
  repeat {
    paths <- c(paths, current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  paths
}

find_repo <- function(explicit = NULL) {
  if (!is.null(explicit)) {
    resolved <- normalizePath(path.expand(explicit), mustWork = FALSE)
    if (!is_bifrost_repo(resolved)) {
      stop_cli("--repo does not identify a bifrost source checkout: ", resolved)
    }
    return(resolved)
  }
  candidates <- ancestor_paths(getwd())
  matches <- candidates[vapply(candidates, is_bifrost_repo, logical(1L))]
  if (length(matches) < 1L) return(NULL)
  matches[[1L]]
}

data_input_paths <- function(data_dir) {
  setNames(
    file.path(data_dir, c(
      "passerine_bodyplan_tree.tre", "passerine_bodyplan_data.RDS"
    )),
    c("tree", "traits")
  )
}

is_data_dir <- function(path) {
  !is.null(path) && dir.exists(path) && all(file.exists(data_input_paths(path)))
}

find_data_dir <- function(explicit = NULL) {
  if (!is.null(explicit)) {
    resolved <- normalizePath(path.expand(explicit), mustWork = FALSE)
    if (!is_data_dir(resolved)) {
      stop_cli(
        "--data-dir must contain passerine_bodyplan_tree.tre and ",
        "passerine_bodyplan_data.RDS: ", resolved
      )
    }
    return(resolved)
  }
  NULL
}

resolve_data_inputs <- function(data_dir = NULL, namespace) {
  inputs <- if (!is.null(data_dir)) {
    data_input_paths(data_dir)
  } else {
    downloader <- get("bifrost_example_file", envir = namespace)
    c(
      tree = downloader("passerine-tree"),
      traits = downloader("passerine-traits")
    )
  }
  if (!identical(names(inputs), c("tree", "traits")) ||
      any(!file.exists(inputs))) {
    stop_cli("Could not resolve the empirical simulation inputs.")
  }
  vapply(inputs, normalizePath, character(1L), mustWork = TRUE)
}

detect_workers <- function(explicit = NULL) {
  if (!is.null(explicit)) return(explicit)
  detected <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (length(detected) != 1L || is.na(detected) || detected < 1L) {
    detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  }
  if (length(detected) != 1L || is.na(detected) || detected < 1L) return(1L)
  as.integer(detected)
}

new_group <- function(id, scenario, tree_tips, shifts, seed,
                      searches_per_replicate, full_replicates = 100L) {
  data.frame(
    id = id, scenario = scenario, tree_tip_count = as.integer(tree_tips),
    num_shifts = as.integer(shifts), seed = as.integer(seed),
    searches_per_replicate = as.integer(searches_per_replicate),
    full_replicates = as.integer(full_replicates), stringsAsFactors = FALSE
  )
}

build_groups <- function(mode = "full") {
  groups <- rbind(
    new_group("null-n100", "null", 100L, 0L, 5L, 4L),
    new_group("null-n200", "null", 200L, 0L, 5L, 4L),
    new_group("null-n300", "null", 300L, 0L, 5L, 4L),
    new_group("proportional-n250-s5", "proportional", 250L, 5L, 5L, 2L),
    new_group("proportional-n350-s10", "proportional", 350L, 10L, 1L, 2L),
    new_group(
      "integration-rate-n250-s5", "integration-rate", 250L, 5L, 5L, 2L
    ),
    new_group(
      "integration-rate-n350-s10", "integration-rate", 350L, 10L, 1L, 2L
    )
  )
  rownames(groups) <- NULL
  groups$index <- seq_len(nrow(groups))
  groups$n_replicates <- if (identical(mode, "smoke")) 1L else groups$full_replicates
  groups$total_searches <- groups$n_replicates * groups$searches_per_replicate
  groups[, c(
    "index", "id", "scenario", "tree_tip_count", "num_shifts", "seed",
    "searches_per_replicate", "full_replicates", "n_replicates",
    "total_searches"
  )]
}

search_configs <- function(group) {
  if (identical(group$scenario, "null")) {
    grid <- expand.grid(
      IC = c("GIC", "BIC"), shift_acceptance_threshold = c(2, 10),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    grid <- grid[order(
      match(grid$IC, c("GIC", "BIC")), grid$shift_acceptance_threshold
    ), , drop = FALSE]
  } else {
    grid <- data.frame(
      IC = c("GIC", "BIC"), shift_acceptance_threshold = c(10, 10),
      stringsAsFactors = FALSE
    )
  }
  grid$config_id <- sprintf(
    "%s-%s-d%s", group$id, tolower(grid$IC),
    format(grid$shift_acceptance_threshold, scientific = FALSE, trim = TRUE)
  )
  grid$weighted <- !identical(group$scenario, "null")
  rownames(grid) <- NULL
  grid[, c("config_id", "IC", "shift_acceptance_threshold", "weighted")]
}

build_jobs <- function(groups) {
  jobs <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    group <- groups[i, , drop = FALSE]
    data.frame(
      group_index = group$index, group_id = group$id,
      replicate = seq_len(group$n_replicates),
      job_id = sprintf("%s-r%03d", group$id, seq_len(group$n_replicates)),
      n_searches = group$searches_per_replicate,
      work_score = (group$tree_tip_count^2) * group$searches_per_replicate,
      stringsAsFactors = FALSE
    )
  }))
  rownames(jobs) <- NULL
  jobs$job_index <- seq_len(nrow(jobs))
  jobs[, c(
    "job_index", "job_id", "group_index", "group_id", "replicate",
    "n_searches", "work_score"
  )]
}

select_jobs <- function(groups, jobs, only = NULL, job_index = NULL) {
  if (!is.null(job_index)) {
    if (job_index > nrow(jobs)) {
      stop_cli("--job-index must be between 1 and ", nrow(jobs), ".")
    }
    return(jobs[job_index, , drop = FALSE])
  }
  if (is.null(only) || identical(only, "all")) return(jobs)
  unknown <- setdiff(only, groups$id)
  if (length(unknown) > 0L) {
    stop_cli("Unknown group ID(s): ", paste(unknown, collapse = ", "))
  }
  jobs[jobs$group_id %in% unique(only), , drop = FALSE]
}

print_group_list <- function(groups) {
  cat(sprintf(
    "%3s  %-31s  %6s  %8s  %8s\n",
    "#", "GROUP_ID", "JOBS", "EACH", "SEARCHES"
  ))
  for (i in seq_len(nrow(groups))) {
    cat(sprintf(
      "%3d  %-31s  %6d  %8d  %8d\n",
      groups$index[i], groups$id[i], groups$n_replicates[i],
      groups$searches_per_replicate[i], groups$total_searches[i]
    ))
  }
}

print_job_list <- function(jobs) {
  cat(sprintf("%4s  %-37s  %8s\n", "#", "JOB_ID", "SEARCHES"))
  for (i in seq_len(nrow(jobs))) {
    cat(sprintf(
      "%4d  %-37s  %8d\n",
      jobs$job_index[i], jobs$job_id[i], jobs$n_searches[i]
    ))
  }
}

print_design <- function(repo, data_dir, mode, workers, groups, jobs, selected,
                         output_dir, executing = FALSE) {
  selected_groups <- groups[groups$id %in% unique(selected$group_id), , drop = FALSE]
  cat("Bifrost manuscript simulation suite\n")
  cat("===================================\n")
  cat("Repository: ", if (is.null(repo)) "not supplied" else repo, "\n", sep = "")
  cat(
    "Data inputs: ",
    if (is.null(data_dir)) "bifrost_example_file() at execution" else data_dir,
    "\n", sep = ""
  )
  cat("Development package: use the active R library\n")
  cat("Output root: ", output_dir, "\n", sep = "")
  cat("Mode: ", mode, "\n", sep = "")
  cat("Response traits: 12\n")
  cat("Workers: ", workers, "\n", sep = "")
  cat("Parallelism: one core per replicate job; workers remain on this machine\n")
  cat("Queue: dynamically scheduled single-core replicate jobs\n")
  cat("Replicate groups: ", nrow(groups), "\n", sep = "")
  cat("Replicate jobs: ", nrow(jobs), "\n", sep = "")
  cat("Search conditions: 20\n")
  cat("Total searches: ", sum(jobs$n_searches), "\n", sep = "")
  cat("Selected replicate jobs: ", nrow(selected), "\n", sep = "")
  cat("Selected searches: ", sum(selected$n_searches), "\n\n", sep = "")
  print_group_list(selected_groups)
  if (nrow(selected) <= 20L) {
    cat("\nSelected jobs\n")
    print_job_list(selected)
  }
  if (!executing) {
    cat("\nNo simulations were run. Supply --execute to start the selected jobs.\n")
  }
}

resolve_output_dir <- function(repo, supplied = NULL) {
  environment_path <- Sys.getenv("BIFROST_SIM_OUTPUT_DIR", unset = "")
  path <- if (is.null(supplied) && nzchar(environment_path)) {
    path.expand(environment_path)
  } else if (is.null(supplied) && !is.null(repo)) {
    file.path(repo, "local-cache", "manuscript-simulation-suite-12d")
  } else if (is.null(supplied)) {
    file.path(getwd(), "manuscript-simulation-suite-12d")
  } else if (grepl("^(/|[A-Za-z]:[/\\\\])", path.expand(supplied))) {
    path.expand(supplied)
  } else if (!is.null(repo)) {
    file.path(repo, supplied)
  } else file.path(getwd(), supplied)
  normalizePath(path, mustWork = FALSE)
}

set_single_threaded_math <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
}

load_development_package <- function() {
  if (!requireNamespace("bifrost", quietly = TRUE)) {
    stop_cli(
      "The development version of bifrost must be available in the active R library."
    )
  }
  namespace <- loadNamespace("bifrost")
  required <- c(
    "createSimulationTemplate", "simulateNullDataset", "simulateShiftedDataset",
    "searchOptimalConfiguration", "evaluateShiftRecovery", "bifrost_example_file"
  )
  missing <- required[!vapply(
    required, exists, logical(1L), envir = namespace, inherits = FALSE
  )]
  if (length(missing) > 0L) {
    stop_cli(
      "The available bifrost development version is missing required functions: ",
      paste(missing, collapse = ", "), "."
    )
  }
  namespace
}

build_template <- function(inputs, namespace) {
  tree_path <- inputs[["tree"]]
  trait_path <- inputs[["traits"]]
  bird_tree <- ape::read.tree(tree_path)
  bodyplan_data <- readRDS(trait_path)
  if (length(setdiff(bird_tree$tip.label, rownames(bodyplan_data))) > 0L) {
    stop_cli("The body-plan data are missing tips present in the phylogeny.")
  }
  bodyplan_data <- as.matrix(bodyplan_data[bird_tree$tip.label, , drop = FALSE])
  if (!"vertnet_mass" %in% colnames(bodyplan_data)) {
    stop_cli("The body-plan data must contain the predictor column vertnet_mass.")
  }
  response_columns <- setdiff(colnames(bodyplan_data), "vertnet_mass")
  if (length(response_columns) != 12L) {
    stop_cli(
      "The empirical simulation template must contain exactly 12 response traits; found ",
      length(response_columns), "."
    )
  }
  bodyplan_data <- bodyplan_data[, c(response_columns, "vertnet_mass"), drop = FALSE]
  formula <- sprintf(
    "trait_data[, 1:%d] ~ trait_data[, %d]",
    length(response_columns), ncol(bodyplan_data)
  )
  create_template <- get("createSimulationTemplate", envir = namespace)
  template <- create_template(
    baseline_tree = bird_tree, trait_data = bodyplan_data, formula = formula,
    response_columns = seq_along(response_columns),
    predictor_columns = ncol(bodyplan_data), method = "LL", error = TRUE
  )
  if (!identical(template$n_response_traits, 12L) ||
      !identical(dim(template$residual_covariance), c(12L, 12L))) {
    stop_cli("The fitted simulation template is not 12-dimensional.")
  }
  template
}

git_value <- function(repo, args) {
  if (is.null(repo)) return(NA_character_)
  result <- tryCatch(
    system2("git", c("-C", shQuote(repo), args), stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)
  )
  if (length(result) < 1L) NA_character_ else paste(result, collapse = "\n")
}

collect_provenance <- function(repo, data_dir, inputs, template, namespace) {
  package_path <- normalizePath(
    getNamespaceInfo(namespace, "path"), mustWork = TRUE
  )
  package_files <- c(
    file.path(package_path, "DESCRIPTION"), file.path(package_path, "NAMESPACE"),
    file.path(package_path, "R", c("bifrost.rdb", "bifrost.rdx")),
    file.path(package_path, "Meta", "package.rds")
  )
  package_files <- package_files[file.exists(package_files)]
  package_md5 <- unname(tools::md5sum(package_files))
  names(package_md5) <- substring(package_files, nchar(package_path) + 2L)
  data_md5 <- unname(tools::md5sum(inputs))
  names(data_md5) <- names(inputs)
  description <- tryCatch(
    read.dcf(file.path(package_path, "DESCRIPTION")),
    error = function(e) matrix(character(0), nrow = 0L, ncol = 0L)
  )
  description_value <- function(field) {
    if (nrow(description) < 1L || !field %in% colnames(description)) {
      return(NA_character_)
    }
    unname(description[1L, field])
  }
  own_path <- script_path()
  source_identity <- list(
    package_version = as.character(utils::packageVersion("bifrost")),
    package_md5 = package_md5, data_md5 = data_md5,
    remote_sha = description_value("RemoteSha"),
    github_sha1 = description_value("GithubSHA1")
  )
  list(
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    git_commit = git_value(repo, c("rev-parse", "HEAD")),
    git_branch = git_value(repo, c("branch", "--show-current")),
    tracked_status = git_value(repo, c("status", "--short", "--untracked-files=no")),
    repository = repo,
    data_source = if (is.null(data_dir)) "bifrost_example_file" else "data_dir",
    data_directory = data_dir, data_files = inputs, package_path = package_path,
    package_version = source_identity$package_version,
    package_remote_sha = source_identity$remote_sha,
    package_github_sha1 = source_identity$github_sha1,
    package_built = description_value("Built"),
    package_md5 = package_md5, data_md5 = data_md5,
    source_fingerprint = digest::digest(source_identity, algo = "sha256"),
    runner_sha256 = if (!is.null(own_path) && file.exists(own_path)) {
      digest::digest(file = own_path, algo = "sha256", serialize = FALSE)
    } else NA_character_,
    template_traits = template$n_response_traits, template_tips = template$n_tips,
    template_covariance_sha256 = digest::digest(
      template$residual_covariance, algo = "sha256"
    ),
    template_residual_df = template$residual_df,
    platform = R.version$platform, R = R.version.string
  )
}

as_scalar_list <- function(row) {
  lapply(row, function(x) if (length(x) == 1L) unname(x) else x)
}

suite_contract <- function(groups, mode) {
  list(
    schema_version = 2L, mode = mode, response_traits = 12L,
    simulation_generator = "empirical",
    groups = groups[, c(
      "id", "scenario", "tree_tip_count", "num_shifts", "seed",
      "n_replicates", "searches_per_replicate"
    )],
    search = list(
      formula = "trait_data ~ 1", min_descendant_tips = 10L,
      ic_uncertainty_threshold = 10, method = "LL", error = FALSE,
      fuzzy_distance = 2L, IC = c("GIC", "BIC"),
      null_thresholds = c(2, 10), recovery_threshold = 10,
      weighted_recovery = TRUE
    ),
    shifted_simulation = list(
      min_shift_tips = 10L, max_shift_tips = 20L, buffer = 3L,
      scale_factor_range = c(0.1, 2.0), exclude_range = c(0.5, 1.5),
      integration_power_range = c(0.5, 1.25),
      integration_exclude_range = c(0.8, 1.1), eigen_floor = 1e-8
    )
  )
}

stable_seed <- function(base_seed, key) {
  hash <- digest::digest(
    paste(base_seed, key, sep = ":"), algo = "xxhash32", serialize = FALSE
  )
  as.integer(strtoi(substr(hash, 1L, 7L), base = 16L)) + 1L
}

job_seeds <- function(group, replicate) {
  sampling_key <- if (identical(group$scenario, "null")) {
    sprintf("null-n%d-r%d", group$tree_tip_count, replicate)
  } else {
    sprintf("shifted-n%d-s%d-r%d", group$tree_tip_count, group$num_shifts, replicate)
  }
  list(
    simulation = stable_seed(group$seed, paste0("simulation:", sampling_key)),
    search = stable_seed(group$seed, paste0("search:", group$id, ":", replicate))
  )
}

simulation_options <- function(group) {
  if (identical(group$scenario, "null")) {
    return(list(simulation_generator = "empirical"))
  }
  options <- list(
    simulation_generator = "empirical", num_shifts = group$num_shifts,
    min_shift_tips = 10L, max_shift_tips = 20L,
    scale_mode = if (identical(group$scenario, "proportional")) {
      "proportional"
    } else "correlation",
    scale_factor_range = c(0.1, 2.0), exclude_range = c(0.5, 1.5), buffer = 3L
  )
  if (identical(group$scenario, "integration-rate")) {
    options$integration_power_range <- c(0.5, 1.25)
    options$integration_exclude_range <- c(0.8, 1.1)
    options$eigen_floor <- 1e-8
  }
  options
}

base_search_options <- function(config) {
  list(
    formula = "trait_data ~ 1", min_descendant_tips = 10L, num_cores = 1L,
    ic_uncertainty_threshold = 10,
    shift_acceptance_threshold = config$shift_acceptance_threshold,
    uncertaintyweights = FALSE,
    uncertaintyweights_par = isTRUE(config$weighted), plot = FALSE,
    IC = config$IC, store_model_fit_history = FALSE, verbose = FALSE,
    progress = FALSE, method = "LL", error = FALSE
  )
}

job_output_path <- function(output_dir, mode, job) {
  file.path(output_dir, mode, "replicates", job$group_id, paste0(job$job_id, ".rds"))
}

summary_output_path <- function(output_dir, mode, config_id) {
  file.path(output_dir, mode, "summaries", paste0(config_id, ".rds"))
}

atomic_save_rds <- function(object, path, overwrite = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !overwrite) stop_cli("Output already exists: ", path)
  temporary <- tempfile(
    pattern = paste0(".", basename(path), "."), tmpdir = dirname(path),
    fileext = ".tmp"
  )
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = "gzip", version = 3)
  check <- readRDS(temporary)
  if (!identical(check$schema_version, object$schema_version)) {
    stop_cli("Temporary RDS failed read-back validation: ", temporary)
  }
  if (!file.rename(temporary, path)) {
    stop_cli("Could not atomically move completed output into place: ", path)
  }
  invisible(path)
}

atomic_write_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(".", basename(path), "."), tmpdir = dirname(path),
    fileext = ".tmp"
  )
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(data, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    stop_cli("Could not atomically move completed CSV into place: ", path)
  }
  invisible(path)
}

result_failed <- function(result) {
  error <- result$error
  !is.null(error) && length(error) > 0L &&
    any(!is.na(error) & nzchar(as.character(error)))
}

compact_search_result <- function(result) {
  keep <- c(
    "shift_nodes_no_uncertainty",
    "num_candidates",
    "candidate_nodes",
    "optimal_ic",
    "baseline_ic",
    "IC_used",
    "ic_weights",
    "warnings",
    "error"
  )
  compact <- result[intersect(keep, names(result))]
  class(compact) <- c("bifrost_compact_search_result", "list")
  compact
}

simulation_content_hash <- function(simulation) {
  content <- simulation
  content$user_input <- NULL
  digest::digest(content, algo = "sha256")
}

validate_search_result <- function(result, true_nodes = integer(0)) {
  candidates <- result$candidate_nodes
  n_candidates <- result$num_candidates
  valid <- is.numeric(candidates) && !anyNA(candidates) &&
    all(is.finite(candidates)) && all(candidates == floor(candidates)) &&
    !anyDuplicated(candidates)
  if (!valid || length(candidates) != n_candidates) {
    stop_cli("A search result is missing a valid candidate_nodes vector.")
  }
  if (length(setdiff(true_nodes, candidates)) > 0L) {
    stop_cli("A planted shift is absent from the eligible candidate-node set.")
  }
  invisible(TRUE)
}

simulate_job_dataset <- function(group, seed, namespace, template) {
  simulator <- if (identical(group$scenario, "null")) {
    get("simulateNullDataset", envir = namespace)
  } else get("simulateShiftedDataset", envir = namespace)
  do.call(
    simulator,
    c(
      list(template = template, tree_tip_count = group$tree_tip_count, seed = seed),
      simulation_options(group)
    )
  )
}

with_local_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  code()
}

eligible_candidate_nodes <- function(tree, min_descendant_tips) {
  tree <- ape::as.phylo(tree)
  n_tips <- ape::Ntip(tree)
  root <- n_tips + 1L
  internal <- setdiff(unique(as.integer(tree$edge[, 1L])), root)
  children <- split(as.integer(tree$edge[, 2L]), as.integer(tree$edge[, 1L]))
  descendant_tip_count <- function(node) {
    stack <- as.character(node)
    tips <- integer(0)
    while (length(stack) > 0L) {
      current <- stack[[length(stack)]]
      stack <- stack[-length(stack)]
      descendants <- children[[current]]
      if (is.null(descendants)) next
      tips <- c(tips, descendants[descendants <= n_tips])
      internal_descendants <- descendants[descendants > n_tips]
      if (length(internal_descendants) > 0L) {
        stack <- c(stack, as.character(internal_descendants))
      }
    }
    length(unique(tips))
  }
  counts <- vapply(internal, descendant_tip_count, integer(1L))
  as.integer(internal[counts >= min_descendant_tips])
}

search_job_dataset <- function(simulation, group, config, seed, namespace) {
  tree <- if (identical(group$scenario, "null")) {
    simulation$tree
  } else simulation$paintedTree
  trait_data <- if (!is.null(simulation$trait_data)) {
    simulation$trait_data
  } else if (!is.null(simulation$data)) {
    simulation$data
  } else simulation$simulatedData
  search_options <- base_search_options(config)
  search <- get("searchOptimalConfiguration", envir = namespace)
  result <- with_local_seed(seed, function() {
    tryCatch(
      withCallingHandlers(
        do.call(
          search,
          c(
            search_options,
            list(baseline_tree = ape::as.phylo(tree), trait_data = trait_data)
          )
        ),
        bifrost_search_settings_warning = function(warning) {
          invokeRestart("muffleWarning")
        }
      ),
      error = function(error) {
        candidates <- eligible_candidate_nodes(
          tree,
          min_descendant_tips = search_options$min_descendant_tips
        )
        list(
          shift_nodes_no_uncertainty = integer(0),
          num_candidates = length(candidates),
          candidate_nodes = candidates,
          ic_weights = data.frame(
            node = integer(0),
            ic_with_shift = numeric(0),
            ic_without_shift = numeric(0),
            delta_ic = numeric(0),
            ic_weight_withshift = numeric(0),
            ic_weight_withoutshift = numeric(0),
            evidence_ratio = numeric(0)
          ),
          error = conditionMessage(error)
        )
      }
    )
  })
  compact <- compact_search_result(result)
  rm(result)
  invisible(gc(verbose = FALSE))
  compact
}

new_job_wrapper <- function(job, group, configs, seeds, simulation, results,
                            provenance, design_fingerprint, mode,
                            started_at, finished_at) {
  list(
    schema_version = 2L, suite = "bifrost-manuscript-simulation-suite-12d",
    mode = mode, job = as_scalar_list(job), group = as_scalar_list(group),
    search_configs = configs, seeds = seeds,
    design_fingerprint = design_fingerprint,
    source_fingerprint = provenance$source_fingerprint,
    runner_sha256 = provenance$runner_sha256, provenance = provenance,
    started_at = started_at, finished_at = finished_at,
    simulation_hash = simulation_content_hash(simulation),
    simulation = simulation, results = results
  )
}

validate_job_wrapper <- function(wrapper, job, group, configs, provenance,
                                 design_fingerprint) {
  if (!is.list(wrapper) || !identical(wrapper$schema_version, 2L) ||
      !identical(wrapper$suite, "bifrost-manuscript-simulation-suite-12d") ||
      !identical(wrapper$job, as_scalar_list(job)) ||
      !identical(wrapper$group, as_scalar_list(group)) ||
      !identical(wrapper$search_configs, configs) ||
      !identical(wrapper$design_fingerprint, design_fingerprint)) {
    stop_cli("A replicate output does not match the current suite design.")
  }
  if (!identical(wrapper$source_fingerprint, provenance$source_fingerprint)) {
    stop_cli("A replicate output used different package sources: ", job$job_id)
  }
  simulation <- wrapper$simulation
  if (!identical(simulation$simulation_generator, "empirical") ||
      ncol(as.matrix(simulation$trait_data)) != 12L ||
      !identical(wrapper$simulation_hash, simulation_content_hash(simulation))) {
    stop_cli("A replicate output failed simulation validation: ", job$job_id)
  }
  if (!identical(names(wrapper$results), configs$config_id)) {
    stop_cli("A replicate output has incomplete search configurations: ", job$job_id)
  }
  true_nodes <- if (identical(group$scenario, "null")) integer(0) else simulation$shiftNodes
  for (result in wrapper$results) validate_search_result(result, true_nodes)
  invisible(TRUE)
}

run_replicate_job <- function(job, group, configs, output_dir, mode,
                              namespace, template, provenance,
                              design_fingerprint, overwrite = FALSE) {
  tryCatch({
    set_single_threaded_math()
    path <- job_output_path(output_dir, mode, job)
    if (file.exists(path) && !overwrite) {
      validate_job_wrapper(
        readRDS(path), job, group, configs, provenance, design_fingerprint
      )
      return(list(status = "skipped", job_id = job$job_id, path = path))
    }
    seeds <- job_seeds(group, job$replicate)
    started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    simulation <- simulate_job_dataset(group, seeds$simulation, namespace, template)
    results <- setNames(vector("list", nrow(configs)), configs$config_id)
    for (i in seq_len(nrow(configs))) {
      results[[i]] <- search_job_dataset(
        simulation, group, configs[i, , drop = FALSE], seeds$search, namespace
      )
    }
    finished_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    wrapper <- new_job_wrapper(
      job, group, configs, seeds, simulation, results, provenance,
      design_fingerprint, mode, started_at, finished_at
    )
    validate_job_wrapper(
      wrapper, job, group, configs, provenance, design_fingerprint
    )
    atomic_save_rds(wrapper, path, overwrite = overwrite)
    list(status = "written", job_id = job$job_id, path = path)
  }, error = function(e) {
    list(status = "error", job_id = job$job_id, error = conditionMessage(e))
  })
}

load_group_wrappers <- function(group, jobs, output_dir, mode, configs,
                                provenance, design_fingerprint) {
  group_jobs <- jobs[jobs$group_id == group$id, , drop = FALSE]
  wrappers <- vector("list", nrow(group_jobs))
  for (i in seq_len(nrow(group_jobs))) {
    job <- group_jobs[i, , drop = FALSE]
    path <- job_output_path(output_dir, mode, job)
    if (!file.exists(path)) return(NULL)
    wrapper <- readRDS(path)
    validate_job_wrapper(
      wrapper, job, group, configs, provenance, design_fingerprint
    )
    wrappers[[i]] <- wrapper
  }
  wrappers
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 1L) NA_real_ else mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 1L) NA_real_ else stats::median(x)
}

compare_candidate_and_legacy <- function(evaluator, simdata, results) {
  candidate <- evaluator(
    simdata, results, fuzzy_distance = 2L, weighted = TRUE, verbose = FALSE
  )
  legacy_results <- lapply(results, function(result) {
    result$candidate_nodes <- NULL
    result
  })
  legacy <- evaluator(
    simdata, legacy_results, fuzzy_distance = 2L, weighted = TRUE, verbose = FALSE
  )
  checked <- c("specificity", "fpr", "balanced_accuracy")
  agrees <- all(vapply(c("strict", "fuzzy"), function(matching) {
    isTRUE(all.equal(
      unlist(candidate[[matching]][checked]), unlist(legacy[[matching]][checked]),
      tolerance = 1e-12, check.attributes = FALSE
    ))
  }, logical(1L)))
  if (!agrees) {
    stop_cli(
      "Candidate-aware and legacy TN calculations disagree despite all planted ",
      "shifts being eligible."
    )
  }
  candidate
}

summarize_config <- function(group, config, wrappers, namespace,
                             provenance, design_fingerprint, mode) {
  results <- lapply(wrappers, function(x) x$results[[config$config_id]])
  completed <- !vapply(results, result_failed, logical(1L))
  n_candidates <- vapply(results, `[[`, numeric(1L), "num_candidates")
  n_inferred <- vapply(seq_along(results), function(i) {
    if (!completed[i]) NA_real_ else length(results[[i]]$shift_nodes_no_uncertainty)
  }, numeric(1L))
  per_replicate <- data.frame(
    replicate = seq_along(results),
    simulation_hash = vapply(wrappers, `[[`, character(1L), "simulation_hash"),
    status = ifelse(completed, "ok", "error"), n_candidates = n_candidates,
    n_inferred_shifts = n_inferred,
    error = vapply(results, function(result) {
      if (!result_failed(result)) NA_character_ else paste(result$error, collapse = "; ")
    }, character(1L)), stringsAsFactors = FALSE
  )
  summary_row <- data.frame(
    config_id = config$config_id, scenario = group$scenario,
    tree_tip_count = group$tree_tip_count, num_shifts = group$num_shifts,
    IC = config$IC, shift_acceptance_threshold = config$shift_acceptance_threshold,
    min_descendant_tips = 10L, n_replicates = length(results),
    n_completed = sum(completed), n_failed = sum(!completed),
    completion_rate = mean(completed), mean_inferred_shifts = safe_mean(n_inferred),
    stringsAsFactors = FALSE
  )
  evaluation <- NULL
  if (identical(group$scenario, "null")) {
    fp_rate <- ifelse(completed & n_candidates > 0, n_inferred / n_candidates, NA)
    per_replicate$false_positive_rate <- fp_rate
    summary_row$mean_false_positive_rate <- safe_mean(fp_rate)
    summary_row$median_false_positive_rate <- safe_median(fp_rate)
    summary_row$fraction_any_false_positive <- safe_mean(n_inferred[completed] > 0)
  } else {
    simdata <- lapply(wrappers, `[[`, "simulation")
    evaluator <- get("evaluateShiftRecovery", envir = namespace)
    evaluation <- compare_candidate_and_legacy(evaluator, simdata, results)
    for (metric in names(evaluation$fuzzy)) {
      summary_row[[paste0("fuzzy_", metric)]] <- evaluation$fuzzy[[metric]]
    }
    for (metric in names(evaluation$strict)) {
      summary_row[[paste0("strict_", metric)]] <- evaluation$strict[[metric]]
    }
    if (!is.null(evaluation$weighted)) {
      for (metric in names(evaluation$weighted$fuzzy)) {
        summary_row[[paste0("weighted_fuzzy_", metric)]] <-
          evaluation$weighted$fuzzy[[metric]]
      }
    }
  }
  list(
    schema_version = 2L, suite = "bifrost-manuscript-simulation-suite-12d",
    mode = mode, group = as_scalar_list(group), config = as_scalar_list(config),
    design_fingerprint = design_fingerprint,
    source_fingerprint = provenance$source_fingerprint,
    per_replicate = per_replicate, evaluation = evaluation,
    summary_row = summary_row
  )
}

aggregate_complete_groups <- function(groups, jobs, output_dir, mode,
                                      namespace, provenance,
                                      design_fingerprint) {
  summary_rows <- list()
  complete_groups <- character(0)
  for (i in seq_len(nrow(groups))) {
    group <- groups[i, , drop = FALSE]
    configs <- search_configs(group)
    wrappers <- load_group_wrappers(
      group, jobs, output_dir, mode, configs, provenance, design_fingerprint
    )
    if (is.null(wrappers)) next
    complete_groups <- c(complete_groups, group$id)
    for (j in seq_len(nrow(configs))) {
      config <- configs[j, , drop = FALSE]
      summary <- summarize_config(
        group, config, wrappers, namespace, provenance, design_fingerprint, mode
      )
      atomic_save_rds(
        summary, summary_output_path(output_dir, mode, config$config_id),
        overwrite = TRUE
      )
      summary_rows[[length(summary_rows) + 1L]] <- summary$summary_row
    }
  }
  if (length(summary_rows) < 1L) {
    return(list(n_complete_groups = 0L, n_summaries = 0L))
  }
  all_columns <- unique(unlist(lapply(summary_rows, names)))
  normalized <- lapply(summary_rows, function(row) {
    for (name in setdiff(all_columns, names(row))) row[[name]] <- NA
    row[, all_columns, drop = FALSE]
  })
  summary_table <- do.call(rbind, normalized)
  rownames(summary_table) <- NULL
  suite_summary <- list(
    schema_version = 2L, suite = "bifrost-manuscript-simulation-suite-12d",
    mode = mode, design_fingerprint = design_fingerprint,
    source_fingerprint = provenance$source_fingerprint,
    complete_groups = complete_groups,
    all_groups_complete = length(complete_groups) == nrow(groups),
    summary_table = summary_table
  )
  summary_root <- file.path(output_dir, mode, "summaries")
  atomic_save_rds(
    suite_summary, file.path(summary_root, "suite-summary.rds"), overwrite = TRUE
  )
  atomic_write_csv(summary_table, file.path(summary_root, "suite-summary.csv"))
  list(n_complete_groups = length(complete_groups), n_summaries = nrow(summary_table))
}

run_suite <- function(cli) {
  repo <- find_repo(cli$repo)
  data_dir <- find_data_dir(cli$data_dir)
  workers <- detect_workers(cli$cores)
  groups <- build_groups(cli$mode)
  jobs <- build_jobs(groups)
  selected <- select_jobs(groups, jobs, cli$only, cli$job_index)
  output_dir <- resolve_output_dir(repo, cli$output_dir)
  if (cli$list) {
    print_group_list(groups)
    return(invisible(NULL))
  }
  if (cli$list_jobs) {
    print_job_list(jobs)
    return(invisible(NULL))
  }
  print_design(
    repo, data_dir, cli$mode, workers, groups, jobs, selected, output_dir,
    cli$execute
  )
  if (!cli$execute) return(invisible(NULL))
  if (.Platform$OS.type == "windows" && workers > 1L) {
    stop_cli("Parallel execution requires a Unix-like machine.")
  }

  set_single_threaded_math()
  namespace <- load_development_package()
  inputs <- resolve_data_inputs(data_dir, namespace)
  template <- build_template(inputs, namespace)
  provenance <- collect_provenance(repo, data_dir, inputs, template, namespace)
  contract <- suite_contract(groups, cli$mode)
  design_fingerprint <- digest::digest(contract, algo = "sha256")

  pending <- list()
  for (i in seq_len(nrow(selected))) {
    job <- selected[i, , drop = FALSE]
    group <- groups[groups$id == job$group_id, , drop = FALSE]
    configs <- search_configs(group)
    path <- job_output_path(output_dir, cli$mode, job)
    if (file.exists(path) && !cli$overwrite) {
      validate_job_wrapper(
        readRDS(path), job, group, configs, provenance, design_fingerprint
      )
    } else {
      pending[[length(pending) + 1L]] <- list(job = job, group = group, configs = configs)
    }
  }

  if (length(pending) > 0L) {
    scores <- vapply(pending, function(x) x$job$work_score, numeric(1L))
    indices <- vapply(pending, function(x) x$job$job_index, integer(1L))
    pending <- pending[order(-scores, indices)]
    worker <- function(item) {
      run_replicate_job(
        item$job, item$group, item$configs, output_dir, cli$mode,
        namespace, template, provenance, design_fingerprint, cli$overwrite
      )
    }
    active_workers <- min(workers, length(pending))
    message(
      "Dispatching ", length(pending), " replicate job(s) to ",
      active_workers, " worker(s) on this machine."
    )
    results <- if (active_workers == 1L) lapply(pending, worker) else {
      parallel::mclapply(
        pending, worker, mc.cores = active_workers, mc.preschedule = FALSE,
        mc.set.seed = FALSE, mc.cleanup = TRUE
      )
    }
    failed <- vapply(results, function(x) identical(x$status, "error"), logical(1L))
    if (any(failed)) {
      for (failure in results[failed]) {
        message("Failed ", failure$job_id, ": ", failure$error)
      }
      stop_cli(
        sum(failed), " replicate job(s) failed. Completed outputs were retained; ",
        "rerun the same command to resume."
      )
    }
    message(
      "Completed ", sum(vapply(results, function(x) x$status == "written", logical(1L))),
      " replicate job(s)."
    )
  } else message("All selected replicate outputs already exist and validated.")

  aggregation <- aggregate_complete_groups(
    groups, jobs, output_dir, cli$mode, namespace, provenance, design_fingerprint
  )
  message(
    "Aggregated ", aggregation$n_summaries, " search condition(s) across ",
    aggregation$n_complete_groups, " complete replicate group(s)."
  )
  invisible(NULL)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cli <- parse_cli(args)
  if (cli$help) {
    usage()
    return(invisible(NULL))
  }
  run_suite(cli)
}

if (sys.nframe() == 0L) {
  tryCatch(
    main(),
    error = function(e) {
      message("Error: ", conditionMessage(e))
      quit(save = "no", status = 1L, runLast = FALSE)
    }
  )
}
