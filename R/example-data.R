.bifrost_example_manifest_schema <- 2L
.bifrost_example_manifest_url <- paste0(
  "https://raw.githubusercontent.com/jakeberv/bifrost/",
  "main/data-remote/empirical-artifacts.json"
)
.bifrost_example_raw_base <-
  "https://raw.githubusercontent.com/jakeberv/bifrost/main/"
.bifrost_example_identifiers <- c(
  "jaw-tree", "jaw-landmarks", "passerine-tree", "passerine-traits",
  "passerine-search", "passerine-sensitivity", "passerine-posthoc",
  "simulation-preview-tables"
)

.bifrost_example_sha256 <- function(path) {
  digest::digest(path, algo = "sha256", serialize = FALSE, file = TRUE)
}

.bifrost_read_example_manifest <- function(path) {
  tryCatch(
    jsonlite::read_json(path, simplifyVector = FALSE),
    error = function(error) {
      stop(
        "Could not read bifrost example-data manifest at ", path,
        ": ", conditionMessage(error),
        call. = FALSE
      )
    }
  )
}

.bifrost_manifest_error <- function(source, message) {
  stop(
    "Invalid bifrost example-data manifest from ", source, ": ", message,
    call. = FALSE
  )
}

.bifrost_is_scalar_string <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

.bifrost_is_safe_repository_path <- function(path) {
  if (!.bifrost_is_scalar_string(path) || startsWith(path, "/") ||
      startsWith(path, "\\") || grepl("^[A-Za-z]:", path) ||
      grepl("\\\\", path)) {
    return(FALSE)
  }

  segments <- strsplit(path, "/", fixed = TRUE)[[1L]]
  length(segments) > 0L &&
    all(nzchar(segments)) &&
    !any(segments %in% c(".", ".."))
}

.bifrost_is_package_version <- function(value) {
  if (!.bifrost_is_scalar_string(value)) {
    return(FALSE)
  }

  !inherits(try(base::package_version(value), silent = TRUE), "try-error")
}

.bifrost_validate_example_manifest <- function(manifest, source) {
  if (!is.list(manifest) || is.null(manifest$schema_version) ||
      length(manifest$schema_version) != 1L ||
      !is.numeric(manifest$schema_version) ||
      is.na(manifest$schema_version) ||
      manifest$schema_version != .bifrost_example_manifest_schema) {
    .bifrost_manifest_error(source, "schema version must be exactly 2")
  }

  if (!is.list(manifest$artifacts)) {
    .bifrost_manifest_error(source, "artifacts must be a list")
  }

  identifiers <- character()
  paths <- character()
  downloader_fields <- c(
    "artifact_id", "path", "sha256", "size_bytes", "minimum_bifrost_version"
  )

  for (index in seq_along(manifest$artifacts)) {
    artifact <- manifest$artifacts[[index]]
    label <- paste0("artifact ", index)
    if (!is.list(artifact)) {
      .bifrost_manifest_error(source, paste0(label, " must be a record"))
    }

    if (!.bifrost_is_safe_repository_path(artifact$path)) {
      .bifrost_manifest_error(
        source,
        paste0(label, " has an unsafe repository-relative path")
      )
    }
    paths <- c(paths, artifact$path)

    if (!is.null(artifact$sha256) &&
        (!.bifrost_is_scalar_string(artifact$sha256) ||
         !grepl("^[0-9a-f]{64}$", artifact$sha256))) {
      .bifrost_manifest_error(source, paste0(label, " has an invalid SHA-256"))
    }

    has_downloader_metadata <- !is.null(artifact$artifact_id) ||
      !is.null(artifact$size_bytes) ||
      !is.null(artifact$minimum_bifrost_version) ||
      startsWith(artifact$path, "data-remote/")
    if (!has_downloader_metadata) {
      next
    }

    missing_fields <- downloader_fields[vapply(
      downloader_fields,
      function(field) is.null(artifact[[field]]),
      logical(1)
    )]
    if (length(missing_fields) > 0L) {
      .bifrost_manifest_error(
        source,
        paste0(label, " is missing downloader field(s): ",
               paste(missing_fields, collapse = ", "))
      )
    }
    if (!startsWith(artifact$path, "data-remote/")) {
      .bifrost_manifest_error(
        source,
        paste0(label, " downloader path must live under data-remote/")
      )
    }
    if (!.bifrost_is_scalar_string(artifact$artifact_id) ||
        !grepl("^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$", artifact$artifact_id)) {
      .bifrost_manifest_error(
        source,
        paste0(label, " has an invalid lowercase hyphenated artifact_id")
      )
    }
    if (length(artifact$size_bytes) != 1L || !is.numeric(artifact$size_bytes) ||
        is.na(artifact$size_bytes) || !is.finite(artifact$size_bytes) ||
        artifact$size_bytes <= 0 || artifact$size_bytes != floor(artifact$size_bytes)) {
      .bifrost_manifest_error(
        source,
        paste0(label, " must have a positive whole-number byte size")
      )
    }
    if (!.bifrost_is_package_version(artifact$minimum_bifrost_version)) {
      .bifrost_manifest_error(
        source,
        paste0(label, " has an invalid minimum_bifrost_version")
      )
    }
    identifiers <- c(identifiers, artifact$artifact_id)
  }

  if (anyDuplicated(paths)) {
    .bifrost_manifest_error(source, "duplicate path entries are not allowed")
  }
  if (anyDuplicated(identifiers)) {
    .bifrost_manifest_error(source, "duplicate artifact_id entries are not allowed")
  }
  if (!setequal(identifiers, .bifrost_example_identifiers)) {
    .bifrost_manifest_error(
      source,
      paste0(
        "downloader artifacts must have exactly the supported identifiers: ",
        paste(.bifrost_example_identifiers, collapse = ", ")
      )
    )
  }

  manifest
}

.bifrost_example_entry <- function(manifest, name) {
  identifiers <- vapply(
    manifest$artifacts,
    function(artifact) {
      if (is.null(artifact$artifact_id)) "" else artifact$artifact_id
    },
    character(1)
  )
  index <- match(name, identifiers)
  manifest$artifacts[[index]]
}

.bifrost_verify_example_file <- function(path, entry, source) {
  identifier <- entry$artifact_id
  if (!utils::file_test("-f", path)) {
    stop(
      "Example-data artifact '", identifier, "' from ", source,
      " is missing at ", path,
      call. = FALSE
    )
  }
  if (!identical(.bifrost_example_sha256(path), entry$sha256)) {
    stop(
      "checksum mismatch for example-data artifact '", identifier,
      "' from ", source,
      call. = FALSE
    )
  }
  if (unname(file.info(path)$size) != entry$size_bytes) {
    stop(
      "Example-data artifact '", identifier, "' from ", source,
      " has a byte-size mismatch",
      call. = FALSE
    )
  }
  invisible(path)
}

.bifrost_installed_version <- function() {
  utils::packageVersion("bifrost")
}

.bifrost_assert_example_version <- function(entry) {
  installed_version <- .bifrost_installed_version()
  required_version <- base::package_version(entry$minimum_bifrost_version)
  if (installed_version < required_version) {
    stop(
      "Example-data artifact '", entry$artifact_id, "' requires bifrost >= ",
      entry$minimum_bifrost_version, "; installed version is ",
      as.character(installed_version),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.bifrost_resolve_local_example <- function(name, root) {
  source <- "BIFROST_ARTIFACT_DIR"
  manifest_path <- file.path(root, "empirical-artifacts.json")
  manifest <- .bifrost_validate_example_manifest(
    .bifrost_read_example_manifest(manifest_path),
    source
  )
  entry <- .bifrost_example_entry(manifest, name)
  .bifrost_assert_example_version(entry)

  relative_path <- sub("^data-remote/", "", entry$path)
  candidate <- file.path(root, relative_path)
  if (!file.exists(candidate)) {
    .bifrost_verify_example_file(candidate, entry, source)
  }
  normalized_root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  normalized_path <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
  if (!startsWith(normalized_path, paste0(normalized_root, "/"))) {
    stop(
      "Example-data artifact '", entry$artifact_id, "' from ", source,
      " resolves outside the configured directory",
      call. = FALSE
    )
  }
  .bifrost_verify_example_file(normalized_path, entry, source)
  normalized_path
}

.bifrost_example_cache_dir <- function() {
  file.path(tools::R_user_dir("bifrost", "cache"), "example-data")
}

.bifrost_prepare_example_cache <- function(cache_dir) {
  root_link_target <- Sys.readlink(cache_dir)
  if (!is.na(root_link_target) && nzchar(root_link_target)) {
    stop(
      "Bifrost example-data cache directory must not be a symlink",
      call. = FALSE
    )
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(cache_dir)) {
    stop("Could not create the bifrost example-data cache directory", call. = FALSE)
  }
  root <- normalizePath(cache_dir, winslash = "/", mustWork = TRUE)
  child_paths <- vapply(c("artifacts", "manifests"), function(child) {
    path <- file.path(root, child)
    link_target <- Sys.readlink(path)
    if (!is.na(link_target) && nzchar(link_target)) {
      stop(
        "Bifrost example-data cache ", child,
        " directory must not be a symlink",
        call. = FALSE
      )
    }
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(path)) {
      stop(
        "Could not create the bifrost example-data cache ", child,
        " directory",
        call. = FALSE
      )
    }
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }, character(1))
  stats::setNames(c(root = root, child_paths), c("root", "artifacts", "manifests"))
}

.bifrost_example_cache_file_is_symlink <- function(path) {
  link_target <- Sys.readlink(path)
  !is.na(link_target) && nzchar(link_target)
}

.bifrost_assert_example_cache_file <- function(path, cache_root, child,
                                                allow_symlink = FALSE) {
  expected_parent <- file.path(cache_root, child)
  if (!identical(dirname(path), expected_parent)) {
    stop(
      "Bifrost example-data ", child,
      " cache file resolves outside the cache root",
      call. = FALSE
    )
  }
  if (.bifrost_example_cache_file_is_symlink(path)) {
    if (allow_symlink) {
      return(FALSE)
    }
    stop(
      "Bifrost example-data ", sub("s$", "", child),
      " cache file must not be a symlink",
      call. = FALSE
    )
  }
  TRUE
}

.bifrost_download_example_file <- function(url, destfile, quiet) {
  if (!.bifrost_is_scalar_string(url) || !startsWith(url, "https://")) {
    stop("Example-data downloads require an HTTPS URL", call. = FALSE)
  }

  # The transport adapter itself requires an external HTTPS service. Resolver
  # behavior is tested with an injected local transport instead. # nocov start
  previous_timeout <- getOption("timeout")
  timeout <- if (is.null(previous_timeout)) 60 else previous_timeout
  options(timeout = max(300, timeout))
  on.exit(options(timeout = previous_timeout), add = TRUE)

  status <- tryCatch(
    utils::download.file(url, destfile, quiet = quiet, mode = "wb"),
    error = function(error) {
      stop(
        "Could not download example-data file from ", url, ": ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
  if (!identical(as.integer(status), 0L)) {
    stop(
      "Could not download example-data file from ", url,
      ": download.file() returned status ", status,
      call. = FALSE
    )
  }
  0L
  # nocov end
}

.bifrost_example_artifact_cache_path <- function(cache_dir, entry) {
  file.path(
    cache_dir,
    "artifacts",
    paste0(entry$artifact_id, "--", entry$sha256, "--", basename(entry$path))
  )
}

.bifrost_cached_example_artifact_is_valid <- function(path, cache_dir, entry) {
  if (!.bifrost_assert_example_cache_file(
    path, cache_dir, "artifacts", allow_symlink = TRUE
  ) || !utils::file_test("-f", path)) {
    return(FALSE)
  }
  isTRUE(tryCatch({
    .bifrost_verify_example_file(path, entry, "example-data cache")
    TRUE
  }, error = function(error) FALSE))
}

.bifrost_example_manifest_cache_path <- function(cache_dir, manifest_path) {
  file.path(
    cache_dir,
    "manifests",
    paste0("manifest--", .bifrost_example_sha256(manifest_path), ".json")
  )
}

.bifrost_read_cached_example_manifests <- function(cache_dir) {
  cache_paths <- .bifrost_prepare_example_cache(cache_dir)
  cache_dir <- cache_paths[["root"]]
  manifest_dir <- cache_paths[["manifests"]]
  paths <- list.files(
    manifest_dir,
    pattern = "^manifest--[0-9a-f]{64}\\.json$",
    full.names = TRUE
  )
  if (length(paths) == 0L) {
    return(list())
  }
  paths <- paths[!vapply(
    paths, .bifrost_example_cache_file_is_symlink, logical(1)
  )]
  if (length(paths) == 0L) {
    return(list())
  }
  paths <- paths[order(file.info(paths)$mtime, decreasing = TRUE)]
  manifests <- lapply(paths, function(path) {
    .bifrost_assert_example_cache_file(path, cache_dir, "manifests")
    expected_path <- .bifrost_example_manifest_cache_path(cache_dir, path)
    if (!identical(path, expected_path)) {
      return(NULL)
    }
    tryCatch(
      .bifrost_validate_example_manifest(
        .bifrost_read_example_manifest(path),
        paste0("cached manifest ", path)
      ),
      error = function(error) NULL
    )
  })
  Filter(Negate(is.null), manifests)
}

.bifrost_resolve_cached_example <- function(name, cache_dir) {
  cache_paths <- .bifrost_prepare_example_cache(cache_dir)
  cache_dir <- cache_paths[["root"]]
  manifests <- .bifrost_read_cached_example_manifests(cache_dir)
  for (manifest in manifests) {
    entry <- .bifrost_example_entry(manifest, name)
    .bifrost_assert_example_version(entry)
    path <- .bifrost_example_artifact_cache_path(cache_dir, entry)
    if (.bifrost_cached_example_artifact_is_valid(path, cache_dir, entry)) {
      .bifrost_assert_example_cache_file(path, cache_dir, "artifacts")
      return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }
  }
  NULL
}

.bifrost_example_remote_failure <- function(name, source_url, cache_location,
                                            error) {
  stop(
    "Could not obtain bifrost example-data artifact '", name,
    "' from ", source_url, ". Expected cache location: ", cache_location,
    ". Set BIFROST_ARTIFACT_DIR to a verified local artifact directory. Details: ",
    conditionMessage(error),
    call. = FALSE
  )
}

.bifrost_prune_example_cache <- function(cache_dir, name, keep_artifact) {
  cache_paths <- .bifrost_prepare_example_cache(cache_dir)
  cache_dir <- cache_paths[["root"]]
  artifact_dir <- cache_paths[["artifacts"]]
  .bifrost_assert_example_cache_file(
    keep_artifact, cache_dir, "artifacts", allow_symlink = TRUE
  )
  artifact_paths <- list.files(artifact_dir, full.names = TRUE)
  matching_artifacts <- artifact_paths[startsWith(
    basename(artifact_paths), paste0(name, "--")
  )]
  lapply(matching_artifacts, .bifrost_assert_example_cache_file,
        cache_root = cache_dir, child = "artifacts", allow_symlink = TRUE)
  unlink(setdiff(matching_artifacts, keep_artifact), force = TRUE)
  invisible(NULL)
}

.bifrost_resolve_remote_example <- function(name, cache_dir, quiet) {
  cache_paths <- .bifrost_prepare_example_cache(cache_dir)
  cache_dir <- cache_paths[["root"]]
  artifact_dir <- cache_paths[["artifacts"]]
  manifest_dir <- cache_paths[["manifests"]]

  manifest_temp <- tempfile("manifest-", tmpdir = manifest_dir, fileext = ".json")
  on.exit(unlink(manifest_temp), add = TRUE)
  manifest <- tryCatch({
    .bifrost_download_example_file(
      .bifrost_example_manifest_url, manifest_temp, quiet = quiet
    )
    .bifrost_validate_example_manifest(
      .bifrost_read_example_manifest(manifest_temp),
      .bifrost_example_manifest_url
    )
  }, error = function(error) {
    .bifrost_example_remote_failure(
      name, .bifrost_example_manifest_url, cache_dir, error
    )
  })

  entry <- .bifrost_example_entry(manifest, name)
  .bifrost_assert_example_version(entry)

  artifact_path <- .bifrost_example_artifact_cache_path(cache_dir, entry)
  artifact_url <- paste0(.bifrost_example_raw_base, entry$path)
  artifact_is_valid <- .bifrost_cached_example_artifact_is_valid(
    artifact_path, cache_dir, entry
  )
  if (!artifact_is_valid) {
    artifact_temp <- tempfile("artifact-", tmpdir = artifact_dir)
    on.exit(unlink(artifact_temp), add = TRUE)
    tryCatch({
      .bifrost_download_example_file(artifact_url, artifact_temp, quiet = quiet)
      .bifrost_verify_example_file(artifact_temp, entry, artifact_url)
    }, error = function(error) {
      .bifrost_example_remote_failure(name, artifact_url, artifact_path, error)
    })
    if (.bifrost_cached_example_artifact_is_valid(
      artifact_path, cache_dir, entry
    )) {
      unlink(artifact_temp, force = TRUE)
    } else {
      if (.bifrost_example_cache_file_is_symlink(artifact_path) ||
          file.exists(artifact_path)) {
        unlink(artifact_path, force = TRUE)
      }
      # An atomic rename can fail for platform/filesystem reasons that cannot be
      # reproduced portably; preserve a clear failure at that boundary. # nocov start
      if (!file.rename(artifact_temp, artifact_path)) {
        stop(
          "Could not move verified example-data artifact '", entry$artifact_id,
          "' into expected cache path ", artifact_path,
          call. = FALSE
        )
      }
      # nocov end
    }
  }

  manifest_path <- .bifrost_example_manifest_cache_path(cache_dir, manifest_temp)
  manifest_is_valid <- .bifrost_assert_example_cache_file(
    manifest_path, cache_dir, "manifests", allow_symlink = TRUE
  ) && file.exists(manifest_path) && isTRUE(tryCatch(
    identical(
      .bifrost_example_sha256(manifest_path),
      .bifrost_example_sha256(manifest_temp)
    ),
    error = function(error) FALSE
  ))
  if (!manifest_is_valid) {
    if (.bifrost_example_cache_file_is_symlink(manifest_path) ||
        file.exists(manifest_path)) {
      unlink(manifest_path, force = TRUE)
    }
    # See the equivalent artifact rename boundary above. # nocov start
    if (!file.rename(manifest_temp, manifest_path)) {
      stop(
        "Could not move verified example-data manifest into expected cache path ",
        manifest_path,
        call. = FALSE
      )
    }
    # nocov end
  }

  .bifrost_prune_example_cache(cache_dir, entry$artifact_id, artifact_path)
  .bifrost_assert_example_cache_file(artifact_path, cache_dir, "artifacts")
  normalizePath(artifact_path, winslash = "/", mustWork = TRUE)
}

#' Locate a verified empirical example-data artifact
#'
#' Resolves a named empirical artifact to a locally verified file. A directory
#' supplied through `BIFROST_ARTIFACT_DIR` always takes precedence, including
#' when `refresh = TRUE`. Ordinary cache hits work offline. Artifacts may change
#' between package releases, so use `refresh = TRUE` to retrieve the latest
#' manifest. The returned value is a path; this function does not deserialize
#' the file.
#'
#' @param name A registered example-data identifier.
#' @param refresh Whether to bypass the cached manifest and check the remote
#'   manifest for a newer artifact.
#' @param quiet Whether downloads should suppress progress output.
#'
#' @return A normalized path to the checksum-verified artifact.
#' @export
#'
#' @examples
#' \dontrun{
#' bifrost_example_file("jaw-tree")
#' bifrost_example_file("jaw-tree", refresh = TRUE)
#' }
bifrost_example_file <- function(name, refresh = FALSE, quiet = FALSE) {
  if (!.bifrost_is_scalar_string(name)) {
    stop("`name` must be a single, non-missing identifier", call. = FALSE)
  }
  if (!(name %in% .bifrost_example_identifiers)) {
    stop(
      "`name` must be a supported identifier: ",
      paste(.bifrost_example_identifiers, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.logical(refresh) || length(refresh) != 1L || is.na(refresh)) {
    stop("`refresh` must be a single, non-missing logical value", call. = FALSE)
  }
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("`quiet` must be a single, non-missing logical value", call. = FALSE)
  }

  local_root <- Sys.getenv("BIFROST_ARTIFACT_DIR", unset = "")
  if (nzchar(local_root)) {
    return(.bifrost_resolve_local_example(name, local_root))
  }

  cache_dir <- .bifrost_example_cache_dir()
  if (!refresh) {
    cached_path <- .bifrost_resolve_cached_example(name, cache_dir)
    if (!is.null(cached_path)) {
      return(cached_path)
    }
  }
  .bifrost_resolve_remote_example(name, cache_dir, quiet)
}
