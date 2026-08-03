example_data_fixture_contents <- function() {
  c(
    "jaw-tree" = "tree bytes",
    "jaw-landmarks" = "landmark bytes",
    "passerine-tree" = "passerine tree bytes",
    "passerine-traits" = "passerine trait bytes",
    "passerine-search" = "passerine search bytes",
    "passerine-sensitivity" = "passerine sensitivity bytes",
    "passerine-posthoc" = "passerine posthoc bytes",
    "simulation-preview-tables" = "simulation preview bytes"
  )
}

write_example_data_fixture <- function(
    root,
    contents = example_data_fixture_contents(),
    minimum_version = "0.2.0") {
  dir.create(file.path(root, "fixtures"), recursive = TRUE, showWarnings = FALSE)
  artifacts <- lapply(names(contents), function(id) {
    filename <- paste0(id, ".rds")
    local_path <- file.path(root, "fixtures", filename)
    writeBin(charToRaw(unname(contents[[id]])), local_path)
    list(
      artifact_id = id,
      path = file.path("data-remote", "fixtures", filename),
      sha256 = digest::digest(
        local_path,
        algo = "sha256",
        serialize = FALSE,
        file = TRUE
      ),
      size_bytes = unname(file.info(local_path)$size),
      minimum_bifrost_version = minimum_version,
      source_id = "fixture-source",
      source_location = filename,
      transformation = list(method = "Test fixture bytes.", script = "Test helper."),
      license_id = "fixture-license"
    )
  })
  manifest <- list(
    schema_version = 2L,
    description = "Downloader unit-test manifest.",
    scope = list(directories = "data-remote/fixtures", extensions = ".rds"),
    sources = stats::setNames(list(list(
      title = "Fixture", doi = "10.0000/fixture",
      url = "https://example.org/fixture", version = "1"
    )), "fixture-source"),
    license_records = stats::setNames(list(list(
      repository_license = "GPL (>= 2)", upstream_status = "Fixture",
      evidence_url = "https://example.org/license", reuse_action = "Test only."
    )), "fixture-license"),
    artifacts = artifacts
  )
  jsonlite::write_json(
    manifest,
    file.path(root, "empirical-artifacts.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )
  invisible(manifest)
}

read_example_data_fixture <- function(root) {
  jsonlite::read_json(
    file.path(root, "empirical-artifacts.json"),
    simplifyVector = FALSE
  )
}

write_example_data_manifest <- function(root, manifest) {
  jsonlite::write_json(
    manifest,
    file.path(root, "empirical-artifacts.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )
}

local_fixture_downloader <- function(remote, calls = NULL) {
  force(remote)
  function(url, destfile, quiet = FALSE) {
    if (!is.null(calls)) {
      calls$urls <- c(calls$urls, url)
      calls$quiet <- c(calls$quiet, quiet)
    }
    source_path <- if (identical(url, .bifrost_example_manifest_url)) {
      file.path(remote, "empirical-artifacts.json")
    } else {
      prefix <- .bifrost_example_raw_base
      if (!startsWith(url, prefix)) {
        stop("Unexpected fixture URL: ", url, call. = FALSE)
      }
      file.path(remote, sub("^data-remote/", "", sub(prefix, "", url, fixed = TRUE)))
    }
    if (!file.exists(source_path)) {
      stop("Fixture source does not exist: ", source_path, call. = FALSE)
    }
    dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(source_path, destfile, overwrite = TRUE)) {
      stop("Could not copy fixture source", call. = FALSE)
    }
    0L
  }
}

local_example_data_bindings <- function(...) {
  testthat::local_mocked_bindings(
    ..., .package = "bifrost", .env = parent.frame()
  )
  invisible(NULL)
}
