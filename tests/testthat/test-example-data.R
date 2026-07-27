read_fixture_bytes <- function(path) {
  rawToChar(readBin(path, "raw", n = file.info(path)$size))
}

testthat::test_that("loading bifrost is network- and cache-free", {
  source_root <- normalizePath(
    testthat::test_path("../.."), winslash = "/", mustWork = TRUE
  )
  source_mode <- file.exists(file.path(source_root, "R", "example-data.R"))
  script_path <- tempfile("bifrost-clean-load-", fileext = ".R")
  cache_root <- tempfile("bifrost-clean-cache-")
  on.exit(unlink(script_path), add = TRUE)
  on.exit(unlink(cache_root, recursive = TRUE), add = TRUE)

  load_command <- if (source_mode) {
    sprintf(
      "pkgload::load_all(%s, export_all = FALSE, helpers = FALSE, quiet = TRUE)",
      encodeString(source_root, quote = "\"")
    )
  } else {
    package_path <- system.file(package = "bifrost")
    c(
      sprintf(
        ".libPaths(c(%s, .libPaths()))",
        encodeString(dirname(package_path), quote = "\"")
      ),
      "loadNamespace(\"bifrost\")"
    )
  }
  writeLines(c(
    "Sys.unsetenv(\"BIFROST_ARTIFACT_DIR\")",
    sprintf(
      "Sys.setenv(R_USER_CACHE_DIR = %s)",
      encodeString(cache_root, quote = "\"")
    ),
    "cache <- file.path(tools::R_user_dir(\"bifrost\", \"cache\"), \"example-data\")",
    "stopifnot(!dir.exists(cache), !\"bifrost\" %in% loadedNamespaces())",
    load_command,
    "stopifnot(\"bifrost\" %in% loadedNamespaces(), !dir.exists(cache))"
  ), script_path)

  rscript <- file.path(R.home("bin"), "Rscript")
  if (.Platform$OS.type == "windows") {
    rscript <- paste0(rscript, ".exe")
  }
  output <- system2(
    rscript, c("--vanilla", shQuote(script_path)), stdout = TRUE, stderr = TRUE
  )
  testthat::expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
})

testthat::test_that("schema-2 manifests expose the exact public data contract", {
  root <- withr::local_tempdir()
  write_example_data_fixture(root)
  manifest <- read_example_data_fixture(root)

  testthat::expect_no_error(
    .bifrost_validate_example_manifest(manifest, "fixture")
  )
  testthat::expect_setequal(
    Filter(nzchar, vapply(
      manifest$artifacts,
      function(artifact) {
        if (is.null(artifact$artifact_id)) "" else artifact$artifact_id
      },
      character(1)
    )),
    .bifrost_example_identifiers
  )

  provenance <- list(
    path = "vignettes/example.png",
    sha256 = paste(rep("a", 64L), collapse = ""),
    source_id = "fixture-source",
    source_location = "example.png",
    transformation = list(method = "Fixture.", script = "Fixture."),
    license_id = "fixture-license"
  )
  manifest$artifacts <- c(manifest$artifacts, list(provenance))
  testthat::expect_no_error(
    .bifrost_validate_example_manifest(manifest, "fixture")
  )
  testthat::expect_identical(
    .bifrost_example_entry(manifest, "jaw-tree")$artifact_id,
    "jaw-tree"
  )
  manifest$artifacts[[length(manifest$artifacts)]]$sha256 <- "bad"
  testthat::expect_error(
    .bifrost_validate_example_manifest(manifest, "fixture"),
    "invalid SHA-256"
  )
})

testthat::test_that("manifest validation rejects malformed downloader records", {
  root <- withr::local_tempdir()
  write_example_data_fixture(root)
  original <- read_example_data_fixture(root)
  cases <- list(
    list("schema version", function(x) { x$schema_version <- 3L; x }),
    list("artifacts must be a list", function(x) { x$artifacts <- "bad"; x }),
    list("must be a record", function(x) { x$artifacts[[1L]] <- "bad"; x }),
    list("unsafe repository-relative path", function(x) {
      x$artifacts[[1L]]$path <- "/tmp/data.rds"; x
    }),
    list("downloader path must live under data-remote", function(x) {
      x$artifacts[[1L]]$path <- "vignettes/jaw-tree.rds"; x
    }),
    list("invalid lowercase hyphenated artifact_id", function(x) {
      x$artifacts[[1L]]$artifact_id <- "jaw_tree"; x
    }),
    list("duplicate artifact_id", function(x) {
      x$artifacts[[2L]]$artifact_id <- x$artifacts[[1L]]$artifact_id; x
    }),
    list("duplicate path", function(x) {
      x$artifacts[[2L]]$path <- x$artifacts[[1L]]$path; x
    }),
    list("invalid SHA-256", function(x) {
      x$artifacts[[1L]]$sha256 <- "bad"; x
    }),
    list("positive whole-number byte size", function(x) {
      x$artifacts[[1L]]$size_bytes <- 0; x
    }),
    list("invalid minimum_bifrost_version", function(x) {
      x$artifacts[[1L]]$minimum_bifrost_version <- "not a version"; x
    }),
    list("invalid minimum_bifrost_version", function(x) {
      x$artifacts[[1L]]$minimum_bifrost_version <- 1; x
    }),
    list("missing downloader field", function(x) {
      x$artifacts[[1L]]$size_bytes <- NULL; x
    }),
    list("supported identifiers: jaw-tree, jaw-landmarks", function(x) {
      x$artifacts[[1L]]$artifact_id <- "extra-artifact"; x
    })
  )
  for (case in cases) {
    testthat::expect_error(
      .bifrost_validate_example_manifest(case[[2L]](original), "fixture"),
      case[[1L]],
      info = case[[1L]]
    )
  }
})

testthat::test_that("the local override verifies all eight artifacts", {
  root <- withr::local_tempdir()
  write_example_data_fixture(root)
  withr::local_envvar(BIFROST_ARTIFACT_DIR = root)

  paths <- vapply(
    .bifrost_example_identifiers, bifrost_example_file, character(1)
  )
  testthat::expect_true(all(file.exists(paths)))
  normalized_root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  testthat::expect_true(all(startsWith(paths, normalized_root)))
  testthat::expect_identical(
    read_fixture_bytes(paths[["jaw-tree"]]), "tree bytes"
  )
})

testthat::test_that("public inputs fail before cache or network access", {
  withr::local_envvar(BIFROST_ARTIFACT_DIR = NA_character_)
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  local_example_data_bindings(
    .bifrost_example_cache_dir = function() stop("cache reached"),
    .bifrost_download_example_file = function(...) {
      calls$n <- calls$n + 1L
      stop("network reached")
    }
  )

  testthat::expect_error(bifrost_example_file("missing"), "supported identifier")
  testthat::expect_error(bifrost_example_file(NA_character_), "`name`")
  testthat::expect_error(bifrost_example_file("jaw-tree", refresh = NA), "`refresh`")
  testthat::expect_error(bifrost_example_file("jaw-tree", quiet = 1), "`quiet`")
  testthat::expect_identical(calls$n, 0L)
})

testthat::test_that("the local override rejects incompatible, changed, and escaped files", {
  root <- withr::local_tempdir()
  write_example_data_fixture(root)
  withr::local_envvar(BIFROST_ARTIFACT_DIR = root)

  jaw <- file.path(root, "fixtures", "jaw-tree.rds")
  writeBin(charToRaw("changed"), jaw)
  testthat::expect_error(bifrost_example_file("jaw-tree"), "checksum mismatch")

  write_example_data_fixture(root)
  unlink(jaw)
  testthat::expect_error(bifrost_example_file("jaw-tree"), "is missing")

  dir.create(jaw)
  testthat::expect_error(bifrost_example_file("jaw-tree"), "is missing")
  unlink(jaw, recursive = TRUE)

  write_example_data_fixture(root)
  manifest <- read_example_data_fixture(root)
  manifest$artifacts[[1L]]$size_bytes <- manifest$artifacts[[1L]]$size_bytes + 1L
  write_example_data_manifest(root, manifest)
  testthat::expect_error(bifrost_example_file("jaw-tree"), "byte-size mismatch")

  write_example_data_fixture(root, minimum_version = "9.9.9")
  local_example_data_bindings(
    .bifrost_installed_version = function() base::package_version("0.1.0")
  )
  testthat::expect_error(
    bifrost_example_file("jaw-tree"),
    "requires bifrost >= 9.9.9; installed version is 0.1.0"
  )
})

testthat::test_that("the local override rejects files escaping through symlinks", {
  testthat::skip_on_os("windows")
  root <- withr::local_tempdir()
  outside <- withr::local_tempdir()
  write_example_data_fixture(root)
  external <- file.path(outside, "jaw-tree.rds")
  writeBin(charToRaw("tree bytes"), external)
  jaw <- file.path(root, "fixtures", "jaw-tree.rds")
  unlink(jaw)
  testthat::skip_if_not(file.symlink(external, jaw), "symlinks unavailable")
  withr::local_envvar(BIFROST_ARTIFACT_DIR = root)

  testthat::expect_error(bifrost_example_file("jaw-tree"), "outside")
})

testthat::test_that("downloads populate a verified cache that works offline", {
  remote <- withr::local_tempdir()
  cache <- withr::local_tempdir()
  write_example_data_fixture(remote)
  calls <- new.env(parent = emptyenv())
  calls$urls <- character()
  calls$quiet <- logical()
  withr::local_envvar(BIFROST_ARTIFACT_DIR = NA_character_)
  local_example_data_bindings(
    .bifrost_example_cache_dir = function() cache,
    .bifrost_download_example_file = local_fixture_downloader(remote, calls)
  )

  path <- bifrost_example_file("jaw-tree", quiet = TRUE)
  testthat::expect_identical(read_fixture_bytes(path), "tree bytes")
  testthat::expect_length(calls$urls, 2L)
  testthat::expect_true(all(calls$quiet))

  calls$urls <- character()
  local_example_data_bindings(
    .bifrost_download_example_file = function(...) stop("network reached")
  )
  testthat::expect_identical(bifrost_example_file("jaw-tree"), path)
  testthat::expect_length(calls$urls, 0L)
})

testthat::test_that("cached artifacts retain the manifest version gate", {
  remote <- withr::local_tempdir()
  cache <- withr::local_tempdir()
  write_example_data_fixture(remote, minimum_version = "9.9.9")
  withr::local_envvar(BIFROST_ARTIFACT_DIR = NA_character_)

  local({
    local_example_data_bindings(
      .bifrost_example_cache_dir = function() cache,
      .bifrost_download_example_file = local_fixture_downloader(remote),
      .bifrost_installed_version = function() base::package_version("10.0.0")
    )
    bifrost_example_file("jaw-tree")
  })

  local({
    local_example_data_bindings(
      .bifrost_example_cache_dir = function() cache,
      .bifrost_download_example_file = function(...) stop("network reached"),
      .bifrost_installed_version = function() base::package_version("0.1.0")
    )
    testthat::expect_error(
      bifrost_example_file("jaw-tree"),
      "requires bifrost >= 9.9.9; installed version is 0.1.0"
    )
  })
})

testthat::test_that("refresh reuses unchanged bytes and updates changed bytes", {
  remote <- withr::local_tempdir()
  cache <- withr::local_tempdir()
  write_example_data_fixture(remote)
  calls <- new.env(parent = emptyenv())
  calls$urls <- character()
  calls$quiet <- logical()
  withr::local_envvar(BIFROST_ARTIFACT_DIR = NA_character_)
  local_example_data_bindings(
    .bifrost_example_cache_dir = function() cache,
    .bifrost_download_example_file = local_fixture_downloader(remote, calls)
  )

  old_path <- bifrost_example_file("jaw-tree")
  calls$urls <- character()
  testthat::expect_identical(
    bifrost_example_file("jaw-tree", refresh = TRUE), old_path
  )
  testthat::expect_identical(calls$urls, .bifrost_example_manifest_url)

  contents <- example_data_fixture_contents()
  contents[["jaw-tree"]] <- "new tree bytes"
  write_example_data_fixture(remote, contents = contents)
  new_path <- bifrost_example_file("jaw-tree", refresh = TRUE)
  testthat::expect_identical(read_fixture_bytes(new_path), "new tree bytes")
  testthat::expect_false(file.exists(old_path))
  testthat::expect_false(identical(old_path, new_path))
})

testthat::test_that("refreshing one artifact preserves unrelated offline cache entries", {
  remote <- withr::local_tempdir()
  cache <- withr::local_tempdir()
  write_example_data_fixture(remote)
  withr::local_envvar(BIFROST_ARTIFACT_DIR = NA_character_)

  old_paths <- local({
    local_example_data_bindings(
      .bifrost_example_cache_dir = function() cache,
      .bifrost_download_example_file = local_fixture_downloader(remote)
    )
    c(
      jaw = bifrost_example_file("jaw-tree"),
      passerine = bifrost_example_file("passerine-tree")
    )
  })
  contents <- example_data_fixture_contents()
  contents[["jaw-tree"]] <- "new tree bytes"
  contents[["passerine-tree"]] <- "new passerine tree bytes"
  write_example_data_fixture(remote, contents = contents)

  local({
    local_example_data_bindings(
      .bifrost_example_cache_dir = function() cache,
      .bifrost_download_example_file = local_fixture_downloader(remote)
    )
    bifrost_example_file("jaw-tree", refresh = TRUE)
  })
  testthat::expect_true(file.exists(old_paths[["passerine"]]))
  testthat::expect_false(file.exists(old_paths[["jaw"]]))
  testthat::expect_gte(
    length(list.files(file.path(cache, "manifests"), pattern = "^manifest--")),
    2L
  )

  offline_path <- local({
    local_example_data_bindings(
      .bifrost_example_cache_dir = function() cache,
      .bifrost_download_example_file = function(...) stop("network reached")
    )
    bifrost_example_file("passerine-tree")
  })
  testthat::expect_identical(offline_path, old_paths[["passerine"]])
  testthat::expect_identical(
    read_fixture_bytes(offline_path), "passerine tree bytes"
  )
})

testthat::test_that("corrupt cache entries are replaced by verified bytes", {
  remote <- withr::local_tempdir()
  cache <- withr::local_tempdir()
  write_example_data_fixture(remote)
  withr::local_envvar(BIFROST_ARTIFACT_DIR = NA_character_)
  local_example_data_bindings(
    .bifrost_example_cache_dir = function() cache,
    .bifrost_download_example_file = local_fixture_downloader(remote)
  )

  path <- bifrost_example_file("jaw-tree")
  writeBin(charToRaw("corrupt"), path)
  replacement <- bifrost_example_file("jaw-tree")
  testthat::expect_identical(replacement, path)
  testthat::expect_identical(read_fixture_bytes(path), "tree bytes")
})

testthat::test_that("failed refreshes preserve an older verified artifact", {
  remote <- withr::local_tempdir()
  cache <- withr::local_tempdir()
  write_example_data_fixture(remote)
  withr::local_envvar(BIFROST_ARTIFACT_DIR = NA_character_)

  old_path <- local({
    local_example_data_bindings(
      .bifrost_example_cache_dir = function() cache,
      .bifrost_download_example_file = local_fixture_downloader(remote)
    )
    bifrost_example_file("jaw-tree")
  })
  manifest <- read_example_data_fixture(remote)
  manifest$artifacts[[1L]]$sha256 <- paste(rep("0", 64L), collapse = "")
  write_example_data_manifest(remote, manifest)

  local({
    local_example_data_bindings(
      .bifrost_example_cache_dir = function() cache,
      .bifrost_download_example_file = local_fixture_downloader(remote)
    )
    testthat::expect_error(
      bifrost_example_file("jaw-tree", refresh = TRUE), "checksum mismatch"
    )
  })
  testthat::expect_true(file.exists(old_path))
  testthat::expect_identical(read_fixture_bytes(old_path), "tree bytes")

  local({
    local_example_data_bindings(
      .bifrost_example_cache_dir = function() cache,
      .bifrost_download_example_file = function(...) stop("network reached")
    )
    testthat::expect_identical(bifrost_example_file("jaw-tree"), old_path)
  })
})

testthat::test_that("cache directory symlinks are rejected before writes", {
  testthat::skip_on_os("windows")
  remote <- withr::local_tempdir()
  write_example_data_fixture(remote)
  withr::local_envvar(BIFROST_ARTIFACT_DIR = NA_character_)

  external <- withr::local_tempdir()
  sentinel <- file.path(external, "sentinel")
  writeLines("keep", sentinel)
  root_link <- tempfile("cache-link-")
  testthat::skip_if_not(file.symlink(external, root_link), "symlinks unavailable")
  local({
    local_example_data_bindings(
      .bifrost_example_cache_dir = function() root_link,
      .bifrost_download_example_file = local_fixture_downloader(remote)
    )
    testthat::expect_error(
      bifrost_example_file("jaw-tree"), "cache directory must not be a symlink"
    )
  })
  testthat::expect_identical(readLines(sentinel), "keep")

  cache <- withr::local_tempdir()
  unlink(file.path(cache, "artifacts"), recursive = TRUE)
  testthat::skip_if_not(
    file.symlink(external, file.path(cache, "artifacts")),
    "child symlinks unavailable"
  )
  local({
    local_example_data_bindings(
      .bifrost_example_cache_dir = function() cache,
      .bifrost_download_example_file = local_fixture_downloader(remote)
    )
    testthat::expect_error(
      bifrost_example_file("jaw-tree"),
      "cache artifacts directory must not be a symlink"
    )
  })
  testthat::expect_identical(readLines(sentinel), "keep")
})

testthat::test_that("cache setup and file guards fail closed", {
  remote <- withr::local_tempdir()
  write_example_data_fixture(remote)
  withr::local_envvar(BIFROST_ARTIFACT_DIR = NA_character_)
  testthat::expect_match(
    .bifrost_example_cache_dir(), "bifrost.*/example-data$"
  )

  cache_file <- tempfile("cache-file-")
  writeLines("occupied", cache_file)
  local({
    local_example_data_bindings(
      .bifrost_example_cache_dir = function() cache_file,
      .bifrost_download_example_file = local_fixture_downloader(remote)
    )
    testthat::expect_error(
      bifrost_example_file("jaw-tree"), "Could not create.*cache directory"
    )
  })

  cache <- withr::local_tempdir()
  writeLines("occupied", file.path(cache, "artifacts"))
  local({
    local_example_data_bindings(
      .bifrost_example_cache_dir = function() cache,
      .bifrost_download_example_file = local_fixture_downloader(remote)
    )
    testthat::expect_error(
      bifrost_example_file("jaw-tree"), "Could not create.*cache artifacts"
    )
  })

  guarded <- .bifrost_prepare_example_cache(withr::local_tempdir())
  testthat::expect_error(
    .bifrost_assert_example_cache_file(
      tempfile(tmpdir = guarded[["root"]]), guarded[["root"]], "artifacts"
    ),
    "outside the cache root"
  )

  if (.Platform$OS.type != "windows") {
    external <- tempfile("external-")
    writeLines("keep", external)
    link <- file.path(guarded[["artifacts"]], "linked.rds")
    if (file.symlink(external, link)) {
      testthat::expect_error(
        .bifrost_assert_example_cache_file(
          link, guarded[["root"]], "artifacts"
        ),
        "must not be a symlink"
      )
    }
  }
})

testthat::test_that("symlinked artifact cache entries are replaced safely", {
  testthat::skip_on_os("windows")
  remote <- withr::local_tempdir()
  cache <- withr::local_tempdir()
  outside <- withr::local_tempdir()
  write_example_data_fixture(remote)
  withr::local_envvar(BIFROST_ARTIFACT_DIR = NA_character_)
  local_example_data_bindings(
    .bifrost_example_cache_dir = function() cache,
    .bifrost_download_example_file = local_fixture_downloader(remote)
  )

  path <- bifrost_example_file("jaw-tree")
  external <- file.path(outside, "artifact.rds")
  file.copy(path, external)
  unlink(path)
  testthat::skip_if_not(file.symlink(external, path), "symlinks unavailable")
  replacement <- bifrost_example_file("jaw-tree")

  testthat::expect_identical(replacement, path)
  testthat::expect_identical(Sys.readlink(path), "")
  testthat::expect_identical(read_fixture_bytes(external), "tree bytes")
})

testthat::test_that("unusable cached manifests are ignored", {
  testthat::skip_on_os("windows")
  remote <- withr::local_tempdir()
  cache <- withr::local_tempdir()
  outside <- withr::local_tempdir()
  write_example_data_fixture(remote)
  cache_paths <- .bifrost_prepare_example_cache(cache)
  source <- file.path(remote, "empirical-artifacts.json")

  wrong_name <- file.path(
    cache_paths[["manifests"]],
    paste0("manifest--", paste(rep("0", 64L), collapse = ""), ".json")
  )
  file.copy(source, wrong_name)
  testthat::expect_length(.bifrost_read_cached_example_manifests(cache), 0L)
  unlink(wrong_name)

  external <- file.path(outside, "manifest.json")
  file.copy(source, external)
  linked <- .bifrost_example_manifest_cache_path(cache, source)
  testthat::skip_if_not(file.symlink(external, linked), "symlinks unavailable")
  testthat::expect_length(.bifrost_read_cached_example_manifests(cache), 0L)
})

testthat::test_that("a concurrent verified artifact fill is reused", {
  remote <- withr::local_tempdir()
  cache <- withr::local_tempdir()
  write_example_data_fixture(remote)
  cache_paths <- .bifrost_prepare_example_cache(cache)
  manifest <- read_example_data_fixture(remote)
  entry <- .bifrost_example_entry(manifest, "jaw-tree")
  artifact_path <- .bifrost_example_artifact_cache_path(
    cache_paths[["root"]], entry
  )
  source <- file.path(remote, sub("^data-remote/", "", entry$path))
  downloader <- local_fixture_downloader(remote)
  interleaved_downloader <- function(url, destfile, quiet = FALSE) {
    status <- downloader(url, destfile, quiet)
    if (!identical(url, .bifrost_example_manifest_url)) {
      file.copy(source, artifact_path, overwrite = TRUE)
    }
    status
  }
  withr::local_envvar(BIFROST_ARTIFACT_DIR = NA_character_)
  local_example_data_bindings(
    .bifrost_example_cache_dir = function() cache,
    .bifrost_download_example_file = interleaved_downloader
  )

  testthat::expect_identical(
    bifrost_example_file("jaw-tree", refresh = TRUE), artifact_path
  )
  testthat::expect_identical(read_fixture_bytes(artifact_path), "tree bytes")
})

testthat::test_that("a symlinked manifest cache entry is replaced safely", {
  testthat::skip_on_os("windows")
  remote <- withr::local_tempdir()
  cache <- withr::local_tempdir()
  outside <- withr::local_tempdir()
  write_example_data_fixture(remote)
  cache_paths <- .bifrost_prepare_example_cache(cache)
  source <- file.path(remote, "empirical-artifacts.json")
  manifest_path <- .bifrost_example_manifest_cache_path(cache, source)
  external <- file.path(outside, "manifest.json")
  writeLines("keep", external)
  testthat::skip_if_not(
    file.symlink(external, manifest_path), "symlinks unavailable"
  )
  withr::local_envvar(BIFROST_ARTIFACT_DIR = NA_character_)
  local_example_data_bindings(
    .bifrost_example_cache_dir = function() cache,
    .bifrost_download_example_file = local_fixture_downloader(remote)
  )

  testthat::expect_identical(
    bifrost_example_file("jaw-tree", refresh = TRUE),
    .bifrost_example_artifact_cache_path(
      cache_paths[["root"]], .bifrost_example_entry(
        read_example_data_fixture(remote), "jaw-tree"
      )
    )
  )
  testthat::expect_identical(Sys.readlink(manifest_path), "")
  testthat::expect_identical(readLines(external), "keep")
})

testthat::test_that("the transport rejects non-HTTPS URLs before download", {
  destination <- withr::local_tempfile()
  testthat::expect_error(
    .bifrost_download_example_file("http://example.org/file", destination, TRUE),
    "require an HTTPS URL"
  )
})

testthat::test_that("invalid remote manifests fail without cached artifacts", {
  remote <- withr::local_tempdir()
  cache <- withr::local_tempdir()
  dir.create(remote, recursive = TRUE, showWarnings = FALSE)
  writeLines("not json", file.path(remote, "empirical-artifacts.json"))
  withr::local_envvar(BIFROST_ARTIFACT_DIR = NA_character_)
  local_example_data_bindings(
    .bifrost_example_cache_dir = function() cache,
    .bifrost_download_example_file = local_fixture_downloader(remote)
  )

  testthat::expect_error(
    bifrost_example_file("jaw-tree"),
    "Could not obtain bifrost example-data artifact.*Expected cache location"
  )
  testthat::expect_length(list.files(file.path(cache, "artifacts")), 0L)
})
