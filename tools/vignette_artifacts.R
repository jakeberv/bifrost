#!/usr/bin/env Rscript

find_repo_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(path, "_pkgdown.yml"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find repository root from ", getwd(), call. = FALSE)
    }
    path <- parent
  }
}

repo_root <- find_repo_root()
setwd(repo_root)

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for vignette artifact tooling.", call. = FALSE)
  }
}

require_pkg("digest")
require_pkg("jsonlite")

normalize_rel <- function(path) {
  path <- gsub("\\\\", "/", path)
  sub("^\\./", "", path)
}

repo_rel <- function(path) {
  if (length(path) == 0L) {
    return(character())
  }
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- paste0(repo_root, "/")
  out <- ifelse(startsWith(path, root), sub(root, "", path, fixed = TRUE), path)
  normalize_rel(out)
}

list_repo_files <- function(path) {
  if (!dir.exists(path)) {
    return(character())
  }
  files <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files[file.exists(files) & !dir.exists(files)]
}

discover_vignette_slugs <- function() {
  slugs <- character()
  parsed <- tryCatch(yaml::read_yaml("_pkgdown.yml"), error = function(e) NULL)
  if (!is.null(parsed)) {
    collect_contents <- function(x) {
      out <- character()
      if (is.list(x)) {
        if (!is.null(x$contents)) {
          out <- c(out, unlist(x$contents, use.names = FALSE))
        }
        out <- c(out, unlist(lapply(x, collect_contents), use.names = FALSE))
      }
      out
    }
    slugs <- collect_contents(parsed$articles)
    slugs <- basename(slugs)
    slugs <- sub("\\.Rmd$", "", slugs, ignore.case = TRUE)
    slugs <- sub("\\.html$", "", slugs, ignore.case = TRUE)
    slugs <- unique(slugs[nzchar(slugs)])
    slugs <- slugs[file.exists(file.path("vignettes", paste0(slugs, ".Rmd")))]
  }

  if (length(slugs) == 0L) {
    files <- list.files("vignettes", pattern = "\\.Rmd$", full.names = FALSE)
    slugs <- sub("\\.Rmd$", "", files, ignore.case = TRUE)
  }

  sort(unique(slugs))
}

extract_call_paths <- function(text, call_pattern) {
  pattern <- paste0(call_pattern, "\\s*\\(\\s*(['\"])([^'\"]+)\\1")
  hits <- gregexpr(pattern, text, perl = TRUE)
  matches <- regmatches(text, hits)
  unlist(lapply(matches, function(x) {
    if (length(x) == 0L || identical(x, character(0))) {
      return(character())
    }
    sub(pattern, "\\2", x, perl = TRUE)
  }), use.names = FALSE)
}

extract_extension_paths <- function(text) {
  pattern <- "['\"]\\.?/?([^'\"]+\\.(?:R|csv|gif|jpe?g|png|rds|RDS|svg|txt))['\"]"
  hits <- gregexpr(pattern, text, perl = TRUE, ignore.case = TRUE)
  matches <- regmatches(text, hits)
  unlist(lapply(matches, function(x) {
    if (length(x) == 0L || identical(x, character(0))) {
      return(character())
    }
    sub(pattern, "\\1", x, perl = TRUE, ignore.case = TRUE)
  }), use.names = FALSE)
}

resolve_vignette_paths <- function(path) {
  path <- normalize_rel(path)
  if (!nzchar(path) ||
      grepl("^(https?:)?//", path) ||
      grepl("^mailto:", path) ||
      grepl("^data:", path) ||
      grepl("^#", path)) {
    return(NA_character_)
  }

  candidates <- c(
    file.path("vignettes", path),
    path
  )
  existing <- candidates[file.exists(candidates) & !dir.exists(candidates)]
  if (length(existing) > 0L) {
    return(repo_rel(existing[1L]))
  }

  unresolved <- repo_rel(candidates)
  unresolved[!startsWith(unresolved, "/") & !startsWith(unresolved, "../")]
}

vignette_dependencies <- function(slug) {
  rmd <- file.path("vignettes", paste0(slug, ".Rmd"))
  if (!file.exists(rmd)) {
    stop("Missing vignette source: ", rmd, call. = FALSE)
  }
  text <- paste(readLines(rmd, warn = FALSE), collapse = "\n")

  direct <- c(
    extract_call_paths(text, "(?:source|readRDS|readLines|read\\.csv|utils::read\\.csv|knitr::include_graphics|include_graphics|knitr::image_uri|image_uri)"),
    extract_extension_paths(text)
  )

  rate_map_files <- extract_call_paths(
    text,
    "(?:rate_map_artifact|rate_map_print_file|rate_map_table|rate_map_table_html|rate_map_table_data|rate_map_save_arc_figure)"
  )
  direct <- c(direct, file.path("rate-map-jaw-shape", rate_map_files))

  deps <- unique(na.omit(unlist(lapply(direct, resolve_vignette_paths),
                               use.names = FALSE)))
  sort(unique(c(rmd, deps)))
}

shared_dependencies <- function() {
  files <- c(
    "DESCRIPTION",
    "_pkgdown.yml",
    "tools/vignette_artifacts.R",
    "tools/render-vignette-pdf.R",
    "tools/build-colab-notebook.py",
    "tools/test-vignette-artifacts.py",
    repo_rel(list_repo_files("R")),
    repo_rel(list_repo_files("inst/extdata"))
  )
  sort(unique(files[file.exists(files) & !dir.exists(files)]))
}

hash_file <- function(path) {
  digest::digest(file = path, algo = "sha256")
}

slug_hashes <- function(slugs = discover_vignette_slugs()) {
  shared <- shared_dependencies()
  stats::setNames(lapply(slugs, function(slug) {
    files <- sort(unique(c(vignette_dependencies(slug), shared)))
    files <- files[file.exists(files) & !dir.exists(files)]
    hashes <- vapply(files, hash_file, character(1))
    digest::digest(paste(c(files, hashes), collapse = "\n"), algo = "sha256")
  }), slugs)
}

manifest_path <- function(cache_dir) {
  file.path(cache_dir, "manifest.json")
}

read_manifest <- function(cache_dir) {
  path <- manifest_path(cache_dir)
  if (!file.exists(path)) {
    return(list(slugs = list()))
  }
  jsonlite::read_json(path, simplifyVector = FALSE)
}

write_manifest <- function(cache_dir, hashes) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  slugs <- lapply(names(hashes), function(slug) {
    list(hash = unname(hashes[[slug]]))
  })
  names(slugs) <- names(hashes)
  jsonlite::write_json(
    list(updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
         slugs = slugs),
    manifest_path(cache_dir),
    auto_unbox = TRUE,
    pretty = TRUE
  )
}

artifact_path <- function(slug, target, artifact_dir) {
  ext <- switch(
    target,
    pdf = ".pdf",
    colab = ".ipynb",
    stop("Unknown target: ", target, call. = FALSE)
  )
  file.path(artifact_dir, paste0(slug, ext))
}

git_changed_files <- function(base) {
  if (is.null(base) || !nzchar(base)) {
    return(character())
  }
  out <- tryCatch(
    suppressWarnings(system2(
      "git",
      c("diff", "--name-only", paste0(base, "...HEAD")),
      stdout = TRUE,
      stderr = TRUE
    )),
    error = function(e) {
      stop("Could not compute changed vignette inputs: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop(
      "Could not compute changed vignette inputs against '", base, "':\n",
      paste(out, collapse = "\n"),
      call. = FALSE
    )
  }
  normalize_rel(out[nzchar(out)])
}

shared_changed <- function(files) {
  shared <- shared_dependencies()
  any(files %in% shared) ||
    any(startsWith(files, "R/")) ||
    any(startsWith(files, "inst/extdata/")) ||
    any(startsWith(files, ".github/workflows/")) ||
    any(files %in% c("DESCRIPTION", "_pkgdown.yml"))
}

changed_by_git <- function(slugs, base, target = NULL, artifact_dir = NULL) {
  files <- git_changed_files(base)
  if (length(files) == 0L) {
    return(character())
  }
  if (shared_changed(files)) {
    return(slugs)
  }
  unique(unlist(lapply(slugs, function(slug) {
    deps <- vignette_dependencies(slug)
    if (!is.null(target) && !is.null(artifact_dir)) {
      deps <- c(deps, repo_rel(artifact_path(slug, target, artifact_dir)))
    }
    if (any(files %in% deps)) {
      slug
    } else {
      character()
    }
  }), use.names = FALSE))
}

changed_by_manifest <- function(slugs, cache_dir) {
  manifest <- read_manifest(cache_dir)
  hashes <- slug_hashes(slugs)
  unique(unlist(lapply(slugs, function(slug) {
    cached <- manifest$slugs[[slug]]$hash
    if (is.null(cached) || !identical(unname(hashes[[slug]]), cached)) {
      slug
    } else {
      character()
    }
  }), use.names = FALSE))
}

missing_artifacts <- function(slugs, target, artifact_dir) {
  slugs[!file.exists(vapply(slugs, artifact_path, character(1),
                            target = target, artifact_dir = artifact_dir))]
}

plan_slugs <- function(target, artifact_dir, cache_dir, base = NULL) {
  slugs <- discover_vignette_slugs()
  changed <- if (!is.null(base) && nzchar(base)) {
    changed_by_git(slugs, base, target, artifact_dir)
  } else {
    changed_by_manifest(slugs, cache_dir)
  }
  sort(unique(c(changed, missing_artifacts(slugs, target, artifact_dir))))
}

plan_changed_slugs <- function(target, artifact_dir, cache_dir, base = NULL,
                               ignore_missing = FALSE) {
  slugs <- plan_slugs(target, artifact_dir, cache_dir, base)
  if (!isTRUE(ignore_missing)) {
    return(slugs)
  }

  all_slugs <- discover_vignette_slugs()
  changed <- if (!is.null(base) && nzchar(base)) {
    changed_by_git(all_slugs, base, target, artifact_dir)
  } else {
    changed_by_manifest(all_slugs, cache_dir)
  }
  sort(unique(changed))
}

copy_cached_pdfs <- function(pdf_dir, cache_dir) {
  slugs <- discover_vignette_slugs()
  cache_pdf_dir <- file.path(cache_dir, "pdfs")
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
  copied <- character()
  for (slug in slugs) {
    src <- file.path(cache_pdf_dir, paste0(slug, ".pdf"))
    dest <- file.path(pdf_dir, paste0(slug, ".pdf"))
    if (file.exists(src)) {
      file.copy(src, dest, overwrite = TRUE)
      copied <- c(copied, slug)
    }
  }
  copied
}

update_pdf_cache <- function(pdf_dir, cache_dir) {
  slugs <- discover_vignette_slugs()
  cache_pdf_dir <- file.path(cache_dir, "pdfs")
  dir.create(cache_pdf_dir, recursive = TRUE, showWarnings = FALSE)
  for (slug in slugs) {
    src <- file.path(pdf_dir, paste0(slug, ".pdf"))
    dest <- file.path(cache_pdf_dir, paste0(slug, ".pdf"))
    if (file.exists(src)) {
      file.copy(src, dest, overwrite = TRUE)
    } else {
      warning("PDF not found for cache update: ", src, call. = FALSE)
    }
  }
  write_manifest(cache_dir, slug_hashes(slugs))
}

opt_value <- function(args, name, default = NULL) {
  hit <- which(args == name)
  if (length(hit) == 0L) {
    return(default)
  }
  if (hit[1L] == length(args)) {
    stop("Missing value for ", name, call. = FALSE)
  }
  args[[hit[1L] + 1L]]
}

has_opt <- function(args, name) {
  any(args == name)
}

emit_slugs <- function(slugs, sep = "\n") {
  if (identical(sep, "space")) {
    cat(paste(slugs, collapse = " "))
    cat("\n")
  } else {
    writeLines(slugs)
  }
}

main <- function(args) {
  if (length(args) == 0L) {
    stop("Usage: vignette_artifacts.R <slugs|plan|restore-pdfs|update-cache> [options]",
         call. = FALSE)
  }
  cmd <- args[[1L]]
  rest <- args[-1L]
  cache_dir <- opt_value(rest, "--cache-dir", ".cache/vignette-artifacts")
  sep <- opt_value(rest, "--sep", "\n")

  if (identical(cmd, "slugs")) {
    emit_slugs(discover_vignette_slugs(), sep)
    return(invisible(NULL))
  }

  if (identical(cmd, "plan")) {
    target <- opt_value(rest, "--target", "pdf")
    artifact_dir <- opt_value(
      rest,
      "--artifact-dir",
      if (identical(target, "colab")) "vignettes/colab" else "docs/articles"
    )
    base <- opt_value(rest, "--base", "")
    emit_slugs(
      plan_changed_slugs(
        target,
        artifact_dir,
        cache_dir,
        base,
        ignore_missing = has_opt(rest, "--ignore-missing")
      ),
      sep
    )
    return(invisible(NULL))
  }

  if (identical(cmd, "restore-pdfs")) {
    pdf_dir <- opt_value(rest, "--pdf-dir", "docs/articles")
    copied <- copy_cached_pdfs(pdf_dir, cache_dir)
    message("Restored cached PDFs for ", length(copied), " vignette(s).")
    return(invisible(NULL))
  }

  if (identical(cmd, "update-cache")) {
    pdf_dir <- opt_value(rest, "--pdf-dir", "docs/articles")
    update_pdf_cache(pdf_dir, cache_dir)
    message("Updated vignette artifact cache.")
    return(invisible(NULL))
  }

  stop("Unknown command: ", cmd, call. = FALSE)
}

if (sys.nframe() == 0L) {
  main(commandArgs(trailingOnly = TRUE))
}
