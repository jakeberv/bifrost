if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the package root.")
}

if (!requireNamespace("xml2", quietly = TRUE) ||
    !requireNamespace("rvest", quietly = TRUE)) {
  stop("Rebuilding the vignette cache requires xml2 and rvest.")
}

html_path <- file.path("vignettes", "simulation-study-vignette.html")
if (!file.exists(html_path)) {
  stop("Render vignettes/simulation-study-vignette.Rmd before rebuilding the cache.")
}

doc <- xml2::read_html(html_path)
tables <- rvest::html_elements(doc, "table.simulation-preview-table")
captions <- trimws(rvest::html_text(rvest::html_element(tables, "caption")))
table_data <- lapply(tables, rvest::html_table, trim = TRUE)
names(table_data) <- captions

if (length(table_data) != 8L) {
  stop("Expected 8 preview tables in the rendered vignette, found ", length(table_data))
}

cache_obj <- list(
  medium_case = list(
    null = table_data[[1]],
    gic_proportional = table_data[[2]],
    bic_proportional = table_data[[3]],
    gic_correlation = table_data[[4]],
    bic_correlation = table_data[[5]]
  ),
  tuning = list(
    gic = table_data[[6]],
    bic = table_data[[7]],
    selected = table_data[[8]]
  )
)

out_dir <- file.path("inst", "extdata", "simulation-study-cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "passerine_preview_tables.rds")

saveRDS(cache_obj, out_path, compress = "xz")
message("Wrote ", out_path)
