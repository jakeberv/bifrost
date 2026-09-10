#!/usr/bin/env Rscript

# The website defaults to light; GitHub's README follows the reader's theme.
# Match the initial page theme before the browser can fetch dark alternatives.
# Keep both sources so the existing page-theme switching continues to work.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 1L)
path <- if (length(args)) args[[1L]] else "docs/index.html"
html <- xml2::read_html(path)
pictures <- xml2::xml_find_all(html, paste0(
  "//picture[",
  "contains(concat(' ', normalize-space(@class), ' '), ' cran-downloads-picture ')",
  " or contains(concat(' ', normalize-space(@class), ' '), ' schmidt-sciences-picture ')",
  "]"
))
stopifnot(length(pictures) == 2L,
          all(lengths(xml2::xml_find_all(pictures, "./img", flatten = FALSE)) == 1L))
sources <- xml2::xml_find_all(pictures, "./source[@data-theme]")
stopifnot(length(sources) == 4L)
theme <- xml2::xml_attr(xml2::xml_find_first(html, "/html"), "data-bs-theme")
if (is.na(theme)) theme <- "light"
stopifnot(theme %in% c("light", "dark"))
xml2::xml_set_attr(sources, "media", ifelse(
  xml2::xml_attr(sources, "data-theme") == theme, "all", "not all"
))
xml2::write_html(html, path)
