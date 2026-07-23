reader_visible_r_chunks_in_vignette <- function(path) {
  text <- readLines(path, warn = FALSE)
  starts <- grep("^```\\{r(?:[ ,}])", text, perl = TRUE)

  chunks <- lapply(starts, function(start) {
    ends <- which(text == "```" & seq_along(text) > start)
    end <- ends[[1L]]
    header <- text[[start]]
    body <- text[seq.int(start + 1L, end - 1L)]
    body <- body[nzchar(trimws(body))]

    list(header = header, body = body)
  })

  chunks[!vapply(chunks, function(chunk) {
    grepl(
      "(?:echo|include)\\s*=\\s*FALSE",
      chunk$header,
      ignore.case = TRUE,
      perl = TRUE
    )
  }, logical(1L))]
}

test_that("reader-visible vignette chunks open with intent comments", {
  vignette_dir <- testthat::test_path("../../vignettes")
  paths <- sort(list.files(
    vignette_dir,
    pattern = "\\.Rmd$",
    full.names = TRUE,
    recursive = TRUE
  ))
  testthat::skip_if(
    length(paths) == 0L,
    "Vignette sources are unavailable in installed tests"
  )

  missing_comments <- unlist(lapply(paths, function(path) {
    chunks <- reader_visible_r_chunks_in_vignette(path)
    missing <- vapply(chunks, function(chunk) {
      length(chunk$body) == 0L || !grepl("^# ", chunk$body[[1L]])
    }, logical(1L))

    sprintf(
      "%s: %s",
      basename(path),
      vapply(chunks[missing], `[[`, character(1L), "header")
    )
  }), use.names = FALSE)

  testthat::expect_equal(missing_comments, character())
})
