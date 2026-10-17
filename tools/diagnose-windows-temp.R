# Run each cancellation mode in a separate R session/temp root on Windows.
root <- normalizePath("output/windows-temp-investigation", mustWork = TRUE)
for (mode in c("interrupt", "cooperative")) {
  work <- file.path(root, mode)
  dir.create(work, showWarnings = FALSE)
  script <- file.path(root, paste0(mode, ".R"))
  writeLines(c(
    'library(future)',
    'plan(multisession, workers = 2L)',
    'a <- future({Sys.sleep(0.05); stop("synthetic worker failure")})',
    'b <- future({Sys.sleep(1); Sys.getpid()})',
    'try(value(a), silent = TRUE)',
    sprintf('cancel(b, interrupt = %s)', if (mode == "interrupt") "TRUE" else "FALSE"),
    'cat("Cancelled state:", b$state, "\\n")',
    'plan(sequential)'
  ), script)
  vars <- c("TMPDIR", "TMP", "TEMP")
  old <- Sys.getenv(vars, unset = NA_character_)
  do.call(Sys.setenv, as.list(setNames(rep(work, 3L), vars)))
  status <- system2(file.path(R.home("bin"), "Rscript.exe"), c("--vanilla", shQuote(script)))
  for (v in vars) if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, setNames(list(old[[v]]), v))
  Sys.sleep(2)
  cat("MODE:", mode, "status:", status, "remaining files:\n")
  print(list.files(work, all.files = TRUE, recursive = TRUE))
}
