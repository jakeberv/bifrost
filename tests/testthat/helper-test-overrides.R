local_rebind <- function(name, value, env) {
  had_local_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- get(name, envir = env, inherits = TRUE)
  was_locked <- had_local_binding && bindingIsLocked(name, env)
  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  withr::defer({
    if (had_local_binding) {
      assign(name, old_value, envir = env)
      if (was_locked) {
        lockBinding(name, env)
      }
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
  }, envir = parent.frame())
}

local_search_rebind <- local_rebind

serial_future_lapply <- function(X, FUN, ...) {
  dots <- list(...)
  future_arguments <- startsWith(names(dots), "future.")
  dots <- dots[!future_arguments]
  do.call(lapply, c(list(X = X, FUN = FUN), dots))
}
