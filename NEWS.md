# bifrost 0.1.1

* Addressed CRAN reviewer feedback following review of 0.1.0:
  - Replaced all uses of shorthand `T`/`F` with `TRUE`/`FALSE`.
  - Ensured all informational output is suppressible via `message()`/`warning()` and controlled by a `verbose` flag.
  - Redirected all on-disk output generated during model fitting to `tempdir()` to comply with CRAN file system policies and avoid writing to the user’s working directory.
  - Ensured graphical parameters and global options are restored using immediate `on.exit()` calls.
  - Refined parallelization behavior to be CRAN-safe and cross-platform:
    - Parallel candidate evaluation uses `future` with `multicore` on Unix outside RStudio and `multisession` otherwise.
    - BLAS/OpenMP threads are capped to one per worker during parallel execution to avoid CPU oversubscription.
    - Sequential execution remains the default when `num_cores = 1`.
  - Improved documentation clarity around parallel execution, verbosity, and model fit history storage.

# bifrost 0.1.0

* Initial CRAN submission.
