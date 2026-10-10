# Compatibility initialization for the historical bifrost vignette baseline.
# Adapted from .startGuess(), .rate_guess(), and BM/BMM input preparation in
# mvMORPH 1.2.1 development code, copyright Julien Clavel, licensed GPL (>= 2).
# Source: https://github.com/JClavel/mvMORPH/tree/023134e993c8b174cf716378503892fdc4d6616d
# These initializers also match the experimental fork used for the cached
# vignettes. They include the random second-tip safeguard added after the CRAN 1.2.1 release.
# Keep this module separate so it can be removed when native initialization is
# adopted. The installed mvMORPH still evaluates the grid and fits the model.

.bifrost_mvgls <- function(formula, data = list(), tree, model,
                          method = "PL-LOOCV", REML = TRUE, ...,
                          start_strategy = c("legacy_1.2.1", "native")) {
  start_strategy <- match.arg(start_strategy)
  args <- list(...)
  used <- start_strategy
  if (!is.null(args$start)) {
    used <- "supplied"
  } else if (!is.null(args$grid.search) && !args$grid.search) {
    used <- "native_no_grid"
  } else if (start_strategy == "legacy_1.2.1") {
    args$start <- .bifrost_mvgls_start_historical(
      formula, data, tree, model, method, REML, args
    )
  }
  fit_args <- list(formula = formula, tree = tree, model = model)
  # Preserve omitted arguments in mvgls's saved call as well as their defaults.
  if (!missing(data)) fit_args$data <- data
  if (!missing(method)) fit_args$method <- method
  if (!missing(REML)) fit_args$REML <- REML
  fit <- do.call(mvMORPH::mvgls, c(fit_args, args))
  attr(fit, "bifrost_initialization") <- list(
    requested = start_strategy, used = used,
    mvMORPH_version = as.character(utils::packageVersion("mvMORPH"))
  )
  fit
}

# This is the only bridge to mvMORPH's private fitting machinery. Fail visibly
# on an incompatible interface rather than silently change the starting policy.
.bifrost_mvgls_start_backend <- function() {
  required <- list(
    .setBounds = c("penalty", "model", "lower", "upper", "tol", "mserr",
                   "penalized", "corrModel", "k"),
    .loocvPhylo = c("par", "cvmethod", "targM", "corrStr", "penalty", "error", "nobs")
  )
  backend <- lapply(names(required), function(name) {
    fun <- get0(name, envir = asNamespace("mvMORPH"), inherits = FALSE)
    if (!is.function(fun) || !all(required[[name]] %in% names(formals(fun)))) {
      stop("This mvMORPH version is incompatible with historical initialization; ",
           "use start_strategy = 'native'.", call. = FALSE)
    }
    fun
  })
  names(backend) <- names(required)
  backend
}

.bifrost_mvgls_start_historical <- function(formula, data, tree, model, method, REML, args) {
  method <- match.arg(method[1L], c("PL-LOOCV", "LOOCV", "LL", "H&L",
                                   "Mahalanobis", "EmpBayes"))
  if (method == "PL-LOOCV") method <- "LOOCV"
  penalty <- if (is.null(args$penalty)) "RidgeArch" else args$penalty
  if (!model %in% c("BM", "BMM") || method == "EmpBayes" ||
      !penalty %in% c("RidgeArch", "RidgeAlt", "LASSO")) {
    stop("This model, method, or penalty is unsupported by historical initialization; ",
         "use start_strategy = 'native' for this fit.",
         call. = FALSE)
  }
  if (!inherits(tree, "phylo") || (model == "BMM" && !inherits(tree, "simmap"))) {
    stop("Historical BM/BMM initialization requires a phylo tree (simmap for BMM).",
         call. = FALSE)
  }
  if (method %in% c("H&L", "Mahalanobis") && penalty != "RidgeArch") {
    stop("H&L and Mahalanobis require RidgeArch penalization.", call. = FALSE)
  }
  data <- data[!vapply(data, inherits, logical(1), "phylo")]
  frame <- stats::model.frame(formula = formula, data = data)
  X <- stats::model.matrix(attr(frame, "terms"), data = frame,
                           contrasts.arg = args$contrasts)
  Y <- if (is.null(args$response)) stats::model.response(frame) else args$response
  if (nrow(frame) != length(tree$tip.label) || !is.matrix(Y) ||
      ncol(Y) < 2L || anyNA(Y)) {
    stop("Historical initialization requires complete multivariate datasets with one row per tip.",
         call. = FALSE)
  }
  if (all(rownames(frame) %in% tree$tip.label)) {
    X <- X[tree$tip.label, , drop = FALSE]
    Y <- Y[tree$tip.label, , drop = FALSE]
  }
  # mvgls assumes tip order when row names are absent; attach those names for
  # the subtree calculations without changing the user's data.
  rownames(X) <- rownames(Y) <- tree$tip.label
  qrx <- qr(X)
  X <- X[, qrx$pivot[seq_len(qrx$rank)], drop = FALSE]
  mserr <- if (isFALSE(args$error)) NULL else args$error
  if (method == "LL" && nrow(Y) < ncol(Y)) {
    stop("There are more variables than observations; use a penalized method.", call. = FALSE)
  }
  root <- if (is.null(args$root)) "stationary" else args$root
  precalc <- list(randomRoot = if (is.null(args$randomRoot)) TRUE else args$randomRoot,
                  root = root, root_std = if (root == "stationary") 1L else 0L)
  if (!is.null(args$scale.height) && args$scale.height) {
    tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  }
  corr <- list(Y = Y, X = X, REML = REML, mserr = mserr, model = model,
               structure = tree, p = ncol(Y), nobs = nrow(Y), m = qrx$rank,
               nloo = seq_len(nrow(Y)), precalc = precalc)
  backend <- .bifrost_mvgls_start_backend()
  corr$bounds <- backend$.setBounds(
    penalty = penalty, model = model, lower = args$lower, upper = args$upper,
    tol = args$tol, mserr = mserr, penalized = method != "LL", corrModel = corr,
    k = if (model == "BMM") ncol(tree$mapped.edge) else 1L
  )
  .bifrost_mvgls_start_grid_historical(
    corr, method, penalty, args$target, args$tol, backend$.loocvPhylo
  )
}

.bifrost_mvgls_rate_historical <- function(tree, Y, X) {
  transform <- mvMORPH::pruning(tree, trans = FALSE)$sqrtM
  X <- crossprod(transform, X)
  Y <- crossprod(transform, Y)
  residuals <- Y - X %*% (corpcor::pseudoinverse(X) %*% Y)
  residuals
}

.bifrost_mvgls_regime_starts_historical <- function(tree, Y, X) {
  terminal <- tree$edge[, 2L] %in% seq_len(ape::Ntip(tree))
  maps <- vapply(tree$maps[terminal], function(x) names(x)[length(x)], character(1))
  k <- ncol(tree$mapped.edge)
  if (length(unique(maps)) < k) {
    # Preserve the original whole-tree PIC rule when a regime has no tips.
    rate <- crossprod(apply(Y, 2L, ape::pic, phy = tree)) / ape::Ntip(tree)
    return(as.list(rep(sqrt(mean(diag(rate))), k)))
  }
  lapply(colnames(tree$mapped.edge), function(regime) {
    tips <- which(maps == regime)
    removed <- tree$tip.label[!tree$tip.label %in% tree$tip.label[tips]]
    # Preserve the historical safeguard and its RNG draw: a singleton needs a
    # second species to estimate a rate. Do this before pruning the tree.
    if (ape::Ntip(tree) - length(removed) <= 1L) {
      removed <- removed[-sample(length(removed), size = 1)]
    }
    subtree <- ape::drop.tip(tree, removed)
    if (ape::Ntip(subtree) <= 1L) subtree <- tree
    sqrt(mean(apply(.bifrost_mvgls_rate_historical(
      subtree, Y[subtree$tip.label, , drop = FALSE],
      X[subtree$tip.label, , drop = FALSE]
    ), 2L, stats::var)))
  })
}

.bifrost_mvgls_start_grid_historical <- function(corr, method, penalty, target, tol, score) {
  tuning <- NULL
  if (method != "LL") {
    tuning <- switch(penalty,
      RidgeArch = c(1e-6, .01, .05, .1, .15, .2, .3, .5, .7, .9),
      RidgeAlt = log(c(1e-11, 1e-9, 1e-6, .01, .1, 1, 10, 100, 1000, 10000)),
      LASSO = log(c(1e-6, .01, .1, 1, 10, 100, 1000))
    )
    if (!is.null(tol)) {
      cutoff <- if (penalty == "RidgeArch") tol else log(tol)
      tuning <- tuning[tuning > cutoff]
    }
    if (!length(tuning)) stop("tol excludes every historical initialization grid point.", call. = FALSE)
  }
  parameters <- list(tuning, 1) # Preserve the original dummy BM parameter.
  if (corr$model == "BMM") {
    rates <- .bifrost_mvgls_regime_starts_historical(corr$structure, corr$Y, corr$X)[-1L]
    parameters <- c(list(tuning), rates)
  }
  if (!is.null(corr$mserr)) {
    error_grid <- c(.001, .01, .1, 1, 10)
    if (corr$model == "BMM") error_grid <- sqrt(error_grid * mean(unlist(rates)^2))
    parameters <- c(parameters, list(error_grid))
  }
  parameters <- parameters[!vapply(parameters, is.null, logical(1))]
  grid <- expand.grid(parameters)
  values <- apply(grid, 1L, score, cvmethod = method,
                  targM = if (is.null(target)) "unitVariance" else target,
                  corrStr = corr, penalty = penalty, error = corr$mserr,
                  nobs = corr$nobs)
  best <- which.min(values)
  if (!length(best)) stop("No usable historical starting values were found.", call. = FALSE)
  as.numeric(grid[best, ])
}
