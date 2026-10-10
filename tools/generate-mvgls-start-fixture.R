# Regenerate against the historical mvMORPH initializer, without installing it.
# From the repository root:
# Rscript tools/generate-mvgls-start-fixture.R /path/to/historical/mvMORPH
# The historical fork's .startGuess and .rate_guess match upstream commit
# 023134e993c8b174cf716378503892fdc4d6616d (2025-08-25).
# Reference source: https://github.com/JClavel/mvMORPH/tree/023134e993c8b174cf716378503892fdc4d6616d
# The source directory supplies the reference functions, not a package dependency.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
source_dir <- normalizePath(args[1L])
suppressPackageStartupMessages(library(mvMORPH))
RhpcBLASctl::blas_set_num_threads(1L); RhpcBLASctl::omp_set_num_threads(1L)
reference <- new.env(parent = asNamespace('mvMORPH'))
for (path in list.files(file.path(source_dir, 'R'), pattern = '\\.[rR]$', full.names = TRUE)) {
  sys.source(path, reference)
}
stopifnot(is.function(reference$mvgls), is.function(reference$.startGuess),
          is.function(reference$.rate_guess))
source_hashes <- vapply(c('.startGuess', '.rate_guess'), function(name) {
  f <- get(name, reference)
  digest::digest(list(formals(f), body(f)), algo = 'sha256')
}, character(1))
expected_hashes <- c(
  .startGuess = 'a9b8e1b89be0205d8174ef944315c4fa10df1af999a75c8571b3220fd2143147',
  .rate_guess = '958ccb623233abf9300df7a40aff500e13504b58a2aef75c0cc9a0380f118bee'
)
if (!identical(source_hashes, expected_hashes)) {
  stop('Reference initializers differ from the pinned upstream commit.')
}
tr <- ape::stree(16, type='balanced'); tr$edge.length <- rep(.25,nrow(tr$edge))
paint <- function(kind) {
  tree <- tr; states <- rep('0', nrow(tree$edge)); node <- 18L
  if(kind != 'BM') states[tree$edge[,2] %in% c(node,phytools::getDescendants(tree,node))] <- '1'
  if(kind == 'singleton') states[tree$edge[,2] == 1L] <- '2'
  if(kind == 'tipless') states[tree$edge[,2] == node] <- '2'
  labs <- sort(unique(states))
  tree$maps <- lapply(seq_along(states),function(i)setNames(tree$edge.length[i],states[i]))
  tree$mapped.edge <- matrix(0,nrow(tree$edge),length(labs),dimnames=list(NULL,labs))
  tree$mapped.edge[cbind(seq_along(states),match(states,labs))] <- tree$edge.length
  class(tree) <- c('simmap','phylo'); tree
}
# Pin the fixture's RNG configuration across R versions (no binomial draws).
suppressWarnings(RNGversion("4.4.2"))
set.seed(620)
Y <- t(chol(ape::vcv(tr))) %*% matrix(rnorm(16*20),16,20) * .2
rownames(Y) <- tr$tip.label
x <- rnorm(16); group <- factor(rep(c('a','b'),8))
data <- data.frame(y1=Y[,1]+x,y2=Y[,2]-x,x=x,group=group,row.names=tr$tip.label)
cases <- list(
  bm_hl=list(kind='BM',method='H&L',error=TRUE),
  bm_ll=list(kind='BM',method='LL',error=TRUE,formula='cbind(y1, y2) ~ x + group'),
  bm_ll_noerror=list(kind='BM',method='LL',error=FALSE,formula='cbind(y1, y2) ~ x'),
  tolerance_hl=list(kind='multi',method='H&L',error=TRUE,tol=.2),
  variance_target=list(kind='multi',method='H&L',error=TRUE,target='Variance'),
  bmm_hl=list(kind='multi',method='H&L',error=TRUE),
  singleton_hl=list(kind='singleton',method='H&L',error=TRUE),
  tipless_hl=list(kind='tipless',method='H&L',error=TRUE),
  formula_ll=list(kind='multi',method='LL',error=TRUE,formula='cbind(y1, y2) ~ x + group'),
  no_intercept_ll=list(kind='multi',method='LL',error=FALSE,formula='cbind(y1, y2) ~ 0 + x'),
  loocv_ridgealt=list(kind='multi',method='LOOCV',error=TRUE,penalty='RidgeAlt',formula='cbind(y1, y2) ~ x'),
  loocv_lasso=list(kind='multi',method='LOOCV',error=FALSE,penalty='LASSO',formula='cbind(y1, y2) ~ x'),
  mahalanobis=list(kind='multi',method='Mahalanobis',error=FALSE),
  ml_hl=list(kind='multi',method='H&L',error=TRUE,REML=FALSE),
  scaled_height=list(kind='multi',method='H&L',error=TRUE,scale.height=TRUE)
)
result <- lapply(names(cases),function(name) {
  z <- cases[[name]]; kind <- z$kind; z$kind <- NULL
  f <- if(is.null(z$formula)) 'Y ~ 1' else z$formula
  z$formula <- stats::as.formula(f); environment(z$formula) <- environment()
  if(f != 'Y ~ 1') z$data <- data
  z$tree <- paint(kind); z$model <- if(kind=='BM') 'BM' else 'BMM'
  if (name == 'scaled_height') {
    z$tree$edge.length <- z$tree$edge.length * 3
    z$tree$maps <- lapply(z$tree$maps, function(x) x * 3)
    z$tree$mapped.edge <- z$tree$mapped.edge * 3
  }
  set.seed(71)
  rng_before <- .Random.seed
  fit <- suppressWarnings(do.call(reference$mvgls,z))
  rng_after <- .Random.seed
  z$formula <- f
  if(f == 'Y ~ 1') z$data <- list(Y=Y)
  cat(name, as.numeric(fit$start_values), '\n')
  list(args=z,start=as.numeric(fit$start_values),objective=fit$opt$value,
       parameters=as.numeric(fit$opt$par),convergence=fit$opt$convergence,
       rng_before=rng_before,rng_after=rng_after)
})
names(result) <- names(cases)
saveRDS(list(source='historical mvMORPH development initializer',
             upstream_commit='023134e993c8b174cf716378503892fdc4d6616d',
             source_hashes=source_hashes,cases=result),
        'tests/testthat/fixtures/mvgls-starts-historical.rds',version=2)
