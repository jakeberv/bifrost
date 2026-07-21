#!/usr/bin/env Rscript

source_repo <- Sys.getenv(
  "BODYPLAN_SOURCE_REPO",
  unset = file.path("..", "passerine-bodyplan-evolution")
)
source_repo_id <- "jakeberv/passerine-bodyplan-evolution"
manuscript_script_id <- "TemporalAnalyses.R"
object_dir <- file.path(source_repo, "my_session_objects")
output_path <- file.path(
  "inst",
  "extdata",
  "avian-skeleton",
  "passerine_bodyplan_posthoc_integration_compact.RDS"
)
pca_min_tips <- as.numeric(Sys.getenv("BODYPLAN_POSTHOC_PCA_MIN_TIPS", unset = "10"))
fit_min_tips <- as.numeric(Sys.getenv("BODYPLAN_POSTHOC_FIT_MIN_TIPS", unset = "2"))
high_corr_threshold <- as.numeric(Sys.getenv("BODYPLAN_POSTHOC_CORR_THRESHOLD", unset = "0.95"))

if (!is.finite(pca_min_tips) || pca_min_tips < 1) {
  stop("BODYPLAN_POSTHOC_PCA_MIN_TIPS must be a finite number >= 1.", call. = FALSE)
}
if (!is.finite(fit_min_tips) || fit_min_tips < 1) {
  stop("BODYPLAN_POSTHOC_FIT_MIN_TIPS must be a finite number >= 1.", call. = FALSE)
}
if (!is.finite(high_corr_threshold) ||
    high_corr_threshold <= 0 ||
    high_corr_threshold >= 1) {
  stop("BODYPLAN_POSTHOC_CORR_THRESHOLD must be between 0 and 1.", call. = FALSE)
}

read_cached <- function(name) {
  path <- file.path(object_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Cached manuscript object not found: ", path, call. = FALSE)
  }
  readRDS(path)
}


upper_pair_names <- function(trait_names) {
  idx <- which(upper.tri(matrix(NA, length(trait_names), length(trait_names))),
               arr.ind = TRUE)
  paste(trait_names[idx[, 1]], trait_names[idx[, 2]], sep = "__")
}

mean_upper <- function(mat) {
  mean(mat[upper.tri(mat, diag = FALSE)], na.rm = TRUE)
}

mean_between <- function(mat, rows, cols) {
  mean(mat[rows, cols, drop = FALSE], na.rm = TRUE)
}

runs_with_posthoc <- read_cached("runs_with_posthoc")
vars_cors <- read_cached("all_vars_cors")
pgls_fit <- read_cached("posthoc_integration_pgls")

coef_se <- sqrt(diag(pgls_fit$vcov))
coef_t <- pgls_fit$coefficients / coef_se
coef_df <- pgls_fit$n - length(pgls_fit$coefficients)
coef_p <- 2 * stats::pt(abs(coef_t), df = coef_df, lower.tail = FALSE)
pgls_summary <- list(
  run = "min20.ic20.gic",
  coefficients = data.frame(
    term = names(pgls_fit$coefficients),
    estimate = unname(pgls_fit$coefficients),
    std_error = unname(coef_se),
    statistic = unname(coef_t),
    p_value = unname(coef_p),
    row.names = NULL,
    check.names = FALSE
  ),
  r_squared = unname(pgls_fit$r.squared),
  adj_r_squared = unname(pgls_fit$adj.r.squared),
  n = unname(pgls_fit$n),
  df_model = length(pgls_fit$coefficients) - 1L,
  df_residual = coef_df,
  model = pgls_fit$model,
  formula = paste(deparse(pgls_fit$formula), collapse = " "),
  settings = list(
    pgls_phylogeny = "collapsePhylogenyByStates(min20.ic20.gic$tree_no_uncertainty_untransformed)",
    remove_high_corr = TRUE,
    corr_threshold = high_corr_threshold
  )
)

posthoc <- runs_with_posthoc$min10.ic20.gic$posthoc
tip_counts <- vapply(
  posthoc,
  function(fit) {
    if (is.null(fit$dims$n)) NA_real_ else as.numeric(fit$dims$n)
  },
  numeric(1)
)
has_matrix <- vapply(
  posthoc,
  function(fit) !is.null(fit$sigma$Pinv) && is.matrix(fit$sigma$Pinv),
  logical(1)
)
regime_ids <- names(posthoc)[!is.na(tip_counts) & tip_counts > pca_min_tips & has_matrix]

correlation_matrices <- lapply(
  posthoc[regime_ids],
  function(fit) stats::cov2cor(fit$sigma$Pinv)
)

regime_ages <- vapply(
  posthoc[regime_ids],
  function(fit) {
    phy <- fit$corrSt$phy
    if (is.null(phy)) {
      return(NA_real_)
    }
    max(phytools::nodeHeights(phy), na.rm = TRUE)
  },
  numeric(1)
)

internal_trait_names <- rownames(correlation_matrices[[1L]])
trait_labels <- c(
  tarsus = "Tibiotarsus",
  metatarsus = "Tarsometatarsus",
  femur = "Femur",
  humerus = "Humerus",
  ulna = "Ulna",
  radius = "Radius",
  carpometacarpus = "Carpometacarpus",
  second_digit = "2nd digit phalanx",
  cv.keel.1 = "Keel",
  cv.furcula.1 = "Furcula",
  sclerotic_ring = "Sclerotic ring",
  cv.skull.1 = "Skull-bill length"
)
trait_labels <- trait_labels[internal_trait_names]

vectorized_matrix <- do.call(
  rbind,
  lapply(correlation_matrices, function(mat) mat[upper.tri(mat, diag = FALSE)])
)
rownames(vectorized_matrix) <- regime_ids
colnames(vectorized_matrix) <- upper_pair_names(unname(trait_labels))

correlation_pca <- stats::prcomp(vectorized_matrix, center = TRUE, scale. = TRUE)
variance_explained <- correlation_pca$sdev^2 / sum(correlation_pca$sdev^2)

modules <- list(
  wing_flight = c(
    "Carpometacarpus",
    "2nd digit phalanx",
    "Radius",
    "Ulna",
    "Humerus",
    "Keel",
    "Furcula"
  ),
  hindlimb = c("Tibiotarsus", "Tarsometatarsus", "Femur"),
  cranial = c("Skull-bill length", "Sclerotic ring")
)
modules$hindlimb_cranial <- c(modules$hindlimb, modules$cranial)

module_scores <- data.frame(
  regime_id = regime_ids,
  within_wing = NA_real_,
  hindcran_wing = NA_real_,
  cranial_hindlimb = NA_real_,
  cranial_wing = NA_real_
)
for (i in seq_along(correlation_matrices)) {
  mat <- correlation_matrices[[i]]
  dimnames(mat) <- list(unname(trait_labels), unname(trait_labels))

  wing <- rownames(mat) %in% modules$wing_flight
  hindcran <- rownames(mat) %in% modules$hindlimb_cranial
  cranial <- rownames(mat) %in% modules$cranial
  hindlimb <- rownames(mat) %in% modules$hindlimb

  module_scores$within_wing[i] <- mean_upper(mat[wing, wing, drop = FALSE])
  module_scores$hindcran_wing[i] <- mean_between(mat, hindcran, wing)
  module_scores$cranial_hindlimb[i] <- mean_between(mat, cranial, hindlimb)
  module_scores$cranial_wing[i] <- mean_between(mat, cranial, wing)
}

module_diagnostics <- list(
  scores = module_scores,
  correlations = data.frame(
    diagnostic = c(
      "cor(PC1, within-wing mean)",
      "cor(PC1, hindlimb+cranial vs wing mean)",
      "cor(PC2, cranial vs hindlimb mean)",
      "cor(PC2, cranial vs wing mean)"
    ),
    value = c(
      stats::cor(correlation_pca$x[, 1], module_scores$within_wing),
      stats::cor(correlation_pca$x[, 1], module_scores$hindcran_wing),
      stats::cor(correlation_pca$x[, 2], module_scores$cranial_hindlimb),
      stats::cor(correlation_pca$x[, 2], module_scores$cranial_wing)
    )
  )
)

compact <- list(
  vars_cors = vars_cors,
  pgls_summary = pgls_summary,
  correlation_matrices = correlation_matrices,
  pca_inputs = list(
    vectorized_matrix = vectorized_matrix,
    regime_ids = regime_ids,
    regime_ages = regime_ages,
    tip_counts = tip_counts[regime_ids],
    min_n = pca_min_tips,
    use_correlation = TRUE
  ),
  pca_summary = list(
    variance_explained = variance_explained,
    scores = correlation_pca$x,
    loadings = correlation_pca$rotation
  ),
  trait_labels = trait_labels,
  modules = modules,
  module_diagnostics = module_diagnostics,
  settings = list(
    fit_min_tips = fit_min_tips,
    pca_min_tips = pca_min_tips,
    pca_tip_filter = "tip_count > pca_min_tips",
    remove_high_corr = TRUE,
    corr_threshold = high_corr_threshold
  ),
  provenance = list(
    generated_on = as.character(Sys.Date()),
    source_repo_path = source_repo_id,
    source_object_names = c(
      "runs_with_posthoc",
      "all_vars_cors",
      "posthoc_integration_pgls"
    ),
    manuscript_script = manuscript_script_id,
    manuscript_lines = list(
      posthoc_block = "3144-3183",
      pca_block = "3184-3444",
      module_diagnostics = "3445-3545"
    ),
    notes = c(
      "Compact object stores independent post-hoc regime matrices, not search$VCVs.",
      "search$VCVs are proportional joint-model matrices from the bifrost shift model.",
      "Spatial local assemblage covariance analyses are intentionally excluded."
    )
  )
)
attr(compact, "provenance") <- compact$provenance

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(compact, output_path, version = 3)

size_mb <- file.info(output_path)$size / 1024^2
message(sprintf("Wrote %s (%.3f MB)", output_path, size_mb))
message(sprintf(
  "PCA PC1-PC4 variance: %s",
  paste(round(100 * variance_explained[1:4], 2), collapse = ", ")
))
message(sprintf(
  "Module diagnostics: %s",
  paste(round(module_diagnostics$correlations$value, 3), collapse = ", ")
))
