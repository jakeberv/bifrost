simulationFormulaToObject <- function(formula) {
  if (inherits(formula, "formula")) {
    return(formula)
  }
  if (is.character(formula) && length(formula) == 1L && !is.na(formula)) {
    return(stats::as.formula(formula))
  }
  stop("formula must be a single character string or a formula object.")
}

simulationFormulaToString <- function(formula) {
  paste(deparse(simulationFormulaToObject(formula)), collapse = "")
}

resolveSimulationColumns <- function(columns, data_names, what) {
  if (is.null(columns)) {
    return(NULL)
  }
  if (length(columns) == 0L) {
    return(integer(0))
  }
  validateSimulationIdentifiers(data_names, "trait_data column names")
  if (is.numeric(columns)) {
    if (any(!is.finite(columns)) || any(columns != floor(columns))) {
      stop(sprintf("%s must contain whole-number column positions.", what))
    }
    if (any(columns < 1) || any(columns > length(data_names))) {
      stop(sprintf("%s contains invalid column positions.", what))
    }
    return(unique(as.integer(columns)))
  }
  if (is.character(columns)) {
    resolve_one <- function(column) {
      by_name <- match(column, data_names)
      if (!is.na(by_name)) {
        return(as.integer(by_name))
      }
      if (!is.na(column) && grepl("^[+]?[0-9]+$", column)) {
        by_position <- suppressWarnings(as.numeric(column))
        if (is.finite(by_position) &&
            by_position >= 1 &&
            by_position <= length(data_names)) {
          return(as.integer(by_position))
        }
      }
      NA_integer_
    }
    idx <- vapply(columns, resolve_one, integer(1))
    if (anyNA(idx)) {
      stop(sprintf("%s contains names not found in trait_data.", what))
    }
    return(unique(idx))
  }
  stop(sprintf("%s must be NULL, numeric positions, or character column names.", what))
}

validateSimulationIdentifiers <- function(identifiers, what) {
  if (is.null(identifiers) ||
      anyNA(identifiers) ||
      any(!nzchar(trimws(as.character(identifiers))))) {
    stop(sprintf("%s must be non-empty and unique.", what))
  }
  if (anyDuplicated(as.character(identifiers))) {
    stop(sprintf("%s must be unique.", what))
  }
  invisible(identifiers)
}

isSimulationInterceptOnly <- function(formula_obj) {
  formula_terms <- stats::terms(formula_obj)
  length(attr(formula_terms, "term.labels")) == 0L &&
    identical(attr(formula_terms, "intercept"), 1L)
}

validateSimulationStudyFormula <- function(formula) {
  formula_obj <- simulationFormulaToObject(formula)
  if (!isSimulationInterceptOnly(formula_obj)) {
    stop(
      "Simulation studies currently support intercept-only search formulas only. ",
      "Use formula = \"trait_data ~ 1\" (or an equivalent intercept-only response formula)."
    )
  }
  formula
}

isLegacyTraitDataReference <- function(expr) {
  is.symbol(expr) && identical(as.character(expr), "trait_data")
}

isLegacyTraitDataSubset <- function(expr) {
  is.call(expr) &&
    identical(expr[[1L]], as.name("[")) &&
    length(expr) >= 2L &&
    is.symbol(expr[[2L]]) &&
    identical(as.character(expr[[2L]]), "trait_data")
}

isSupportedNamedResponse <- function(expr) {
  if (is.symbol(expr) && !identical(as.character(expr), "trait_data")) {
    return(TRUE)
  }
  if (is.call(expr) && identical(expr[[1L]], as.name("cbind"))) {
    args <- as.list(expr)[-1L]
    return(length(args) > 0L && all(vapply(args, is.symbol, logical(1))))
  }
  FALSE
}

containsSimulationDot <- function(expr) {
  if (is.symbol(expr) && identical(as.character(expr), ".")) {
    return(TRUE)
  }
  if (!is.call(expr)) {
    return(FALSE)
  }
  any(vapply(as.list(expr), containsSimulationDot, logical(1)))
}

extractNamedResponseNames <- function(lhs_expr, data_names) {
  if (is.symbol(lhs_expr)) {
    nm <- as.character(lhs_expr)
    if (!nm %in% data_names) {
      stop("Response variable referenced in formula was not found in trait_data.")
    }
    return(nm)
  }

  if (is.call(lhs_expr) && identical(lhs_expr[[1L]], as.name("cbind"))) {
    response_names <- vapply(as.list(lhs_expr)[-1L], as.character, character(1))
    if (!all(response_names %in% data_names)) {
      stop("One or more response variables referenced in formula were not found in trait_data.")
    }
    if (anyDuplicated(response_names)) {
      stop("Response variables in formula must be unique.")
    }
    return(response_names)
  }

  stop("Unsupported response form in formula.")
}

resolveLegacyTraitDataColumns <- function(expr, data_names) {
  placeholder <- as.data.frame(
    setNames(
      lapply(seq_along(data_names), function(i) rep.int(i, 3L)),
      data_names
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  value <- tryCatch(
    eval(expr, envir = list(trait_data = placeholder)),
    error = function(e) e
  )
  if (inherits(value, "error")) {
    stop("Legacy trait_data indexing could not be resolved from the supplied formula.")
  }

  if (is.data.frame(value) || is.matrix(value)) {
    cols <- colnames(value)
    if (is.null(cols) || !all(cols %in% data_names)) {
      stop("Legacy trait_data indexing could not be resolved to named columns.")
    }
    return(match(cols, data_names))
  }

  unique_vals <- unique(as.integer(value))
  unique_vals <- unique_vals[!is.na(unique_vals)]
  if (length(unique_vals) == 0L || any(unique_vals < 1L) || any(unique_vals > length(data_names))) {
    stop("Legacy trait_data indexing could not be resolved to valid columns.")
  }
  unique_vals
}

collectLegacyTraitDataColumns <- function(expr, data_names) {
  if (isLegacyTraitDataSubset(expr)) {
    return(resolveLegacyTraitDataColumns(expr, data_names))
  }
  if (!is.call(expr)) {
    return(integer(0))
  }
  unique(unlist(lapply(as.list(expr)[-1L], collectLegacyTraitDataColumns, data_names = data_names)))
}

rewriteLegacyTraitDataExpr <- function(expr, data_names, side = c("lhs", "rhs")) {
  side <- match.arg(side)

  if (isLegacyTraitDataReference(expr)) {
    if (side == "lhs") {
      cols <- seq_along(data_names)
      names <- data_names[cols]
      if (length(names) == 1L) {
        return(as.name(names))
      }
      return(as.call(c(as.name("cbind"), lapply(names, as.name))))
    }
    return(expr)
  }

  if (isLegacyTraitDataSubset(expr)) {
    cols <- resolveLegacyTraitDataColumns(expr, data_names)
    names <- data_names[cols]
    if (side == "rhs" && length(names) != 1L) {
      stop("Legacy indexed predictors must resolve to single raw columns.")
    }
    if (length(names) == 1L) {
      return(as.name(names))
    }
    return(as.call(c(as.name("cbind"), lapply(names, as.name))))
  }

  if (!is.call(expr)) {
    return(expr)
  }

  as.call(lapply(as.list(expr), rewriteLegacyTraitDataExpr, data_names = data_names, side = side))
}

predictorSchemaForData <- function(data, predictor_names) {
  if (length(predictor_names) == 0L) {
    return(list())
  }

  setNames(lapply(predictor_names, function(name) {
    column <- data[[name]]
    if (is.character(column)) {
      stop(
        "Character predictors are not supported in simulation templates. ",
        "Convert them to factors explicitly before fitting."
      )
    }
    if (is.factor(column)) {
      return(list(
        type = if (is.ordered(column)) "ordered" else "factor",
        levels = levels(column)
      ))
    }
    if (is.logical(column)) {
      return(list(type = "logical", levels = NULL))
    }
    if (is.numeric(column) || is.integer(column)) {
      return(list(type = "numeric", levels = NULL))
    }
    stop(
      "Predictors in simulation templates must be numeric, logical, factor, or ordered factor."
    )
  }), predictor_names)
}

normalizeSimulationFormulaSpec <- function(formula,
                                           trait_data,
                                           response_columns = NULL,
                                           predictor_columns = NULL) {
  formula_obj <- simulationFormulaToObject(formula)
  formula_original_chr <- paste(deparse(formula_obj), collapse = "")
  data_names <- colnames(trait_data)

  response_override <- resolveSimulationColumns(response_columns, data_names, "response_columns")
  predictor_override <- resolveSimulationColumns(predictor_columns, data_names, "predictor_columns")

  if (!is.null(response_override) && length(response_override) == 0L) {
    stop("response_columns must identify at least one response column.")
  }
  if (!is.null(predictor_override) && length(predictor_override) == 0L) {
    stop("predictor_columns must identify at least one predictor column when supplied.")
  }

  if (containsSimulationDot(formula_obj[[3L]])) {
    stop("'.' shorthand is not supported in simulation formulas.")
  }

  lhs_expr <- formula_obj[[2L]]
  rhs_expr <- formula_obj[[3L]]
  intercept_only <- isSimulationInterceptOnly(formula_obj)

  if (!intercept_only && is.call(lhs_expr) &&
      !isLegacyTraitDataSubset(lhs_expr) &&
      !isSupportedNamedResponse(lhs_expr)) {
    stop("Transformed responses are not supported in simulation formulas.")
  }

  if (intercept_only) {
    response_idx <- response_override
    if (is.null(response_idx)) {
      if (isLegacyTraitDataReference(lhs_expr)) {
        response_idx <- seq_along(data_names)
      } else if (isLegacyTraitDataSubset(lhs_expr)) {
        response_idx <- resolveLegacyTraitDataColumns(lhs_expr, data_names)
      } else if (isSupportedNamedResponse(lhs_expr)) {
        response_idx <- match(extractNamedResponseNames(lhs_expr, data_names), data_names)
      } else {
        stop("Unsupported response form in intercept-only simulation formula.")
      }
    }
    if (length(response_idx) == 0L) {
      stop("response_columns must identify at least one response column.")
    }
    if (!is.null(predictor_override) && length(predictor_override) > 0L) {
      stop("predictor_columns should not be supplied for intercept-only formulas.")
    }

    response_names <- data_names[response_idx]
    return(list(
      formula_original = formula_obj,
      formula_original_chr = formula_original_chr,
      formula_normalized_obj = stats::as.formula("trait_data ~ 1"),
      formula_normalized_chr = "trait_data ~ 1",
      formula_mode = "intercept_only",
      response_columns = as.integer(response_idx),
      response_column_names = response_names,
      predictor_columns = integer(0),
      predictor_column_names = character(0),
      predictor_schema = list(),
      intercept_only = TRUE,
      data_prototype = NULL
    ))
  }

  if (isLegacyTraitDataReference(lhs_expr)) {
    if (is.null(response_override)) {
      stop("response_columns must be supplied when the non-intercept formula uses 'trait_data' on the left-hand side.")
    }
    response_idx <- response_override
    predictor_idx <- predictor_override
    rhs_rewritten <- rewriteLegacyTraitDataExpr(rhs_expr, data_names, side = "rhs")
    formula_mode <- "legacy_indexed"
  } else if (isLegacyTraitDataSubset(lhs_expr)) {
    response_idx <- resolveLegacyTraitDataColumns(lhs_expr, data_names)
    if (!is.null(response_override) && !identical(sort(response_idx), sort(response_override))) {
      stop("response_columns conflict with the response block implied by formula.")
    }
    predictor_idx <- predictor_override
    if (is.null(predictor_idx)) {
      predictor_idx <- collectLegacyTraitDataColumns(rhs_expr, data_names)
      predictor_idx <- setdiff(predictor_idx, response_idx)
    }
    rhs_rewritten <- rewriteLegacyTraitDataExpr(rhs_expr, data_names, side = "rhs")
    formula_mode <- "legacy_indexed"
  } else if (isSupportedNamedResponse(lhs_expr)) {
    response_names <- extractNamedResponseNames(lhs_expr, data_names)
    response_idx <- match(response_names, data_names)
    if (!is.null(response_override) && !identical(sort(response_idx), sort(response_override))) {
      stop("response_columns conflict with the response block implied by formula.")
    }

    predictor_vars <- unique(setdiff(all.vars(rhs_expr), response_names))
    if (!all(predictor_vars %in% data_names)) {
      stop("Simulation formulas may only reference raw columns present in trait_data.")
    }
    predictor_idx <- match(predictor_vars, data_names)
    if (!is.null(predictor_override) && !identical(sort(predictor_idx), sort(predictor_override))) {
      stop("predictor_columns conflict with the predictors implied by formula.")
    }
    rhs_rewritten <- rhs_expr
    formula_mode <- "named_formula"
  } else {
    stop("Unsupported response form in simulation formula.")
  }

  if (length(response_idx) == 0L) {
    stop("Simulation formula must identify at least one response column.")
  }
  if (is.null(predictor_idx)) {
    predictor_idx <- setdiff(seq_along(data_names), response_idx)
  }
  predictor_idx <- unique(as.integer(predictor_idx))
  if (length(predictor_idx) == 0L) {
    stop("Simulation formula must identify at least one predictor column for non-intercept workflows.")
  }
  if (length(intersect(response_idx, predictor_idx)) > 0L) {
    stop("response_columns and predictor_columns must not overlap.")
  }

  response_names <- data_names[response_idx]
  predictor_names <- data_names[predictor_idx]

  lhs_rewritten <- if (length(response_names) == 1L) {
    as.name(response_names)
  } else {
    as.call(c(as.name("cbind"), lapply(response_names, as.name)))
  }
  normalized_formula_obj <- stats::as.formula(
    call("~", lhs_rewritten, rhs_rewritten),
    env = parent.frame()
  )
  normalized_formula_chr <- paste(deparse(normalized_formula_obj), collapse = "")

  predictor_schema <- predictorSchemaForData(
    if (is.data.frame(trait_data)) trait_data else as.data.frame(trait_data, check.names = FALSE),
    predictor_names
  )

  list(
    formula_original = formula_obj,
    formula_original_chr = formula_original_chr,
    formula_normalized_obj = normalized_formula_obj,
    formula_normalized_chr = normalized_formula_chr,
    formula_mode = formula_mode,
    response_columns = as.integer(response_idx),
    response_column_names = response_names,
    predictor_columns = as.integer(predictor_idx),
    predictor_column_names = predictor_names,
    predictor_schema = predictor_schema,
    intercept_only = FALSE,
    data_prototype = NULL
  )
}

normalizeMvglsFormulaCall <- function(formula, trait_data, args_list, allow_single_response = FALSE) {
  if (missing(formula)) {
    stop("A character formula or formula object must be provided.")
  }

  if (inherits(formula, "formula")) {
    formula_obj <- formula
  } else if (is.character(formula) && length(formula) == 1L && !is.na(formula)) {
    formula_obj <- as.formula(formula)
  } else {
    stop("A character formula or formula object must be provided.")
  }

  formula_chr <- paste(deparse(formula_obj, width.cutoff = 500), collapse = " ")
  data_obj <- if ("data" %in% names(args_list)) args_list$data else trait_data
  if (is.matrix(data_obj)) {
    data_obj <- as.data.frame(data_obj, stringsAsFactors = FALSE)
  }
  if (!is.null(rownames(trait_data)) && nrow(data_obj) == nrow(trait_data)) {
    rownames(data_obj) <- rownames(trait_data)
  }

  if (grepl("\\btrait_data\\b", formula_chr)) {
    formula_terms <- stats::terms(formula_obj)
    intercept_only <- length(attr(formula_terms, "term.labels")) == 0L &&
      identical(attr(formula_terms, "intercept"), 1L)
    lhs_is_trait_data <- is.symbol(formula_obj[[2L]]) &&
      identical(as.character(formula_obj[[2L]]), "trait_data")
    response_only_df <- is.data.frame(data_obj) &&
      ncol(data_obj) > 0L &&
      all(vapply(data_obj, is.numeric, logical(1)))

    if (intercept_only && lhs_is_trait_data && response_only_df) {
      response_names <- colnames(data_obj)
      if (is.null(response_names) || any(!nzchar(response_names))) {
        response_names <- paste0("Y", seq_len(ncol(data_obj)))
        colnames(data_obj) <- response_names
      }
      response_names <- make.unique(make.names(response_names))
      colnames(data_obj) <- response_names
      formula_obj <- if (length(response_names) == 1L) {
        stats::as.formula(paste0(response_names, " ~ 1"))
      } else {
        stats::as.formula(
          paste0("cbind(", paste(response_names, collapse = ", "), ") ~ 1")
        )
      }
      args_list$data <- data_obj
      formula_chr <- paste(deparse(formula_obj, width.cutoff = 500), collapse = " ")
    } else if (intercept_only && lhs_is_trait_data) {
      return(list(
        formula   = formula_obj,
        args_list = args_list
      ))
    } else {
      data_names <- colnames(data_obj)
      has_missing_names <- is.null(data_names) || any(!nzchar(data_names))

      if (has_missing_names) {
        if (is.null(data_names)) {
          data_names <- rep("", ncol(data_obj))
        }
        missing_idx <- which(!nzchar(data_names))
        data_names[missing_idx] <- paste0("V", seq_along(missing_idx))
        data_names <- make.unique(make.names(data_names))
        colnames(data_obj) <- data_names
      }

      lhs_expr <- formula_obj[[2L]]
      rhs_expr <- formula_obj[[3L]]
      has_legacy_refs <- isLegacyTraitDataReference(lhs_expr) ||
        isLegacyTraitDataSubset(lhs_expr) ||
        length(collectLegacyTraitDataColumns(rhs_expr, colnames(data_obj))) > 0L

      if (!has_legacy_refs) {
        return(list(
          formula   = formula_obj,
          args_list = args_list
        ))
      }

      formula_spec <- normalizeSimulationFormulaSpec(
        formula = formula_obj,
        trait_data = data_obj
      )

      if (isTRUE(formula_spec$intercept_only)) {
        data_obj <- data_obj[, formula_spec$response_column_names, drop = FALSE]
        rownames(data_obj) <- rownames(trait_data)
        args_list$data <- data_obj
        formula_obj <- stats::as.formula("trait_data ~ 1")
        formula_chr <- paste(deparse(formula_obj, width.cutoff = 500), collapse = " ")

        response_only_df <- is.data.frame(data_obj) &&
          ncol(data_obj) > 0L &&
          all(vapply(data_obj, is.numeric, logical(1)))

        if (response_only_df) {
          response_names <- colnames(data_obj)
          if (is.null(response_names) || any(!nzchar(response_names))) {
            response_names <- paste0("Y", seq_len(ncol(data_obj)))
            colnames(data_obj) <- response_names
          }
          response_names <- make.unique(make.names(response_names))
          colnames(data_obj) <- response_names
          formula_obj <- if (length(response_names) == 1L) {
            stats::as.formula(paste0(response_names, " ~ 1"))
          } else {
            stats::as.formula(
              paste0("cbind(", paste(response_names, collapse = ", "), ") ~ 1")
            )
          }
          args_list$data <- data_obj
          formula_chr <- paste(deparse(formula_obj, width.cutoff = 500), collapse = " ")
        }
      } else {
        formula_obj <- formula_spec$formula_normalized_obj
        formula_chr <- paste(deparse(formula_obj, width.cutoff = 500), collapse = " ")
        args_list$data <- data_obj
      }
    }
  }

  mf <- stats::model.frame(formula_obj, data = data_obj, na.action = stats::na.pass)
  response_block <- stats::model.response(mf)
  response_dim <- if (is.null(dim(response_block))) 1L else ncol(response_block)

  if (response_dim < 2L && !allow_single_response) {
    stop("Named formulas require a multivariate response for mvgls.")
  }

  if (response_dim >= 2L) {
    response_names <- colnames(response_block)
    if (is.null(response_names) || any(!nzchar(response_names))) {
      response_names <- paste0("Y", seq_len(response_dim))
    }
    response_names <- make.unique(make.names(response_names))
    colnames(response_block) <- response_names
  } else {
    lhs_expr <- formula_obj[[2L]]
    response_names <- if (is.symbol(lhs_expr)) {
      as.character(lhs_expr)
    } else {
      NULL
    }
    if (is.null(response_names) || any(!nzchar(response_names))) {
      response_names <- "Y1"
    }
    response_names <- make.unique(make.names(response_names))
    response_block <- as.matrix(response_block)
    colnames(response_block) <- response_names
  }

  rhs_terms <- stats::delete.response(attr(mf, "terms"))
  xmat <- stats::model.matrix(rhs_terms, mf)
  if ("(Intercept)" %in% colnames(xmat)) {
    xmat <- xmat[, colnames(xmat) != "(Intercept)", drop = FALSE]
  }

  if (ncol(xmat) > 0) {
    predictor_names <- make.unique(make.names(colnames(xmat)))
    colnames(xmat) <- predictor_names
    working_data <- data.frame(
      response_block,
      xmat,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    formula_fit <- as.formula(
      paste0(
        "cbind(",
        paste(response_names, collapse = ", "),
        ") ~ ",
        paste(predictor_names, collapse = " + ")
      )
    )
  } else {
    working_data <- data.frame(
      response_block,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    formula_fit <- as.formula(
      paste0("cbind(", paste(response_names, collapse = ", "), ") ~ 1")
    )
  }

  rownames(working_data) <- rownames(data_obj)

  list(
    formula   = formula_fit,
    args_list = {
      args_list$data <- working_data
      args_list
    }
  )
}
