simulationStudyFormulaCases <- list(
  string = "trait_data ~ 1",
  formula = stats::as.formula("trait_data ~ 1")
)

simulationStudyFormulaErrors <- list(
  non_intercept = paste0(
    "^Simulation studies currently support intercept-only search formulas only\\. ",
    "Use formula = \\\"trait_data ~ 1\\\" ",
    "\\(or an equivalent intercept-only response formula\\)\\.$"
  ),
  invalid_type = "^formula must be a single character string or a formula object\\.$"
)
