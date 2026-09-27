test_that("covariance validity does not depend on an overall change of units", {
  invalid <- matrix(c(1, 2, 2, 1), 2)
  asymmetric <- matrix(c(1, 0.2, 0.8, 1), 2)
  for (scale in c(1e-200, 1e-10, 1, 1e200)) {
    expect_error(
      summarize_regime_covariances(list(invalid = scale * invalid)),
      "Regime `invalid`.*positive semidefinite"
    )
    expect_error(
      summarize_regime_covariances(list(asymmetric = scale * asymmetric)),
      "Regime `asymmetric`.*symmetric"
    )
    expect_error(
      regime_correlation_pca(list(invalid = scale * invalid, valid = scale * diag(2))),
      "Regime `invalid`.*positive semidefinite"
    )
  }
})

test_that("valid covariance summaries retain scale and allow numerical roundoff", {
  valid <- matrix(c(1, 0.5, 0.5, 4), 2)
  near_symmetric <- matrix(c(1, 0.2 + 1e-10, 0.2, 1), 2)
  singular <- matrix(1, 2, 2)
  near_psd <- matrix(c(1, 1 + 1e-10, 1 + 1e-10, 1), 2)
  for (scale in c(1e-200, 1e-10, 1, 1e200)) {
    out <- summarize_regime_covariances(list(valid = scale * valid))
    expect_equal(out$status, "ok")
    expect_equal(out$mean_variance / scale, 2.5)
    expect_equal(out$mean_abs_correlation, 0.25)
    expect_equal(out$fisher_z_mean_abs_correlation, atanh(0.25))
    expect_no_error(summarize_regime_covariances(list(near = scale * near_symmetric)))
    expect_no_error(summarize_regime_covariances(list(singular = scale * singular)))
    expect_no_error(summarize_regime_covariances(list(near_psd = scale * near_psd)))
    expect_error(summarize_regime_covariances(list(zero = matrix(0, 2, 2))),
                 "positive diagonal variances")
  }
})
