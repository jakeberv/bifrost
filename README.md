# bifrost <img src="man/figures/logo.png" align="right" height="138" alt="" />
<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/jakeberv/bifrost/graph/badge.svg)](https://app.codecov.io/gh/jakeberv/bifrost)
[![R-CMD-check](https://github.com/jakeberv/bifrost/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jakeberv/bifrost/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![GitHub version](https://img.shields.io/github/r-package/v/jakeberv/bifrost?label=dev%20version)](https://github.com/jakeberv/bifrost)
[![pkgdown](https://github.com/jakeberv/bifrost/actions/workflows/pkgdown.yaml/badge.svg)](https://jakeberv.github.io/bifrost/)
<!-- badges: end -->

# Branch-level Inference Framework for Recognizing Optimal Shifts in Traits

# bifrost

Branch-level inference of multiregime trait evolution on phylogenies

## Overview

The primary function in `bifrost` is `searchOptimalConfiguration()`. It performs a greedy stepwise search for evolutionary regime shifts on a SIMMAP-style phylogeny, using multivariate fits from `mvMORPH::mvgls`. Models are fit in trait space without PCA, assuming a multi-rate Brownian Motion with proportional VCV scaling across regimes. The function can search, score, and accept shifts by an information criterion, then optionally compute per-shift IC weights.

This repository is under active development.

## Installation

Install the development version from GitHub.

```r
# install.packages("remotes")
remotes::install_github("jakeberv/bifrost")
```

On Windows, install Rtools for your R version and ensure it is on the PATH.
