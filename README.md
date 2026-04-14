# bifrost <img src="man/figures/logo.png" align="right" height="140" alt="bifrost hex sticker" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/jakeberv/bifrost/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jakeberv/bifrost/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/jakeberv/bifrost/graph/badge.svg)](https://app.codecov.io/gh/jakeberv/bifrost)
[![CRAN status](https://www.r-pkg.org/badges/version/bifrost)](https://CRAN.R-project.org/package=bifrost)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/bifrost)](https://cran.r-project.org/package=bifrost)
[![bioRxiv preprint](https://img.shields.io/endpoint?url=https%3A%2F%2Fjakeberv.github.io%2Fbiorxiv-badge%2Fbadges%2F10.64898__2026.04.12.718036.json)](https://doi.org/10.64898/2026.04.12.718036)
[![License: GPL (>= 2)](https://img.shields.io/badge/license-GPL%20(%E2%89%A5%202)-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- badges: end -->

<span style="font-size:1.1em;">
  <strong>
    <span style="color:#c40000;">B</span>ranch-level
    <span style="color:#c40000;">I</span>nference
    <span style="color:#c40000;">F</span>ramework for
    <span style="color:#c40000;">R</span>ecognizing
    <span style="color:#c40000;">O</span>ptimal
    <span style="color:#c40000;">S</span>hifts in
    <span style="color:#c40000;">T</span>raits
  </strong>
</span>

`bifrost` performs branch-level inference of multi-regime, multivariate trait evolution on a phylogeny using [penalized-likelihood multivariate GLS fits](https://doi.org/10.1093/sysbio/syy045). The current version searches for evolutionary model shifts under a multi-rate Brownian Motion (BMM) model with proportional regime VCV scaling, operates directly in trait space (for example, without PCA), and is designed for high-dimensional datasets (`p > n`) and large trees (`> 1000` tips). The method works with fossil tip-dated trees and with a wide range of multivariate comparative data, including GPA-aligned morphometric coordinates, linear dimensions, and related trait matrices. A future major release will add support for the [multivariate scalar Ornstein-Uhlenbeck process](https://doi.org/10.1093/sysbio/syy005).

## Installation

### Stable release

```r
install.packages("bifrost")
```

### Development version

```r
# install.packages("remotes")
remotes::install_github("jakeberv/bifrost")
```

**Windows users:**  
Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for your R version and ensure that it is added to your system `PATH`.

**macOS users:**  
You may need to install [XQuartz](https://www.xquartz.org/) to build or run packages that depend on certain graphical or system libraries.

## Overview

- **Primary goal.** Infer *where*, *when*, and *how* patterns of phenotypic evolution change across a tree using many traits simultaneously.
- **Model.** Multi-rate Brownian Motion with regime-specific VCVs estimated via penalized-likelihood (`mvMORPH::mvgls`), supporting `p >= n`.
- **Search.** Greedy, step-wise acceptance of shifts guided by information criteria (`GIC` or `BIC`), with optional per-shift IC weights.
- **Scale.** Parallel candidate scoring using the `future` ecosystem; practical on thousands of taxa x traits.

## Key features

- Joint multivariate modeling without information loss or [distortion due to PCA](https://doi.org/10.1093/sysbio/syv019).
- Under BMM, [proportional VCV scaling](https://doi.org/10.1111/j.1558-5646.1999.tb05414.x) across regimes for tractability at high `p`.
- Provides a multivariate phylogenetic GLS (mvPGLS)-like framework in which hidden branch-specific rate regimes are inferred and incorporated when estimating predictor effects.
- Candidate shift nodes are determined by a minimum clade size specified by the user.
- Greedy [step-wise heuristic search](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.19099) using GIC/BIC delta-IC thresholds, with optional IC-weight support for inferred shifts.
- Output includes estimated VCV per regime, shift weights, and SIMMAP-style mappings for downstream visualization and analysis.
- Parallelization via `future` and `future.apply`.

## Documentation

### Background and theory

- [Brownian Motion and Multivariate Shifts](https://jakeberv.com/bifrost/articles/theoretical-background-vignette.html)  
  Conceptual background on Brownian motion, multivariate covariance, and what branch-specific shifts mean in `bifrost`.
- [Whole-Tree Models, PCA, and bifrost](https://jakeberv.com/bifrost/articles/pca-model-selection-and-bifrost-vignette.html)  
  Explains why whole-tree homogeneous models and PCA truncation can mislead inference when evolutionary processes vary across the tree.

### Using bifrost

- [Quick Start with bifrost](https://jakeberv.com/bifrost/articles/quick-start-vignette.html)  
  A practical introduction to the core `bifrost` workflow using a minimal simulated example, including setup, key arguments, outputs, and interpretation.
- [Detecting Evolutionary Shifts in Paleozoic Fish Jaw Shape with bifrost](https://jakeberv.com/bifrost/articles/jaw-shape-vignette.html)  
  A full empirical case study using the packaged fossil jaw-shape dataset, showing how to run, inspect, and interpret a real `bifrost` analysis end to end.

## Additional note

Though `bifrost` was initially developed as a framework for inferring macroevolutionary regime shifts in multivariate trait data, it can also be applied to perform multivariate phylogenetic generalized least squares (pGLS) analyses with factors or continuous predictors (e.g., `cbind(trait1, trait2, ...) ~ predictor`, or `"trait_data[, 1:5] ~ trait_data[, 6]"` when working directly with a matrix). In this context, `bifrost` identifies branch-specific rate variation under a multi-rate Brownian Motion model and fits the pGLS conditional on the resulting residual (phylogenetic) covariance structure, so estimated effect sizes and uncertainties account for "hidden" rate variation not explained by the predictors. This is conceptually similar to hidden-state approaches (e.g., [Boyko et al. 2023](https://doi.org/10.1093/evolut/qpad002)), except that here the regimes influence variance and evolutionary rate rather than introducing regime-specific means. This use case is an active area of ongoing methodological development.

## Citation

If you use `bifrost`, please cite the package and methods references below. The same set is also available from:

```r
citation("bifrost")
```

### Recommended citations

1. `bifrost` methods / application paper  
   Berv JS, Probst CM, Claramunt S, Shipley JR, Friedman M, Smith SA, Fouhey DF, Weeks BC (2026). *Rates of passerine body plan evolution in time and space*. *Nature Ecology & Evolution*. In press.

2. `bifrost` preprint  
   Berv JS, Fox N, Thorstensen MJ, Lloyd-Laney H, Troyer EM, Rivero-Vega RA, Smith SA, Friedman M, Fouhey DF, Weeks BC (2026). *bifrost: an R package for scalable inference of phylogenetic shifts in multivariate evolutionary dynamics*. *bioRxiv*. [https://doi.org/10.64898/2026.04.12.718036](https://doi.org/10.64898/2026.04.12.718036)

3. `bifrost` software citation  
   Berv JS, Fox N, Thorstensen MJ, Lloyd-Laney H, Troyer EM, Rivero-Vega RA, Smith SA, Friedman M, Fouhey DF, Weeks BC (2026). *Branch-Level Inference Framework for Recognizing Optimal Shifts in Traits*. R package version 0.1.3. [https://CRAN.R-project.org/package=bifrost](https://CRAN.R-project.org/package=bifrost)

4. `mvMORPH` package paper  
   Clavel J, Escarguel G, Merceron G (2015). *mvmorph: an R package for fitting multivariate evolutionary models to morphometric data*. *Methods in Ecology and Evolution*, 6(11), 1311-1319. [https://doi.org/10.1111/2041-210X.12420](https://doi.org/10.1111/2041-210X.12420)

5. Penalized-likelihood framework paper  
   Clavel J, Aristide L, Morlon H (2019). *A Penalized Likelihood Framework for High-Dimensional Phylogenetic Comparative Methods and an Application to New-World Monkeys Brain Evolution*. *Systematic Biology*, 68(1), 93-116. [https://doi.org/10.1093/sysbio/syy045](https://doi.org/10.1093/sysbio/syy045)

## Contributing

Bug reports, feature requests, and pull requests are welcome. Please open an issue at [https://github.com/jakeberv/bifrost/issues](https://github.com/jakeberv/bifrost/issues).

## License

This project is released under the GPL (>= 2) License. See the `LICENSE` file for details.

## Acknowledgements and dependencies

`bifrost` builds on substantial work from `mvMORPH`, `phytools`, `ape`, `future`, and `future.apply`. The greedy search algorithm is adapted from [Mitov et al. 2019](https://www.pnas.org/doi/10.1073/pnas.1813823116) and [Smith et al. 2023](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.19099). See the `DESCRIPTION` file for complete dependency and version information.

The name of our R package is inspired by the Bifrost, the rainbow bridge of Norse mythology that connects Earth (Midgard) and Asgard within the cosmic structure of Yggdrasil, the Tree of Life, echoing how this framework links observable data to hidden evolutionary shifts across the history of life.

Development of the `bifrost` R package was supported by the [Oxford Research Software Engineering Group](https://www.rse.ox.ac.uk/schmidt-ai-science), with support from [Schmidt Sciences, LLC.](https://www.schmidtsciences.org/ai-in-science/) and the [Michigan Institute for Data Science and AI in Society](https://midas.umich.edu/).

<p align="center" style="display:flex; justify-content:center; align-items:center; gap:50px; padding:30px 0;">
  <img src="https://jakeberv.com/images/SchmidtSciencesLogo.png"
       alt="Schmidt Sciences logo"
       style="height:90px !important; width:auto !important; max-width:100%;" />
  <img src="https://www.rse.ox.ac.uk/sites/default/files/rse/site-logo/2024_oxrse_next_to_oxford.svg"
       alt="Oxford RSE logo"
       style="height:90px !important; width:auto !important; max-width:100%;" />
</p>
