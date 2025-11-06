# bifrost <img src="man/figures/logo.png" align="right" height="140" alt="bifrost hex sticker" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/jakeberv/bifrost/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jakeberv/bifrost/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/jakeberv/bifrost/graph/badge.svg)](https://app.codecov.io/gh/jakeberv/bifrost)
[![CRAN status](https://www.r-pkg.org/badges/version/bifrost)](https://CRAN.R-project.org/package=bifrost)
[![License: GPL (≥ 2)](https://img.shields.io/badge/license-GPL%20(%E2%89%A5%202)-blue.svg)](LICENSE)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- badges: end -->

**<span style="color:#c40000; font-size:1.2em;">B</span>ranch-level 
<span style="color:#c40000; font-size:1.2em;">I</span>nference 
<span style="color:#c40000; font-size:1.2em;">F</span>ramework for 
<span style="color:#c40000; font-size:1.2em;">R</span>ecognizing 
<span style="color:#c40000; font-size:1.2em;">O</span>ptimal 
<span style="color:#c40000; font-size:1.2em;">S</span>hifts in 
<span style="color:#c40000; font-size:1.2em;">T</span>raits**

`bifrost` performs branch-level inference of multi-regime, multivariate trait evolution on a phylogeny using [penalized-likelihood multivariate GLS fits](https://academic.oup.com/sysbio/article/67/4/662/4827615). The current version searches for evolutionary model shifts under a multi-rate Brownian Motion (BMM) model with proportional regime VCV scaling, operating directly in trait space (e.g., no PCA), and is designed for high-dimensional datasets (p > n) and large trees (> 1000 tips). The method will work with fossil tip-dated trees, and will accept most forms of multivariate comparative data (e.g., GPA aligned morphometric coordinates, linear dimensions, and others). The next major release will enable usage of the [multivariate scalar Ornstein–Uhlenbeck process](https://academic.oup.com/sysbio/article/67/4/662/4827615).

---

## Overview

- **Primary goal.** Infer *where*, *when*, and *how* patterns of phenotypic evolution change across a tree using many traits simultaneously.
- **Model.** Multi-rate Brownian Motion with regime-specific VCVs estimated via penalized-likelihood (`mvMORPH::mvgls`), supporting p ≳ n.
- **Search.** Greedy, step-wise acceptance of shifts guided by information criteria (**GIC** or **BIC**), with optional post-hoc pruning and per-shift IC weights.
- **Scale.** Parallel candidate scoring using the `future` ecosystem; practical on thousands of taxa × traits.

---

## Key features

- Joint multivariate modeling without information loss or [distortion due to PCA](https://academic.oup.com/sysbio/article/64/4/677/1649888).
- Under BMM, [proportional VCV scaling](https://doi.org/10.1111/j.1558-5646.1999.tb05414.x) across regimes for tractability at high p.
- Candidate shift nodes are determined by a minimum clade size specified by the user.
- Greedy [step-wise heuristic search](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.19099) using GIC/BIC ΔIC threshold set by the user; uncertainty estimation with IC weights.
- Output includes estimated VCV per regime, shift weights, SIMMAP-style output for cross-compatibility.
- Parallelization steps via `future` / `future.apply`.


---

📄 **Vignette 1:** [Getting Started with bifrost](https://jakeberv.com/bifrost/articles/jaw-shape-vignette.html)

## Installation (development version)

```r
# install.packages("remotes")
remotes::install_github("jakeberv/bifrost")
```

**Windows users:**  
Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for your R version and ensure it is added to your system PATH.

**macOS users:**  
You may need to install [XQuartz](https://www.xquartz.org/) to build or run packages that depend on certain graphical or system libraries.

---

### Quick start

```r
library(bifrost)
library(ape)

set.seed(1)
tree   <- rtree(50)
traits <- matrix(rnorm(50 * 5), ncol = 5)
rownames(traits) <- tree$tip.label   # critical: rownames must match tip labels

res <- searchOptimalConfiguration(
  baseline_tree = tree,
  trait_data = traits,
  IC = "GIC",
  min_descendant_tips = 5,
  num_cores = 2,
  shift_acceptance_threshold = 10,
  plot = FALSE
)

# This example doesn’t detect any shift nodes, so the next two lines will have no visible output.

res$shift_nodes           # no detected shifts in this toy dataset
plotSimmap(res$tree_no_uncertainty)  # no output produced here
```

### Data requirements

- **Tree and data alignment.** `rownames(trait_data)` must match `tree$tip.label` (same order and names).  
- **Branch lengths.** Interpreted in units of time; ultrametric not required.  
- **SIMMAP style.** Internally, regimes are stored using SIMMAP conventions.  
- **Multi-dimensional traits.** Works directly in trait space; tune penalties/methods in `mvgls` options for your data.  
- **Thresholds.** Use conservative `shift_acceptance_threshold` and `ic_uncertainty_threshold` to limit false positives; explore sensitivity.

### Core workflow

```mermaid
flowchart TD
A([Start])
A --> B[Check IC]
B --> C[Validate inputs]
C --> D[Paint baseline]
D --> E[Generate candidates]
E --> F[Fit baseline]
F --> G[Set parallel plan]
G --> H[Fit candidates in parallel]
H --> I[Restore sequential]
I --> J[Compute delta IC and sort]
J --> K[Init search state]
K --> L{More candidates?}
L -->|Yes| M[Next candidate]
M --> N[Add shift]
N --> O[Fit shifted model]
O --> P{Improvement >= threshold?}
P -->|Yes| Q[Accept and update]
P -->|No| R[Reject]
O --> S{Store history?}
S -->|Yes| T[Save iteration]
S -->|No| U[Continue]
O --> V{Warnings?}
V -->|Yes| W[Collect warnings]
V -->|No| X[Continue]
Q --> L
R --> L
L -->|No| Y[Finalize best model]
Y --> Z{Compute weights?}
Z -->|Yes| ZA{Parallel weights?}
ZA -->|Yes| ZB[Parallel drop-one refits]
ZA -->|No| ZC[Serial drop-one refits]
Z -->|No| ZD[Skip weights]
ZB --> ZE[Build ic weights]
ZC --> ZE
ZD --> ZE
ZE --> ZF[Build output trees]
ZF --> ZG[Assemble result list]
ZG --> ZH[Extract VCVs]
ZH --> ZI([Return])
```

-----

### Primary functions

  - `searchOptimalConfiguration()`: The main function for end-to-end greedy search: candidate generation → parallel fitting → iterative acceptance → optional pruning/IC weights.
  - `plot_ic_acceptance_matrix()`: Visualize shift acceptance and information-criterion (IC) differences across search iterations.
  
### Helper functions (not exported)
  
  - **Candidate generation**: `generatePaintedTrees()`
  - **Model fitting helpers**: `fitMvglsAndExtractGIC()`, `fitMvglsAndExtractBIC()`, and formula variants.
  - **IC utilities**: `calculateAllDeltaGIC()`
  - **Tree painting utilities**: `paintSubTree_mod()`, `addShiftToModel()`, `removeShiftFromTree()`, `paintSubTree_removeShift()`, `whichShifts()`
  - **Regime VCVs**: `extractRegimeVCVs()`

-----

### Outputs

`searchOptimalConfiguration()` returns a comprehensive list containing:

- **`user_input`** — All arguments passed to `searchOptimalConfiguration()`, including the tree, trait data, IC choice, thresholds, and other parameters for reproducibility.
- **`tree_no_uncertainty_transformed`** — Optimal SIMMAP tree with accepted shifts, using transformed branch lengths (if a branch-length transformation was applied).
- **`tree_no_uncertainty_untransformed`** — Same optimal SIMMAP tree but with original, untransformed branch lengths.
- **`model_no_uncertainty`** — Final fitted `mvgls` model object (BM or multi-rate BMM), including estimated parameters, log-likelihood, and variance–covariance matrices.
- **`shift_nodes_no_uncertainty`** — Node numbers corresponding to accepted evolutionary shifts.
- **`optimal_ic`** — Information criterion (IC) value for the optimal model.
- **`baseline_ic`** — IC value for the null (single-rate) baseline model.
- **`IC_used`** — The information criterion applied (`"GIC"` or `"BIC"`).
- **`num_candidates`** — Total number of candidate models evaluated during the search.
- **`model_fit_history`** — Per-iteration log of model fits, IC values, and acceptance decisions; useful for visualizing or debugging search behavior.
- **`VCVs`** — Regime-specific penalized-likelihood variance–covariance matrices.
- **`ic_weights`** — Data frame of per-shift IC weights and evidence ratios (if `uncertaintyweights_par = TRUE`), providing support values for individual shifts.

-----


### Performance and scalability

Enable parallel processing using the `future` package:

```r
library(future)
plan(multisession)   # or multicore on Linux/macOS
```

  - **Reduce plotting** (`plot = FALSE`) for large trees.
  - **Increase memory** for heavy runs, especially with high `p`.
  - Consider **larger `min_descendant_tips`** and stricter IC thresholds on very large problems.
  - **Repeat searches** with different seeds and thresholds to check for robustness.

-----

### Reproducibility

  - **Record `sessionInfo()`** and the `mvMORPH` version.
  - For projects, consider using `renv` to lock package versions.

-----

### Additional note

Though `bifrost` was initially developed as a framework for inferring macroevolutionary regime shifts in multivariate trait data, it can also be applied to perform multivariate phylogenetic generalized least squares (pGLS) analyses with factors or continuous predictors (e.g., `Y ~ trophic_niche`). In this context, `bifrost` identifies branch-specific rate variation under a multi-rate Brownian Motion model and fits the pGLS conditional on that residual (phylogenetic) covariance, so estimated effect sizes and uncertainties account for “hidden” rate variation not explained by the predictors. This is conceptually similar to hidden-state approaches (e.g., [Boyko et al. 2023](https://doi.org/10.1093/evolut/qpad002)), except that here the regimes influence variance and evolutionary rate rather than introducing regime-specific means. This use case has not yet been explored in detail.

### Citation

If you use `bifrost`, please cite:

> TBD

Also, run the following to obtain a BibTeX entry when available:

```r
citation("bifrost")
```

-----

### Contributing

Bug reports, feature requests, and pull requests are welcome. Please open an issue at [https://github.com/jakeberv/bifrost/issues](https://github.com/jakeberv/bifrost/issues).

-----

### License

This project is released under the GPL >= 2 License. See the `LICENSE` file for details.

-----

### Acknowledgements and dependencies

`bifrost` builds on substantial work from `mvMORPH`, `phytools`, `ape`, `future`, and `future.apply`. The greedy search algorithim is adapted from [Mitov et al 2019](https://www.pnas.org/doi/10.1073/pnas.1813823116) and [Smith et al 2023](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.19099). See the `DESCRIPTION` file for complete dependency and version information.

Development of the `bifrost` R package was supported by the [Oxford Research Software Engineering Group](https://www.rse.ox.ac.uk/schmidt-ai-science), with support from [Schmidt Sciences, LLC.](https://www.schmidtsciences.org/ai-in-science/) and the [Michigan Institute for Data Science and AI in Society](https://midas.umich.edu/).

<p align="center" style="display:flex; justify-content:center; align-items:center; gap:50px; padding:30px 0;">
  <img src="https://jakeberv.com/images/SchmidtSciencesLogo.png"
       alt="Schmidt Sciences logo"
       style="height:90px !important; width:auto !important; max-width:100%;" />
  <img src="https://www.rse.ox.ac.uk/sites/default/files/rse/site-logo/2024_oxrse_next_to_oxford.svg"
       alt="Oxford RSE logo"
       style="height:90px !important; width:auto !important; max-width:100%;" />
</p>


