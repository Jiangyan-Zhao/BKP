
# BKP: An R Package for Beta Kernel Process Modeling <img src="man/figures/logo.png" align="right" height="140"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/BKP)](https://cran.r-project.org/package=BKP)
![Total downloads](https://cranlogs.r-pkg.org/badges/grand-total/BKP)
[![R-CMD-check](https://github.com/Jiangyan-Zhao/BKP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Jiangyan-Zhao/BKP/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/Jiangyan-Zhao/BKP/graph/badge.svg)](https://app.codecov.io/gh/Jiangyan-Zhao/BKP)
<!-- badges: end -->

We present **BKP**, a user-friendly and extensible **R** package that
implements the **Beta Kernel Process (BKP)**—a fully nonparametric and
computationally efficient framework for modeling spatially varying
binomial probabilities. The BKP model combines localized kernel-weighted
likelihoods with conjugate beta priors, resulting in closed-form
posterior inference without requiring latent variable augmentation or
intensive MCMC sampling. The package supports binary and aggregated
binomial responses, allows flexible choices of kernel functions and
prior specification, and provides loss-based kernel hyperparameter
tuning procedures. In addition, BKP extends naturally to the **Dirichlet
Kernel Process (DKP)** for modeling spatially varying multinomial data.

## Features

- ✅ Bayesian modeling for binomial and multinomial count data
- ✅ Kernel-based local information sharing
- ✅ Posterior prediction and uncertainty quantification
- ✅ Class label prediction using threshold or MAP rule
- ✅ Simulation from posterior (Beta or Dirichlet) distributions

## Installation

You can install the stable version of **BKP** from
[CRAN](https://CRAN.R-project.org/package=BKP) with:

``` r
install.packages("BKP")
```

Or install the development version from
[GitHub](https://github.com/Jiangyan-Zhao/BKP) with:

``` r
# install.packages("pak")
pak::pak("Jiangyan-Zhao/BKP")
```

## Documentation

The statistical foundations and example applications are described in
the following vignette:

- [**BKP User Guide
  (PDF)**](https://github.com/Jiangyan-Zhao/BKP/blob/master/doc/vignettes.pdf)

## Citing

If you use **BKP** in your work, please cite both the methodology paper
and the R package:

- **Methodology paper**  
  Zhao, J., Qing, K., and Xu, J. (2025). *BKP: An R Package for Beta
  Kernel Process Modeling.*  
  arXiv:2508.10447. <https://arxiv.org/abs/2508.10447>

- **R package**  
  Zhao, J., Qing, K., and Xu, J. (2025). *BKP: Beta Kernel Process
  Modeling.*  
  R package version 0.2.1. <https://cran.r-project.org/package=BKP>

You can also obtain the citation information directly within R:

``` r
citation("BKP")
```

## Development

The BKP package is under active development. Contributions and
suggestions are welcome via GitHub issues or pull requests.
