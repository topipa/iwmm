
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src='man/figures/logo.png' align="right" height="160" />

# iwmm

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/iwmm)](https://CRAN.R-project.org/package=iwmm)
[![R-CMD-check](https://github.com/topipa/iwmm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/topipa/iwmm/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/topipa/iwmm/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/topipa/iwmm/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

## Overview

**IWMM** is an R package for adaptive importance sampling. The package
implements the importance weighted moment matching (IWMM) algorithm.
Given draws from a probability distribution, IWMM adapts the draws based
on a given target distribution and an optional expectation function. The
draws can, but do not have to be from the distribution over which the
expectation is defined.

The method is described in detail in the [Implicitly Adaptive Importance
Sampling](https://doi.org/10.1007/s11222-020-09982-2) paper.

## Installation

<!-- You can install the released version of iwmm from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("iwmm") -->
<!-- ``` -->
<!-- And the development version from [GitHub](https://github.com/) with: -->
<!-- ``` r -->
<!-- # install.packages("devtools") -->
<!-- devtools::install_github("topipa/iwmm") -->
<!-- ``` -->

You can install the the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("topipa/iwmm")
```

<!-- ## Example -->
<!-- This is a basic example which shows you how to solve a common problem: -->
<!-- ```{r example} -->
<!-- library(iwmm) -->
<!-- ## basic example code -->
<!-- ``` -->
<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub! -->

## Usage

IWMM takes a sample of draws, and functions that define the proposal and
target distributions, as well as an optional expectation function.
Without an expectation function, the draws are adapted to match the
proposal distribution. If an expectation function is given, the
expectation is computed after adapting the draws specifically for this
expectation.

IWMM can also work with models fitted with
[rstan](https://github.com/stan-dev/rstan) or
[cmdstanr](https://github.com/stan-dev/cmdstanr) where the fitted model
posterior is treated as the proposal distribution.

### Example

## References

Paananen, T., Piironen, J., Bürkner, P.-C., and Vehtari, A. Implicitly
Adaptive Importance Sampling. *Stat Comput* **31**, 16 (2021).
([paper](https://doi.org/10.1007/s11222-020-09982-2))([code](https://github.com/topipa/iter-mm-paper))
