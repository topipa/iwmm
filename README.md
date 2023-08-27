
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
install.packages("remotes")
remotes::install_github("topipa/iwmm")
```

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

Consider a simple univariate normal model (available via
`example_iwmm_model("normal_model")`):

``` stan
data {
int<lower=0> N;
vector[N] x;
}
parameters {
  real mu;
  real log_sigma;
}
transformed parameters {
  real<lower=0> sigma = exp(log_sigma);
}
model {
  target += normal_lpdf(x | mu, sigma);
}
```

We first fit the model using Stan:

``` r
library("iwmm")

normal_model <- example_iwmm_model("normal_model")

fit <- rstan::stan(
  model_code = normal_model$model_code,
  data = normal_model$data,
  refresh = FALSE,
  seed = 1234
)
```

After fitting the model, let us define an expectation function that we
are interested in, and compute the expectation:

``` r
expectation_fun_first_moment <- function(draws, ...) {
  draws
}

first_moment <- moment_match(
  fit,
  expectation_fun = expectation_fun_first_moment
)
```

We can check that the expectation matches the posterior mean

``` r
first_moment$expectation
#> [1] 4.002850 2.000641
colMeans(as.matrix(fit))[c("mu", "log_sigma")]
#>        mu log_sigma 
#>  4.002850  2.000641
```

The package can also compute expectations over distributions that are
different than the one where the draws are from. Let us try to evaluate
the posterior mean when the last observation is removed. We achieve this
by defining a `target_observation_weights` which we can use to indicate
the target distribution as a vector of weights for each observation in
the data used to fit the model. A vector of ones means using the
existing data. A zero value for some observation means that the target
distribution is the posterior without that specific observation.

``` r
first_moment_loo <- moment_match(
  fit,
  target_observation_weights=append(rep(1,9), 0),
  expectation_fun = expectation_fun_first_moment
)
```

Let us compare this to the posterior mean we get by actually fitting the
model again without the last observation.

``` r
loo_data <- normal_model$data
loo_data$x <- loo_data$x[-loo_data$N]
loo_data$N <- loo_data$N - 1

fit_loo <- rstan::stan(
  model_code = normal_model$model_code,
  data = loo_data,
  refresh = FALSE,
  seed = 1234
)
```

``` r
first_moment_loo$expectation
#> [1]  1.8939496 -0.1475459
colMeans(as.matrix(fit_loo))[c("mu", "log_sigma")]
#>        mu log_sigma 
#>  1.895802 -0.151749
```

## References

Paananen, T., Piironen, J., Bürkner, P.-C., and Vehtari, A. Implicitly
Adaptive Importance Sampling. *Stat Comput* **31**, 16 (2021).
([paper](https://doi.org/10.1007/s11222-020-09982-2))([code](https://github.com/topipa/iter-mm-paper))
