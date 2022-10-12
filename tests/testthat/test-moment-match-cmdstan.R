
library(cmdstanr)

# Very simple example model
# Gaussian 1-dimensional data with unknown location and scale

stancode <- "data {
  int<lower=0> N;
  vector[N] x;
  }
  parameters {
    real mu;
    real<lower=0> sigma;
  }
  model {
    x ~ normal(mu, sigma);
  }
  generated quantities {
    vector[N] log_lik;
    for (n in 1:N)
      log_lik[n] = normal_lpdf(x[n] | mu, sigma);
  }"

stanmodel <- cmdstan_model(stan_file = write_stan_file(stancode))

# Let us generate data from the true data generating mechanism
SEED <- 24
set.seed(SEED)
data_var = 1
n = as.integer(30)
x = rnorm(n = n, mean = 0,sd = data_var)
standata = list(N = n, x = x)
fit <- stanmodel$sample(data = standata, chains = 4, iter_sampling = 1000, iter_warmup = 1000,
                        refresh = 0, seed = SEED)

fit$init_model_methods()

xloo <- x
xloo[1] <- 9
standata_loo = list(N = n, x = xloo)
fit_loo <- stanmodel$sample(data = standata_loo, chains = 4, iter_sampling = 1000, iter_warmup = 1000, refresh = 0, seed = SEED)

fit_loo$init_model_methods()


test_that("moment_match.CmdStanFit works", {

  # define target
  target_density <- function(draws, ...) {
    dnorm(draws[,1], 2, 0.136, log = TRUE) + dnorm(draws[,2], 0.754, 0.104,
                                                   log = TRUE)
  }
  ratio_density <- function(draws, fit, ...) {
    target_density(draws, ...) - apply(draws, 1, fit$log_prob, ...)
  }

  iw1 <- moment_match.CmdStanFit(fit,
                                 log_prob_target_fun = target_density)

  
  iw2 <- moment_match.CmdStanFit(fit,
                             log_ratio_fun = ratio_density)

  expect_equal(matrixStats::colWeightedMeans(iw1$draws,w = exp(iw1$log_weights)), c(1.99, 0.75), tolerance = 1e-2)
  expect_equal(matrixStats::colWeightedMeans(iw1$draws^2,w = exp(iw1$log_weights)) - matrixStats::colWeightedMeans(iw1$draws,w = exp(iw1$log_weights))^2, c(0.017, 0.011), tolerance = 1e-2)
  expect_equal(iw1$pareto_k, 0.3, tolerance = 1e-2)

  expect_equal(iw1, iw2)

  iw_ex1 <- moment_match.CmdStanFit(fit,
                                log_prob_target_fun = target_density,
                                expectation_fun =  function(draws, ...) {
                                  draws})

  iw_ex2 <- moment_match.CmdStanFit(fit,
                                log_ratio_fun = ratio_density,
                                expectation_fun =  function(draws, ...) {
                                  draws})

  expect_equal(matrixStats::colWeightedMeans(iw1$draws,w = exp(iw1$log_weights)), iw_ex1$expectation, tolerance = 1e-6)
  expect_equal(iw1$pareto_k, iw_ex1$pareto_k, tolerance = 1e-6)
  expect_equal(iw_ex1$pareto_kf, c(0.33, 0.27), tolerance = 1e-2)

  expect_equal(iw_ex1, iw_ex2)

})
