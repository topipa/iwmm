
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

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
    x ~ normal(mu,sigma);
  }
  generated quantities {
    vector[N] log_lik;
    for (n in 1:N)
      log_lik[n] = normal_lpdf(x[n] | mu, sigma);
  }"

stanmodel <- stan_model(model_code = stancode)

# Let us generate data from the true data generating mechanism
SEED <- 24
set.seed(SEED)
data_var = 1
n = as.integer(30)
x = rnorm(n = n, mean = 0,sd = data_var)
standata = list(N = n, x = x)
fit <- sampling(stanmodel, data = standata, chains = 4, iter = 2000,
                refresh = 0, seed = SEED)

xloo <- x
xloo[1] <- 9
standata_loo = list(N = n, x = xloo)
fit_loo <- sampling(stanmodel, data = standata_loo, chains = 4, iter = 2000, refresh = 0, seed = SEED)




test_that("moment_match.stanfit works", {

  # define target
  target_density <- function(draws, ...) {
    dnorm(draws[,1], 2, 0.136, log = TRUE) + dnorm(draws[,2], 0.754, 0.104,
                                                   log = TRUE)
  }
  ratio_density <- function(draws, ...) {
    target_density(draws, ...) - log_prob_upars.stanfit(draws, ...)
  }

  iw1 <- moment_match.stanfit(fit,
                             log_prob_target_fun = target_density)
  iw2 <- moment_match.stanfit(fit,
                             log_ratio_fun = ratio_density)

  expect_equal(matrixStats::colWeightedMeans(iw1$draws,w = exp(iw1$log_weights)), c(1.9984141, 0.7520985), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw1$draws^2,w = exp(iw1$log_weights)) - matrixStats::colWeightedMeans(iw1$draws,w = exp(iw1$log_weights))^2, c(0.01772950, 0.01102333), tolerance = 1e-6)
  expect_equal(iw1$pareto_k, 0.377199, tolerance = 1e-6)

  expect_equal(iw1,iw2)

  iw_ex1 <- moment_match.stanfit(fit,
                                log_prob_target_fun = target_density,
                                expectation_fun =  function(draws, ...) {
                                  draws})

  iw_ex2 <- moment_match.stanfit(fit,
                                log_ratio_fun = ratio_density,
                                expectation_fun =  function(draws, ...) {
                                  draws})

  expect_equal(matrixStats::colWeightedMeans(iw1$draws,w = exp(iw1$log_weights)), iw_ex1$expectation, tolerance = 1e-6)
  expect_equal(iw1$pareto_k, iw_ex1$pareto_k, tolerance = 1e-6)
  expect_equal(iw_ex1$pareto_kf, c(0.2773750, 0.4282701), tolerance = 1e-6)

  expect_equal(iw_ex1, iw_ex2)

})




test_that("moment_match.stanfit works with obs_weights formulation", {


  log_lik_stanfit <- function(stanfit, upars, parameter_name = "log_lik",
                              ...) {
    ll <- loo::extract_log_lik(stanfit, parameter_name, merge_chains = TRUE)
    S <- nrow(upars)
    n <- ncol(ll)
    out <- matrix(0,S,n)
    for (s in seq_len(S)) {
      out[s,] <- rstan::constrain_pars(stanfit, upars = upars[s, ])[[parameter_name]]
    }
    out
  }

  ratio_density <- function(draws, stanfit, obs_weights, ...) {
    log_lik <- log_lik_stanfit(stanfit, draws)
    colSums((obs_weights - 1) * t(log_lik))
  }
  target_density <- function(draws, ...) {
    ratio_density(draws, ...) + log_prob_upars.stanfit(draws, ...)
  }


  iw1 <- moment_match.stanfit(fit_loo,
                             log_ratio_fun = ratio_density,
                             k_threshold = 0.2,
                             obs_weights = c(0, rep(1,n - 1)))

  iw2 <- moment_match.stanfit(fit_loo,
                              log_prob_target_fun = target_density,
                              k_threshold = 0.2,
                              obs_weights = c(0, rep(1,n - 1)))

  expect_equal(matrixStats::colWeightedMeans(iw1$draws,w = exp(iw1$log_weights)), c(-0.1544952, -0.2800384), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw1$draws^2,w = exp(iw1$log_weights)) - matrixStats::colWeightedMeans(iw1$draws,w = exp(iw1$log_weights))^2, c(0.01990700, 0.01874187), tolerance = 1e-6)
  expect_equal(iw1$pareto_k, 0.1153987, tolerance = 1e-6)

  expect_equal(iw1,iw2)


  iw_ex1 <- moment_match.stanfit(fit_loo,
                                log_ratio_fun = ratio_density,
                                expectation_fun =  function(draws, ...) {
                                  draws},
                                k_threshold = 0.2,
                                obs_weights = c(0, rep(1,n - 1)))

  iw_ex2 <- moment_match.stanfit(fit_loo,
                                log_prob_target_fun = target_density,
                                expectation_fun =  function(draws, ...) {
                                  draws},
                                k_threshold = 0.2,
                                obs_weights = c(0, rep(1,n - 1)))

  expect_equal(matrixStats::colWeightedMeans(iw1$draws,w = exp(iw1$log_weights)), iw_ex1$expectation, tolerance = 1e-6)
  expect_equal(iw1$pareto_k, iw_ex1$pareto_k, tolerance = 1e-6)
  expect_equal(iw_ex1$pareto_kf, c(0.1758149, 0.1814282), tolerance = 1e-6)


  expect_equal(iw_ex1,iw_ex2)
})



