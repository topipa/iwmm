
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



test_that("moment_match.stanfit works", {




  # Let us generate data from the true data generating mechanism

  SEED <- 24
  set.seed(SEED)
  data_var = 1
  n = as.integer(30)
  x = rnorm(n = n, mean = 0,sd = data_var)

  # fit the model
  # We will use this model's posterior as the basis for importance sampling
  standata = list(N = n, x = x)
  fit <- sampling(stanmodel, data = standata, chains = 4, iter = 2000,
                  refresh = 0, seed = SEED)

  # define target
  target_density <- function(draws, x, ...) {
    dnorm(draws[,1], 2, 0.136, log = TRUE) + dnorm(draws[,2], 0.754, 0.104,
                                                   log = TRUE)
  }


  iw <- moment_match.stanfit(fit,
                             log_prob_target_fun = target_density)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), c(1.9984141, 0.7520985), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw$draws^2,w = exp(iw$log_weights)) - matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights))^2, c(0.01772950, 0.01102333), tolerance = 1e-6)
  expect_equal(iw$pareto_k, 0.377199, tolerance = 1e-6)




  # different formulation

  ratio_density <- function(draws, x, ...) {
    target_density(draws, x, ...) - log_prob_upars.stanfit(draws, x, ...)
  }

  iw <- moment_match.stanfit(fit,
                             log_ratio_fun = ratio_density)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), c(1.9984141, 0.7520985), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw$draws^2,w = exp(iw$log_weights)) - matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights))^2, c(0.01772950, 0.01102333), tolerance = 1e-6)
  expect_equal(iw$pareto_k, 0.377199, tolerance = 1e-6)


})



test_that("moment_match.stanfit works with obs_weights formulation", {






  # Let us generate data from the true data generating mechanism
  # but set one outlier

  SEED <- 24
  set.seed(SEED)
  data_var = 1
  n = as.integer(30)
  x = rnorm(n = n, mean = 0,sd = data_var)
  x[1] <- 9


  # fit the model
  # We will use this model's posterior as the basis for importance sampling
  standata = list(N = n, x = x)
  fit <- sampling(stanmodel, data = standata, chains = 4, iter = 2000, refresh = 0, seed = SEED)

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

  iw <- moment_match.stanfit(fit,
                             log_ratio_fun = ratio_density,
                             k_threshold = 0.0,
                             obs_weights = c(0, rep(1,n - 1)))


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), c(-0.1609051, -0.2804901), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw$draws^2,w = exp(iw$log_weights)) - matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights))^2, c(0.02082176, 0.01873438), tolerance = 1e-6)
  expect_equal(iw$pareto_k, 0.04992905, tolerance = 1e-6)



  # another way


  target_density <- function(draws, x, ...) {
    ratio_density(draws, x, ...) + log_prob_upars.stanfit(draws, x, ...)
  }

  iw <- moment_match.stanfit(fit,
                             log_prob_target_fun = target_density,
                             k_threshold = 0.0,
                             obs_weights = c(0, rep(1,n - 1)))


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), c(-0.1609051, -0.2804901), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw$draws^2,w = exp(iw$log_weights)) - matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights))^2, c(0.02082176, 0.01873438), tolerance = 1e-6)
  expect_equal(iw$pareto_k, 0.04992905, tolerance = 1e-6)


})




test_that("moment_match.stanfit with expectation works", {



  # Let us generate data from the true data generating mechanism

  SEED <- 24
  set.seed(SEED)
  data_var = 1
  n = as.integer(30)
  x = rnorm(n = n, mean = 0,sd = data_var)

  # fit the model
  # We will use this model's posterior as the basis for importance sampling
  standata = list(N = n, x = x)
  fit <- sampling(stanmodel, data = standata, chains = 4, iter = 2000, refresh = 0, seed = SEED)

  # define target
  target_density <- function(draws, x, ...) {
    dnorm(draws[,1], 2, 0.136, log = TRUE) + dnorm(draws[,2], 0.754, 0.104, log = TRUE)
  }


  ex_mm <- moment_match.stanfit(fit,
                                log_prob_target_fun = target_density,
                                expectation_fun =  function(draws, ...) {
                                  draws})

  iw <- moment_match.stanfit(fit,
                             log_prob_target_fun = target_density)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), ex_mm$expectation, tolerance = 1e-6)
  expect_equal(iw$pareto_k, ex_mm$pareto_k, tolerance = 1e-6)
  expect_equal(ex_mm$pareto_kf, c(0.2773750, 0.4282701), tolerance = 1e-6)




  # different formulation

  ratio_density <- function(draws, x, ...) {
    target_density(draws, x, ...) - log_prob_upars.stanfit(draws, x, ...)
  }

  ex_mm <- moment_match.stanfit(fit,
                                log_ratio_fun = ratio_density,
                                expectation_fun =  function(draws, ...) {
                                  draws})

  iw <- moment_match.stanfit(fit,
                             log_ratio_fun = ratio_density)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), ex_mm$expectation, tolerance = 1e-6)
  expect_equal(iw$pareto_k, ex_mm$pareto_k, tolerance = 1e-6)
  expect_equal(ex_mm$pareto_kf, c(0.2773750, 0.4282701), tolerance = 1e-6)


})






test_that("moment_match.stanfit with expectation works with obs_weights formulation", {


  # Let us generate data from the true data generating mechanism
  # but set one outlier

  SEED <- 24
  set.seed(SEED)
  data_var = 1
  n = as.integer(30)
  x = rnorm(n = n, mean = 0,sd = data_var)
  x[1] <- 9


  # fit the model
  # We will use this model's posterior as the basis for importance sampling
  standata = list(N = n, x = x)
  fit <- sampling(stanmodel, data = standata, chains = 4, iter = 2000, refresh = 0, seed = SEED)

  obs_weights <- c(0, rep(1,n - 1))

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



  ratio_density <- function(draws, stanfit, ...) {
    log_lik <- log_lik_stanfit(stanfit, draws)
    colSums((obs_weights - 1) * t(log_lik))
  }

  ex_mm <- moment_match.stanfit(fit,
                                log_ratio_fun = ratio_density,
                                expectation_fun =  function(draws, ...) {
                                  draws},
                                k_threshold = 0.2)

  iw <- moment_match.stanfit(fit,
                             log_ratio_fun = ratio_density,
                             k_threshold = 0.2)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), ex_mm$expectation, tolerance = 1e-6)
  expect_equal(iw$pareto_k, ex_mm$pareto_k, tolerance = 1e-6)
  expect_equal(ex_mm$pareto_kf, c(0.1758149, 0.1814282), tolerance = 1e-6)




  # another way


  target_density <- function(draws, x, ...) {
    ratio_density(draws, x, ...) + log_prob_upars.stanfit(draws, x, ...)
  }

  ex_mm <- moment_match.stanfit(fit,
                                log_prob_target_fun = target_density,
                                expectation_fun =  function(draws, ...) {
                                  draws},
                                k_threshold = 0.2)

  iw <- moment_match.stanfit(fit,
                             log_prob_target_fun = target_density,
                             k_threshold = 0.2)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), ex_mm$expectation, tolerance = 1e-6)
  expect_equal(iw$pareto_k, ex_mm$pareto_k, tolerance = 1e-6)
  expect_equal(ex_mm$pareto_kf, c(0.1758149, 0.1814282), tolerance = 1e-6)



})



