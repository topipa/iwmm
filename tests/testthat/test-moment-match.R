


test_that("moment_match works", {

  set.seed(7)

  S <- 4000

  prop_mean <- 0
  prop_var <- sqrt(0.1)

  target_mean <- 5
  target_var <- 1

  prop_sample <- matrix(rnorm(2 * S, prop_mean, prop_var), S, 2)
  prop_density <- function(draws, ...) {
    dnorm(draws[,1], prop_mean, prop_var, log = TRUE) + dnorm(draws[,2], prop_mean, prop_var, log = TRUE)
  }

  target_density <- function(draws ,...) {
    dnorm(draws[,1], target_mean, target_var, log = TRUE) + dnorm(draws[,2], target_mean, target_var, log = TRUE)
  }

  ratio_density <- function(draws ,...) {
    target_density(draws, ...) - prop_density(draws, ...)
  }

  iw <- moment_match(prop_sample,
              log_prob_prop_draws_fun = prop_density,
              log_ratio_draws_fun = ratio_density)

  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), c(4.988783, 4.996672), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw$draws^2,w = exp(iw$log_weights)) - matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights))^2, c(0.9243542, 1.0270310), tolerance = 1e-6)
  expect_equal(iw$pareto_k, 0.3605297, tolerance = 1e-6)


  # another definition

  iw <- moment_match(prop_sample,
             log_prob_prop_draws_fun = prop_density,
             log_prob_target_draws_fun = target_density)

  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), c(4.988783, 4.996672), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw$draws^2,w = exp(iw$log_weights)) - matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights))^2, c(0.9243542, 1.0270310), tolerance = 1e-6)
  expect_equal(iw$pareto_k, 0.3605297, tolerance = 1e-6)



  iw <- moment_match(prop_sample,
                     log_prob_prop_draws_fun = prop_density,
                     log_prob_target_draws_fun = target_density,
                     dummy_arg = 123)

})





test_that("moment_match with model works", {

  set.seed(7)

  S <- 4000

  prop_mean <- 0
  prop_var <- sqrt(0.1)

  target_mean <- 5
  target_var <- 1

  prop_sample <- matrix(rnorm(2 * S, prop_mean, prop_var), S, 2)
  prop_density <- function(draws, x, ...) {
    # print(length(x$draws))
    dnorm(draws[,1], prop_mean, prop_var, log = TRUE) + dnorm(draws[,2], prop_mean, prop_var, log = TRUE)
  }

  target_density <- function(draws, x, ...) {
    dnorm(draws[,1], target_mean, target_var, log = TRUE) + dnorm(draws[,2], target_mean, target_var, log = TRUE)
  }

  ratio_density <- function(draws, x, ...) {
    target_density(draws, x, ...) - prop_density(draws, x, ...)
  }

  test_model <- list()
  test_model$draws <- prop_sample
  test_post_draws_fun <- function(x, ...) {
    x$draws
  }
  unconstrain_pars_fun <- function(x, pars, ...) {
    pars
  }

  prop_sample_pars <- test_post_draws_fun(test_model)
  prop_sample <- unconstrain_pars_fun(test_model, prop_sample_pars)

  iw <- moment_match(prop_sample,
                     log_prob_prop_draws_fun = prop_density,
                     log_ratio_draws_fun = ratio_density,
                     x = test_model)



  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), c(4.988783, 4.996672), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw$draws^2,w = exp(iw$log_weights)) - matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights))^2, c(0.9243542, 1.0270310), tolerance = 1e-6)
  expect_equal(iw$pareto_k, 0.3605297, tolerance = 1e-6)


  # another definition

  iw <- moment_match(prop_sample,
                     log_prob_prop_draws_fun = prop_density,
                     log_prob_target_draws_fun = target_density,
                     x = test_model)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), c(4.988783, 4.996672), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw$draws^2,w = exp(iw$log_weights)) - matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights))^2, c(0.9243542, 1.0270310), tolerance = 1e-6)
  expect_equal(iw$pareto_k, 0.3605297, tolerance = 1e-6)





})


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
  fit <- sampling(stanmodel, data = standata, chains = 4, iter = 2000, refresh = 0)

  # and compute the posterior mean
  # postmean <- colMeans(as.matrix(fit)[,1:2])
  # postsd <- matrixStats::colSds(as.matrix(fit)[,1:2])

  # define target as gaussian with mean of mu = 2
  target_density <- function(draws, x, ...) {
    dnorm(draws[,1], 2, 0.136, log = TRUE) + dnorm(draws[,2], 0.754, 0.104, log = TRUE)
  }


  iw <- moment_match.stanfit(fit,
                     log_prob_target_draws_fun = target_density)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), c(1.9916292, 0.7592228), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw$draws^2,w = exp(iw$log_weights)) - matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights))^2, c(0.01779679, 0.01091765), tolerance = 1e-6)
  expect_equal(iw$pareto_k, 0.1848982, tolerance = 1e-6)




  # different formulation

  ratio_density <- function(draws, x, ...) {
    target_density(draws, x, ...) - log_prob_upars.stanfit(draws, x, ...)
  }

  iw <- moment_match.stanfit(fit,
                             log_ratio_draws_fun = ratio_density)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), c(1.9916292, 0.7592228), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw$draws^2,w = exp(iw$log_weights)) - matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights))^2, c(0.01779679, 0.01091765), tolerance = 1e-6)
  expect_equal(iw$pareto_k, 0.1848982, tolerance = 1e-6)


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
  fit <- sampling(stanmodel, data = standata, chains = 4, iter = 2000, refresh = 0)

  # and compute the posterior mean
  # postmean <- colMeans(unconstrain_pars.stanfit(fit, as.matrix(fit)))

  # fit_loo <- sampling(stanmodel, data = list(N = n - 1, x = x[-c(1)]), chains = 4, iter = 2000, refresh = 0)
  # postmean_loo <- colMeans(unconstrain_pars.stanfit(fit_loo, as.matrix(fit_loo)))


  obs_weights <- c(0, rep(1,n - 1))

  log_lik_stanfit <- function(x, upars, parameter_name = "log_lik",
                                      ...) {
    ll <- loo::extract_log_lik(x, parameter_name, merge_chains = TRUE)
    S <- nrow(upars)
    n <- ncol(ll)
    out <- matrix(0,S,n)
    for (s in seq_len(S)) {
      out[s,] <- rstan::constrain_pars(x, upars = upars[s, ])[[parameter_name]]
    }
    out
  }



  ratio_density <- function(draws, x, ...) {
    log_lik <- log_lik_stanfit(x, draws)
    colSums((obs_weights - 1) * t(log_lik))
  }

  iw <- moment_match.stanfit(fit,
                             log_ratio_draws_fun = ratio_density,
                             k_threshold = 0.0)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), c(-0.1572105, -0.2797750), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw$draws^2,w = exp(iw$log_weights)) - matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights))^2, c(0.01987301, 0.01801787), tolerance = 1e-6)
  expect_equal(iw$pareto_k, -0.03150002, tolerance = 1e-6)



  # another way


  target_density <- function(draws, x, ...) {
    ratio_density(draws, x, ...) + log_prob_upars.stanfit(draws, x, ...)
  }

  iw <- moment_match.stanfit(fit,
                             log_prob_target_draws_fun = target_density,
                             k_threshold = 0.0)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), c(-0.1572105, -0.2797750), tolerance = 1e-6)
  expect_equal(matrixStats::colWeightedMeans(iw$draws^2,w = exp(iw$log_weights)) - matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights))^2, c(0.01987301, 0.01801787), tolerance = 1e-6)
  expect_equal(iw$pareto_k, -0.03150002, tolerance = 1e-6)


})








test_that("moment_match with expectation works for target and ratio", {


  set.seed(6)
  S <- 4000

  prop_mean <- 0
  prop_var <- sqrt(0.1)

  target_mean <- 5
  target_var <- 1

  prop_sample <- matrix(rnorm(2 * S, prop_mean, prop_var), S, 2)
  prop_density <- function(draws, ...) {
    dnorm(draws[,1], prop_mean, prop_var, log = TRUE) + dnorm(draws[,2], prop_mean, prop_var, log = TRUE)
  }

  target_density <- function(draws ,...) {
    dnorm(draws[,1], target_mean, target_var, log = TRUE) + dnorm(draws[,2], target_mean, target_var, log = TRUE)
  }

  ratio_density <- function(draws ,...) {
    target_density(draws, ...) - prop_density(draws, ...)
  }

  target_sample <- matrix(rnorm(2 * S, target_mean, target_var), S, 2)


  iw_mean <- moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    draws[,1:2]},
    log_prob_prop_draws_fun = prop_density,
    log_prob_target_draws_fun = target_density)


  iw2_mean <- moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    draws[,1:2]},
    log_prob_prop_draws_fun = prop_density,
    log_prob_target_draws_fun = target_density,
    k_threshold = 0.0)



  iw3_mean <- moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,1])},
    log_prob_prop_draws_fun = prop_density,
    log_prob_target_draws_fun = target_density,
    k_threshold = -10.0, split = TRUE, cov_transform = TRUE)

  iw3b_mean <- moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,2])},
    log_prob_prop_draws_fun = prop_density,
    log_prob_target_draws_fun = target_density,
    k_threshold = -10.0, split = TRUE, cov_transform = TRUE)

  iw4_mean <- moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,1])},
    log_prob_prop_draws_fun = prop_density,
    log_prob_target_draws_fun = target_density,
    k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)

  iw4b_mean <- moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,2])},
    log_prob_prop_draws_fun = prop_density,
    log_prob_target_draws_fun = target_density,
    k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)



  expect_equal(iw_mean$expectation, c(5.005028, 4.998967), tolerance = 1e-6)
  expect_equal(iw2_mean$expectation, c(5.005028, 4.998967), tolerance = 1e-6)
  expect_equal(c(iw3_mean$expectation,iw3b_mean$expectation), c(5.005665, 5.000551), tolerance = 1e-6)








  # using ratio



  iw_mean <- moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    draws[,1:2]},
    log_prob_prop_draws_fun = prop_density,
    log_ratio_draws_fun = ratio_density)


  iw2_mean <- moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    draws[,1:2]},
    log_prob_prop_draws_fun = prop_density,
    log_ratio_draws_fun = ratio_density,
    k_threshold = 0.0, cov_transform = TRUE)



  iw3_mean <- moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,1])},
    log_prob_prop_draws_fun = prop_density,
    log_ratio_draws_fun = ratio_density,
    k_threshold = -10.0, split = TRUE, cov_transform = TRUE)

  iw3b_mean <- moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,2])},
    log_prob_prop_draws_fun = prop_density,
    log_ratio_draws_fun = ratio_density,
    k_threshold = -10.0, split = TRUE, cov_transform = TRUE)

  iw4_mean <- moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,1])},
    log_prob_prop_draws_fun = prop_density,
    log_ratio_draws_fun = ratio_density,
    k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)

  iw4b_mean <- moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,2])},
    log_prob_prop_draws_fun = prop_density,
    log_ratio_draws_fun = ratio_density,
    k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)




  expect_equal(iw_mean$expectation, c(5.005028, 4.998967), tolerance = 1e-6)
  expect_equal(iw2_mean$expectation, c(5.005028, 4.998967), tolerance = 1e-6)
  expect_equal(c(iw3_mean$expectation,iw3b_mean$expectation), c(5.005665, 5.000551), tolerance = 1e-6)









})


test_that("moment_match with expectation works for simple Monte Carlo case", {


  set.seed(6)
  S <- 4000

  target_mean <- 5
  target_var <- 1

  target_density <- function(draws ,...) {
    dnorm(draws[,1], target_mean, target_var, log = TRUE) + dnorm(draws[,2], target_mean, target_var, log = TRUE)
  }

  target_sample <- matrix(rnorm(2 * S, target_mean, target_var), S, 2)


  iw_mean <- moment_match(target_sample, expectation_fun =  function(draws, ...) {
    draws[,1:2]},
    log_prob_prop_draws_fun = target_density)


  iw2_mean <- moment_match(target_sample, expectation_fun =  function(draws, ...) {
    draws[,1:2]},
    log_prob_prop_draws_fun = target_density,
    k_threshold = 0.0, cov_transform = TRUE)



  iw3_mean <- moment_match(target_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,1])},
    log_prob_prop_draws_fun = target_density,
    k_threshold = -10.0, split = TRUE, cov_transform = TRUE)

  iw3b_mean <- moment_match(target_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,2])},
    log_prob_prop_draws_fun = target_density,
    k_threshold = -10.0, split = TRUE, cov_transform = TRUE)

  iw4_mean <- moment_match(target_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,1])},
    log_prob_prop_draws_fun = target_density,
    k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)

  iw4b_mean <- moment_match(target_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,2])},
    log_prob_prop_draws_fun = target_density,
    k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)



  expect_equal(iw_mean$expectation, c(5.005361,4.999969), tolerance = 1e-6)
  expect_equal(iw2_mean$expectation, c(5.005361, 4.999969), tolerance = 1e-6)
  expect_equal(c(iw3_mean$expectation,iw3b_mean$expectation), c(5.005361, 4.998519), tolerance = 1e-6)






})





test_that("moment_match with expectation with model works", {

  set.seed(7)

  S <- 4000

  prop_mean <- 0
  prop_var <- sqrt(0.1)

  target_mean <- 5
  target_var <- 1

  prop_sample <- matrix(rnorm(2 * S, prop_mean, prop_var), S, 2)
  prop_density <- function(draws, x, ...) {
    # print(length(x$draws))
    dnorm(draws[,1], prop_mean, prop_var, log = TRUE) + dnorm(draws[,2], prop_mean, prop_var, log = TRUE)
  }

  target_density <- function(draws, x, ...) {
    dnorm(draws[,1], target_mean, target_var, log = TRUE) + dnorm(draws[,2], target_mean, target_var, log = TRUE)
  }

  ratio_density <- function(draws, x, ...) {
    target_density(draws, x, ...) - prop_density(draws, x, ...)
  }

  test_model <- list()
  test_model$draws <- prop_sample
  test_post_draws_fun <- function(x, ...) {
    x$draws
  }
  unconstrain_pars_fun <- function(x, pars, ...) {
    pars
  }

  prop_sample_pars <- test_post_draws_fun(test_model)
  prop_sample <- unconstrain_pars_fun(test_model, prop_sample_pars)




  ex_mm <- moment_match(prop_sample,
                        expectation_fun =  function(draws, ...) {
                          draws},
                        log_prob_prop_draws_fun = prop_density,
                        log_ratio_draws_fun = ratio_density,
                        x = test_model)

  iw <- moment_match(prop_sample,
                     log_prob_prop_draws_fun = prop_density,
                     log_ratio_draws_fun = ratio_density,
                     x = test_model)



  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), ex_mm$expectation, tolerance = 1e-6)
  expect_equal(iw$pareto_k, ex_mm$pareto_k, tolerance = 1e-6)
  expect_equal(ex_mm$pareto_kf, c(0.3875001, 0.3379889), tolerance = 1e-6)




  # another definition

  ex_mm <- moment_match(prop_sample,
                        expectation_fun =  function(draws, ...) {
                          draws},
                        log_prob_prop_draws_fun = prop_density,
                        log_prob_target_draws_fun = target_density,
                        x = test_model)

  iw <- moment_match(prop_sample,
                     log_prob_prop_draws_fun = prop_density,
                     log_prob_target_draws_fun = target_density,
                     x = test_model)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), ex_mm$expectation, tolerance = 1e-6)
  expect_equal(iw$pareto_k, ex_mm$pareto_k, tolerance = 1e-6)
  expect_equal(ex_mm$pareto_kf, c(0.3875001, 0.3379889), tolerance = 1e-6)





})




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
  fit <- sampling(stanmodel, data = standata, chains = 4, iter = 2000, refresh = 0)

  # and compute the posterior mean
  # postmean <- colMeans(as.matrix(fit)[,1:2])
  # postsd <- matrixStats::colSds(as.matrix(fit)[,1:2])

  # define target as gaussian with mean of mu = 2
  target_density <- function(draws, x, ...) {
    dnorm(draws[,1], 2, 0.136, log = TRUE) + dnorm(draws[,2], 0.754, 0.104, log = TRUE)
  }


  ex_mm <- moment_match.stanfit(fit,
                                log_prob_target_draws_fun = target_density,
                                expectation_fun =  function(draws, ...) {
                                  draws})

  iw <- moment_match.stanfit(fit,
                             log_prob_target_draws_fun = target_density)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), ex_mm$expectation, tolerance = 1e-6)
  expect_equal(iw$pareto_k, ex_mm$pareto_k, tolerance = 1e-6)
  expect_equal(ex_mm$pareto_kf, c(0.2083834, 0.1870787), tolerance = 1e-6)




  # different formulation

  ratio_density <- function(draws, x, ...) {
    target_density(draws, x, ...) - log_prob_upars.stanfit(draws, x, ...)
  }

  ex_mm <- moment_match.stanfit(fit,
                                log_ratio_draws_fun = ratio_density,
                                expectation_fun =  function(draws, ...) {
                                  draws})

  iw <- moment_match.stanfit(fit,
                             log_ratio_draws_fun = ratio_density)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), ex_mm$expectation, tolerance = 1e-6)
  expect_equal(iw$pareto_k, ex_mm$pareto_k, tolerance = 1e-6)
  expect_equal(ex_mm$pareto_kf, c(0.2083834, 0.1870787), tolerance = 1e-6)


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
  fit <- sampling(stanmodel, data = standata, chains = 4, iter = 2000, refresh = 0)

  # and compute the posterior mean
  # postmean <- colMeans(unconstrain_pars.stanfit(fit, as.matrix(fit)))

  # fit_loo <- sampling(stanmodel, data = list(N = n - 1, x = x[-c(1)]), chains = 4, iter = 2000, refresh = 0)
  # postmean_loo <- colMeans(unconstrain_pars.stanfit(fit_loo, as.matrix(fit_loo)))


  obs_weights <- c(0, rep(1,n - 1))

  log_lik_stanfit <- function(x, upars, parameter_name = "log_lik",
                              ...) {
    ll <- loo::extract_log_lik(x, parameter_name, merge_chains = TRUE)
    S <- nrow(upars)
    n <- ncol(ll)
    out <- matrix(0,S,n)
    for (s in seq_len(S)) {
      out[s,] <- rstan::constrain_pars(x, upars = upars[s, ])[[parameter_name]]
    }
    out
  }



  ratio_density <- function(draws, x, ...) {
    log_lik <- log_lik_stanfit(x, draws)
    colSums((obs_weights - 1) * t(log_lik))
  }

  ex_mm <- moment_match.stanfit(fit,
                                log_ratio_draws_fun = ratio_density,
                                expectation_fun =  function(draws, ...) {
                                  draws},
                                k_threshold = 0.0)

  iw <- moment_match.stanfit(fit,
                             log_ratio_draws_fun = ratio_density,
                             k_threshold = 0.0)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), ex_mm$expectation, tolerance = 1e-6)
  expect_equal(iw$pareto_k, ex_mm$pareto_k, tolerance = 1e-6)
  expect_equal(ex_mm$pareto_kf, c(-0.0577350, -0.2304642), tolerance = 1e-6)




  # another way


  target_density <- function(draws, x, ...) {
    ratio_density(draws, x, ...) + log_prob_upars.stanfit(draws, x, ...)
  }

  ex_mm <- moment_match.stanfit(fit,
                                log_prob_target_draws_fun = target_density,
                                expectation_fun =  function(draws, ...) {
                                  draws},
                                k_threshold = 0.0)

  iw <- moment_match.stanfit(fit,
                             log_prob_target_draws_fun = target_density,
                             k_threshold = 0.0)


  expect_equal(matrixStats::colWeightedMeans(iw$draws,w = exp(iw$log_weights)), ex_mm$expectation, tolerance = 1e-6)
  expect_equal(iw$pareto_k, ex_mm$pareto_k, tolerance = 1e-6)
  expect_equal(ex_mm$pareto_kf, c(-0.0577350, -0.2304642), tolerance = 1e-6)


})












