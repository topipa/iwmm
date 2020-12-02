
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

  target_sample <- matrix(rnorm(2 * S, target_mean, target_var), S, 2)

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


})
