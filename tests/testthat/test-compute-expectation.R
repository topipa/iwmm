test_that("compute_expectation.matrix works", {
  set.seed(7)

  S <- 4000

  prop_mean <- 0
  prop_var <- sqrt(0.1)

  target_mean <- 5
  target_var <- 1

  prop_sample <- matrix(rnorm(2 * S, prop_mean, prop_var), S, 2)
  prop_density <- function(draws, ...) {
    dnorm(draws[, 1], prop_mean, prop_var, log = TRUE) +
      dnorm(draws[, 2], prop_mean, prop_var, log = TRUE)
  }

  target_density <- function(draws, ...) {
    dnorm(draws[, 1], target_mean, target_var, log = TRUE) +
      dnorm(draws[, 2], target_mean, target_var, log = TRUE)
  }

  ratio_density <- function(draws, ...) {
    target_density(draws, ...) - prop_density(draws, ...)
  }

  iw <- moment_match(prop_sample,
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density
  )

  iw_expectation_first_moment <- compute_expectation(iw$draws, iw$log_weights, expectation_fun = function(draws, ...) {
    draws
  })
  iw_expectation_second_moment <- compute_expectation(iw$draws, iw$log_weights, expectation_fun = function(draws, ...) {
    draws^2
  })

  expect_equal(
    iw_expectation_first_moment$expectation,
    c(5, 5),
    tolerance = 1e-2
  )
  expect_equal(
    iw_expectation_second_moment$expectation -
      iw_expectation_first_moment$expectation^2,
    c(1, 1),
    tolerance = 1e-1
  )
  expect_equal(iw_expectation_first_moment$diagnostics$pareto_kf, 0.4, tolerance = 1e-1)
  expect_equal(iw_expectation_second_moment$diagnostics$pareto_kf, 0.4, tolerance = 1e-1)
})

test_that("compute_expectation.adapted_importance_sampling works", {
  set.seed(7)

  S <- 4000

  prop_mean <- 0
  prop_var <- sqrt(0.1)

  target_mean <- 5
  target_var <- 1

  prop_sample <- matrix(rnorm(2 * S, prop_mean, prop_var), S, 2)
  prop_density <- function(draws, ...) {
    dnorm(draws[, 1], prop_mean, prop_var, log = TRUE) +
      dnorm(draws[, 2], prop_mean, prop_var, log = TRUE)
  }

  target_density <- function(draws, ...) {
    dnorm(draws[, 1], target_mean, target_var, log = TRUE) +
      dnorm(draws[, 2], target_mean, target_var, log = TRUE)
  }

  ratio_density <- function(draws, ...) {
    target_density(draws, ...) - prop_density(draws, ...)
  }

  iw <- moment_match(prop_sample,
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density
  )

  iw_expectation_first_moment <- compute_expectation(iw, expectation_fun = function(draws, ...) {
    draws
  })
  iw_expectation_second_moment <- compute_expectation(iw, expectation_fun = function(draws, ...) {
    draws^2
  })

  expect_equal(
    iw_expectation_first_moment$expectation,
    c(5, 5),
    tolerance = 1e-2
  )
  expect_equal(
    iw_expectation_second_moment$expectation -
      iw_expectation_first_moment$expectation^2,
    c(1, 1),
    tolerance = 1e-1
  )
  expect_equal(iw_expectation_first_moment$diagnostics$pareto_kf, 0.4, tolerance = 1e-1)
  expect_equal(iw_expectation_second_moment$diagnostics$pareto_kf, 0.4, tolerance = 1e-1)
})
