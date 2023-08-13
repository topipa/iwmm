test_that("moment_match works", {
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

  expect_equal(
    matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights)),
    c(4.988783, 4.996672),
    tolerance = 1e-6
  )
  expect_equal(
    matrixStats::colWeightedMeans(iw$draws^2, w = exp(iw$log_weights)) -
      matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights))^2,
    c(0.9243542, 1.0270310),
    tolerance = 1e-6
  )
  expect_equal(iw$pareto_k, 0.3605297, tolerance = 1e-6)


  # another definition

  iw <- moment_match(
    prop_sample,
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density
  )

  expect_equal(matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights)),
    c(4.988783, 4.996672),
    tolerance = 1e-6
  )
  expect_equal(
    matrixStats::colWeightedMeans(iw$draws^2, w = exp(iw$log_weights)) -
      matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights))^2,
    c(0.9243542, 1.0270310),
    tolerance = 1e-6
  )
  expect_equal(iw$pareto_k, 0.3605297, tolerance = 1e-6)



  iw <- moment_match(
    prop_sample,
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density,
    dummy_arg = 123
  )
})

test_that("moment_match with model works", {
  set.seed(7)

  S <- 4000

  prop_mean <- 0
  prop_var <- sqrt(0.1)

  target_mean <- 5
  target_var <- 1

  prop_sample <- matrix(rnorm(2 * S, prop_mean, prop_var), S, 2)
  prop_density <- function(draws, model, ...) {
    dnorm(draws[, 1], prop_mean, prop_var, log = TRUE) +
      dnorm(draws[, 2], prop_mean, prop_var, log = TRUE)
  }

  target_density <- function(draws, model, ...) {
    dnorm(draws[, 1], target_mean, target_var, log = TRUE) +
      dnorm(draws[, 2], target_mean, target_var, log = TRUE)
  }

  ratio_density <- function(draws, model, ...) {
    target_density(draws, model, ...) - prop_density(draws, model, ...)
  }

  test_model <- list()
  test_model$draws <- prop_sample
  test_post_fun <- function(model, ...) {
    model$draws
  }
  unconstrain_pars_fun <- function(model, pars, ...) {
    pars
  }

  prop_sample_pars <- test_post_fun(test_model)
  prop_sample <- unconstrain_pars_fun(test_model, prop_sample_pars)

  iw <- moment_match(
    prop_sample,
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density,
    model = test_model
  )

  expect_equal(matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights)),
    c(4.988783, 4.996672),
    tolerance = 1e-6
  )
  expect_equal(
    matrixStats::colWeightedMeans(iw$draws^2, w = exp(iw$log_weights)) -
      matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights))^2,
    c(0.9243542, 1.0270310),
    tolerance = 1e-6
  )
  expect_equal(iw$pareto_k, 0.3605297, tolerance = 1e-6)

  # another definition

  iw <- moment_match(
    prop_sample,
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density,
    model = test_model
  )


  expect_equal(matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights)),
    c(4.988783, 4.996672),
    tolerance = 1e-6
  )
  expect_equal(
    matrixStats::colWeightedMeans(iw$draws^2, w = exp(iw$log_weights)) -
      matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights))^2,
    c(0.9243542, 1.0270310),
    tolerance = 1e-6
  )
  expect_equal(iw$pareto_k, 0.3605297, tolerance = 1e-6)
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

  target_sample <- matrix(rnorm(2 * S, target_mean, target_var), S, 2)

  iw_mean <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) {
      draws[, 1:2]
    },
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density
  )

  iw2_mean <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) {
      draws[, 1:2]
    },
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density,
    k_threshold = 0.0
  )

  iw3_mean <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) {
      matrix(draws[, 1])
    },
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density,
    k_threshold = -10.0,
    split = TRUE,
    cov_transform = TRUE
  )

  iw3b_mean <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) {
      matrix(draws[, 2])
    },
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density,
    k_threshold = -10.0,
    split = TRUE,
    cov_transform = TRUE
  )

  iw4_mean <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) {
      matrix(draws[, 1])
    },
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density,
    k_threshold = -10.0,
    split = TRUE,
    cov_transform = TRUE,
    restart_transform = TRUE
  )

  iw4b_mean <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) {
      matrix(draws[, 2])
    },
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density,
    k_threshold = -10.0,
    split = TRUE,
    cov_transform = TRUE,
    restart_transform = TRUE
  )

  expect_equal(iw_mean$expectation, c(5.005028, 4.998967), tolerance = 1e-6)
  expect_equal(iw2_mean$expectation, c(5.005028, 4.998967), tolerance = 1e-6)
  expect_equal(c(iw3_mean$expectation, iw3b_mean$expectation),
    c(5.005665, 5.000551),
    tolerance = 1e-6
  )

  # using ratio

  iw_mean <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) {
      draws[, 1:2]
    },
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density
  )

  iw2_mean <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) {
      draws[, 1:2]
    },
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density,
    k_threshold = 0.0,
    cov_transform = TRUE
  )

  iw3_mean <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) {
      matrix(draws[, 1])
    },
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density,
    k_threshold = -10.0,
    split = TRUE,
    cov_transform = TRUE
  )

  iw3b_mean <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) {
      matrix(draws[, 2])
    },
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density,
    k_threshold = -10.0,
    split = TRUE,
    cov_transform = TRUE
  )

  iw4_mean <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) {
      matrix(draws[, 1])
    },
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density,
    k_threshold = -10.0,
    split = TRUE,
    cov_transform = TRUE,
    restart_transform = TRUE
  )

  iw4b_mean <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) {
      matrix(draws[, 2])
    },
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density,
    k_threshold = -10.0,
    split = TRUE,
    cov_transform = TRUE,
    restart_transform = TRUE
  )

  expect_equal(iw_mean$expectation, c(5.005028, 4.998967), tolerance = 1e-6)
  expect_equal(iw2_mean$expectation, c(5.005028, 4.998967), tolerance = 1e-6)
  expect_equal(c(iw3_mean$expectation, iw3b_mean$expectation),
    c(5.005665, 5.000551),
    tolerance = 1e-6
  )
})


test_that("moment_match with expectation works for simple Monte Carlo case", {
  set.seed(6)
  S <- 4000

  target_mean <- 5
  target_var <- 1

  target_density <- function(draws, ...) {
    dnorm(draws[, 1], target_mean, target_var, log = TRUE) +
      dnorm(draws[, 2], target_mean, target_var, log = TRUE)
  }

  target_sample <- matrix(rnorm(2 * S, target_mean, target_var), S, 2)

  iw_mean <- moment_match(
    target_sample,
    expectation_fun = function(draws, ...) {
      draws[, 1:2]
    },
    log_prob_prop_fun = target_density
  )


  iw2_mean <- moment_match(
    target_sample,
    expectation_fun = function(draws, ...) {
      draws[, 1:2]
    },
    log_prob_prop_fun = target_density,
    k_threshold = 0.0,
    cov_transform = TRUE
  )

  iw3_mean <- moment_match(
    target_sample,
    expectation_fun = function(draws, ...) {
      matrix(draws[, 1])
    },
    log_prob_prop_fun = target_density,
    k_threshold = -10.0,
    split = TRUE,
    cov_transform = TRUE
  )

  iw3b_mean <- moment_match(
    target_sample,
    expectation_fun = function(draws, ...) {
      matrix(draws[, 2])
    },
    log_prob_prop_fun = target_density,
    k_threshold = -10.0,
    split = TRUE,
    cov_transform = TRUE
  )

  iw4_mean <- moment_match(
    target_sample,
    expectation_fun = function(draws, ...) {
      matrix(draws[, 1])
    },
    log_prob_prop_fun = target_density,
    k_threshold = -10.0,
    split = TRUE,
    cov_transform = TRUE,
    restart_transform = TRUE
  )

  iw4b_mean <- moment_match(
    target_sample,
    expectation_fun = function(draws, ...) {
      matrix(draws[, 2])
    },
    log_prob_prop_fun = target_density,
    k_threshold = -10.0,
    split = TRUE,
    cov_transform = TRUE,
    restart_transform = TRUE
  )



  expect_equal(iw_mean$expectation, c(5.005361, 4.999969), tolerance = 1e-6)
  expect_equal(iw2_mean$expectation, c(5.005361, 4.999969), tolerance = 1e-6)
  expect_equal(c(iw3_mean$expectation, iw3b_mean$expectation),
    c(5.005361, 4.998519),
    tolerance = 1e-6
  )
})


test_that("moment_match with expectation with model works", {
  set.seed(7)

  S <- 4000

  prop_mean <- 0
  prop_var <- sqrt(0.1)

  target_mean <- 5
  target_var <- 1

  prop_sample <- matrix(rnorm(2 * S, prop_mean, prop_var), S, 2)
  prop_density <- function(draws, model, ...) {
    dnorm(draws[, 1], prop_mean, prop_var, log = TRUE) +
      dnorm(draws[, 2], prop_mean, prop_var, log = TRUE)
  }

  target_density <- function(draws, model, ...) {
    dnorm(draws[, 1], target_mean, target_var, log = TRUE) +
      dnorm(draws[, 2], target_mean, target_var, log = TRUE)
  }

  ratio_density <- function(draws, model, ...) {
    target_density(draws, model, ...) - prop_density(draws, model, ...)
  }

  test_model <- list()
  test_model$draws <- prop_sample
  test_post_fun <- function(model, ...) {
    model$draws
  }
  unconstrain_pars_fun <- function(model, pars, ...) {
    pars
  }

  prop_sample_pars <- test_post_fun(test_model)
  prop_sample <- unconstrain_pars_fun(test_model, prop_sample_pars)

  ex_mm <- moment_match(prop_sample,
    expectation_fun = function(draws, ...) {
      draws
    },
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density,
    model = test_model
  )

  iw <- moment_match(prop_sample,
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density,
    model = test_model
  )

  expect_equal(matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights)),
    ex_mm$expectation,
    tolerance = 1e-6
  )
  expect_equal(iw$pareto_k, ex_mm$pareto_k, tolerance = 1e-6)
  expect_equal(ex_mm$pareto_kf, c(0.3875001, 0.3379889), tolerance = 1e-6)

  # another definition

  ex_mm <- moment_match(prop_sample,
    expectation_fun = function(draws, ...) {
      draws
    },
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density,
    model = test_model
  )

  iw <- moment_match(prop_sample,
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density,
    model = test_model
  )


  expect_equal(matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights)),
    ex_mm$expectation,
    tolerance = 1e-6
  )
  expect_equal(iw$pareto_k, ex_mm$pareto_k, tolerance = 1e-6)
  expect_equal(ex_mm$pareto_kf, c(0.3875001, 0.3379889), tolerance = 1e-6)
})
