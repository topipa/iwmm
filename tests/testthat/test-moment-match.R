test_that("moment_match.matrix works", {

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
    c(5, 5),
    tolerance = 1e-2
  )
  expect_equal(
    matrixStats::colWeightedMeans(iw$draws^2, w = exp(iw$log_weights)) -
      matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights))^2,
    c(1, 1),
    tolerance = 1e-1
  )
  expect_equal(iw$diagnostics$pareto_k, 0.35, tolerance = 1e-1)


  # another definition

  iw <- moment_match(
    prop_sample,
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density
  )

  expect_equal(matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights)),
    c(5, 5),
    tolerance = 1e-2
  )
  expect_equal(
    matrixStats::colWeightedMeans(iw$draws^2, w = exp(iw$log_weights)) -
      matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights))^2,
    c(1, 1),
    tolerance = 1e-1
  )
  expect_equal(iw$diagnostics$pareto_k, 0.35, tolerance = 1e-1)
})

test_that("moment_match.draws_matrix works", {
  set.seed(7)

  S <- 4000

  prop_mean <- 0
  prop_var <- sqrt(0.1)

  target_mean <- 5
  target_var <- 1

  prop_sample <- matrix(rnorm(2 * S, prop_mean, prop_var), S, 2)
  colnames(prop_sample) <- paste0("V", 1:2)
  prop_sample <- posterior::as_draws_matrix(prop_sample)

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

  iw_means <- posterior::summarise_draws(
    iw$draws,
    matrixStats::weightedMean,
    .args = list(w = exp(iw$log_weights))
  )
  iw_vars <- posterior::summarise_draws(
    iw$draws,
    matrixStats::weightedVar,
    .args = list(w = 1000 * exp(iw$log_weights))
  )

  expect_equal(
    iw_means$`matrixStats::weightedMean`,
    c(5, 5),
    tolerance = 1e-1
  )
  expect_equal(
    iw_vars$`matrixStats::weightedVar`,
    c(1, 1),
    tolerance = 1e-1
  )
  expect_equal(iw$diagnostics$pareto_k, 0.4, tolerance = 1e-1)


  # another definition

  iw <- moment_match(
    prop_sample,
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density
  )

  iw_means <- posterior::summarise_draws(
    iw$draws,
    matrixStats::weightedMean,
    .args = list(w = exp(iw$log_weights))
  )
  # TODO: why does thia fail when weights are small and sum to one?
  iw_vars <- posterior::summarise_draws(
    iw$draws,
    matrixStats::weightedVar,
    .args = list(w = 1000 * exp(iw$log_weights))
  )

  expect_equal(
    iw_means$`matrixStats::weightedMean`,
    c(5, 5),
    tolerance = 1e-1
  )
  expect_equal(
    iw_vars$`matrixStats::weightedVar`,
    c(1, 1),
    tolerance = 1e-1
  )
  expect_equal(iw$diagnostics$pareto_k, 0.4, tolerance = 1e-1)
})


test_that("moment_match with other posterior formats works", {
  set.seed(7)

  S <- 400

  prop_mean <- 0
  prop_var <- sqrt(0.1)

  target_mean <- 5
  target_var <- 1

  prop_sample <- matrix(rnorm(2 * S, prop_mean, prop_var), S, 2)
  colnames(prop_sample) <- paste0("V", 1:2)
  prop_sample <- posterior::as_draws_matrix(prop_sample)

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

  iw_draws_matrix <- moment_match(prop_sample,
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density
  )

  iw_draws_array <- moment_match(posterior::as_draws_array(prop_sample),
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density
  )
  iw_draws_df <- moment_match(posterior::as_draws_df(prop_sample),
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density
  )
  iw_draws_list <- moment_match(posterior::as_draws_list(prop_sample),
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density
  )
  iw_draws_rvars <- moment_match(posterior::as_draws_rvars(prop_sample),
    log_prob_prop_fun = prop_density,
    log_ratio_fun = ratio_density
  )

  expect_equal(
    iw_draws_matrix, iw_draws_array
  )
  expect_equal(
    iw_draws_matrix, iw_draws_df
  )
  expect_equal(
    iw_draws_matrix, iw_draws_list
  )
  expect_equal(
    iw_draws_matrix, iw_draws_rvars
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
    c(5, 5),
    tolerance = 1e-2
  )
  expect_equal(
    matrixStats::colWeightedMeans(iw$draws^2, w = exp(iw$log_weights)) -
      matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights))^2,
    c(1, 1),
    tolerance = 1e-1
  )
  expect_equal(iw$diagnostics$pareto_k, 0.35, tolerance = 1e-1)

  # another definition

  iw <- moment_match(
    prop_sample,
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density,
    model = test_model
  )


  expect_equal(matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights)),
    c(5, 5),
    tolerance = 1e-2
  )
  expect_equal(
    matrixStats::colWeightedMeans(iw$draws^2, w = exp(iw$log_weights)) -
      matrixStats::colWeightedMeans(iw$draws, w = exp(iw$log_weights))^2,
    c(1, 1),
    tolerance = 1e-1
  )
  expect_equal(iw$diagnostics$pareto_k, 0.35, tolerance = 1e-1)
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

  expect_equal(iw_mean$expectation, c(5, 5), tolerance = 1e-2)
  expect_equal(iw2_mean$expectation, c(5, 5), tolerance = 1e-2)
  expect_equal(c(iw3_mean$expectation, iw3b_mean$expectation),
    c(5, 5),
    tolerance = 1e-2
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

  expect_equal(iw_mean$expectation, c(5, 5), tolerance = 1e-2)
  expect_equal(iw2_mean$expectation, c(5, 5), tolerance = 1e-2)
  expect_equal(c(iw3_mean$expectation, iw3b_mean$expectation),
    c(5, 5),
    tolerance = 1e-2
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


  iw2_mean <- suppressWarnings(moment_match(
    target_sample,
    expectation_fun = function(draws, ...) {
      draws[, 1:2]
    },
    log_prob_prop_fun = target_density,
    k_threshold = 0.0,
    cov_transform = TRUE
  ))

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



  expect_equal(iw_mean$expectation, c(5, 5), tolerance = 1e-2)
  expect_equal(iw2_mean$expectation, c(5, 5), tolerance = 1e-2)
  expect_equal(c(iw3_mean$expectation, iw3b_mean$expectation),
    c(5, 5),
    tolerance = 1e-2
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
    tolerance = 1e-2
  )
  expect_equal(iw$diagnostics$pareto_k, ex_mm$diagnostics$pareto_k, tolerance = 1e-2)
  expect_equal(ex_mm$diagnostics$pareto_kf, c(0.4, 0.4), tolerance = 1e-1)

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
    tolerance = 1e-2
  )
  expect_equal(iw$diagnostics$pareto_k, ex_mm$diagnostics$pareto_k, tolerance = 1e-6)
  expect_equal(ex_mm$diagnostics$pareto_kf, c(0.4, 0.4), tolerance = 1e-1)
})

test_that("moment_match with expectation works with transformation", {
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
      exp(draws[, 1:2])
    },
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density,
    draws_transformation_fun = NULL
  )

  draws_transformation_fun <- function(draws, ...) {
    exp(draws[, 1:2])
  }

  iw2_mean <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) {
      draws[, 1:2]
    },
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density,
    draws_transformation_fun = draws_transformation_fun
  )


  expect_equal(iw_mean$expectation, c(245, 245), tolerance = 1e-2)
  expect_equal(iw2_mean$expectation, iw_mean$expectation, tolerance = 1e-6)
})
