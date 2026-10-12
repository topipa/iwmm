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

test_that("split = TRUE is exact when the scaling step is accepted", {
  # Proposal N(0, I). The expectation is a likelihood that is narrow in
  # both coordinates, so moment matching shifts and scales the draws.
  # The split estimate then uses the inverse of the accumulated affine
  # transformation, which was exact only for pure shifts.
  y <- 1.5
  s <- 0.05
  exact <- 2 * dnorm(y, 0, sqrt(1 + s^2), log = TRUE)
  err <- vapply(1:3, function(seed) {
    set.seed(seed)
    draws <- matrix(rnorm(8000), 4000, 2)
    mm <- moment_match(
      draws,
      log_prob_prop_fun = function(draws, ...) -0.5 * rowSums(draws^2),
      expectation_fun = function(draws, ...) {
        matrix(dnorm(y, draws[, 1], s, log = TRUE) +
          dnorm(y, draws[, 2], s, log = TRUE))
      },
      log_expectation_fun = TRUE,
      k_threshold = 0.5,
      cov_transform = TRUE,
      split = TRUE
    )
    log(mm$expectation) - exact
  }, numeric(1))
  expect_lt(max(abs(err)), 0.05)
})

test_that("split = TRUE evaluates the proposal at the inverse transformations", {
  # The returned draws are T1(theta) in the first half and T2(theta) in
  # the second half. The last two calls of log_prob_prop_fun evaluate
  # T2^-1 of the first half and T1^-1 of the second half. With a target
  # density, T1 shifts and scales the draws; with restart_transform =
  # TRUE, T2 is found again from the original draws, so T2 differs from
  # T1 and both inverses are tested.
  set.seed(6)
  S <- 4000
  prop_sample <- matrix(rnorm(2 * S, 0, sqrt(0.1)), S, 2)
  calls <- list()
  prop_density <- function(draws, ...) {
    calls[[length(calls) + 1]] <<- draws
    dnorm(draws[, 1], 0, sqrt(0.1), log = TRUE) +
      dnorm(draws[, 2], 0, sqrt(0.1), log = TRUE)
  }
  target_density <- function(draws, ...) {
    dnorm(draws[, 1], 5, 1, log = TRUE) + dnorm(draws[, 2], 5, 1, log = TRUE)
  }
  mm <- moment_match(
    prop_sample,
    expectation_fun = function(draws, ...) matrix(draws[, 1]),
    log_prob_prop_fun = prop_density,
    log_prob_target_fun = target_density,
    k_threshold = -10,
    split = TRUE,
    cov_transform = TRUE,
    restart_transform = TRUE
  )
  first <- seq_len(S / 2)
  second <- -first
  # the affine maps T1 and T2 from the (theta, T(theta)) pairs
  affine <- function(from, to) {
    x <- cbind(1, from)
    coef <- solve(crossprod(x), crossprod(x, to))
    function(z) cbind(1, z) %*% coef
  }
  T1 <- affine(prop_sample[first, ], mm$draws[first, ])
  T2 <- affine(prop_sample[second, ], mm$draws[second, ])
  # the two transformations must differ, else the test cannot detect a
  # mix-up of their shifts
  expect_gt(max(abs(T1(prop_sample[1:10, ]) - T2(prop_sample[1:10, ]))), 0.01)
  n_calls <- length(calls)
  inv1 <- calls[[n_calls - 1]][first, , drop = FALSE]
  inv2 <- calls[[n_calls]][second, , drop = FALSE]
  expect_equal(unname(T2(inv1)), unname(mm$draws[first, ]), tolerance = 1e-8)
  expect_equal(unname(T1(inv2)), unname(mm$draws[second, ]), tolerance = 1e-8)
})
