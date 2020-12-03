

# TODO
# test variance
# make deterministic instead of rnorm


test_that("expectation_moment_match works", {


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


  iw_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    draws[,1:2]},
    log_prob_prop_draws_fun = prop_density,
    log_prob_target_draws_fun = target_density)


  iw2_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    draws[,1:2]},
    log_prob_prop_draws_fun = prop_density,
    log_prob_target_draws_fun = target_density,
    moment_match = TRUE, k_threshold = 0.0, cov_transform = TRUE)



  iw3_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,1])},
    log_prob_prop_draws_fun = prop_density,
    log_prob_target_draws_fun = target_density,
    moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE)

  iw3b_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,2])},
    log_prob_prop_draws_fun = prop_density,
    log_prob_target_draws_fun = target_density,
    moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE)

  iw4_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,1])},
    log_prob_prop_draws_fun = prop_density,
    log_prob_target_draws_fun = target_density,
    moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)

  iw4b_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,2])},
    log_prob_prop_draws_fun = prop_density,
    log_prob_target_draws_fun = target_density,
    moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)


  # true means are (5,5)
  iw_mean$expectation
  iw2_mean$expectation
  c(iw3_mean$expectation,iw3b_mean$expectation)



expect_equal(iw_mean$expectation, c(5.005028, 4.998967), tolerance = 1e-6)
expect_equal(iw2_mean$expectation, c(5.005028, 4.998967), tolerance = 1e-6)
expect_equal(c(iw3_mean$expectation,iw3b_mean$expectation), c(5.005665, 5.000551), tolerance = 1e-6)








# using ratio



iw_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  draws[,1:2]},
  log_prob_prop_draws_fun = prop_density,
  log_ratio_draws_fun = ratio_density)


iw2_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  draws[,1:2]},
  log_prob_prop_draws_fun = prop_density,
  log_ratio_draws_fun = ratio_density,
  moment_match = TRUE, k_threshold = 0.0, cov_transform = TRUE)



iw3_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  matrix(draws[,1])},
  log_prob_prop_draws_fun = prop_density,
  log_ratio_draws_fun = ratio_density,
  moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE)

iw3b_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  matrix(draws[,2])},
  log_prob_prop_draws_fun = prop_density,
  log_ratio_draws_fun = ratio_density,
  moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE)

iw4_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  matrix(draws[,1])},
  log_prob_prop_draws_fun = prop_density,
  log_ratio_draws_fun = ratio_density,
  moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)

iw4b_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  matrix(draws[,2])},
  log_prob_prop_draws_fun = prop_density,
  log_ratio_draws_fun = ratio_density,
  moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)


# true means are (5,5)
iw_mean$expectation
iw2_mean$expectation
c(iw3_mean$expectation,iw3b_mean$expectation)



expect_equal(iw_mean$expectation, c(5.005028, 4.998967), tolerance = 1e-6)
expect_equal(iw2_mean$expectation, c(5.005028, 4.998967), tolerance = 1e-6)
expect_equal(c(iw3_mean$expectation,iw3b_mean$expectation), c(5.005665, 5.000551), tolerance = 1e-6)









































  # set.seed(6)
  #
  # S <- 40000
  #
  # prop_mean <- 0
  # prop_var <- sqrt(0.1)
  #
  # target_mean <- 5
  # target_var <- 1
  #
  # prop_sample <- matrix(rnorm(2 * S, prop_mean, prop_var), S, 2)
  # prop_density <- function(draws, ...) {
  #   dnorm(draws[,1], prop_mean, prop_var, log = TRUE) + dnorm(draws[,2], prop_mean, prop_var, log = TRUE)
  # }
  #
  # target_density <- function(draws ,...) {
  #   dnorm(draws[,1], target_mean, target_var, log = TRUE) + dnorm(draws[,2], target_mean, target_var, log = TRUE)
  # }
  #
  # ratio_density <- function(draws ,...) {
  #   target_density(draws, ...) - prop_density(draws, ...)
  # }
  #
  # target_sample <- matrix(rnorm(2 * S, target_mean, target_var), S, 2)
  #
  #
  # iw_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  #   draws[,1:2]},
  #   log_prob_prop_draws_fun = prop_density,
  #   log_prob_target_draws_fun = target_density)
  #
  #
  # iw2_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  #   draws[,1:2]},
  #   log_prob_prop_draws_fun = prop_density,
  #   log_prob_target_draws_fun = target_density,
  #   moment_match = TRUE, k_threshold = 0.0, cov_transform = TRUE)
  #
  #
  #
  # iw3_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  #   matrix(draws[,1])},
  #   log_prob_prop_draws_fun = prop_density,
  #   log_prob_target_draws_fun = target_density,
  #   moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE)
  #
  # iw3b_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  #   matrix(draws[,2])},
  #   log_prob_prop_draws_fun = prop_density,
  #   log_prob_target_draws_fun = target_density,
  #   moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE)
  #
  # iw4_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  #   matrix(draws[,1])},
  #   log_prob_prop_draws_fun = prop_density,
  #   log_prob_target_draws_fun = target_density,
  #   moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)
  #
  # iw4b_mean <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  #   matrix(draws[,2])},
  #   log_prob_prop_draws_fun = prop_density,
  #   log_prob_target_draws_fun = target_density,
  #   moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)
  #
  #
  # # true means are (5,5)
  # iw_mean$expectation
  # iw2_mean$expectation
  # c(iw3_mean$expectation,iw3b_mean$expectation)
  #
  #
  #
  # iw_moment2 <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  #   draws[,1:2]^2},
  #   log_prob_prop_draws_fun = prop_density,
  #   log_prob_target_draws_fun = target_density, cov_transform = TRUE)
  #
  #
  #
  #
  # iw2_moment2 <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  #   draws[,1:2]^2},
  #   log_prob_prop_draws_fun = prop_density,
  #   log_prob_target_draws_fun = target_density,
  #   moment_match = TRUE, k_threshold = 0.0, cov_transform = TRUE)
  #
  #
  # iw3_moment2 <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  #   matrix(draws[,1]^2)},
  #   log_prob_prop_draws_fun = prop_density,
  #   log_prob_target_draws_fun = target_density,
  #   moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = F)
  #
  # iw3b_moment2 <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  #   matrix(draws[,2]^2)},
  #   log_prob_prop_draws_fun = prop_density,
  #   log_prob_target_draws_fun = target_density,
  #   moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = F)
  #
  # iw4_moment2 <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  #   matrix(draws[,1]^2)},
  #   log_prob_prop_draws_fun = prop_density,
  #   log_prob_target_draws_fun = target_density,
  #   moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = F, restart_transform = TRUE)
  #
  # iw4b_moment2 <- expectation_moment_match(prop_sample, expectation_fun =  function(draws, ...) {
  #   matrix(draws[,2]^2)},
  #   log_prob_prop_draws_fun = prop_density,
  #   log_prob_target_draws_fun = target_density,
  #   moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = F, restart_transform = TRUE)
  #
  #
  # iw_var <- iw_moment2$expectation - iw_mean$expectation^2
  # iw2_var <- iw2_moment2$expectation - iw2_mean$expectation^2
  # iw3_var <- c(iw3_moment2$expectation,iw3b_moment2$expectation) - c(iw3_mean$expectation,iw3b_mean$expectation)^2
  # iw4_var <- c(iw4_moment2$expectation,iw4b_moment2$expectation) - c(iw4_mean$expectation,iw4b_mean$expectation)^2
  #
  # iw_moment2$expectation
  # iw2_moment2$expectation
  # c(iw3_moment2$expectation,iw3b_moment2$expectation)
  # c(iw4_moment2$expectation,iw4b_moment2$expectation)
  # # true variances are (1,1)
  # iw_var
  # iw2_var
  # iw3_var
  # iw4_var







})


test_that("expectation_moment_match works for simple Monte Carlo case", {


  set.seed(6)
  S <- 4000

  target_mean <- 5
  target_var <- 1

  target_density <- function(draws ,...) {
    dnorm(draws[,1], target_mean, target_var, log = TRUE) + dnorm(draws[,2], target_mean, target_var, log = TRUE)
  }

  target_sample <- matrix(rnorm(2 * S, target_mean, target_var), S, 2)


  iw_mean <- expectation_moment_match(target_sample, expectation_fun =  function(draws, ...) {
    draws[,1:2]},
    log_prob_prop_draws_fun = target_density)


  iw2_mean <- expectation_moment_match(target_sample, expectation_fun =  function(draws, ...) {
    draws[,1:2]},
    log_prob_prop_draws_fun = target_density,
    moment_match = TRUE, k_threshold = 0.0, cov_transform = TRUE)



  iw3_mean <- expectation_moment_match(target_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,1])},
    log_prob_prop_draws_fun = target_density,
    moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE)

  iw3b_mean <- expectation_moment_match(target_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,2])},
    log_prob_prop_draws_fun = target_density,
    moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE)

  iw4_mean <- expectation_moment_match(target_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,1])},
    log_prob_prop_draws_fun = target_density,
    moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)

  iw4b_mean <- expectation_moment_match(target_sample, expectation_fun =  function(draws, ...) {
    matrix(draws[,2])},
    log_prob_prop_draws_fun = target_density,
    moment_match = TRUE, k_threshold = -10.0, split = TRUE, cov_transform = TRUE, restart_transform = TRUE)


  # true means are (5,5)
  iw_mean$expectation
  iw2_mean$expectation
  c(iw3_mean$expectation,iw3b_mean$expectation)



  expect_equal(iw_mean$expectation, c(5.005361,4.999969), tolerance = 1e-6)
  expect_equal(iw2_mean$expectation, c(5.005361, 4.999969), tolerance = 1e-6)
  expect_equal(c(iw3_mean$expectation,iw3b_mean$expectation), c(5.005361, 4.998519), tolerance = 1e-6)













})
