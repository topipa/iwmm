library(rstan)
# rstan_options(auto_write = TRUE)
options(mc.cores = 1)

stancode <- "data {
  int<lower=0> N;
  vector[N] x;
  int<lower=0> prior_only;
  real mu0;
  real<lower=0> nu0;
  real<lower=0> kappa0;
  real<lower=0> sigma0;
  }
  parameters {
    real mu;
    real<lower=0> sigma_sq;
  }
  transformed parameters {
    real<lower=0> sigma = sqrt(sigma_sq);
  }
  model {
    target += normal_lpdf(mu | mu0, sigma / sqrt(kappa0));
    target += scaled_inv_chi_square_lpdf(sigma_sq | nu0, sigma0);

    if (prior_only == 0) {
      target += normal_lpdf(x | mu, sigma);
    }

  }
  generated quantities {
    vector[N] log_lik;
    for (n in 1:N) log_lik[n] = normal_lpdf(x[n] | mu, sigma);
  }"

stanmodel <- stan_model(model_code = stancode)


# Proposal = prior, target = posterior, ratio = likelihood
# gaussian model, known sigma

# generate data
SEED <- 123
set.seed(SEED)
n <- as.integer(10)
data_var <- 1
x <- rnorm(n = n, mean = 2, sd = data_var)

mu0 <- 0
nu0 <- 10
kappa0 <- 1
sigma0 <- 1


standata_full <- list(
  N = n,
  x = x,
  prior_only = 0,
  mu0 = mu0,
  nu0 = nu0,
  kappa0 = kappa0,
  sigma0 = sigma0
)

standata_prior <- list(
  N = n,
  x = x,
  prior_only = 1,
  mu0 = mu0,
  nu0 = nu0,
  kappa0 = kappa0,
  sigma0 = sigma0
)

# fit model (just sampling from prior and calculating log_lik)
fit_full <- sampling(
  stanmodel,
  data = standata_full,
  chains = 4,
  iter = 1000,
  refresh = 0,
  seed = SEED
)

fit_prior <- sampling(
  stanmodel,
  data = standata_prior,
  chains = 4,
  iter = 1000,
  refresh = 0,
  seed = SEED
)

# analytical posterior for unknown sigma
# from BDA3 p
# mu ~ student_t(mu_n, sigma_n / kappa_n)
# sigma^2 ~ scaled_inv_chi_square(nu_n, sigma_n)

# mu_n = kappa_0 / (kappa_0 + n) * mu_0 + n / (kappa_0 + n) * ybar
# kappa_n = kappa_0 + n
# nu_n = nu_0 + n
# nu_n*sigma_sq_n = nu_0 * sigma_sq_0 + (n - 1) * s^2 +
# (kappa_0 * n) / (kappa_0 + n) * (ybar - mu_0)^2

ybar <- mean(x)
s_sq <- var(x)
nu_n <- nu0 + n
kappa_n <- kappa0 + n

mu_n <- kappa0 / kappa_n * mu0 + n / kappa_n * ybar

sigma_sq_n <- (nu0 * sigma0^2 + (n - 1) * s_sq + (kappa0 * n) /
                 kappa_n * (ybar - mu0)^2) / nu_n # styler: off

sigma_sq_post_mean <- nu_n * sigma_sq_n / (nu_n - 2)

sigma_sq_post_var <- 2 * nu_n^2 * sigma_sq_n^2 / ((nu_n - 2)^2 * (nu_n - 4))

sigma_sq_post_sd <- sqrt(sigma_sq_post_var)

mu_post_mean <- mu_n
mu_post_var <- (sigma_sq_n / kappa_n) * (nu_n / (nu_n - 2))
mu_post_sd <- sqrt(mu_post_var)

mean_analytical <- c(mu = mu_n, sigma_sq = sigma_sq_post_mean)
sd_analytical <- c(mu = mu_post_sd, sigma_sq = sqrt(sigma_sq_post_var))


test_that("moment_match.stanfit matches analytical results", {
  # TODO: implement this test with expectation_fun
  # ratio = jointlikelihood
  joint_log_lik <- function(draws, fit, ...) {
    cdraws <- constrain_draws.stanfit(fit, draws)
    ll <- posterior::merge_chains(
      posterior::subset_draws(cdraws, variable = "log_lik")
    )
    apply(ll, 2, rowSums)
  }

  iw_prior <- moment_match.stanfit(
    fit_prior,
    log_ratio_fun = joint_log_lik,
    k_threshold = -Inf # ensure moment-matching is used
  )

  draws_mm_prior <- posterior::subset_draws(
    posterior::as_draws_matrix(iw_prior$draws),
    variable = c("mu", "sigma_sq")
  )

  weights_mm_prior <- exp(iw_prior$log_weights)
  mean_mm_prior <- matrixStats::colWeightedMeans(
    draws_mm_prior,
    w = weights_mm_prior
  )

  var_weighted <- function(x, w) {
    stats::cov.wt(cbind(x), wt = w)$cov
  }

  sd_mm_prior <- apply(
    posterior::as_draws_matrix(draws_mm_prior),
    2,
    function(x) sqrt(var_weighted(x = x, w = weights_mm_prior))
  )

  expect_equal(
    mean_mm_prior[c(1, 2)],
    mean_analytical,
    tolerance = 0.1
  )

  expect_equal(
    sd_mm_prior[c(1, 2)],
    sd_analytical,
    tolerance = 0.1
  )

})


test_that("moment_match.stanfit works with expectation using constrained parameters", {
  expectation_fun_first_moment <- function(draws, ...) {
    draws
  }
  expectation_fun_second_moment <- function(draws, ...) {
    draws^2
  }

  iw_first_moment <- suppressWarnings(moment_match.stanfit(
    fit_full,
    expectation_fun = expectation_fun_first_moment,
    k_threshold = -Inf # ensure moment-matching is used
  ))

  iw_second_moment <- suppressWarnings(moment_match.stanfit(
    fit_full,
    expectation_fun = expectation_fun_second_moment,
    k_threshold = -Inf # ensure moment-matching is used
  ))


  draws_fit_full <- posterior::subset_draws(
    posterior::as_draws_matrix(fit_full),
    variable = c("mu", "sigma_sq")
  )

  fit_means <- colMeans(draws_fit_full[, c("mu", "sigma_sq")])
  fit_sds <- matrixStats::colSds(draws_fit_full[, c("mu", "sigma_sq")])

  iwmm_mean <- iw_first_moment$expectation[, 1:2]
  iwmm_sd <- sqrt(iw_second_moment$expectation[, 1:2] - iw_first_moment$expectation[, 1:2]^2)

  expect_equal(
    iwmm_mean,
    fit_means,
    tolerance = 0.1
  )
  expect_equal(
    iwmm_sd,
    fit_sds,
    tolerance = 0.1
  )
})


test_that("moment_match.stanfit works with expectation on unconstrained parameters", {
  expectation_fun_first_moment <- function(draws, ...) {
    draws
  }
  expectation_fun_second_moment <- function(draws, ...) {
    draws^2
  }

  iw_first_moment <- suppressWarnings(moment_match.stanfit(
    fit_full,
    expectation_fun = expectation_fun_first_moment,
    k_threshold = -Inf, # ensure moment-matching is used
    constrain_draws = FALSE
  ))

  iw_second_moment <- suppressWarnings(moment_match.stanfit(
    fit_full,
    expectation_fun = expectation_fun_second_moment,
    k_threshold = -Inf, # ensure moment-matching is used
    constrain_draws = FALSE
  ))


  draws_fit_full <- posterior::subset_draws(
    posterior::as_draws_matrix(fit_full),
    variable = c("mu", "sigma_sq")
  )

  draws_fit_full <- posterior::mutate_variables(draws_fit_full, log_sigma_sq = log(sigma_sq))

  fit_means <- colMeans(draws_fit_full[, c("mu", "log_sigma_sq")])
  fit_sds <- matrixStats::colSds(draws_fit_full[, c("mu", "log_sigma_sq")])

  iwmm_mean <- iw_first_moment$expectation
  names(iwmm_mean) <- c("mu", "log_sigma_sq")

  iwmm_sd <- sqrt(iw_second_moment$expectation - iw_first_moment$expectation^2)
  names(iwmm_sd) <- c("mu", "log_sigma_sq")

  expect_equal(
    iwmm_mean,
    fit_means,
    tolerance = 0.1
  )
  expect_equal(
    iwmm_sd,
    fit_sds,
    tolerance = 0.1
  )
})
