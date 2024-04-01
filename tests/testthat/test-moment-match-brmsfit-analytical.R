brms_available <- require(brms)

if (brms_available) {
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


  standata_nine_obs <- list(
    N = n - 1,
    x = x[-1],
    mu0 = mu0,
    nu0 = nu0,
    kappa0 = kappa0,
    sigma0 = sigma0
  )

  standata_single_obs <- list(
    N = 1,
    x = x[1],
    mu0 = mu0,
    nu0 = nu0,
    kappa0 = kappa0,
    sigma0 = sigma0
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

  mean_analytical <- c(b_Intercept = mu_n, sigma_sq = sigma_sq_post_mean)
  sd_analytical <- c(b_Intercept = mu_post_sd, sigma_sq = sqrt(sigma_sq_post_var))




  bprior <- c(
    prior("", class = "sigma"),
    prior("", class = "Intercept"),
    set_prior("target += normal_lpdf(Intercept | mu0, sigma / sqrt(kappa0))", check = FALSE),
    set_prior("target += scaled_inv_chi_square_lpdf(sigma_sq | nu0, sigma0)", check = FALSE)
  )




  stanvars <- c(
    stanvar(mu0,
      scode = "real mu0;",
      block = "data"
    ),
    stanvar(nu0,
      scode = "real<lower=0> nu0;",
      block = "data"
    ),
    stanvar(kappa0,
      scode = "real<lower=0> kappa0;",
      block = "data"
    ),
    stanvar(sigma0,
      scode = "real<lower=0> sigma0;",
      block = "data"
    ),
    stanvar(
      scode = "real<lower=0> sigma_sq = square(sigma);",
      block = "tparameters"
    )
  )


  fit_single_obs <- brm(x ~ 1,
    data = standata_single_obs,
    prior = bprior,
    stanvars = stanvars,
    save_pars = save_pars(all = TRUE),
    chains = 4,
    iter = 1000,
    refresh = 0,
    seed = SEED
  )



  test_that("moment_match.brmsfit matches analytical results", {
    # TODO: implement this test with expectation_fun

    joint_log_lik_extra_data <- function(draws, fit, extra_data, ...) {
      fit <- brms:::.update_pars(x = fit, upars = draws)
      ll <- log_lik(fit, newdata = extra_data)
      rowSums(ll)
    }

    iw_single_obs <- moment_match.brmsfit(
      fit_single_obs,
      log_ratio_fun = joint_log_lik_extra_data,
      extra_data = standata_nine_obs,
      k_threshold = -Inf # ensure moment-matching is used
    )

    draws_mm_single_obs <- posterior::subset_draws(
      posterior::as_draws_matrix(iw_single_obs$draws),
      variable = c("b_Intercept", "sigma_sq")
    )

    weights_mm_single_obs <- exp(iw_single_obs$log_weights)
    mean_mm_single_obs <- matrixStats::colWeightedMeans(
      draws_mm_single_obs,
      w = weights_mm_single_obs
    )

    var_weighted <- function(x, w) {
      stats::cov.wt(cbind(x), wt = w)$cov
    }

    sd_mm_single_obs <- apply(
      posterior::as_draws_matrix(draws_mm_single_obs),
      2,
      function(x) sqrt(var_weighted(x = x, w = weights_mm_single_obs))
    )

    expect_equal(
      mean_mm_single_obs[c(1, 2)],
      mean_analytical,
      tolerance = 0.1
    )

    expect_equal(
      sd_mm_single_obs[c(1, 2)],
      sd_analytical,
      tolerance = 0.1
    )
  })

  test_that("moment_match.brmsfit works with target_observation_weights", {
    normal_model <- example_iwmm_model("normal_model")

    bprior <- c(
      prior("", class = "sigma"),
      prior("", class = "Intercept")
    )

    fit <- brm(x ~ 1,
                          data = normal_model$data,
                          prior=bprior,
                          save_pars = save_pars(all = TRUE),
                          chains = 4,
                          iter = 1000,
                          refresh = 0,
                          seed = 1234
    )

    expectation_fun_first_moment <- function(draws, ...) {
      draws
    }

    first_moment <- moment_match(
      fit,
      expectation_fun = expectation_fun_first_moment,
    )

    expect_equal(
      first_moment$expectation,
      posterior_summary(fit)[,1],
      tolerance = 0.001
    )
    #TODO: why lprior and lp__ are zeros?
    expect_equal(
      posterior_summary(fit, variable = c("b_Intercept", "sigma", "Intercept")),
      posterior_summary(first_moment$fit, variable = c("b_Intercept", "sigma", "Intercept")),
      tolerance = 0.001
    )

    # leave-one-out expectation

    first_moment_loo <- moment_match(
      fit,
      target_observation_weights=append(rep(1,9), 0),
      expectation_fun = expectation_fun_first_moment,
    )

    loo_data <- normal_model$data
    loo_data$x <- loo_data$x[-loo_data$N]
    loo_data$N <- loo_data$N - 1

    fit_loo <- brm(x ~ 1,
                   data = loo_data,
                   prior=bprior,
                   save_pars = save_pars(all = TRUE),
                   chains = 4,
                   iter = 1000,
                   refresh = 0,
                   seed = 1234
    )

    #TODO: why lprior and lp__ are zeros?
    expect_equal(
      first_moment_loo$expectation,
      posterior_summary(fit_loo)[,1],
      tolerance = 0.001
    )

    expect_equal(
      first_moment_loo$expectation[c("b_Intercept", "sigma", "Intercept")],
      posterior_summary(fit_loo, variable = c("b_Intercept", "sigma", "Intercept"))[,1],
      tolerance = 0.1
    )
    #TODO: why lprior and lp__ are zeros?
    #TODO: check other estimares than mean
    expect_equal(
      posterior_summary(fit_loo, variable = c("b_Intercept", "sigma", "Intercept"))[,1],
      posterior_summary(first_moment_loo$fit, variable = c("b_Intercept", "sigma", "Intercept"))[,1],
      tolerance = 0.1
    )

  })

}
