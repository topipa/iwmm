library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

# Proposal = prior, target = posterior, ratio = likelihood
# gaussian model, known sigma

stancode <- "data {
  int<lower=0> N;
  vector[N] x;
  real<lower=0> sigma;
  }
  parameters {
    real mu;
  }
  model {
    mu ~ std_normal();
  }
  generated quantities {
    real log_lik = normal_lpdf(x | mu, sigma);
  }"

stanmodel <- stan_model(model_code = stancode)

# generate data
SEED <- 24
set.seed(SEED)
n <- as.integer(10)
data_var <- 1
x <- rnorm(n = n, mean = 10, sd = data_var)
standata <- list(N = n, x = x, sigma = data_var)

# fit model (just sampling from prior and calculating log_lik)
fit <- sampling(
  stanmodel,
  data = standata,
  chains = 4,
  iter = 2000,
  refresh = 0,
  seed = SEED
)

# analytical posterior
# from BDA3 p 41-42
sigma0 <- 1
mu0 <- 1
sigma <- 1
ybar <- mean(x)
mu_post <- (1 / sigma0^2 * mu0 + n / sigma^2 * ybar) /
  (1 / sigma0^2 + n / sigma^2)
sigma_post <- sqrt(1 / (1 / sigma0^2 + n / sigma^2))

test_that("moment_match.stanfit matches analytical results", {

  # ratio = likelihood
  ratio_density <- function(draws, stanfit, ...) {
    udraws <- unconstrain_draws_stanfit(stanfit, draws)
    cdraws <- constrain_draws_stanfit(stanfit, udraws)
    posterior::extract_variable(cdraws, "log_lik")
  }

  target_density <- function(draws, ...) {
    dnorm(draws, mu_post, sigma_post, log = TRUE)
  }

    iw_ratio <- moment_match.stanfit(fit, log_ratio_fun = ratio_density)
    draws_mm_ratio <- iw_ratio$draws
    weights_mm_ratio <- exp(iw_ratio$log_weights)
    mu_ratio <- weighted.mean(draws_mm_ratio, w = weights_mm_ratio)

    expect_equal(round(mu_ratio), round(mu_post))
    
    iw_target <- moment_match.stanfit(fit, log_prob_target_fun = target_density)
    draws_mm_target <- iw_target$draws
    weights_mm_target <- exp(iw_target$log_weights)
    mu_target <- weighted.mean(draws_mm_target, w = weights_mm_target)

    expect_equal(round(mu_target), round(mu_post))
})
