##' Example Stan model for importance weighted moment matching
##'
##' Provides example models (with data) that are ready for use with
##' IWMM.
##' @param model Character specifying which model code to
##'   return. Currently "normal_model" is implemented.
##' @return List containing model code and corresponding data.
##' @export
example_iwmm_model <- function(model = "normal_model") {
  examples <- iwmm_examples()

  return(list(
    model_code = examples[[model]][["model_code"]],
    data = examples[[model]][["data"]]
  ))
}

iwmm_examples <- function() {
  list(
    normal_model =
      list(
        model_code = "data {
  int<lower=0> N;
  vector[N] x;
  }
  parameters {
    real mu;
    real log_sigma;
  }
  transformed parameters {
    real<lower=0> sigma = exp(log_sigma);
  }
  model {
    target += normal_lpdf(x | mu, sigma);
  }
  generated quantities {
  vector[N] log_lik;
  for (n in 1:N) log_lik[n] =  normal_lpdf(x[n] | mu, sigma);
}",
        data = list(
          N = as.integer(10),
          x = c(1.4395244, 1.7698225, 3.5587083, 2.0705084, 2.1292877, 2.4609162, 0.7349388, 1.3131471, 1.5543380, 23.7150650)
        )
      )
  )
}
