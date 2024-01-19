##' Return log probability of posterior
##'
##' @param x model fit object
##' @param ... arguments passed to methods
##' @return TODO
log_prob_draws <- function(x, ...) {
  UseMethod("log_prob_draws")
}

##' @export
log_prob_draws.CmdStanFit <- function(fit, draws, ...) {
  apply(
    draws,
    1,
    fit$log_prob,
    jacobian = TRUE
  )
}

##' @export
log_prob_draws.stanfit <- function(fit, draws, ...) {
  apply(
    draws,
    1,
    rstan::log_prob,
    object = fit,
    adjust_transform = TRUE,
    gradient = FALSE
  )
}
