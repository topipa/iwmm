##' Return log probability of posterior
##'
##' @param x model fit object
##' @param ... arguments passed to methods
##' @return TODO
##' @export
log_prob_draws <- function(x, ...) {
  UseMethod("log_prob_draws")
}
