##' Unconstrain all draws from a fitted model
##'
##' @param x model fit object
##' @param ... arguments passed to methods
##' @return unconstrained draws
##' @export
unconstrain_draws <- function(x, ...) {
  UseMethod("unconstrain_draws")
}
