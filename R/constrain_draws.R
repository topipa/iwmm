##' Constrain all draws from a fitted model
##'
##' @param x model fit object
##' @param ... arguments passed to methods
##' @return constrained draws
##' @export
constrain_draws <- function(x, ...) {
  UseMethod("constrain_draws")
}
