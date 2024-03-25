##' Unconstrain all draws from a fitted model
##'
##' @param x model fit object
##' @param ... arguments passed to methods
##' @return unconstrained draws
##' @export
unconstrain_draws <- function(x, ...) {
  UseMethod("unconstrain_draws")
}

##' @export
unconstrain_draws.stanfit <- function(x, draws, ...) {
  skeleton <- .create_skeleton(x@sim$pars_oi, x@par_dims[x@sim$pars_oi])
  udraws <- apply(posterior::as_draws_matrix(draws), 1, FUN = function(draw) {
    rstan::unconstrain_pars(x, pars = .relist(draw, skeleton))
  })
  # for one parameter models
  if (is.null(dim(udraws))) {
    dim(udraws) <- c(1, length(udraws))
  }
  out <- posterior::as_draws_matrix(t(udraws))

  out
}
