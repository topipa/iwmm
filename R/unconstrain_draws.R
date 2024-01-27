##' Unconstrain all draws from a fitted model
##'
##' @param x model fit object
##' @param ... arguments passed to methods
##' @return unconstrained draws
##' @export
unconstrain_draws <- function(x, ...) {
  UseMethod("unconstrain_draws")
}

unconstrain_draws.stanfit <- function(x, draws, ...) {
  skeleton <- .create_skeleton(x@sim$pars_oi, x@par_dims[x@sim$pars_oi])
  udraws <- apply(posterior::as_draws_matrix(draws), 1, FUN = function(draw) {
    rstan::unconstrain_pars(x, pars = .relist(draw, skeleton))
  })
  # for one parameter models
  if (is.null(dim(udraws))) {
    dim(udraws) <- c(1, length(udraws))
  }
  t(udraws)
}


# unconstrain_draws.CmdStanFit <- function(x, draws, ...) {
#   if (!is.matrix(draws)) {
#     draws <- posterior::as_draws_matrix(x)
#   }
#
#   skeleton <- x$variable_skeleton()
#
#   # handle lp__ in draws object correctly
#   if ("lp__" %in% posterior::variables(draws)) {
#     skeleton <- c(list("lp__" = 0), skeleton)
#   }
#
#   udraws <- apply(draws, 1, FUN = function(draw) {
#     x$unconstrain_variables(.relist(draw, skeleton))
#   })
#   # for one parameter models
#   if (is.null(dim(udraws))) {
#     dim(udraws) <- c(1, length(udraws))
#   }
#   t(udraws)
# }
