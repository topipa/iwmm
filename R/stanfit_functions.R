

#' Generic importance weighted moment matching algorithm for `stanfit` objects.
#' See additional arguments from `moment_match.matrix`
#'
#' @param x A fitted `stanfit` object.
#' @param log_prob_target_fun Log density of the target.
#' The function takes argument `draws`, which are the unconstrained parameters.
#' @param log_ratio_fun Log of the density ratio (target/proposal).
#' The function takes argument `draws`, which are the unconstrained parameters.
#' @param ... Further arguments passed to `moment_match.matrix`.
#'
#' @return Returns a list with 3 elements: transformed draws, updated
#' importance weights, and the pareto k diagnostic value.
#'
#' @export
moment_match.stanfit <- function(x,
                                 log_prob_target_fun = NULL,
                                 log_ratio_fun = NULL,
                                 ...) {

  pars <- as.matrix(x)
  # transform the model parameters to unconstrained space
  upars <- unconstrain_pars.stanfit(x, pars = pars, ...)

  out <- moment_match.matrix(
    upars,
    log_prob_prop_fun = log_prob_upars.stanfit,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    stanfit = x,
    ...
  )
  out
}


log_prob_upars.stanfit <- function(draws, stanfit, ...) {
  apply(draws, 1, rstan::log_prob,
        object = stanfit,
        adjust_transform = TRUE, gradient = FALSE
  )
}


unconstrain_pars.stanfit <- function(stanfit, pars, ...) {
  skeleton <- .create_skeleton(stanfit@sim$pars_oi, stanfit@par_dims[stanfit@sim$pars_oi])
  upars <- apply(pars, 1, FUN = function(theta) {
    rstan::unconstrain_pars(stanfit, pars = .rstan_relist(theta, skeleton))
  })
  # for one parameter models
  if (is.null(dim(upars))) {
    dim(upars) <- c(1, length(upars))
  }
  t(upars)
}

# -------- will be imported from rstan at some point -------
# create a named list of draws for use with rstan methods
.rstan_relist <- function(stanfit, skeleton) {
  out <- utils::relist(stanfit, skeleton)
  for (i in seq_along(skeleton)) {
    dim(out[[i]]) <- dim(skeleton[[i]])
  }
  out
}

# rstan helper function to get dims of parameters right
.create_skeleton <- function(pars, dims) {
  out <- lapply(seq_along(pars), function(i) {
    len_dims <- length(dims[[i]])
    if (len_dims < 1) {
      return(0)
    }
    return(array(0, dim = dims[[i]]))
  })
  names(out) <- pars
  out
}
# -------- will be imported from rstan at some point -------

