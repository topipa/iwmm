

#' Generic importance weighted moment matching algorithm for `stanfit` objects.
#' See additional arguments from `moment_match.matrix`
#'
#' @param x A fitted `stanfit` object.
#' @param log_prob_target_fun Log density of the target.
#' The function takes argument `draws`, which are the unconstrained parameters.
#' @param log_ratio_fun Log of the density ratio (target/proposal).
#' The function takes argument `draws`, which are the unconstrained parameters.
#' @param constrain Logical specifying whether the returned draws be returned on the constrained space? Default is TRUE.
#' @param ... Further arguments passed to `moment_match.matrix`.
#'
#' @return Returns a list with 3 elements: transformed draws, updated
#' importance weights, and the pareto k diagnostic value.
#'
#' @export
moment_match.stanfit <- function(x,
                                 log_prob_target_fun = NULL,
                                 log_ratio_fun = NULL,
                                 constrain = TRUE,
                                 ...) {

  draws <- posterior::as_draws_matrix(x)
  # transform the model parameters to unconstrained space
  udraws <- unconstrain_draws_stanfit(x, draws = draws, ...)
  
  out <- moment_match.matrix(
    udraws,
    log_prob_prop_fun = log_prob_upars.stanfit,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    stanfit = x,
    ...
  )
  
  if (constrain) {
    out$draws <- constrain_draws_stanfit(x, out$draws, ...)

  }
  
  out
}


log_prob_upars.stanfit <- function(draws, stanfit, ...) {
  apply(draws, 1, rstan::log_prob,
        object = stanfit,
        adjust_transform = TRUE, gradient = FALSE
  )
}


unconstrain_draws_stanfit <- function(x, draws, ...) {
  skeleton <- .create_skeleton(x@sim$pars_oi, x@par_dims[x@sim$pars_oi])
  udraws <- apply(posterior::as_draws_matrix(draws), 1, FUN = function(theta) {
    rstan::unconstrain_pars(x, pars = .rstan_relist(theta, skeleton))
  })
  # for one parameter models
  if (is.null(dim(udraws))) {
    dim(udraws) <- c(1, length(udraws))
  }
  t(udraws)
}

#' @export
constrain_draws_stanfit <- function(x, udraws, ...) {

  # list with one element per posterior draw
  draws <- apply(udraws, 1, rstan::constrain_pars, object = x)
  varnames <- rep(names(draws[[1]]), lengths(draws[[1]]))
  # transform samples
  nsamples <- length(draws)
  draws <- unlist(draws)
  nvars <- length(draws) / nsamples
  dim(draws) <- c(nvars, nsamples)
  rownames(draws) <- varnames
  # lp__ is not computed automatically
#  lp__ <- log_prob_upars.stanfit(x, upars = upars, ...)
 # pars <- rbind(pars, lp__ = lp__)
  # bring samples into the right structure

  new_samples <- named_list(x@sim$fnames_oi[-length(x@sim$fnames_oi)], list(numeric(nsamples)))
  
  new_varnames <- sub("\\[.+", "", names(new_samples))
  new_varnames_unique <- unique(new_varnames)
  for (v in new_varnames_unique) {
    sub_vars <- draws[rownames(draws) == v, , drop = FALSE]
    sel <- which(new_varnames == v)
    for (i in seq_along(sel)) {
      new_samples[[sel[i]]] <- sub_vars[i, ]
    }
  }

  new_samples <- posterior::as_draws_array(new_samples)
    
  new_samples
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

