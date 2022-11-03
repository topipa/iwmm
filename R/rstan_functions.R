#' Generic importance weighted moment matching algorithm for `stanfit` objects.
#' See additional arguments from `moment_match.matrix`
#'
#' @param x A fitted `stanfit` object.
#' @param log_prob_target_fun Log density of the target.  The function
#'   takes argument `draws`, which are the unconstrained draws.
#' @param log_ratio_fun Log of the density ratio (target/proposal).
#'   The function takes argument `draws`, which are the unconstrained
#'   draws.
#' @param constrain Logical specifying whether to return draws on the
#'   constrained space? Default is TRUE.
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

  # ensure draws are in matrix form
  draws <- posterior::as_draws_matrix(x)

  # transform the draws to unconstrained space
  udraws <- unconstrain_draws_stanfit(x, draws = draws, ...)

  out <- moment_match.matrix(
    udraws,
    log_prob_prop_fun = log_prob_draws_stanfit,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    fit = x,
    ...
  )

  if (constrain) {
    out$draws <- constrain_draws_stanfit(x, out$draws, ...)
  }

  out
}



log_prob_draws_stanfit <- function(fit, draws, ...) {
  apply(
    draws,
    1,
    rstan::log_prob,
    object = fit,
    adjust_transform = TRUE,
    gradient = FALSE
  )
}

unconstrain_draws_stanfit <- function(x, draws, ...) {
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

constrain_draws_stanfit <- function(x, udraws, ...) {

  # list with one element per posterior draw
  draws <- apply(udraws, 1, rstan::constrain_pars, object = x)
  varnames <- rep(names(draws[[1]]), lengths(draws[[1]]))
  # transform draws
  ndraws <- length(draws)
  draws <- unlist(draws)
  nvars <- length(draws) / ndraws
  dim(draws) <- c(nvars, ndraws)
  rownames(draws) <- varnames
  # lp__ is not computed automatically
  lp__ <- log_prob_draws_stanfit(x, draws = udraws, ...)
  draws <- rbind(draws, lp__ = lp__)

  # bring draws into the right structure
  new_draws <- named_list(
    x@sim$fnames_oi[-length(x@sim$fnames_oi)],
    list(numeric(ndraws))
  )
  new_varnames <- sub("\\[.+", "", names(new_draws))
  new_varnames_unique <- unique(new_varnames)
  for (v in new_varnames_unique) {
    sub_vars <- draws[rownames(draws) == v, , drop = FALSE]
    sel <- which(new_varnames == v)
    for (i in seq_along(sel)) {
      new_draws[[sel[i]]] <- sub_vars[i, ]
    }
  }

  new_draws <- posterior::as_draws_array(new_draws)

  new_draws
}
