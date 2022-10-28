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
moment_match_CmdStanFit <- function(x,
                                    log_prob_target_fun = NULL,
                                    log_ratio_fun = NULL,
                                    constrain = TRUE,
                                    ...) {

  var_names <- names(x$constrain_pars(skeleton_only = TRUE))
  
  draws <- posterior::as_draws_matrix(x)
  draws <- posterior::subset_draws(draws, variable = var_names)
  # transform the model parameters to unconstrained space
  udraws <- unconstrain_draws_CmdStanFit(x, draws = draws, ...)
  
  out <- moment_match.matrix(
    udraws,
    log_prob_prop_fun = log_prob_upars_cmdstan,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    cmdfit = x,
    ...
  )

  if (constrain) {
    out$draws <- constrain_draws_CmdStanFit(x, udraws = out$draws, ...)
  }
  
  out
}

log_prob_upars_cmdstan <- function(draws, cmdfit, ...) {
  apply(draws, 1, cmdfit$log_prob)
}

unconstrain_draws_CmdStanFit <- function(x, draws, ...) {
  skeleton <- x$constrain_pars(skeleton_only = TRUE)
  udraws <- apply(draws, 1, FUN = function(theta) {
    x$unconstrain_pars(pars = .cmdstan_relist(theta, skeleton))
  })
  # for one parameter models
  if (is.null(dim(udraws))) {
    dim(udraws) <- c(1, length(udraws))
  }
  t(udraws)
}


.cmdstan_relist <- function(x, skeleton) {
  out <- utils::relist(x, skeleton)
  for (i in seq_along(skeleton)) {
    dim(out[[i]]) <- dim(skeleton[[i]])
  }
  out
  }

#' @export
constrain_draws_CmdStanFit <- function(x, udraws, ...) {

  # list with one element per posterior draw
  draws <- apply(udraws, 1, x$constrain_pars)
  varnames <- posterior::variables(x$draws())[-1]
  # transform samples
  nsamples <- length(draws)
  draws <- unlist(draws)
  nvars <- length(draws) / nsamples
  dim(draws) <- c(nvars, nsamples)
  rownames(draws) <- varnames
  # lp__ is not computed automatically
#  lp__ <- log_prob_upars.CmdStanFit(x, upars = upars, ...)
#  pars <- rbind(pars, lp__ = lp__)

  # bring samples into the right structure
  new_samples <- named_list(varnames)
  new_varnames <- sub("\\[.+", "", names(new_samples))
  new_varnames <- names(new_samples)
  new_varnames_unique <- unique(new_varnames)
  for (v in new_varnames_unique) {
    sub_vars <- draws[rownames(draws) == v, , drop = FALSE]
    sel <- which(new_varnames == v)
    for (i in seq_along(sel)) {
      new_samples[[sel[i]]] <- sub_vars[i, ]
    }
  }

  posterior::as_draws_array(new_samples)
}



# initialize a named list
# @param names names of the elements
# @param values optional values of the elements
named_list <- function(names, values = NULL) {
  if (!is.null(values)) {
    if (length(values) <= 1L) {
      values <- replicate(length(names), values)
    }
    values <- as.list(values)
    stopifnot(length(values) == length(names))
  } else {
    values <- vector("list", length(names))
  }
  stats::setNames(values, names)
}
