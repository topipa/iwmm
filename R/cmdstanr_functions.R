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
moment_match.CmdStanFit <- function(x,
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
    x = udraws,
    log_prob_prop_fun = log_prob_draws_CmdStanFit,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    fit = x,
    ...
  )

  if (constrain) {
    out$draws <- constrain_draws.CmdStanFit(x, udraws = out$draws, ...)
  }

  out
}


log_prob_draws_CmdStanFit <- function(fit, draws, ...) {
  apply(
    draws,
    1,
    fit$log_prob,
    jacobian_adjustment = TRUE
  )
}

unconstrain_draws_CmdStanFit <- function(x, draws, ...) {

  if (!is.matrix(draws)) {
    draws <- posterior::as_draws_matrix(x)
  }

  skeleton <- x$constrain_pars(skeleton_only = TRUE)

  # handle lp__ in draws object correctly
  if ("lp__" %in% posterior::variables(draws)) {
    skeleton <- c(list("lp__" = 0), skeleton)
  }

  udraws <- apply(draws, 1, FUN = function(draw) {
    x$unconstrain_pars(pars = .relist(draw, skeleton))
  })
  # for one parameter models
  if (is.null(dim(udraws))) {
    dim(udraws) <- c(1, length(udraws))
  }
  t(udraws)
}

constrain_draws.CmdStanFit <- function(x, udraws, ...) {

  # list with one element per posterior draw
  draws <- apply(udraws, 1, x$constrain_pars)
  varnames <- posterior::variables(x$draws())[-1]
  # transform draws
  ndraws <- length(draws)
  draws <- unlist(draws)
  nvars <- length(draws) / ndraws
  dim(draws) <- c(nvars, ndraws)
  rownames(draws) <- varnames
  # lp__ is not computed automatically
  lp__ <- log_prob_draws_CmdStanFit(x, draws = udraws, ...)
  draws <- rbind(draws, lp__ = lp__)

  # bring draws into the right structure
  new_draws <- named_list(c("lp__", varnames))
  new_varnames <- sub("\\[.+", "", names(new_draws))
  new_varnames <- names(new_draws)
  new_varnames_unique <- unique(new_varnames)
  for (v in new_varnames_unique) {
    sub_vars <- draws[rownames(draws) == v, , drop = FALSE]
    sel <- which(new_varnames == v)
    for (i in seq_along(sel)) {
      new_draws[[sel[i]]] <- sub_vars[i, ]
    }
  }

  posterior::as_draws_array(new_draws)
}
