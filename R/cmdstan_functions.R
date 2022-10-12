#' Generic importance weighted moment matching algorithm for `CmdStanFit` objects.
#' See additional arguments from `moment_match.matrix`
#'
#' @param x A fitted `CmdStanFit` object.
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
moment_match.CmdStanFit <- function(x,
                                    log_prob_target_fun = NULL,
                                    log_ratio_fun = NULL,
                                    ...) {

  par_names <- names(.create_skeleton_cmdstan(x$runset$args$model_variables))
  
  pars <- posterior::as_draws_matrix(x)
  pars <- posterior::subset_draws(pars, variable = par_names)
  # transform the model parameters to unconstrained space
  upars <- unconstrain_pars.CmdStanFit(x, pars = pars, ...)
  
  out <- moment_match.matrix(
    upars,
    log_prob_prop_fun = log_prob_upars.CmdStanFit,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    cmdfit = x,
    ...
  )
  out
}

log_prob_upars.CmdStanFit <- function(draws, cmdfit, ...) {
  apply(draws, 1, cmdfit$log_prob)
}

unconstrain_pars.CmdStanFit <- function(x, pars, ...) {
  skeleton <- .create_skeleton_cmdstan(x$runset$args$model_variables)
  upars <- apply(pars, 1, FUN = function(theta) {
    x$unconstrain_pars(pars = .cmdstan_relist(theta, skeleton))
  })
  # for one parameter models
  if (is.null(dim(upars))) {
    dim(upars) <- c(1, length(upars))
  }
  t(upars)
}


# cmdstan helper function to get dims of parameters right
.create_skeleton_cmdstan <- function(model_variables) {
    model_pars <- model_variables$parameters
    skeleton <- lapply(model_pars, function(par) {
        dims <- par$dimensions
        dims <- ifelse(dims == 0, 1, dims)
        array(0, dim = dims)
    })
    stats::setNames(skeleton, names(model_pars))
    }

.cmdstan_relist <- function(x, skeleton) {
  out <- utils::relist(x, skeleton)
  for (i in seq_along(skeleton)) {
    dim(out[[i]]) <- dim(skeleton[[i]])
  }
  out
  }
