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
moment_match.CmdStanFit <- function(x,
                                    log_prob_target_fun = NULL,
                                    log_ratio_fun = NULL,
                                    constrain_pars = TRUE,
                                    ...) {

  par_names <- names(.create_skeleton_cmdstan(x$runset$args$model_variables))
  
  pars <- posterior::as_draws_matrix(x)
  pars <- posterior::subset_draws(pars, variable = par_names)
  # transform the model parameters to unconstrained space
  upars <- unconstrain_pars.CmdStanFit(x, pars = pars, ...)
  
  out <- moment_match.matrix(
    upars,
    log_prob_prop_fun = log_prob_upars_cmdstan,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    cmdfit = x,
    ...
  )

  if (constrain_pars) {
    out$draws <- constrain_all_pars(x, out$draws, ...)

  }
  
  out
}

log_prob_upars_cmdstan <- function(draws, cmdfit, ...) {
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

#' @export
constrain_all_pars <- function(x, ...) {
  UseMethod("constrain_all_pars")
}

#' @export
constrain_all_pars.CmdStanFit <- function(x, upars, ...) {

  # list with one element per posterior draw
  pars <- apply(upars, 1, x$constrain_pars)
  parnames <- rep(names(pars[[1]]), lengths(pars[[1]]))
  # transform samples
  nsamples <- length(pars)
  pars <- unlist(pars)
  npars <- length(pars) / nsamples
  dim(pars) <- c(npars, nsamples)
  rownames(pars) <- parnames
  # lp__ is not computed automatically
#  lp__ <- log_prob_upars.CmdStanFit(x, upars = upars, ...)
#  pars <- rbind(pars, lp__ = lp__)
  # bring samples into the right structure
  new_samples <- named_list(names(x$runset$args$model_variables$parameters), list(numeric(nsamples)))
  new_parnames <- sub("\\[.+", "", names(new_samples))
  new_parnames_unique <- unique(new_parnames)
  for (p in new_parnames_unique) {
    sub_pars <- pars[rownames(pars) == p, , drop = FALSE]
    sel <- which(new_parnames == p)
    for (i in seq_along(sel)) {
      new_samples[[sel[i]]] <- sub_pars[i, ]
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
