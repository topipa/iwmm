##' Constrain all draws from a fitted model
##'
##' @param x model fit object
##' @param ... arguments passed to methods
##' @return constrained draws
##' @export
constrain_draws <- function(x, ...) {
  UseMethod("constrain_draws")
}

##' @export
constrain_draws.CmdStanFit <- function(x, udraws, ...) {
  # list with one element per posterior draw
  draws <- apply(udraws, 1, x$constrain_variables)
  varnames <- posterior::variables(x$draws())[-1]
  # transform draws
  ndraws <- length(draws)
  draws <- unlist(draws)
  nvars <- length(draws) / ndraws
  dim(draws) <- c(nvars, ndraws)
  rownames(draws) <- varnames
  # lp__ is not computed automatically
  lp__ <- log_prob_draws.CmdStanFit(x, draws = udraws, ...)
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



##' @export
constrain_draws.stanfit <- function(x, udraws, ...) {
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
  lp__ <- log_prob_draws.stanfit(x, draws = udraws, ...)
  draws <- rbind(draws, lp__ = lp__)

  # bring draws into the right structure
  new_draws <- named_list(
    x@sim$fnames_oi,
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
