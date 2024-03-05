#' Generic importance weighted moment matching algorithm for `brmsfit` objects.
#' See additional arguments from `moment_match.matrix`
#'
#' @param x A fitted `brmsfit` object.
#' @param log_prob_target_fun Log density of the target.  The function
#'   takes argument `draws`, which are the unconstrained draws.
#'   Can also take the argument `fit` which is the stan model fit.
#' @param log_ratio_fun Log of the density ratio (target/proposal).
#'   The function takes argument `draws`, which are the unconstrained
#'   draws. Can also take the argument `fit` which is the stan model fit.
#' @param target_observation_weights A vector of weights for observations for
#' defining the target distribution. A value 0 means dropping the observation,
#' a value 1 means including the observation similarly as in the current data,
#' and a value 2 means including the observation twice.
#' @param expectation_fun Optional argument, NULL by default. A
#'   function whose expectation is being computed. The function takes
#'   arguments `draws`.
#' @param log_expectation_fun Logical indicating whether the
#'   expectation_fun returns its values as logarithms or not. Defaults
#'   to FALSE. If set to TRUE, the expectation function must be
#'   nonnegative (before taking the logarithm).  Ignored if
#'   `expectation_fun` is NULL.
#' @param constrain Logical specifying whether to return draws on the
#'   constrained space? Default is TRUE.
#' @param ... Further arguments passed to `moment_match.matrix`.
#'
#' @return Returns a list with 3 elements: transformed draws, updated
#'   importance weights, and the pareto k diagnostic value. If expectation_fun
#'   is given, also returns the expectation.
#'
#' @export
moment_match.brmsfit <- function(x,
                                 log_prob_target_fun = NULL,
                                 log_ratio_fun = NULL,
                                 target_observation_weights = NULL,
                                 expectation_fun = NULL,
                                 log_expectation_fun = FALSE,
                                 constrain = TRUE,
                                 ...) {
  if (!is.null(target_observation_weights) && (!is.null(log_prob_target_fun) || !is.null(log_ratio_fun))) {
    stop("You must give only one of target_observation_weights, log_prob_target_fun, or log_ratio_fun.")
  }

  # ensure draws are in matrix form
  draws <- posterior::as_draws_matrix(x)
  # draws <- as.matrix(draws)

  if (!is.null(target_observation_weights)) {
    out <- tryCatch(brms::log_lik(x),
                    error = function(cond) {
                      message(cond)
                      message("\nYour brmsfit does not include a parameter called log_lik.")
                      message("This should not happen. Perhaps you are using an unsupported observation model?")
                      return(NA)
                    }
    )

    function(draws, fit, extra_data, ...) {
      fit <- .update_pars(x = fit, upars = draws)
      ll <- brms::log_lik(fit, newdata = extra_data)
      rowSums(ll)
    }

    log_ratio_fun <- function(draws, fit, ...) {
      fit <- .update_pars(x = fit, upars = draws)
      ll <- brms::log_lik(fit)
      colSums(t(drop(ll)) * (target_observation_weights - 1))
    }
  }


  # transform the draws to unconstrained space
  udraws <- unconstrain_draws.brmsfit(x, draws = draws, ...)

  out <- moment_match.matrix(
    udraws,
    log_prob_prop_fun = log_prob_draws.brmsfit,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    expectation_fun = expectation_fun,
    log_expectation_fun = log_expectation_fun,
    fit = x,
    ...
  )

  x <- .update_pars(x = x, upars = out$draws)

   if (constrain) {
     out$draws <- posterior::as_draws(x)
   }

  out$fit <- x

  out
}


#' @export
constrain_draws.brmsfit <- function(x, draws, ...) {
  # list with one element per posterior draw
  x <- x$fit
  udraws <- draws
  draws <- apply(draws, 1, rstan::constrain_pars, object = x)
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
    x@sim$fnames_oi_old,
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

  posterior::as_draws_array(new_draws)

}

#' @export
log_prob_draws.brmsfit <- function(fit, draws, ...) {
  # x <- update_misc_env(x, only_windows = TRUE)
  log_prob_draws.stanfit(fit$fit, draws = draws, ...)
}

#' @export
unconstrain_draws.brmsfit <- function(x, draws, ...) {
  unconstrain_draws.stanfit(x$fit, draws = draws, ...)
}


# the following functions are copied from brms

# wrapper around rstan::constrain_pars
# ensures that the right posterior draws are excluded
.constrain_pars <- function(upars, x) {
  out <- rstan::constrain_pars(upars, object = x$fit)
  out[x$exclude] <- NULL
  out
  }

.update_pars <- function(x, upars, ...) {
  # list with one element per posterior draw
  pars <- apply(upars, 1, .constrain_pars, x = x)
  # select required parameters only
  pars <- lapply(pars, "[", x$fit@sim$pars_oi_old)
  # transform draws
  ndraws <- length(pars)
  pars <- unlist(pars)
  npars <- length(pars) / ndraws
  dim(pars) <- c(npars, ndraws)
  # add dummy 'lp__' draws
  pars <- rbind(pars, rep(0, ndraws))
  # bring draws into the right structure
  new_draws <- named_list(x$fit@sim$fnames_oi_old, list(numeric(ndraws)))
  if (length(new_draws) != nrow(pars)) {
    stop("Updating parameters in 'brmsfit' failed. ")
  }
  for (i in seq_len(npars)) {
    new_draws[[i]] <- pars[i, ]
  }

  # create dummy sampler_params for new sim object
  newsamples <- list(new_draws)
  attr(newsamples[[1]], "sampler_params") <- .create_dummy_sampler_params(x)

  # create new sim object to overwrite x$fit@sim
  x$fit@sim <- list(
    samples = newsamples,
    iter = ndraws,
    thin = 1,
    warmup = 0,
    chains = 1,
    n_save = ndraws,
    warmup2 = 0,
    permutation = list(seq_len(ndraws)),
    pars_oi = x$fit@sim$pars_oi_old,
    dims_oi = x$fit@sim$dims_oi_old,
    fnames_oi = x$fit@sim$fnames_oi_old,
    n_flatnames = length(x$fit@sim$fnames_oi_old)
  )
  x$fit@stan_args <- list(
    list(chain_id = 1, iter = ndraws, thin = 1, warmup = 0)
  )

  brms::rename_pars(x)
}

.create_dummy_sampler_params <- function(x) {

  params <- attr(x$fit@sim$samples[[1]], "sampler_params")
  newparams <- params
  for (i in seq_along(params)) {
    newparams[[i]] <- numeric()
  }
  newparams
}


# update .MISC environment of the stanfit object
# allows to call log_prob and other C++ using methods
# on objects not created in the current R session
# or objects created via another backend
# update_misc_env <- function(x, recompile = FALSE, only_windows = FALSE) {
#   stopifnot(is.brmsfit(x))
#   recompile <- as_one_logical(recompile)
#   only_windows <- as_one_logical(only_windows)
#   if (recompile || !has_rstan_model(x)) {
#     x <- add_rstan_model(x, overwrite = TRUE)
#   } else if (os_is_windows() || !only_windows) {
#     # TODO: detect when updating .MISC is not required
#     # TODO: find a more efficient way to update .MISC
#     old_backend <- x$backend
#     x$backend <- "rstan"
#     x$fit@.MISC <- suppressMessages(brm(fit = x, chains = 0))$fit@.MISC
#     x$backend <- old_backend
#   }
#   x
# }
