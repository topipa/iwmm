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

  # if (!is.null(target_observation_weights)) {
  #   out <- tryCatch(posterior::subset_draws(draws, variable = "log_lik"),
  #                   error = function(cond) {
  #                     message(cond)
  #                     message("\nYour stan fit does not include a parameter called log_lik.")
  #                     message("To use target_observation_weights, you must define log_lik in the generated quantities block.")
  #                     return(NA)
  #                   }
  #   )
  #
  #   log_ratio_fun <- function(draws, fit, ...) {
  #     cdraws <- constrain_draws(fit, draws)
  #     ll <- posterior::merge_chains(
  #       posterior::subset_draws(cdraws, variable = "log_lik")
  #     )
  #     colSums(t(drop(ll)) * (target_observation_weights - 1))
  #   }
  # }


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

  if (constrain) {
    out$draws <- constrain_draws.stanfit(x$fit, out$draws, ...)
  }

  out
}


log_prob_draws.brmsfit <- function(fit, draws, ...) {
  # x <- update_misc_env(x, only_windows = TRUE)
  log_prob_draws.stanfit(fit$fit, draws = draws, ...)
}

unconstrain_draws.brmsfit <- function(x, draws, ...) {
  unconstrain_draws.stanfit(x$fit, draws = draws, ...)
}

constrain_draws.brmsfit <- function(x, udraws, ...) {
  out <- rstan::constrain_pars(udraws, object = x$fit)
  out[x$exclude] <- NULL
  out
}

# # transform parameters to the constraint space
# update_pars_brmsfit <- function(x, draws, ...) {
#   # list with one element per posterior draw
#   pars <- apply(draws, 1, constrain_draws.brmsfit, x = x)
#   # select required parameters only
#   pars <- lapply(pars, "[", x$fit@sim$pars_oi_old)
#   # transform draws
#   ndraws <- length(pars)
#   pars <- unlist(pars)
#   npars <- length(pars) / ndraws
#   dim(pars) <- c(npars, ndraws)
#   # add dummy 'lp__' draws
#   pars <- rbind(pars, rep(0, ndraws))
#   # bring draws into the right structure
#   new_draws <- named_list(x$fit@sim$fnames_oi_old, list(numeric(ndraws)))
#   if (length(new_draws) != nrow(pars)) {
#     stop2("Updating parameters in `update_pars_brmsfit' failed. ",
#           "Please report a bug at https://github.com/paul-buerkner/brms.")
#   }
#   for (i in seq_len(npars)) {
#     new_draws[[i]] <- pars[i, ]
#   }
#   # create new sim object to overwrite x$fit@sim
#   x$fit@sim <- list(
#     samples = list(new_draws),
#     iter = ndraws,
#     thin = 1,
#     warmup = 0,
#     chains = 1,
#     n_save = ndraws,
#     warmup2 = 0,
#     permutation = list(seq_len(ndraws)),
#     pars_oi = x$fit@sim$pars_oi_old,
#     dims_oi = x$fit@sim$dims_oi_old,
#     fnames_oi = x$fit@sim$fnames_oi_old,
#     n_flatnames = length(x$fit@sim$fnames_oi_old)
#   )
#   x$fit@stan_args <- list(
#     list(chain_id = 1, iter = ndraws, thin = 1, warmup = 0)
#   )
#   rename_pars(x)
# }

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
