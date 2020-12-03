#' Generic importance weighted moment matching algorithm.
#'
#' Matches a matrix of draws to its importance weighted moments, but
#' does not compute any expectation.
#'
#' @param draws A matrix of draws.
#' @param ... Further arguments passed.
moment_match <- function(draws, ...) {
  UseMethod("moment_match")
}

#' Generic importance weighted moment matching algorithm for matrices.
#'
#' @param draws A matrix of draws.
#' @param log_prob_prop_draws_fun Log density of the proposal.
#' The function takes argument `draws`.
#' @param log_prob_target_draws_fun Log density of the target.
#' The function takes argument `draws`.
#' @param log_ratio_draws_fun Log of the density ratio (target/proposal).
#' The function takes argument `draws`.
#' @param k_threshold Threshold value for Pareto k values above which the moment
#'   matching algorithm is used. The default value is 0.5.
#' @param cov_transform Logical; Indicates whether to match the covariance of
#' the samples or not. If `FALSE`, only the mean and marginal variances are
#'   matched. Default is `TRUE`.
#' @param ... Further arguments passed to `log_prob_prop_draws_fun`,
#' `log_prob_target_draws_fun` and `log_ratio_draws_fun`.
#'
#' @return Returns a list with 3 elements: transformed draws, updated
#' importance weights, and the pareto k diagnostic value.
#'
#' @export
#' @importFrom stats weights
moment_match.matrix <- function(draws,
                                log_prob_prop_draws_fun,
                                log_prob_target_draws_fun = NULL,
                                log_ratio_draws_fun = NULL,
                                k_threshold = 0.5,
                                cov_transform = TRUE, ...) {

  checkmate::assertMatrix(draws, any.missing = FALSE)
  checkmate::assertFunction(log_prob_prop_draws_fun)
  checkmate::assertNumber(k_threshold)
  checkmate::assertLogical(cov_transform)


  if (is.null(log_prob_target_draws_fun) && is.null(log_ratio_draws_fun)) {
    stop("You must give either log_prob_target_draws_fun or log_ratio_draws_fun.")
  }

  orig_log_prob_prop <- log_prob_prop_draws_fun(draws = draws, ...)

  if (!is.null(log_prob_target_draws_fun)) {
    update_quantities <- update_quantities_target
    density_function_list <- list(log_prob_target_draws_fun = log_prob_target_draws_fun)
    lw <- log_prob_target_draws_fun(draws, ...) - orig_log_prob_prop
  }
  if (!is.null(log_ratio_draws_fun)) {
    update_quantities <- update_quantities_ratio
    density_function_list <- list(log_ratio_draws_fun = log_ratio_draws_fun, log_prob_prop_draws_fun = log_prob_prop_draws_fun)
    lw <- log_ratio_draws_fun(draws, ...)
  }



  lw_psis <- suppressWarnings(loo::psis(lw))
  lw <- as.vector(weights(lw_psis))
  k <- lw_psis$diagnostics$pareto_k

  npars <- ncol(draws)
  S <- nrow(draws)
  cov_transform <- cov_transform && S >= 10 * npars

  while (k > k_threshold) {


    # 1. match means
    trans <- shift(draws, lw)
    quantities <- update_quantities(
      draws = trans$draws,
      orig_log_prob_prop = orig_log_prob_prop,
      density_function_list,
      ...
    )
    if (quantities$k < k) {
      draws <- trans$draws

      lw <- quantities$lw
      k <- quantities$k
      next
    }

    # 2. match means and marginal variances
    trans <- shift_and_scale(draws, lw)
    quantities <- update_quantities(
      draws = trans$draws,
      orig_log_prob_prop = orig_log_prob_prop,
      density_function_list,
      ...
    )
    if (quantities$k < k) {
      draws <- trans$draws

      lw <- quantities$lw
      k <- quantities$k
      next
    }

    if (cov_transform) {
      # 3. match means and covariances
      trans <- shift_and_cov(draws, lw)
      quantities <- update_quantities(
        draws = trans$draws,
        orig_log_prob_prop = orig_log_prob_prop,
        density_function_list,
        ...
      )
      if (quantities$k < k) {
        draws <- trans$draws

        lw <- quantities$lw
        k <- quantities$k
        next
      }
    }


    break
  }

  list("draws" = draws, "log_weights" = lw, "pareto_k" = k)
}

#' Function for updating importance weights and pareto k diagnostic
#'
#' @param draws A matrix of draws.
#' @param orig_log_prob Log density of the proposal before moment matching.
#' @param density_function_list List of functions for computing the log
#' importance weights.
#' @return List with the updated log importance weights and the Pareto k.
#'
#' @noRd
update_quantities_ratio <- function(draws, orig_log_prob_prop,
                                    density_function_list,
                                    ...) {

  log_ratio_draws_fun <- density_function_list$log_ratio_draws_fun
  log_prob_prop_draws_fun <- density_function_list$log_prob_prop_draws_fun

  log_ratio_new <- log_ratio_draws_fun(draws = draws, ...)
  log_prop_new <- log_prob_prop_draws_fun(draws = draws, ...)

  lw_new <- log_ratio_new + log_prop_new - orig_log_prob_prop

  psis_new <- suppressWarnings(loo::psis(lw_new))
  k_new <- psis_new$diagnostics$pareto_k
  lw_new <- as.vector(weights(psis_new))

  # gather results
  list(
    lw = lw_new,
    k = k_new
  )
}

#' Function for updating importance weights and pareto k diagnostic
#'
#' @param draws A matrix of draws.
#' @param orig_log_prob Log density of the proposal before moment matching.
#' @param density_function_list List of functions for computing the log
#' importance weights.
#' @return List with the updated log importance weights and the Pareto k.
#'
#' @noRd
update_quantities_target <- function(draws, orig_log_prob_prop,
                                     density_function_list,
                                     ...) {

  log_prob_target_draws_fun <- density_function_list$log_prob_target_draws_fun

  log_prob_target_new <- log_prob_target_draws_fun(draws = draws, ...)

  lw_new <- log_prob_target_new - orig_log_prob_prop

  psis_new <- suppressWarnings(loo::psis(lw_new))
  k_new <- psis_new$diagnostics$pareto_k
  lw_new <- as.vector(weights(psis_new))

  # gather results
  list(
    lw = lw_new,
    k = k_new
  )
}


















# TODO obs_weights is just one possibility along with ratio and target

# TODO powersense takes in psis object. is that necessary? I think not

# TODO r_eff ?

# TODO max_iters?















moment_match.stanfit <- function(x, ...) {
  # TODO: ensure compatibility with objects not created in the current R session
  post_draws <- function(x, ...) {
    # ensure additional arguments are not passed further
    as.matrix(x)
  }
  out <- moment_match_modelfit(
    x,
    post_draws = post_draws,
    unconstrain_pars = unconstrain_pars.stanfit,
    log_prob_prop_draws_fun = log_prob_upars.stanfit,
    # log_ratio_upars = log_ratio_upars.stanfit,
    ...
  )
  out
}

log_prob_upars.stanfit <- function(x, upars, ...) {
  apply(upars, 1, rstan::log_prob,
        object = x,
        adjust_transform = TRUE, gradient = FALSE
  )
}

unconstrain_pars.stanfit <- function(x, pars, ...) {
  skeleton <- .create_skeleton(x@sim$pars_oi, x@par_dims[x@sim$pars_oi])
  upars <- apply(pars, 1, FUN = function(theta) {
    rstan::unconstrain_pars(x, pars = .rstan_relist(theta, skeleton))
  })
  # for one parameter models
  if (is.null(dim(upars))) {
    dim(upars) <- c(1, length(upars))
  }
  t(upars)
}

# -------- will be imported from rstan at some point -------
# create a named list of draws for use with rstan methods
.rstan_relist <- function(x, skeleton) {
  out <- utils::relist(x, skeleton)
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



#' Generic importance weighted moment matching algorithm for matrices.
#'
#' @param x A fitted model object.
#' @param post_draws A function the takes `x` as the first argument and returns
#'   a matrix of posterior draws of the model parameters (`pars`).
#' @param unconstrain_pars A function that takes arguments `x`, and `pars` and
#'   returns posterior draws on the unconstrained space based on the posterior
#'   draws on the constrained space passed via `pars`.
#' @param log_prob_prop_draws_fun Log density of the proposal.
#' The function takes argument `draws`.
#' @param log_prob_target_draws_fun Log density of the target.
#' The function takes argument `draws`.
#' @param log_ratio_draws_fun Log of the density ratio (target/proposal).
#' The function takes argument `draws`.
#' @param k_threshold Threshold value for Pareto k values above which the moment
#'   matching algorithm is used. The default value is 0.5.
#' @param cov_transform Logical; Indicates whether to match the covariance of
#' the samples or not. If `FALSE`, only the mean and marginal variances are
#'   matched. Default is `TRUE`.
#' @param ... Further arguments passed to `log_prob_prop_draws_fun`,
#' `log_prob_target_draws_fun` and `log_ratio_draws_fun`.
#'
#' @return Returns a list with 3 elements: transformed draws, updated
#' importance weights, and the pareto k diagnostic value.
#'
#' @export
#' @importFrom stats weights
moment_match_modelfit <- function(x,
                                  post_draws,
                                  unconstrain_pars,
                                  log_prob_prop_draws_fun,
                                  log_prob_target_draws_fun = NULL,
                                  log_ratio_draws_fun = NULL,
                                  k_threshold = 0.5,
                                  cov_transform = TRUE, ...) {

  checkmate::assertFunction(post_draws)
  checkmate::assertFunction(unconstrain_pars)
  checkmate::assertFunction(log_prob_prop_draws_fun)
  checkmate::assertNumber(k_threshold)
  checkmate::assertLogical(cov_transform)


  if (is.null(log_prob_target_draws_fun) && is.null(log_ratio_draws_fun)) {
    stop("You must give either log_prob_target_draws_fun or log_ratio_draws_fun.")
  }

  pars <- post_draws(x, ...)
  # transform the model parameters to unconstrained space
  draws <- unconstrain_pars(x, pars = pars, ...)

  orig_log_prob_prop <- log_prob_prop_draws_fun(x, draws = draws, ...)

  if (!is.null(log_prob_target_draws_fun)) {
    update_quantities <- update_quantities_target_modelfit
    density_function_list <- list(log_prob_target_draws_fun = log_prob_target_draws_fun)
    lw <- log_prob_target_draws_fun(x, draws, ...) - orig_log_prob_prop
  }
  if (!is.null(log_ratio_draws_fun)) {
    update_quantities <- update_quantities_ratio_modelfit
    density_function_list <- list(log_ratio_draws_fun = log_ratio_draws_fun, log_prob_prop_draws_fun = log_prob_prop_draws_fun)
    lw <- log_ratio_draws_fun(x, draws, ...)
  }



  lw_psis <- suppressWarnings(loo::psis(lw))
  lw <- as.vector(weights(lw_psis))
  k <- lw_psis$diagnostics$pareto_k

  npars <- ncol(draws)
  S <- nrow(draws)
  cov_transform <- cov_transform && S >= 10 * npars

  while (k > k_threshold) {


    # 1. match means
    trans <- shift(draws, lw)
    quantities <- update_quantities(
      x = x,
      draws = trans$draws,
      orig_log_prob_prop = orig_log_prob_prop,
      density_function_list,
      ...
    )
    if (quantities$k < k) {
      draws <- trans$draws

      lw <- quantities$lw
      k <- quantities$k
      next
    }

    # 2. match means and marginal variances
    trans <- shift_and_scale(draws, lw)
    quantities <- update_quantities(
      x = x,
      draws = trans$draws,
      orig_log_prob_prop = orig_log_prob_prop,
      density_function_list,
      ...
    )
    if (quantities$k < k) {
      draws <- trans$draws

      lw <- quantities$lw
      k <- quantities$k
      next
    }

    if (cov_transform) {
      # 3. match means and covariances
      trans <- shift_and_cov(draws, lw)
      quantities <- update_quantities(
        x = x,
        draws = trans$draws,
        orig_log_prob_prop = orig_log_prob_prop,
        density_function_list,
        ...
      )
      if (quantities$k < k) {
        draws <- trans$draws

        lw <- quantities$lw
        k <- quantities$k
        next
      }
    }


    break
  }

  list("draws" = draws, "log_weights" = lw, "pareto_k" = k)
}





#' Function for updating importance weights and pareto k diagnostic
#'
#' @param x A fitted model object.
#' @param draws A matrix of draws.
#' @param orig_log_prob Log density of the proposal before moment matching.
#' @param density_function_list List of functions for computing the log
#' importance weights.
#' @return List with the updated log importance weights and the Pareto k.
#'
#' @noRd
update_quantities_ratio_modelfit <- function(x, draws, orig_log_prob_prop,
                                             density_function_list,
                                             ...) {

  log_ratio_draws_fun <- density_function_list$log_ratio_draws_fun
  log_prob_prop_draws_fun <- density_function_list$log_prob_prop_draws_fun

  log_ratio_new <- log_ratio_draws_fun(x, draws = draws, ...)
  log_prop_new <- log_prob_prop_draws_fun(x, draws = draws, ...)

  lw_new <- log_ratio_new + log_prop_new - orig_log_prob_prop

  psis_new <- suppressWarnings(loo::psis(lw_new))
  k_new <- psis_new$diagnostics$pareto_k
  lw_new <- as.vector(weights(psis_new))

  # gather results
  list(
    lw = lw_new,
    k = k_new
  )
}


#' Function for updating importance weights and pareto k diagnostic
#'
#' @param x A fitted model object.
#' @param draws A matrix of draws.
#' @param orig_log_prob Log density of the proposal before moment matching.
#' @param density_function_list List of functions for computing the log
#' importance weights.
#' @return List with the updated log importance weights and the Pareto k.
#'
#' @noRd
update_quantities_target_modelfit <- function(x, draws, orig_log_prob_prop,
                                              density_function_list,
                                              ...) {

  log_prob_target_draws_fun <- density_function_list$log_prob_target_draws_fun

  log_prob_target_new <- log_prob_target_draws_fun(x, draws = draws, ...)

  lw_new <- log_prob_target_new - orig_log_prob_prop

  psis_new <- suppressWarnings(loo::psis(lw_new))
  k_new <- psis_new$diagnostics$pareto_k
  lw_new <- as.vector(weights(psis_new))

  # gather results
  list(
    lw = lw_new,
    k = k_new
  )
}









