
#' Generic importance weighted moment matching (IWMM) importance sampling
#' algorithm. Matches a matrix of draws to its importance weighted moments, but
#' does not compute any expectation
#'
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
#'   matched.
#' @param ... Further arguments passed to `log_ratio_draws_fun`
#' and `log_prob_prop_draws_fun`.
#'
#' @return Returns a list with 3 elements: transformed draws, updated
#' importance weights, and the pareto k diagnostic value.
#'
#' @export
#' @importFrom stats weights
IWMM <- function(draws,
                 log_prob_prop_draws_fun,
                 log_prob_target_draws_fun = NULL,
                 log_ratio_draws_fun = NULL,
                 k_threshold = 0.5, cov_transform = TRUE, ...) {

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



