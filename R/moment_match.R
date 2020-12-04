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
  if (!is.null(log_prob_target_draws_fun) && !is.null(log_ratio_draws_fun)) {
    stop("You cannot give both log_prob_target_draws_fun and log_ratio_draws_fun.")
  }

  orig_log_prob_prop <- log_prob_prop_draws_fun(draws = draws, ...)

  if (!is.null(log_prob_target_draws_fun)) {
    update_quantities <- update_quantities_target
    density_function_list <- list(expectation = FALSE,
                                  log_prob_target_draws_fun = log_prob_target_draws_fun)
    lw <- log_prob_target_draws_fun(draws, ...) - orig_log_prob_prop
  }
  if (!is.null(log_ratio_draws_fun)) {
    update_quantities <- update_quantities_ratio
    density_function_list <- list(expectation = FALSE,
                                  log_ratio_draws_fun = log_ratio_draws_fun,
                                  log_prob_prop_draws_fun = log_prob_prop_draws_fun)
    lw <- log_ratio_draws_fun(draws, ...)
  }


  lw_psis <- suppressWarnings(loo::psis(lw))
  lw <- as.vector(weights(lw_psis))
  k <- lw_psis$diagnostics$pareto_k

  npars <- ncol(draws)
  S <- nrow(draws)
  cov_transform <- cov_transform && S >= 10 * npars

  trans_loop <- transform_loop(draws,
                               lw,
                               k,
                               update_quantities,
                               density_function_list,
                               orig_log_prob_prop,
                               k_threshold,
                               cov_transform,
                               total_shift = rep(0,npars),
                               total_scaling = rep(1,npars),
                               total_mapping = diag(npars),
                               ...)
  draws <- trans_loop$draws
  lw <- trans_loop$log_weights
  k <- trans_loop$pareto_k


  list("draws" = draws, "log_weights" = lw, "pareto_k" = k)
}

transform_loop <- function(draws,
                      lw,
                      k,
                      update_quantities,
                      density_function_list,
                      orig_log_prob_prop,
                      k_threshold,
                      cov_transform,
                      total_shift,
                      total_scaling,
                      total_mapping,
                      ...) {
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
      total_shift <- total_shift + trans$shift

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
      total_shift <- total_shift + trans$shift
      total_scaling <- total_scaling * trans$scaling

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
        total_shift <- total_shift + trans$shift
        total_mapping <- trans$mapping %*% total_mapping

        lw <- quantities$lw
        k <- quantities$k
        next
      }
    }


    break
  }

  list("draws" = draws, "log_weights" = lw, "pareto_k" = k,
       "total_shift" = total_shift, "total_scaling" = total_scaling,
       "total_mapping" = total_mapping)
}
