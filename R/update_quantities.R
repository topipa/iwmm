#' Function for updating importance weights and pareto k diagnostic.
#'
#' @param draws A matrix of draws.
#' @param orig_log_prob Log density of the proposal before moment matching.
#' @param update_properties List of properties to define how quantities
#' are updated.
#' @param expectation_fun A function whose expectation is being computed.
#' The function takes arguments `draws`.
#' @param log_expectation_fun Logical indicating whether the expectation_fun
#' returns its values as logarithms or not. Defaults to FALSE. If set to TRUE,
#' the expectation function must be nonnegative (before taking the logarithm).
#' @return List with the updated log importance weights and the Pareto k.
#'
#' @noRd
update_quantities <- function(draws, orig_log_prob_prop,
                              update_properties,
                              ...) {
  if (update_properties$target_type == "simple") {
    log_prob_prop_fun <- update_properties$log_prob_prop_fun
    log_prop_new <- log_prob_prop_fun(draws = draws, ...)
    lw_new <- log_prop_new - orig_log_prob_prop
  } else if (update_properties$target_type == "ratio") {
    log_ratio_fun <- update_properties$log_ratio_fun
    log_prob_prop_fun <- update_properties$log_prob_prop_fun
    log_ratio_new <- log_ratio_fun(draws = draws, ...)
    log_prop_new <- log_prob_prop_fun(draws = draws, ...)
    lw_new <- log_ratio_new + log_prop_new - orig_log_prob_prop
  } else if (update_properties$target_type == "target") {
    log_prob_target_fun <- update_properties$log_prob_target_fun
    log_prob_target_new <- log_prob_target_fun(draws = draws, ...)
    lw_new <- log_prob_target_new - orig_log_prob_prop
  }

  if (update_properties$expectation) {
    lw_new <- compute_lwf(
      draws,
      lw_new,
      update_properties$expectation_fun,
      update_properties$log_expectation_fun,
      update_properties$draws_transformation_fun,
      ...
    )
  }

  pareto_smoothed_w_new <- posterior::pareto_smooth(
    lw_new - matrixStats::logSumExp(lw_new),
    are_log_weights = TRUE,
    tail = "right", r_eff = 1, return_k = TRUE,
    verbose = FALSE
  )
  k <- pareto_smoothed_w_new$diagnostics$khat
  lw <- as.vector(pareto_smoothed_w_new$x)

  # normalize log weights
  lw <- lw - matrixStats::logSumExp(lw)

  # gather results
  list(
    lw = lw,
    k = k
  )
}

#' Function for computing expectation-specific log importance weights
#' from common log importance weights.
#'
#' @param draws A matrix of draws.
#' @param lw common log importance weights.
#' @param expectation_fun A function whose expectation is being computed.
#' The function takes arguments `draws`.
#' @param log_expectation_fun Logical indicating whether the expectation_fun
#' returns its values as logarithms or not. Defaults to FALSE. If set to TRUE,
#' the expectation function must be nonnegative (before taking the logarithm).
#' @param draws_transformation_fun Optional argument, NULL by default. A
#'   function that transforms draws before computing expectation. The function takes
#'   arguments `draws`.
#' @return expectation-specific log importance weights.
#'
#' @noRd
compute_lwf <- function(draws, lw,
                        expectation_fun, log_expectation_fun,
                        draws_transformation_fun,
                        ...) {
  if (!is.null(draws_transformation_fun)) {
    draws <- draws_transformation_fun(draws)
  }
  unweighted_expectation <- expectation_fun(draws, ...)

  if (log_expectation_fun) {
    lwf <- lw + unweighted_expectation
  } else {
    lwf <- lw + log(abs(unweighted_expectation))
  }
  lwf
}
