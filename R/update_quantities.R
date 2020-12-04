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
  lw_new <- update_ratio(draws, orig_log_prob_prop,
                         density_function_list,
                         ...)

  update_quantities_core(lw_new)
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

  lw_new <- update_target(draws, orig_log_prob_prop,
                          density_function_list,
                          ...)

  update_quantities_core(lw_new)
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
update_ratio <- function(draws, orig_log_prob_prop,
                         density_function_list,
                         ...) {
  log_ratio_draws_fun <- density_function_list$log_ratio_draws_fun
  log_prob_prop_draws_fun <- density_function_list$log_prob_prop_draws_fun

  log_ratio_new <- log_ratio_draws_fun(draws = draws, ...)
  log_prop_new <- log_prob_prop_draws_fun(draws = draws, ...)
  lw_new <- log_ratio_new + log_prop_new - orig_log_prob_prop

  lw_new
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
update_target <- function(draws, orig_log_prob_prop,
                          density_function_list,
                          ...) {

  log_prob_target_draws_fun <- density_function_list$log_prob_target_draws_fun

  log_prob_target_new <- log_prob_target_draws_fun(draws = draws, ...)
  lw_new <- log_prob_target_new - orig_log_prob_prop

  lw_new
}

update_proposal <- function(draws, orig_log_prob_prop,
                            density_function_list,
                            ...) {

  log_prob_prop_draws_fun <- density_function_list$log_prob_prop_draws_fun

  log_prop_new <- log_prob_prop_draws_fun(draws = draws, ...)
  lw_new <- log_prop_new - orig_log_prob_prop

  lw_new
}

#' Function for updating importance weights and pareto k diagnostic
#'
#' @param lw log importance weights
#' @return List with the updated log importance weights and the Pareto k.
#'
#' @noRd
update_quantities_core <- function(lw, ...) {

  psis <- suppressWarnings(loo::psis(lw))
  k <- psis$diagnostics$pareto_k
  lw <- as.vector(weights(psis))

  # gather results
  list(
    lw = lw,
    k = k
  )
}

#' Function for updating importance weights and pareto k diagnostic for
#' expectation-specific weights. For importance sampling expectations
#' where we have the `log_prob_target_draws_fun`.
#'
#' @param draws A matrix of draws.
#' @param orig_log_prob Log density of the proposal before moment matching.
#' @param density_function_list List of functions for computing the log
#' importance weights.
#' @param expectation_fun A function whose expectation is being computed.
#' The function takes arguments `draws`.
#' @param log_expectation_fun Logical indicating whether the expectation_fun
#' returns its values as logarithms or not. Defaults to FALSE. If set to TRUE,
#' the expectation function must be nonnegative (before taking the logarithm).
#' @return List with the updated log importance weights and the Pareto k.
#'
#' @noRd
update_quantities_target_expectation <- function(draws, orig_log_prob_prop,
                                                 density_function_list,
                                                 expectation_fun,
                                                 log_expectation_fun,
                                                 ...) {

  lw_new <- update_target(draws, orig_log_prob_prop,
                          density_function_list,
                          ...)
  lwf_new <- compute_lwf(draws, lw_new,
                         expectation_fun, log_expectation_fun,
                         ...)

  update_quantities_core(lwf_new)
}

#' Function for updating importance weights and pareto k diagnostic for
#' expectation-specific weights. For importance sampling expectations
#' where we have the `log_ratio_draws_fun`.
#'
#' @param draws A matrix of draws.
#' @param orig_log_prob Log density of the proposal before moment matching.
#' @param density_function_list List of functions for computing the log
#' importance weights.
#' @param expectation_fun A function whose expectation is being computed.
#' The function takes arguments `draws`.
#' @param log_expectation_fun Logical indicating whether the expectation_fun
#' returns its values as logarithms or not. Defaults to FALSE. If set to TRUE,
#' the expectation function must be nonnegative (before taking the logarithm).
#' @return List with the updated log importance weights and the Pareto k.
#'
#' @noRd
update_quantities_ratio_expectation <- function(draws, orig_log_prob_prop,
                                                density_function_list,
                                                expectation_fun,
                                                log_expectation_fun,
                                                ...) {

  lw_new <- update_ratio(draws, orig_log_prob_prop,
                         density_function_list,
                         ...)
  lwf_new <- compute_lwf(draws, lw_new, expectation_fun, log_expectation_fun, ...)

  update_quantities_core(lwf_new)
}

#' Function for updating importance weights and pareto k diagnostic for
#' expectation-specific weights. For simple Monte Carlo expectations.
#'
#' @param draws A matrix of draws.
#' @param orig_log_prob Log density of the proposal before moment matching.
#' @param density_function_list List of functions for computing the log
#' importance weights.
#' @param expectation_fun A function whose expectation is being computed.
#' The function takes arguments `draws`.
#' @param log_expectation_fun Logical indicating whether the expectation_fun
#' returns its values as logarithms or not. Defaults to FALSE. If set to TRUE,
#' the expectation function must be nonnegative (before taking the logarithm).
#' @return List with the updated log importance weights and the Pareto k.
#'
#' @noRd
update_quantities_expectation <- function(draws, orig_log_prob_prop,
                                          density_function_list,
                                          expectation_fun, log_expectation_fun,
                                          ...) {

  lw_new <- update_proposal(draws, orig_log_prob_prop,
                            density_function_list,
                            ...)

  lwf_new <- compute_lwf(draws, lw_new, expectation_fun, log_expectation_fun, ...)

  update_quantities_core(lwf_new)
}

#' Function for computing expectation-specific importance weights
#' from common importance weights.
#'
#' @param draws A matrix of draws.
#' @param lw common importance weights.
#' @param expectation_fun A function whose expectation is being computed.
#' The function takes arguments `draws`.
#' @param log_expectation_fun Logical indicating whether the expectation_fun
#' returns its values as logarithms or not. Defaults to FALSE. If set to TRUE,
#' the expectation function must be nonnegative (before taking the logarithm).
#' @return List with the updated log importance weights and the Pareto k.
#'
#' @noRd
compute_lwf <- function(draws, lw,
                        expectation_fun, log_expectation_fun,
                        ...) {
  if (log_expectation_fun) {
    lwf <- lw + expectation_fun(draws, ...)
  }
  else {
    lwf <- lw + log(abs(expectation_fun(draws, ...)))
  }
  lwf
}
