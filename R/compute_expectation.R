#' Compute an expectation.
#'
#' @export
compute_expectation <- function(x, ...) {
  UseMethod("compute_expectation")
}

#' @rdname compute_expectation
#' @export
compute_expectation.matrix <- function(x,
                                       lw,
                                       expectation_fun,
                                       log_expectation_fun = FALSE,
                                       draws_transformation_for_expectation_fun = NULL,
                                       ...) {
  checkmate::assertVector(lw, any.missing = FALSE)
  checkmate::assertFunction(expectation_fun)
  checkmate::assertLogical(log_expectation_fun)

  if(!is.null(draws_transformation_for_expectation_fun)) {
    transformed_draws_for_expectation <- draws_transformation_for_expectation_fun(x)
    unweighted_expectation <- expectation_fun(transformed_draws_for_expectation, ...)
  } else {
    unweighted_expectation <- expectation_fun(x, ...)
  }

  if (log_expectation_fun) {
    expectation <- exp(matrixStats::colLogSumExps(
      lw + unweighted_expectation
    ))
  } else {
    w <- exp(lw)
    expectation <- colSums(w * unweighted_expectation)
  }

  lwf <- compute_lwf(
    x,
    lw,
    expectation_fun,
    log_expectation_fun,
    draws_transformation_for_expectation_fun,
    ...
  )
  pareto_smoothed_wf <- posterior::pareto_smooth(exp(lwf),
    tail = "right", extra_diags = TRUE, r_eff = 1
  )
  kf <- pareto_smoothed_wf$diagnostics$khat

  ret <- list(
    expectation = expectation,
    lwf = lwf,
    diagnostics = list(
      pareto_k = NA,
      pareto_kf = kf
    )
  )
  class(ret) <- c("importance_sampling_expectation", class(ret))

  return(ret)
}

#' @export
compute_expectation.adapted_importance_sampling <- function(x,
                                                            expectation_fun,
                                                            log_expectation_fun = FALSE,
                                                            draws_transformation_for_expectation_fun = NULL,
                                                            ...) {
  return(compute_expectation.matrix(
    x = x$draws,
    lw = x$log_weights,
    expectation_fun = expectation_fun,
    log_expectation_fun = log_expectation_fun,
    draws_transformation_for_expectation_fun=draws_transformation_for_expectation_fun,
    ...
  ))
}

#' @export
compute_expectation.draws_matrix <- function(x,
                                             lw,
                                             expectation_fun,
                                             log_expectation_fun = FALSE,
                                             draws_transformation_for_expectation_fun = NULL,
                                             ...) {
  return(compute_expectation.matrix(
    x = x,
    lw = lw,
    expectation_fun = expectation_fun,
    log_expectation_fun = log_expectation_fun,
    draws_transformation_for_expectation_fun=draws_transformation_for_expectation_fun,
    ...
  ))
}


#' @export
compute_expectation.draws_array <- function(x,
                                            lw,
                                            expectation_fun,
                                            log_expectation_fun = FALSE,
                                            draws_transformation_for_expectation_fun = NULL,
                                            ...) {
  return(compute_expectation.draws_matrix(
    x = posterior::as_draws_matrix(x),
    lw = lw,
    expectation_fun = expectation_fun,
    log_expectation_fun = log_expectation_fun,
    draws_transformation_for_expectation_fun=draws_transformation_for_expectation_fun,
    ...
  ))
}
