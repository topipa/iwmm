#' Compute an expectation.
#'
#' @export
compute_expectation <- function(x, ...) {
  UseMethod("compute_expectation")
}

#' @rdname compute_expectation
#' @export
compute_expectation.matrix <- function(draws,
                                       lw,
                                       expectation_fun,
                                       transformation_fun = NULL,
                                log_expectation_fun = FALSE,
                                ...) {

  checkmate::assertMatrix(draws, any.missing = FALSE)
  checkmate::assertVector(lw, any.missing = FALSE)
  checkmate::assertFunction(expectation_fun)
  checkmate::assertLogical(log_expectation_fun)

  if (!is.null(transformation_fun)) {
    draws <- transformation_fun(draws)
  }

  if (log_expectation_fun) {
      expectation <- exp(matrixStats::colLogSumExps(
        lw + expectation_fun(draws, ...)
      ))
  } else {
    w <- exp(lw)
    expectation <- colSums(w * expectation_fun(draws, ...))
  }

  lwf <- compute_lwf(
    draws,
    lw,
    expectation_fun,
    log_expectation_fun,
    ...
  )
  pareto_smoothed_wf <- posterior::pareto_smooth(exp(lwf),
                                                tail = "right", extra_diags = TRUE, r_eff = 1
  )
  kf <- pareto_smoothed_wf$diagnostics$khat

  ret <- list(expectation = expectation,
              lwf=lwf,
              diagnostics = list(
                pareto_k = NA,
                pareto_kf = kf
              ))
  class(ret) <- c("importance_sampling_expectation", class(ret))

  return(ret)

}


