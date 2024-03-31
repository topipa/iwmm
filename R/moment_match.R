#' Generic importance weighted moment matching algorithm.
#'
#' @export
moment_match <- function(x, ...) {
  UseMethod("moment_match")
}

#' @export
moment_match.draws_matrix <- function(x,
                                      log_prob_prop_fun,
                                      log_prob_target_fun = NULL,
                                      log_ratio_fun = NULL,
                                      expectation_fun = NULL,
                                      log_expectation_fun = FALSE,
                                      is_method = "psis",
                                      adaptation_method = "iwmm",
                                      k_threshold = 0.5,
                                      cov_transform = TRUE,
                                      split = FALSE,
                                      restart_transform = FALSE,
                                      ...) {
  return(moment_match.matrix(
    x = x,
    log_prob_prop_fun = log_prob_prop_fun,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    expectation_fun = expectation_fun,
    log_expectation_fun = log_expectation_fun,
    is_method = is_method,
    adaptation_method = adaptation_method,
    k_threshold = k_threshold,
    cov_transform = cov_transform,
    split = split,
    restart_transform = restart_transform,
    ...
  ))
}

#' @export
moment_match.draws_array <- function(x,
                                     log_prob_prop_fun,
                                     log_prob_target_fun = NULL,
                                     log_ratio_fun = NULL,
                                     expectation_fun = NULL,
                                     log_expectation_fun = FALSE,
                                     is_method = "psis",
                                     adaptation_method = "iwmm",
                                     k_threshold = 0.5,
                                     cov_transform = TRUE,
                                     split = FALSE,
                                     restart_transform = FALSE,
                                     ...) {
  return(moment_match.draws_matrix(
    x = posterior::as_draws_matrix(x),
    log_prob_prop_fun = log_prob_prop_fun,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    expectation_fun = expectation_fun,
    log_expectation_fun = log_expectation_fun,
    is_method = is_method,
    adaptation_method = adaptation_method,
    k_threshold = k_threshold,
    cov_transform = cov_transform,
    split = split,
    restart_transform = restart_transform,
    ...
  ))
}

#' @export
moment_match.draws_df <- function(x,
                                  log_prob_prop_fun,
                                  log_prob_target_fun = NULL,
                                  log_ratio_fun = NULL,
                                  expectation_fun = NULL,
                                  log_expectation_fun = FALSE,
                                  is_method = "psis",
                                  adaptation_method = "iwmm",
                                  k_threshold = 0.5,
                                  cov_transform = TRUE,
                                  split = FALSE,
                                  restart_transform = FALSE,
                                  ...) {
  return(moment_match.matrix(
    x = posterior::as_draws_matrix(x),
    log_prob_prop_fun = log_prob_prop_fun,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    expectation_fun = expectation_fun,
    log_expectation_fun = log_expectation_fun,
    is_method = is_method,
    adaptation_method = adaptation_method,
    k_threshold = k_threshold,
    cov_transform = cov_transform,
    split = split,
    restart_transform = restart_transform,
    ...
  ))
}

#' @export
moment_match.draws_list <- function(x,
                                    log_prob_prop_fun,
                                    log_prob_target_fun = NULL,
                                    log_ratio_fun = NULL,
                                    expectation_fun = NULL,
                                    log_expectation_fun = FALSE,
                                    is_method = "psis",
                                    adaptation_method = "iwmm",
                                    k_threshold = 0.5,
                                    cov_transform = TRUE,
                                    split = FALSE,
                                    restart_transform = FALSE,
                                    ...) {
  return(moment_match.matrix(
    x = posterior::as_draws_matrix(x),
    log_prob_prop_fun = log_prob_prop_fun,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    expectation_fun = expectation_fun,
    log_expectation_fun = log_expectation_fun,
    is_method = is_method,
    adaptation_method = adaptation_method,
    k_threshold = k_threshold,
    cov_transform = cov_transform,
    split = split,
    restart_transform = restart_transform,
    ...
  ))
}

#' @export
moment_match.draws_rvars <- function(x,
                                     log_prob_prop_fun,
                                     log_prob_target_fun = NULL,
                                     log_ratio_fun = NULL,
                                     expectation_fun = NULL,
                                     log_expectation_fun = FALSE,
                                     is_method = "psis",
                                     adaptation_method = "iwmm",
                                     k_threshold = 0.5,
                                     cov_transform = TRUE,
                                     split = FALSE,
                                     restart_transform = FALSE,
                                     ...) {
  return(moment_match.matrix(
    x = posterior::as_draws_matrix(x),
    log_prob_prop_fun = log_prob_prop_fun,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    expectation_fun = expectation_fun,
    log_expectation_fun = log_expectation_fun,
    is_method = is_method,
    adaptation_method = adaptation_method,
    k_threshold = k_threshold,
    cov_transform = cov_transform,
    split = split,
    restart_transform = restart_transform,
    ...
  ))
}

#' Generic importance weighted moment matching algorithm for matrices.
#'
#'
#' @param x A matrix of draws. Must be unconstrained.
#' @param log_prob_prop_fun Log density of the proposal.  The function
#'   takes argument `draws`.
#' @param log_prob_target_fun Log density of the target for importance
#'   sampling. The function takes argument `draws`.
#' @param log_ratio_fun Log of the density ratio for importance
#'   sampling (target/proposal). The function takes argument `draws`.
#' @param expectation_fun Optional argument, NULL by default. A
#'   function whose expectation is being computed. The function takes
#'   arguments `draws`.
#' @param log_expectation_fun Logical indicating whether the
#'   expectation_fun returns its values as logarithms or not. Defaults
#'   to FALSE. If set to TRUE, the expectation function must be
#'   nonnegative (before taking the logarithm).  Ignored if
#'   `expectation_fun` is NULL.
#' @param draws_transformation_fun Optional argument, NULL by default. A
#'   function that transforms draws before computing expectation. The function takes
#'   arguments `draws`.
#' @param is_method Which importance sampling method to use. Currently only `psis` is supported.
#' @param adaptation_method Which adaptation method to use. Currently only `iwmm` is supported.
#' @param k_threshold Threshold value for Pareto k values above which
#'   the moment matching algorithm is used. The default value is 0.5.
#' @param cov_transform Logical; Indicates whether to match the
#'   covariance of the samples or not. If `FALSE`, only the mean and
#'   marginal variances are matched. Default is `TRUE`.
#' @param split Logical; Indicate whether to do the split
#'   transformation or not at the end of moment matching. FALSE by
#'   default.
#' @param restart_transform Logical; When split is TRUE, indicates
#'   whether to start the second transformation from the original
#'   model parameters or the transformed parameters. If split is
#'   FALSE, this is ignored.
#' @param ... Further arguments passed to `log_prob_prop_fun`,
#'   `log_prob_target_fun` and `log_ratio_fun`.
#'
#' @return Returns a list with: transformed draws, updated importance
#'   weights, and the pareto k diagnostic value.  If expectation_fun
#'   is given, also returns the expectation.
#'
#' @rdname moment_match
#' @export
moment_match.matrix <- function(x,
                                log_prob_prop_fun,
                                log_prob_target_fun = NULL,
                                log_ratio_fun = NULL,
                                expectation_fun = NULL,
                                log_expectation_fun = FALSE,
                                draws_transformation_fun = NULL,
                                is_method = "psis",
                                adaptation_method = "iwmm",
                                k_threshold = 0.5,
                                cov_transform = TRUE,
                                split = FALSE,
                                restart_transform = FALSE,
                                ...) {
  draws <- x

  checkmate::assertMatrix(draws, any.missing = FALSE)
  checkmate::assertFunction(log_prob_prop_fun)
  checkmate::assertNumber(k_threshold)
  checkmate::assertLogical(cov_transform)
  checkmate::assertLogical(split)
  checkmate::assertLogical(restart_transform)
  checkmate::assertLogical(log_expectation_fun)


  if (is.null(log_prob_target_fun) && is.null(log_ratio_fun) && is.null(expectation_fun)) {
    stop("You must give either log_prob_target_fun or log_ratio_fun,
         or give an expectation_fun.")
  }
  if (!is.null(log_prob_target_fun) && !is.null(log_ratio_fun)) {
    stop("You cannot give both log_prob_target_fun and log_ratio_fun.")
  }

  if (is_method != "psis") {
    stop("Currently psis is the only supported is_method.")
  }
  if (adaptation_method != "iwmm") {
    stop("Currently iwmm is the only supported adaptation_method.")
  }

  orig_log_prob_prop <- log_prob_prop_fun(draws = draws, ...)

  npars <- ncol(draws)
  S <- nrow(draws)
  cov_transform <- cov_transform && S >= 10 * npars

  # initialize objects that keep track of the total transformation
  total_shift <- rep(0, npars)
  total_scaling <- rep(1, npars)
  total_mapping <- diag(npars)

  draws_orig <- draws

  if (is.null(log_prob_target_fun) && is.null(log_ratio_fun)) {
    lw <- rep(0, S)
    lw <- lw - matrixStats::logSumExp(lw)
    k <- 0

    lw_orig <- lw
  } else {
    if (!is.null(log_prob_target_fun)) {
      update_properties <- list(
        target_type = "target",
        expectation = FALSE,
        log_prob_target_fun = log_prob_target_fun
      )
      lw <- log_prob_target_fun(draws, ...) - orig_log_prob_prop
    }
    if (!is.null(log_ratio_fun)) {
      update_properties <- list(
        target_type = "ratio",
        expectation = FALSE,
        log_ratio_fun = log_ratio_fun,
        log_prob_prop_fun = log_prob_prop_fun
      )
      lw <- log_ratio_fun(draws, ...)
    }
    if (length(unique(lw)) == 1) {
      stop("All of the importance weights are equal. This indicates that your
           target density is equal to your proposal density.")
    }

    pareto_smoothed_lw <- posterior::pareto_smooth(
      lw - matrixStats::logSumExp(lw),
      are_log_weights = TRUE,
      tail = "right", extra_diags = TRUE, r_eff = 1,
      return_k = TRUE,
      verbose = FALSE
    )
    k <- pareto_smoothed_lw$diagnostics$khat
    lw <- as.vector(pareto_smoothed_lw$x)

    if (any(is.infinite(k))) {
      stop("Something went wrong, and encountered infinite Pareto k values..")
    }

    lw_orig <- lw

    trans_loop <- transform_loop(
      draws,
      lw,
      k,
      update_properties,
      orig_log_prob_prop,
      k_threshold,
      cov_transform,
      total_shift,
      total_scaling,
      total_mapping,
      ...
    )

    draws <- trans_loop$draws
    lw <- trans_loop$log_weights
    k <- trans_loop$pareto_k

    total_shift <- trans_loop$total_shift
    total_scaling <- trans_loop$total_scaling
    total_mapping <- trans_loop$total_mapping
  }

  if (is.null(expectation_fun)) {
    if (!is.null(draws_transformation_fun)) {
      draws <- draws_transformation_fun(draws)
    }
    adapted_draws <- list(
      draws = draws,
      log_weights = lw,
      expectation = NA,
      diagnostics = list(
        pareto_k = k,
        pareto_kf = NA
      )
    )
  } else {
    lwf <- compute_lwf(draws, lw, expectation_fun, log_expectation_fun, draws_transformation_fun, ...)

    pareto_smoothed_lwf <- apply(lwf, 2, function(x) {
      posterior::pareto_smooth(
        x,
        are_log_weights = TRUE,
        tail = "right", extra_diags = TRUE, r_eff = 1,
        return_k = TRUE, verbose = FALSE
      )
    })
    pareto_smoothed_lwf <- do.call(mapply, c(cbind, pareto_smoothed_lwf))
    kf <- as.numeric(pareto_smoothed_lwf$diagnostics["khat", ])

    if (split) {
      # prepare for split and check kfs
      if (restart_transform) {
        draws2 <- draws_orig
        # initialize objects that keep track of the total transformation
        total_shift2 <- rep(0, npars)
        total_scaling2 <- rep(1, npars)
        total_mapping2 <- diag(npars)
        lw <- lw_orig
      } else {
        draws2 <- draws
        # initialize objects that keep track of the total transformation
        total_shift2 <- total_shift
        total_scaling2 <- total_scaling
        total_mapping2 <- total_mapping
      }

      lwf_check <- compute_lwf(draws2, lw, expectation_fun, log_expectation_fun, draws_transformation_fun, ...)
      if (ncol(lwf_check) > 1) {
        stop("Using split = TRUE is not yet supported for expectation functions
              that return a matrix. As a workaround, you can wrap your function
              call using apply.")
      }
      lwf <- as.vector(pareto_smoothed_lwf$x)

      if (is.null(log_prob_target_fun) && is.null(log_ratio_fun)) {
        update_properties <- list(
          target_type = "simple",
          expectation = TRUE,
          expectation_fun = expectation_fun,
          log_expectation_fun = log_expectation_fun,
          log_prob_prop_fun = log_prob_prop_fun,
          draws_transformation_fun = draws_transformation_fun
        )
      } else if (!is.null(log_prob_target_fun)) {
        update_properties <- list(
          target_type = "target",
          expectation = TRUE,
          expectation_fun = expectation_fun,
          log_expectation_fun = log_expectation_fun,
          log_prob_target_fun = log_prob_target_fun,
          draws_transformation_fun = draws_transformation_fun
        )
      } else if (!is.null(log_ratio_fun)) {
        update_properties <- list(
          target_type = "ratio",
          expectation = TRUE,
          expectation_fun = expectation_fun,
          log_expectation_fun = log_expectation_fun,
          log_ratio_fun = log_ratio_fun,
          log_prob_prop_fun = log_prob_prop_fun,
          draws_transformation_fun = draws_transformation_fun
        )
      }

      trans_loop <- transform_loop(
        draws2,
        lwf,
        kf,
        update_properties,
        orig_log_prob_prop,
        k_threshold,
        cov_transform,
        total_shift2,
        total_scaling2,
        total_mapping2,
        ...
      )

      draws2 <- trans_loop$draws
      lwf <- trans_loop$log_weights
      kf <- trans_loop$pareto_k

      total_shift2 <- trans_loop$total_shift
      total_scaling2 <- trans_loop$total_scaling
      total_mapping2 <- trans_loop$total_mapping

      # second trasnsformations are done
      # we have updated draws2, lwf and kf
      # now compute split weights and split pars

      S_half <- as.integer(0.5 * S)
      take <- seq_len(S_half)


      # original parameters
      mean_original <- colMeans(draws_orig)

      # accumulated affine transformation
      draws_T1 <- sweep(draws_orig, 2, mean_original, "-")
      draws_T1 <- sweep(draws_T1, 2, total_scaling, "*")
      if (cov_transform) {
        draws_T1 <- tcrossprod(draws_T1, total_mapping)
      }
      draws_T1 <- sweep(draws_T1, 2, total_shift + mean_original, "+")

      draws_T2 <- sweep(draws_orig, 2, mean_original, "-")
      draws_T2 <- sweep(draws_T2, 2, total_scaling2, "*")
      if (cov_transform) {
        draws_T2 <- tcrossprod(draws_T2, total_mapping2)
      }
      draws_T2 <- sweep(draws_T2, 2, total_shift2 + mean_original, "+")

      # inverse accumulated affine transformation
      draws_T2_T1inv <- sweep(draws_T2, 2, mean_original + total_shift2, "-")
      if (cov_transform) {
        draws_T2_T1inv <- tcrossprod(draws_T2_T1inv, solve(total_mapping))
      }
      draws_T2_T1inv <- sweep(draws_T2_T1inv, 2, total_scaling, "/")
      draws_T2_T1inv <- sweep(
        draws_T2_T1inv, 2,
        mean_original + total_shift2 - total_shift, "+"
      )

      draws_T1_T2inv <- sweep(draws_T1, 2, mean_original + total_shift, "-")
      if (cov_transform) {
        draws_T1_T2inv <- tcrossprod(draws_T1_T2inv, solve(total_mapping2))
      }
      draws_T1_T2inv <- sweep(draws_T1_T2inv, 2, total_scaling2, "/")
      draws_T1_T2inv <- sweep(
        draws_T1_T2inv, 2,
        mean_original + total_shift - total_shift2, "+"
      )

      # these are the real used draws
      # first half of draws_trans are T1(theta)
      # second half are T2(theta)
      draws_trans <- draws_T2
      draws_trans[take, ] <- draws_T1[take, , drop = FALSE]

      # then we need two sets of pseudo draws
      draws_trans_inv1 <- draws_orig
      draws_trans_inv1[take, ] <- draws_T1_T2inv[take, , drop = FALSE]

      draws_trans_inv2 <- draws_T2_T1inv
      draws_trans_inv2[take, ] <- draws_orig[take, , drop = FALSE]

      log_prob_trans_inv1 <- log_prob_prop_fun(draws_trans_inv1, ...)
      log_prob_trans_inv2 <- log_prob_prop_fun(draws_trans_inv2, ...)

      if (is.null(log_prob_target_fun) && is.null(log_ratio_fun)) {
        log_prob_prop_trans <- log_prob_prop_fun(draws_trans, ...)
        lw_trans <- log_prob_prop_trans -
          log(
            exp(log_prob_trans_inv1 - log(prod(total_scaling2)) -
                  log(det(total_mapping2))) + # styler: off
              exp(log_prob_trans_inv2 -
                    log(prod(total_scaling)) - # styler: off
                    log(det(total_mapping))) # styler: off
          )
      } else if (!is.null(log_prob_target_fun)) {
        log_prob_target_trans <- log_prob_target_fun(draws_trans, ...)
        lw_trans <- log_prob_target_trans -
          log(
            exp(log_prob_trans_inv1 - log(prod(total_scaling2)) -
                  log(det(total_mapping2))) + # styler: off
              exp(log_prob_trans_inv2 - log(prod(total_scaling)) -
                    log(det(total_mapping))) # styler: off
          )
      } else if (!is.null(log_ratio_fun)) {
        log_prob_ratio_trans <- log_ratio_fun(draws_trans, ...)
        log_prob_prop_trans <- log_prob_prop_fun(draws_trans, ...)
        lw_trans <- log_prob_ratio_trans + log_prob_prop_trans -
          log(
            exp(log_prob_trans_inv1 - log(prod(total_scaling2))
                - log(det(total_mapping2))) + # styler: off
              exp(log_prob_trans_inv2 - log(prod(total_scaling))
                  - log(det(total_mapping))) # styler: off
          )
      }

      # log_prob_target_trans <- log_prob_target_fun(draws_trans, ...)
      #
      # lw_trans <-  log_prob_target_trans -
      #   log(
      #     exp(log_prob_trans_inv1 - log(prod(total_scaling2)) -
      # log(det(total_mapping2))) +
      #       exp(log_prob_trans_inv2 - log(prod(total_scaling)) -
      # log(det(total_mapping)))
      #   )


      lw_trans <- lw_trans - matrixStats::logSumExp(lw_trans)


      # replace draws and lw
      lw <- lw_trans
      draws <- draws_trans
    } else {
      # if not splitting, warn about high pareto ks
      if (any(kf > k_threshold)) {
        warning("Importance sampling may be unreliable.
                 Consider setting split to TRUE.")
      }
    }

    if (!is.null(draws_transformation_fun)) {
      draws <- draws_transformation_fun(draws)
    }
    unweighted_expectation <- expectation_fun(draws, ...)

    if (log_expectation_fun) {
      expectation <- exp(matrixStats::colLogSumExps(
        lw + unweighted_expectation
      ))
    } else {
      w <- exp(lw)
      expectation <- colSums(w * unweighted_expectation)
    }

    adapted_draws <- list(
      draws = draws,
      log_weights = lw,
      expectation = expectation,
      diagnostics = list(
        pareto_k = k,
        pareto_kf = kf
      )
    )
  }
  class(adapted_draws) <- c("adapted_importance_sampling", class(adapted_draws))

  return(adapted_draws)
}

#' Generic importance weighted moment matching algorithm for
#' `CmdStanFit` objects.  See additional arguments from
#' `moment_match.matrix`
#'
#' @param x A fitted `CmdStanFit` object.
#' @param log_prob_target_fun Log density of the target.  The function
#'   takes argument `draws`, which are the unconstrained draws.
#' @param log_ratio_fun Log of the density ratio (target/proposal).
#'   The function takes argument `draws`, which are the unconstrained
#'   draws.
#' @param constrain_draws Logical specifying whether to return draws on the
#'   constrained space. Draws are also constrained for computing expectations. Default is TRUE.
#' @param ... Further arguments passed to `moment_match.matrix`.
#'
#' @return Returns a list with 3 elements: transformed draws, updated
#' importance weights, and the pareto k diagnostic value.
#'
#' @export
moment_match.CmdStanFit <- function(x,
                                    log_prob_target_fun = NULL,
                                    log_ratio_fun = NULL,
                                    constrain_draws = TRUE,
                                    ...) {
  # TODO: support expectation fun?
  # and add tests for that

  # transform the model parameters to unconstrained space
  udraws <- x$unconstrain_draws(format = "draws_matrix")

  if (constrain_draws) {
    draws_transformation_fun <- function(draws, ...) {
      return(constrain_draws(x, draws, ...))
    }
  } else {
    draws_transformation_fun <- NULL
  }

  out <- moment_match(
    x = udraws,
    log_prob_prop_fun = log_prob_draws.CmdStanFit,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    draws_transformation_fun = draws_transformation_fun,
    fit = x,
    ...
  )

  out
}

#' Generic importance weighted moment matching algorithm for `stanfit` objects.
#' See additional arguments from `moment_match.matrix`
#'
#' @param x A fitted `stanfit` object.
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
#' @param constrain_draws Logical specifying whether to return draws on the
#'   constrained space. Draws are also constrained for computing expectations. Default is TRUE.
#' @param ... Further arguments passed to `moment_match.matrix`.
#'
#' @return Returns a list with 3 elements: transformed draws, updated
#'   importance weights, and the pareto k diagnostic value. If expectation_fun
#'   is given, also returns the expectation.
#'
#' @export
moment_match.stanfit <- function(x,
                                 log_prob_target_fun = NULL,
                                 log_ratio_fun = NULL,
                                 target_observation_weights = NULL,
                                 expectation_fun = NULL,
                                 log_expectation_fun = FALSE,
                                 constrain_draws = TRUE,
                                 ...) {
  if (!is.null(target_observation_weights) && (!is.null(log_prob_target_fun) || !is.null(log_ratio_fun))) {
    stop("You must give only one of target_observation_weights, log_prob_target_fun, or log_ratio_fun.")
  }

  # TODO: should we give option to return only subset of draws?
  # TODO: should draws_transformation_fun be allowed here?
  # TODO: should it be possible to set different values for returning constrained draws and
  # whether or not computing expectation with constrained draws??

  # ensure draws are in matrix form
  draws <- posterior::as_draws_matrix(x)

  if (!is.null(target_observation_weights)) {
    out <- tryCatch(posterior::subset_draws(draws, variable = "log_lik"),
      error = function(cond) {
        message(cond)
        message("\nYour stan fit does not include a parameter called log_lik.")
        message("To use target_observation_weights, you must define log_lik in the generated quantities block.")
        return(NA)
      }
    )

    log_ratio_fun <- function(draws, fit, ...) {
      cdraws <- constrain_draws(fit, draws)
      ll <- posterior::merge_chains(
        posterior::subset_draws(cdraws, variable = "log_lik")
      )
      colSums(t(drop(ll)) * (target_observation_weights - 1))
    }
  }


  # transform the draws to unconstrained space
  udraws <- unconstrain_draws(x, draws = draws, ...)

  if (constrain_draws) {
    draws_transformation_fun <- function(draws, ...) {
      return(constrain_draws(x, draws, ...))
    }
  } else {
    draws_transformation_fun <- NULL
  }

  out <- moment_match.matrix(
    udraws,
    log_prob_prop_fun = log_prob_draws.stanfit,
    log_prob_target_fun = log_prob_target_fun,
    log_ratio_fun = log_ratio_fun,
    expectation_fun = expectation_fun,
    log_expectation_fun = log_expectation_fun,
    draws_transformation_fun = draws_transformation_fun,
    fit = x,
    ...
  )

  # TODO: should this function update the parameters of the stanfit and return it?
  out
}
