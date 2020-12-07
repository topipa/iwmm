#' Generic importance weighted moment matching algorithm.
#'
#'
#' @param draws A matrix of draws.
#' @param ... Further arguments passed.
#'
#' @export
moment_match <- function(draws, ...) {
  UseMethod("moment_match")
}

#' Generic importance weighted moment matching algorithm for matrices.
#'
#'
#' @param draws A matrix of draws. Must be unconstrained.
#' @param log_prob_prop_fun Log density of the proposal.
#' The function takes argument `draws`.
#' @param log_prob_target_fun Log density of the target for
#' importance sampling. The function takes argument `draws`.
#' @param log_ratio_fun Log of the density ratio for importance sampling
#' (target/proposal). The function takes argument `draws`.
#' @param expectation_fun Optional argument, NULL by default. A function whose expectation is
#' being computed. The function takes arguments `draws`.
#' @param log_expectation_fun Logical indicating whether the expectation_fun
#' returns its values as logarithms or not. Defaults to FALSE. If set to TRUE,
#' the expectation function must be nonnegative (before taking the logarithm).
#' Ignored if `expectation_fun` is NULL.
#' @param k_threshold Threshold value for Pareto k values above which the moment
#'   matching algorithm is used. The default value is 0.5.
#' @param cov_transform Logical; Indicates whether to match the covariance of
#' the samples or not. If `FALSE`, only the mean and marginal variances are
#'   matched. Default is `TRUE`.
#' @param split Logical; Indicate whether to do the split transformation or not
#' at the end of moment matching. FALSE by default.
#' @param restart_transform Logical; When split is TRUE, indicates whether to
#' start the second transformation from the original model parameters
#' or the transformed parameters. If split is FALSE, this is ignored.
#' @param ... Further arguments passed to `log_prob_prop_fun`,
#' `log_prob_target_fun` and `log_ratio_fun`.
#'
#' @return Returns a list with: transformed draws, updated
#' importance weights, and the pareto k diagnostic value.
#' If expectation_fun is given, also returns the expectation.
#'
#' @rdname moment_match
#' @export
#' @importFrom stats weights
moment_match.matrix <- function(draws,
                         log_prob_prop_fun,
                         log_prob_target_fun = NULL,
                         log_ratio_fun = NULL,
                         expectation_fun = NULL,
                         log_expectation_fun = FALSE,
                         k_threshold = 0.5,
                         cov_transform = TRUE,
                         split = FALSE,
                         restart_transform = FALSE,
                         ...) {

  checkmate::assertMatrix(draws, any.missing = FALSE)
  checkmate::assertFunction(log_prob_prop_fun)
  checkmate::assertNumber(k_threshold)
  checkmate::assertLogical(cov_transform)
  checkmate::assertLogical(split)
  checkmate::assertLogical(restart_transform)
  checkmate::assertLogical(log_expectation_fun)


  if (is.null(log_prob_target_fun) && is.null(log_ratio_fun)
      && is.null(expectation_fun) ) {
    stop("You must give either log_prob_target_fun or log_ratio_fun,
         or give an expectation_fun.")
  }
  if (!is.null(log_prob_target_fun) && !is.null(log_ratio_fun)) {
    stop("You cannot give both log_prob_target_fun and log_ratio_fun.")
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
    lw <- - matrixStats::logSumExp(rep(0,nrow(draws)))
    k <- 0

    lw_orig <- lw
  } else {
    if (!is.null(log_prob_target_fun)) {
      update_properties <- list(target_type = "target",
                                    expectation = FALSE,
                                    log_prob_target_fun = log_prob_target_fun)
      lw <- log_prob_target_fun(draws, ...) - orig_log_prob_prop
    }
    if (!is.null(log_ratio_fun)) {
      update_properties <- list(target_type = "ratio",
                                    expectation = FALSE,
                                    log_ratio_fun = log_ratio_fun,
                                    log_prob_prop_fun = log_prob_prop_fun)
      lw <- log_ratio_fun(draws, ...)
    }

    lw_psis <- suppressWarnings(loo::psis(lw))
    lw <- as.vector(weights(lw_psis))
    k <- lw_psis$diagnostics$pareto_k

    lw_orig <- lw

    trans_loop <- transform_loop(draws,
                                 lw,
                                 k,
                                 update_quantities,
                                 update_properties,
                                 orig_log_prob_prop,
                                 k_threshold,
                                 cov_transform,
                                 total_shift,
                                 total_scaling,
                                 total_mapping,
                                 ...)

    draws <- trans_loop$draws
    lw <- trans_loop$log_weights
    k <- trans_loop$pareto_k

    total_shift <- trans_loop$total_shift
    total_scaling <- trans_loop$total_scaling
    total_mapping <- trans_loop$total_mapping



  }




  if (is.null(expectation_fun)) {
    list("draws" = draws, "log_weights" = lw, "pareto_k" = k)
  } else {


    # prepare for split and check kfs
    if (restart_transform) {
      draws2 <- draws_orig
      total_shift2 <- rep(0, npars)
      total_scaling2 <- rep(1, npars)
      total_mapping2 <- diag(npars)
      lw <- lw_orig

    }
    else {
      draws2 <- draws
      # initialize objects that keep track of the total transformation
      total_shift2 <- total_shift
      total_scaling2 <- total_scaling
      total_mapping2 <- total_mapping
    }

    lwf <- compute_lwf(draws2, lw, expectation_fun, log_expectation_fun, ...)


    psisf <- suppressWarnings(loo::psis(lwf))
    kf <- psisf$diagnostics$pareto_k

    if (split) {


      if (ncol(lwf) > 1) {
        stop('Using split = TRUE is not yet supported for expectation functions that return a matrix.
           As a workaround, you can wrap your function call using apply.')
      }

      if (is.null(log_prob_target_fun) && is.null(log_ratio_fun)) {
        update_properties <- list(target_type = "simple",
                                  expectation = TRUE,
                                  expectation_fun = expectation_fun,
                                  log_expectation_fun = log_expectation_fun,
                                  log_prob_prop_fun = log_prob_prop_fun)
      } else if (!is.null(log_prob_target_fun)) {
        update_properties <- list(target_type = "target",
                                  expectation = TRUE,
                                  expectation_fun = expectation_fun,
                                  log_expectation_fun = log_expectation_fun,
                                  log_prob_target_fun = log_prob_target_fun)
      } else if (!is.null(log_ratio_fun)) {
        update_properties <- list(target_type = "ratio",
                                  expectation = TRUE,
                                  expectation_fun = expectation_fun,
                                  log_expectation_fun = log_expectation_fun,
                                  log_ratio_fun = log_ratio_fun,
                                  log_prob_prop_fun = log_prob_prop_fun)
      }

      lwf <- as.vector(weights(psisf))


      trans_loop <- transform_loop(draws2,
                                   lwf,
                                   kf,
                                   update_quantities,
                                   update_properties,
                                   orig_log_prob_prop,
                                   k_threshold,
                                   cov_transform,
                                   total_shift2,
                                   total_scaling2,
                                   total_mapping2,
                                   ...)

      draws2 <- trans_loop$draws
      lwf <- trans_loop$log_weights
      kf <- trans_loop$pareto_k

      total_shift2 <- trans_loop$total_shift
      total_scaling2 <- trans_loop$total_scaling
      total_mapping2 <- trans_loop$total_mapping



      # second trasnsformations are done
      # we have updated draws2, lwf and kf
      # now compute split weights and split pars

      S <- nrow(draws2)
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
      draws_T2_T1inv <- sweep(draws_T2_T1inv, 2, mean_original + total_shift2 - total_shift, "+")

      draws_T1_T2inv <- sweep(draws_T1, 2, mean_original + total_shift, "-")
      if (cov_transform) {
        draws_T1_T2inv <- tcrossprod(draws_T1_T2inv, solve(total_mapping2))
      }
      draws_T1_T2inv <- sweep(draws_T1_T2inv, 2, total_scaling2, "/")
      draws_T1_T2inv <- sweep(draws_T1_T2inv, 2, mean_original + total_shift - total_shift2, "+")




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




      log_prob_trans_inv1 <- log_prob_prop_fun(draws = draws_trans_inv1, ...)
      log_prob_trans_inv2 <- log_prob_prop_fun(draws = draws_trans_inv2, ...)




      if (is.null(log_prob_target_fun) && is.null(log_ratio_fun)) {
        log_prob_prop_trans <- log_prob_prop_fun(draws = draws_trans, ...)
        lw_trans <-  log_prob_prop_trans -
          log(
            exp(log_prob_trans_inv1 - log(prod(total_scaling2))  - log(det(total_mapping2))) +
              exp(log_prob_trans_inv2 - log(prod(total_scaling))  - log(det(total_mapping)))
          )
      } else if (!is.null(log_prob_target_fun)) {
        log_prob_target_trans <- log_prob_target_fun(draws = draws_trans, ...)
        lw_trans <-  log_prob_target_trans -
          log(
            exp(log_prob_trans_inv1 - log(prod(total_scaling2))  - log(det(total_mapping2))) +
              exp(log_prob_trans_inv2 - log(prod(total_scaling))  - log(det(total_mapping)))
          )
      } else if (!is.null(log_ratio_fun)) {
        log_prob_ratio_trans <- log_ratio_fun(draws = draws_trans, ...)
        log_prob_prop_trans <- log_prob_prop_fun(draws = draws_trans, ...)
        lw_trans <-  log_prob_ratio_trans + log_prob_prop_trans -
          log(
            exp(log_prob_trans_inv1 - log(prod(total_scaling2))  - log(det(total_mapping2))) +
              exp(log_prob_trans_inv2 - log(prod(total_scaling))  - log(det(total_mapping)))
          )
      }










      # log_prob_target_trans <- log_prob_target_fun(draws = draws_trans, ...)
      #
      # lw_trans <-  log_prob_target_trans -
      #   log(
      #     exp(log_prob_trans_inv1 - log(prod(total_scaling2))  - log(det(total_mapping2))) +
      #       exp(log_prob_trans_inv2 - log(prod(total_scaling))  - log(det(total_mapping)))
      #   )


      # lw_trans_psis <- suppressWarnings(loo::psis(lw_trans))
      # lw_trans <- as.vector(weights(lw_trans_psis))
      lw_trans <- lw_trans - matrixStats::logSumExp(lw_trans)


      # replace draws and lw
      lw <- lw_trans
      draws <- draws_trans



    }
    else {
      # if not splitting, warn about high pareto ks
      if (any(kf > k_threshold)) {
        warning('Importance sampling may be unreliable. Consider setting split to TRUE.')
      }
    }

    if (log_expectation_fun) {
      expectation <- exp(matrixStats::colLogSumExps(lw + expectation_fun(draws, ...)))
    }
    else {
      w <- exp(lw)
      expectation <- colSums(w * expectation_fun(draws, ...))
    }

    list("expectation" = expectation, "pareto_k" = k, "pareto_kf" = kf, "draws" = draws, "log_weights" = lw)
  }







}


