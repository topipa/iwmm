#' Generic importance weighted moment matching algorithm.
#'
#' Matches a matrix of draws to its importance weighted moments, but
#' does not compute any expectation.
#'
#'
#' @param draws A matrix of draws.
#' @param log_prob_prop_draws_fun Log density of the proposal.
#' The function takes argument `draws`.
#' @param expectation_fun A function whose expectation is being computed.
#' The function takes arguments `draws`.
#' @param log_expectation_fun Logical indicating whether the expectation_fun
#' returns its values as logarithms or not. Defaults to FALSE. If set to TRUE,
#' the expectation function must be nonnegative (before taking the logarithm).
#' @param log_prob_target_draws_fun Log density of the target.
#' The function takes argument `draws`.
#' @param log_ratio_draws_fun Log of the density ratio (target/proposal).
#' The function takes argument `draws`.
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
#' @param ... Further arguments passed to `log_prob_prop_draws_fun`,
#' `log_prob_target_draws_fun` and `log_ratio_draws_fun`.
#'
#' @return Returns a list with 4 elements: the expectation, transformed draws, updated
#' importance weights, and the pareto k diagnostic value.
#'
#' @export
#' @importFrom stats weights
expectation_moment_match <- function(draws,
                         log_prob_prop_draws_fun,
                         expectation_fun,
                         log_expectation_fun = FALSE,
                         log_prob_target_draws_fun = NULL,
                         log_ratio_draws_fun = NULL,
                         k_threshold = 0.5,
                         cov_transform = TRUE,
                         split = FALSE,
                         restart_transform = FALSE,
                         ...) {

  checkmate::assertMatrix(draws, any.missing = FALSE)
  checkmate::assertFunction(log_prob_prop_draws_fun)
  checkmate::assertFunction(expectation_fun)
  checkmate::assertNumber(k_threshold)
  checkmate::assertLogical(cov_transform)
  checkmate::assertLogical(split)
  checkmate::assertLogical(restart_transform)
  checkmate::assertLogical(log_expectation_fun)


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

  # initialize objects that keep track of the total transformation
  total_shift <- rep(0, npars)
  total_scaling <- rep(1, npars)
  total_mapping <- diag(npars)

  draws_orig <- draws
  lw_orig <- lw

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

  if (log_expectation_fun) {
    lwf <- lw + expectation_fun(draws2, ...)
    # Here we expect that expectation_fun is nonnegative (before log)')
  }
  else {
    lwf <- lw + log(abs(expectation_fun(draws2, ...)))
  }
  psisf <- suppressWarnings(loo::psis(lwf))
  kf <- psisf$diagnostics$pareto_k

  if (split) {


    if (ncol(lwf) > 1) {
      stop('Using split = TRUE is not yet supported for expectation functions that return a matrix.
           As a workaround, you can wrap your function call using apply.')
    }

    if (!is.null(log_prob_target_draws_fun)) {
      update_quantities_expectation <- update_quantities_target_expectation
      density_function_list <- list(log_prob_target_draws_fun = log_prob_target_draws_fun)
    }
    if (!is.null(log_ratio_draws_fun)) {
      update_quantities_expectation <- update_quantities_ratio_expectation
      density_function_list <- list(log_ratio_draws_fun = log_ratio_draws_fun, log_prob_prop_draws_fun = log_prob_prop_draws_fun)
    }

    lwf <- as.vector(weights(psisf))

    while (kf > k_threshold) {

      # 1. match means
      trans <- shift(draws2, lwf)
      quantities <- update_quantities_expectation(
        draws = trans$draws,
        orig_log_prob_prop = orig_log_prob_prop,
        density_function_list = density_function_list,
        expectation_fun = expectation_fun,
        log_expectation_fun = log_expectation_fun,
        ...
      )
      if (quantities$kf < kf) {
        draws2 <- trans$draws
        total_shift2 <- total_shift2 + trans$shift

        lwf <- quantities$lwf
        kf <- quantities$kf
        next
      }

      # 2. match means and variances
      trans <- shift_and_scale(draws2, lwf)
      quantities <- update_quantities_expectation(
        draws = trans$draws,
        orig_log_prob_prop = orig_log_prob_prop,
        density_function_list = density_function_list,
        expectation_fun = expectation_fun,
        log_expectation_fun = log_expectation_fun,
        ...
      )
      if (quantities$kf < kf) {
        draws2 <- trans$draws
        total_shift2 <- total_shift2 + trans$shift
        total_scaling2 <- total_scaling2 * trans$scaling

        lwf <- quantities$lwf
        kf <- quantities$kf
        next
      }

      if (cov_transform) {
        # 3. match means and covariances
        trans <- shift_and_cov(draws2, lwf)
        quantities <- update_quantities_expectation(
          draws = trans$draws,
          orig_log_prob_prop = orig_log_prob_prop,
          density_function_list = density_function_list,
          expectation_fun = expectation_fun,
          log_expectation_fun = log_expectation_fun,
          ...
        )
        if (quantities$kf < kf) {
          draws2 <- trans$draws
          total_shift2 <- total_shift2 + trans$shift
          total_mapping2 <- trans$mapping %*% total_mapping2

          lwf <- quantities$lwf
          kf <- quantities$kf
          next
        }
      }



      break
    }

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

    log_prob_target_trans <- log_prob_target_draws_fun(draws = draws_trans, ...)

    log_prob_trans_inv1 <- log_prob_prop_draws_fun(draws = draws_trans_inv1, ...)
    log_prob_trans_inv2 <- log_prob_prop_draws_fun(draws = draws_trans_inv2, ...)

    lw_trans <-  log_prob_target_trans -
      log(
        exp(log_prob_trans_inv1 - log(prod(total_scaling2))  - log(det(total_mapping2))) +
          exp(log_prob_trans_inv2 - log(prod(total_scaling))  - log(det(total_mapping)))
      )


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



#' Function for updating importance weights and pareto k diagnostic for
#' expectation-specific weights.
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
                               expectation_fun, log_expectation_fun,
                               ...) {

  log_prob_target_draws_fun <- density_function_list$log_prob_target_draws_fun

  log_prob_target_new <- log_prob_target_draws_fun(draws = draws, ...)
  lw_new <- log_prob_target_new - orig_log_prob_prop

  if (log_expectation_fun) {
    lwf_new <- lw_new + expectation_fun(draws, ...)
  }
  else {
    lwf_new <- lw_new + log(abs(expectation_fun(draws, ...)))
  }

  psisf_new <- suppressWarnings(loo::psis(lwf_new))
  kf_new <- psisf_new$diagnostics$pareto_k
  lwf_new <- as.vector(weights(psisf_new))

  # gather results
  list(
    lwf = lwf_new,
    kf = kf_new
  )
}

update_quantities_ratio_expectation <- function(draws, orig_log_prob_prop,
                                                 density_function_list,
                                                 expectation_fun, log_expectation_fun,
                                                 ...) {

  log_ratio_draws_fun <- density_function_list$log_ratio_draws_fun
  log_prob_prop_draws_fun <- density_function_list$log_prob_prop_draws_fun

  log_ratio_new <- log_ratio_draws_fun(draws = draws, ...)
  log_prop_new <- log_prob_prop_draws_fun(draws = draws, ...)

  lw_new <- log_ratio_new + log_prop_new - orig_log_prob_prop

  if (log_expectation_fun) {
    lwf_new <- lw_new + expectation_fun(draws, ...)
  }
  else {
    lwf_new <- lw_new + log(abs(expectation_fun(draws, ...)))
  }

  psisf_new <- suppressWarnings(loo::psis(lwf_new))
  kf_new <- psisf_new$diagnostics$pareto_k
  lwf_new <- as.vector(weights(psisf_new))

  # gather results
  list(
    lwf = lwf_new,
    kf = kf_new
  )
}

