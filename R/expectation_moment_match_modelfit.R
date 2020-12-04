

expectation_moment_match_modelfit <- function(x,
                                              post_draws_fun,
                                              unconstrain_pars_fun,
                                              log_prob_prop_draws_fun,
                                              expectation_fun,
                                              log_expectation_fun = FALSE,
                                              log_prob_target_draws_fun = NULL,
                                              log_ratio_draws_fun = NULL,
                                              obs_weights = NULL,
                                              log_lik_fun = NULL,
                                              k_threshold = 0.5,
                                              cov_transform = TRUE,
                                              split = FALSE,
                                              restart_transform = FALSE,
                                              ...) {

  checkmate::assertFunction(post_draws_fun)
  checkmate::assertFunction(unconstrain_pars_fun)
  checkmate::assertFunction(log_prob_prop_draws_fun)
  checkmate::assertFunction(expectation_fun)
  checkmate::assertNumber(k_threshold)
  checkmate::assertLogical(cov_transform)
  checkmate::assertLogical(split)
  checkmate::assertLogical(restart_transform)
  checkmate::assertLogical(log_expectation_fun)



  if (!is.null(log_prob_target_draws_fun) + !is.null(log_ratio_draws_fun) + !is.null(obs_weights) > 1) {
    stop("You canonly give one of the three options: log_prob_target_draws_fun,
         log_ratio_draws_fun and obs_weights.")
  }

  pars <- post_draws_fun(x, ...)
  # transform the model parameters to unconstrained space
  draws <- unconstrain_pars_fun(x, pars = pars, ...)

  orig_log_prob_prop <- log_prob_prop_draws_fun(draws = draws, ...)

  npars <- ncol(draws)
  S <- nrow(draws)
  cov_transform <- cov_transform && S >= 10 * npars

  # initialize objects that keep track of the total transformation
  total_shift <- rep(0, npars)
  total_scaling <- rep(1, npars)
  total_mapping <- diag(npars)

  draws_orig <- draws

  if (is.null(log_prob_target_draws_fun) && is.null(log_ratio_draws_fun)) {
    lw <- - matrixStats::logSumExp(rep(0,nrow(draws)))
    k <- 0

    lw_orig <- lw
  } else {
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
    if (!is.null(obs_weights)) {
      if (is.null(log_lik_fun)) {
        stop("You must give log_lik_fun when using obs_weights.")
      }
      log_lik <- log_lik_fun(x, ...)
      S <- nrow(log_lik)
      N <- ncol(log_lik)
      checkmate::assertNumeric(obs_weights, len = N)
      obs_weights <- matrix(c(obs_weights), S, N, byrow = TRUE)
      lw <- rowSums((obs_weights - 1) * log_lik)

      update_quantities <- update_quantities_obs_weights_modelfit
      density_function_list <- list(log_lik_fun = log_lik_fun, log_prob_prop_draws_fun = log_prob_prop_draws_fun, obs_weights = obs_weights)
    }

    lw_psis <- suppressWarnings(loo::psis(lw))
    lw <- as.vector(weights(lw_psis))
    k <- lw_psis$diagnostics$pareto_k

    lw_orig <- lw



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
        total_shift <- total_shift + trans$shift

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
          x = x,
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

    if (is.null(log_prob_target_draws_fun) && is.null(log_ratio_draws_fun) && is.null(obs_weights)) {
      update_quantities_expectation <- update_quantities_expectation_modelfit
      density_function_list <- list(log_prob_prop_draws_fun = log_prob_prop_draws_fun)
    } else if (!is.null(log_prob_target_draws_fun)) {
      update_quantities_expectation <- update_quantities_target_expectation_modelfit
      density_function_list <- list(log_prob_target_draws_fun = log_prob_target_draws_fun)
    } else if (!is.null(log_ratio_draws_fun)) {
      update_quantities_expectation <- update_quantities_ratio_expectation_modelfit
      density_function_list <- list(log_ratio_draws_fun = log_ratio_draws_fun, log_prob_prop_draws_fun = log_prob_prop_draws_fun)
    } else if (!is.null(obs_weights)) {
      update_quantities_expectation <- update_quantities_obs_weights_expectation_modelfit
      density_function_list <- list(log_prob_prop_draws_fun = log_prob_prop_draws_fun, log_lik_fun = log_lik_fun, obs_weights = obs_weights)
    }

    lwf <- as.vector(weights(psisf))

    while (kf > k_threshold) {

      # 1. match means
      trans <- shift(draws2, lwf)
      quantities <- update_quantities_expectation(
        x = x,
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
        x = x,
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
          x = x,
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




    log_prob_trans_inv1 <- log_prob_prop_draws_fun(x, draws = draws_trans_inv1, ...)
    log_prob_trans_inv2 <- log_prob_prop_draws_fun(x, draws = draws_trans_inv2, ...)




    if (is.null(log_prob_target_draws_fun) && is.null(log_ratio_draws_fun) && is.null(obs_weights)) {
      log_prob_prop_trans <- log_prob_prop_draws_fun(x, draws = draws_trans, ...)
      lw_trans <-  log_prob_prop_trans -
        log(
          exp(log_prob_trans_inv1 - log(prod(total_scaling2))  - log(det(total_mapping2))) +
            exp(log_prob_trans_inv2 - log(prod(total_scaling))  - log(det(total_mapping)))
        )
    } else if (!is.null(log_prob_target_draws_fun)) {
      log_prob_target_trans <- log_prob_target_draws_fun(x, draws = draws_trans, ...)
      lw_trans <-  log_prob_target_trans -
        log(
          exp(log_prob_trans_inv1 - log(prod(total_scaling2))  - log(det(total_mapping2))) +
            exp(log_prob_trans_inv2 - log(prod(total_scaling))  - log(det(total_mapping)))
        )
    } else if (!is.null(log_ratio_draws_fun)) {
      log_prob_ratio_trans <- log_ratio_draws_fun(x, draws = draws_trans, ...)
      log_prob_prop_trans <- log_prob_prop_draws_fun(x, draws = draws_trans, ...)
      lw_trans <-  log_prob_ratio_trans + log_prob_prop_trans -
        log(
          exp(log_prob_trans_inv1 - log(prod(total_scaling2))  - log(det(total_mapping2))) +
            exp(log_prob_trans_inv2 - log(prod(total_scaling))  - log(det(total_mapping)))
        )
    } else if (!is.null(obs_weights)) {
      log_lik_trans <- log_lik_fun(x, draws = draws_trans, ...)
      lw_trans <-  rowSums((obs_weights - 1) * log_lik_trans) -
        log(
          exp(log_prob_trans_inv1 - log(prod(total_scaling2))  - log(det(total_mapping2))) +
            exp(log_prob_trans_inv2 - log(prod(total_scaling))  - log(det(total_mapping)))
        )
    }









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


























update_quantities_target_expectation_modelfit <- function(x, draws, orig_log_prob_prop,
                                                          density_function_list,
                                                          expectation_fun, log_expectation_fun,
                                                          ...) {

  log_prob_target_draws_fun <- density_function_list$log_prob_target_draws_fun

  log_prob_target_new <- log_prob_target_draws_fun(x, draws = draws, ...)
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

update_quantities_ratio_expectation_modelfit <- function(x, draws, orig_log_prob_prop,
                                                         density_function_list,
                                                         expectation_fun, log_expectation_fun,
                                                         ...) {

  log_ratio_draws_fun <- density_function_list$log_ratio_draws_fun
  log_prob_prop_draws_fun <- density_function_list$log_prob_prop_draws_fun

  log_ratio_new <- log_ratio_draws_fun(x, draws = draws, ...)
  log_prop_new <- log_prob_prop_draws_fun(x, draws = draws, ...)

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


update_quantities_expectation_modelfit <- function(x, draws, orig_log_prob_prop,
                                                   density_function_list,
                                                   expectation_fun, log_expectation_fun,
                                                   ...) {

  log_prob_prop_draws_fun <- density_function_list$log_prob_prop_draws_fun

  log_prop_new <- log_prob_prop_draws_fun(x, draws = draws, ...)
  lw_new <- log_prop_new - orig_log_prob_prop

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















# this is here because it is not used in expectation_moment_match
update_quantities_obs_weights_expectation_modelfit <- function(x, draws, orig_log_prob_prop,
                                                density_function_list,
                                                expectation_fun, log_expectation_fun,
                                                ...) {

  log_prob_prop_draws_fun <- density_function_list$log_prob_prop_draws_fun
  log_lik_fun <- density_function_list$log_lik_fun
  obs_weights <- density_function_list$obs_weights

  log_prob_prop_new <- log_prob_prop_draws_fun(x, draws = draws, ...)
  log_lik_new <- log_lik_fun(x, draws = draws, ...)

  lw_new <- rowSums((obs_weights - 1) * log_lik_new) + log_prob_prop_new - orig_log_prob_prop

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


