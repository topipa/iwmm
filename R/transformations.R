#' Shift a matrix of draws to their weighted mean.
#'
#' @param draws A matrix of draws.
#' @param lw A vector representing the log-weight of each draw.
#' @return List with the shift that was performed, and the new draws matrix.
#'
#' @export
shift <- function(draws, lw) {
  checkmate::assertMatrix(draws, any.missing = FALSE)
  checkmate::assertNumeric(lw)
  checkmate::assertSetEqual(nrow(draws), length(lw))

  w <- exp(lw) * length(lw)
  # compute moments using log weights
  mean_original <- colMeans(draws)
  mean_weighted <- matrixStats::colWeightedMeans(draws, w = w)

  shift <- mean_weighted - mean_original
  # transform posterior draws
  draws_new <- sweep(draws, 2, shift, "+")
  list(
    draws = draws_new,
    shift = shift
  )
}


#' Shift a matrix of draws to their weighted mean and scale the marginal
#' variances to match the weighted marginal variances.
#'
#' @param draws A matrix of draws.
#' @param lw A vector representing the log-weight of each draw.
#' @return List with the shift and scaling that were performed, and the new
#' draws matrix.
#'
#' @export
shift_and_scale <- function(draws, lw) {
  checkmate::assertMatrix(draws, any.missing = FALSE)
  checkmate::assertNumeric(lw)
  checkmate::assertSetEqual(nrow(draws), length(lw))

  w <- exp(lw) * length(lw)
  # compute moments using log weights
  vars_original <- matrixStats::colVars(draws)
  vars_weighted <- matrixStats::colWeightedVars(draws, w = w)
  if (all(vars_original > 1e-12) && all(vars_weighted > 1e-12)) {
    scaling <- sqrt(vars_weighted / vars_original)

    mean_original <- colMeans(draws)
    mean_weighted <- matrixStats::colWeightedMeans(draws, w = w)
    shift <- mean_weighted - mean_original

    draws_new <- sweep(draws, 2, mean_original, "-")
    draws_new <- sweep(draws_new, 2, scaling, "*")
    draws_new <- sweep(draws_new, 2, mean_weighted, "+")
  } else {
    draws_new <- draws
    scaling <- NA
  }

  list(
    draws = draws_new,
    shift = shift,
    scaling = scaling
  )
}

#' Shift a matrix of draws to their weighted mean and scale the covariance
#' to match the weighted covariance.
#'
#' @param draws A matrix of draws.
#' @param lw A vector representing the log-weight of each draw.
#' @return List with the shift and mapping that were performed, and the new
#' draws matrix.
#'
#' @export
shift_and_cov <- function(draws, lw) {
  checkmate::assertMatrix(draws, any.missing = FALSE)
  checkmate::assertNumeric(lw)
  checkmate::assertSetEqual(nrow(draws), length(lw))

  w <- exp(lw) * length(lw)
  # compute moments using log weights
  covar_original <- stats::cov(draws)
  covar_weighted <- stats::cov.wt(draws, wt = w)$cov
  chol1 <- tryCatch(
    {
      chol(covar_weighted)
    },
    error = function(cond) {
      return(NULL)
    }
  )
  chol2 <- tryCatch(
    {
      chol(covar_original)
    },
    error = function(cond) {
      return(NULL)
    }
  )
  if (is.null(chol1) || is.null(chol2)) {
    draws_new <- draws
    mapping <- diag(ncol(draws))
    shift <- rep(0, ncol(draws))
  } else {
    mapping <- t(chol1) %*% solve(t(chol2))

    mean_original <- colMeans(draws)
    mean_weighted <- matrixStats::colWeightedMeans(draws, w = w)
    shift <- mean_weighted - mean_original

    # transform posterior draws
    draws_new <- sweep(draws, 2, mean_original, "-")
    draws_new <- tcrossprod(draws_new, mapping)
    draws_new <- sweep(draws_new, 2, mean_weighted, "+")
    colnames(draws_new) <- colnames(draws)
  }

  list(
    draws = draws_new,
    shift = shift,
    mapping = mapping
  )
}

transform_loop <- function(draws,
                           lw,
                           k,
                           update_properties,
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
      update_properties,
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
      update_properties,
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
        update_properties,
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

  list(
    "draws" = draws,
    "log_weights" = lw,
    "pareto_k" = k,
    "total_shift" = total_shift,
    "total_scaling" = total_scaling,
    "total_mapping" = total_mapping
  )
}
