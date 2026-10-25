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
  total <- list(
    shift = total_shift,
    scaling = total_scaling,
    mapping = total_mapping
  )
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
      total <- add_step_to_total(total, trans, cov_transform)

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
      total <- add_step_to_total(total, trans, cov_transform)

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
        total <- add_step_to_total(total, trans, cov_transform)

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
    "total_shift" = total$shift,
    "total_scaling" = total$scaling,
    "total_mapping" = total$mapping
  )
}

#' Add one accepted moment matching step to the total transformation
#'
#' Each step maps the current draws `u` to `L (u - mean(u)) + mean(u) + shift`,
#' where `L` is the identity, `diag(scaling)` or `mapping`. The shifts add up,
#' and the linear parts compose in the order in which the steps were accepted.
#' After any accepted step `transform_loop()` starts again with the shift
#' step, so a scaling step can follow a covariance step. With
#' `cov_transform = TRUE` the scaling is therefore applied to the rows of
#' `mapping` (`diag(scaling) %*% mapping`), and `scaling` stays one. With
#' `cov_transform = FALSE` there is no covariance step, and the scalings are
#' kept in `scaling`.
#'
#' @noRd
#' @param total A list with `shift`, `scaling` and `mapping`.
#' @param trans The output of `shift()`, `shift_and_scale()` or
#'   `shift_and_cov()`.
#' @param cov_transform Logical; whether covariance steps are used.
#' @return The updated `total`.
add_step_to_total <- function(total, trans, cov_transform) {
  total$shift <- total$shift + trans$shift
  if (!is.null(trans$scaling)) {
    if (cov_transform) {
      total$mapping <- trans$scaling * total$mapping
    } else {
      total$scaling <- total$scaling * trans$scaling
    }
  }
  if (!is.null(trans$mapping)) {
    total$mapping <- trans$mapping %*% total$mapping
  }
  total
}
