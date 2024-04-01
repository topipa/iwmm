

# the following functions are copied from brms

# wrapper around rstan::constrain_pars
# ensures that the right posterior draws are excluded
.constrain_pars <- function(upars, x) {
  out <- rstan::constrain_pars(upars, object = x$fit)
  out[x$exclude] <- NULL
  out
}

.update_pars <- function(x, upars, ...) {
  # list with one element per posterior draw
  pars <- apply(upars, 1, .constrain_pars, x = x)
  # select required parameters only
  pars <- lapply(pars, "[", x$fit@sim$pars_oi_old)
  # transform draws
  ndraws <- length(pars)
  pars <- unlist(pars)
  npars <- length(pars) / ndraws
  dim(pars) <- c(npars, ndraws)
  # add dummy 'lp__' draws
  pars <- rbind(pars, rep(0, ndraws))
  # bring draws into the right structure
  new_draws <- named_list(x$fit@sim$fnames_oi_old, list(numeric(ndraws)))
  if (length(new_draws) != nrow(pars)) {
    stop("Updating parameters in 'brmsfit' failed. ")
  }
  for (i in seq_len(npars)) {
    new_draws[[i]] <- pars[i, ]
  }

  # create dummy sampler_params for new sim object
  newsamples <- list(new_draws)
  attr(newsamples[[1]], "sampler_params") <- .create_dummy_sampler_params(x)

  # create new sim object to overwrite x$fit@sim
  x$fit@sim <- list(
    samples = newsamples,
    iter = ndraws,
    thin = 1,
    warmup = 0,
    chains = 1,
    n_save = ndraws,
    warmup2 = 0,
    permutation = list(seq_len(ndraws)),
    pars_oi = x$fit@sim$pars_oi_old,
    dims_oi = x$fit@sim$dims_oi_old,
    fnames_oi = x$fit@sim$fnames_oi_old,
    n_flatnames = length(x$fit@sim$fnames_oi_old)
  )
  x$fit@stan_args <- list(
    list(chain_id = 1, iter = ndraws, thin = 1, warmup = 0)
  )

  brms::rename_pars(x)
}

.create_dummy_sampler_params <- function(x) {
  params <- attr(x$fit@sim$samples[[1]], "sampler_params")
  newparams <- params
  for (i in seq_along(params)) {
    newparams[[i]] <- numeric()
  }
  newparams
}


# update .MISC environment of the stanfit object
# allows to call log_prob and other C++ using methods
# on objects not created in the current R session
# or objects created via another backend
# update_misc_env <- function(x, recompile = FALSE, only_windows = FALSE) {
#   stopifnot(is.brmsfit(x))
#   recompile <- as_one_logical(recompile)
#   only_windows <- as_one_logical(only_windows)
#   if (recompile || !has_rstan_model(x)) {
#     x <- add_rstan_model(x, overwrite = TRUE)
#   } else if (os_is_windows() || !only_windows) {
#     # TODO: detect when updating .MISC is not required
#     # TODO: find a more efficient way to update .MISC
#     old_backend <- x$backend
#     x$backend <- "rstan"
#     x$fit@.MISC <- suppressMessages(brm(fit = x, chains = 0))$fit@.MISC
#     x$backend <- old_backend
#   }
#   x
# }
