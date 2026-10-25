test_that("transformations work", {
  draws <- matrix(c(
    1.75, -1.07, -1.27, 0.43, 1.91,
    2.26, -0.13, -0.14
  ), 4, 2)

  lw <- dnorm(draws[, 1], 1, 1, log = TRUE) +
    dnorm(draws[, 2], 1, 1, log = TRUE)
  lw <- lw - matrixStats::logSumExp(lw)


  trans1 <- shift(draws, lw)
  expect_equal(trans1$draws,
    matrix(c(
      2.7129633, -0.1070367, -0.3070367, 1.3929633,
      1.905564, 2.255564, -0.134436, -0.144436
    ), 4, 2),
    tolerance = 1e-7
  )
  expect_equal(trans1$shift, c(0.962963269, -0.004435991),
    tolerance = 1e-7
  )

  trans2 <- shift_and_scale(draws, lw)
  expect_equal(trans2$draws,
    matrix(
      c(
        2.24928499, 0.15977256, 0.01158019, 1.27121534,
        1.84278570, 2.16928580, -0.06024344, -0.06957202
      ),
      4, 2
    ),
    tolerance = 1e-7
  )
  expect_equal(trans2$shift, c(0.962963269, -0.004435991), tolerance = 1e-7)
  expect_equal(trans2$scaling, c(0.7409619, 0.9328574), tolerance = 1e-7)

  trans3 <- shift_and_cov(draws, lw)
  expect_equal(trans3$draws, matrix(c(
    2.43067094, 0.05539964,
    -0.11305932, 1.31884182,
    2.5189668, 1.6583828,
    -0.4939763, 0.1988827
  ), 4, 2),
  tolerance = 1e-7
  )
  expect_equal(trans3$shift, c(0.962963269, -0.004435991), tolerance = 1e-7)
  expect_equal(trans3$mapping,
    matrix(c(0.8422948, 0.4126585, 0.0000000, 0.8660365), 2, 2),
    tolerance = 1e-7
  )
})

test_that("the total transformation composes the accepted steps in order", {
  # transform_loop() starts again with the shift step after every accepted
  # step, so a scaling step can follow a covariance step. The total
  # transformation T(x) = M diag(s) (x - m) + m + shift, which the split step
  # rebuilds from the totals (M only with cov_transform = TRUE), must map the
  # original draws to the draws that the steps produced one after another.
  set.seed(7)
  S <- 400
  draws <- cbind(rnorm(S), 0.5 * rnorm(S) + 0.3 * rnorm(S))
  lw <- dnorm(draws[, 1], 0.5, 0.7, log = TRUE) +
    dnorm(draws[, 2], -0.2, 0.6, log = TRUE)
  lw <- lw - matrixStats::logSumExp(lw)
  m <- colMeans(draws)
  sequences <- list(
    list(cov = TRUE, steps = c("cov", "scale", "shift", "cov", "scale")),
    list(cov = FALSE, steps = c("scale", "shift", "scale"))
  )
  for (sq in sequences) {
    total <- list(shift = c(0, 0), scaling = c(1, 1), mapping = diag(2))
    u <- draws
    for (step in sq$steps) {
      trans <- switch(step,
        shift = shift(u, lw),
        scale = shift_and_scale(u, lw),
        cov = shift_and_cov(u, lw)
      )
      u <- trans$draws
      total <- add_step_to_total(total, trans, sq$cov)
    }
    rebuilt <- sweep(draws, 2, m, "-")
    rebuilt <- sweep(rebuilt, 2, total$scaling, "*")
    if (sq$cov) {
      rebuilt <- tcrossprod(rebuilt, total$mapping)
    }
    rebuilt <- sweep(rebuilt, 2, total$shift + m, "+")
    expect_equal(unname(rebuilt), unname(u))
  }
})

test_that("transformation assertions work", {
  # only 1 argument
  expect_error(shift(matrix(0, 2, 2)))
  # not matrix
  expect_error(shift(c(2, 3), c(1, 2)))
  # dimension mismatch
  expect_error(shift(matrix(0, 2, 2), c(1, 2, 3)))
})
