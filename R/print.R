##' @export
print.adapted_importance_sampling <- function(x, ...) {
  cat("Adapted importance sampling\n")
  if (!any(is.na(x$expectation))) {
    cat("Expectation:\n")
    print(round(x$expectation, 2))
  }
  cat("Draws:\n")
  print(x$draws)

  invisible(x)
}
