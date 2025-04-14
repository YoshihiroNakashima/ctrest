#' Custom print method for ResultDensity objects
#'
#' @param x An object of class 'ResultDensity'
#' @param ... Additional arguments passed to print
#' @keywords internal
#' @export
#' @noRd
print.ResultDensity <- function(x, ...) {
  cat("\n")
  cat("===== WAIC =====\n")
  print(x$WAIC)
  cat("===== Summary of the mean staying time estimates =====\n")
  print(x$summary_result)
  cat("\n(Note) MCMC samples for each parameter can be accessed via `$samples`.\n")
  cat("You can calculate Rhat values and generate trace plots using the MCMCvis package.\n")
}
