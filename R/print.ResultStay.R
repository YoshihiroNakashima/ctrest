#' Custom print method for ResultStay objects
#'
#' @param x An object of class 'ResultStay'
#' @param ... Additional arguments passed to print
#' @keywords internal
#' @export
#' @noRd
print.ResultStay <- function(x, ...) {
  cat("\n")
  cat("===== WAIC =====\n")
  print(x$WAIC)
  cat("===== Bayesian p-value for the best model =====\n")
  print(x$Bayesian_p_value)
  cat("===== Summary of the mean staying time estimates =====\n")
  print(x$summary_result)
  cat("\n(Note) MCMC samples for each parameter can be accessed via `$samples`.\n")
  cat("You can calculate Rhat values and generate trace plots using the MCMCvis package.\n")
}
