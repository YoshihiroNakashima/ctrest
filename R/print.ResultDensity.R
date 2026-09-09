#' Print method for ResultDensity objects
#'
#' Displays a formatted summary of the output from \code{bayes_rest}, including
#' model comparison (WAIC), posterior density estimates, and convergence diagnostics.
#'
#' @param x An object of class \code{"ResultDensity"}, as returned by \code{bayes_rest}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#' @export
print.ResultDensity <- function(x, ...) {
  species_label <- if (!is.null(x$target_species)) x$target_species else "Unknown"
  model_label   <- if (!is.null(x$model))          x$model          else "Unknown"
  family_label  <- if (!is.null(x$stay_family))    x$stay_family    else "Unknown"

  cat("\n=== ctrest: Density estimation ===\n")
  cat(sprintf("Species    : %s\n", species_label))
  cat(sprintf("Model      : %s\n", model_label))
  cat(sprintf("Stay family: %s\n", family_label))

  # --- Model comparison (WAIC) -----------------------------------------------

  .section("Model comparison (WAIC)")
  waic_df  <- x$WAIC
  best_row <- which.min(waic_df$WAIC)
  mark     <- rep("", nrow(waic_df))
  mark[best_row] <- "<- best"
  print_df <- cbind(waic_df, Note = mark)
  print(print_df, row.names = FALSE, right = FALSE)

  # --- Posterior estimates ----------------------------------------------------

  .section("Posterior estimates (best model)")
  sr       <- x$summary_result
  max_rows <- 15L
  if (!is.null(sr) && nrow(sr) > 0) {
    if (nrow(sr) > max_rows) {
      print(utils::head(sr, max_rows), row.names = FALSE)
      cat(sprintf(
        "  ... [%d more row(s)]. Use 'print(x$summary_result)' to see all.\n",
        nrow(sr) - max_rows
      ))
    } else {
      print(sr, row.names = FALSE)
    }
  }

  # --- Convergence -----------------------------------------------------------

  .section("Convergence")
  .check_convergence(x$summary_result)

  # --- Notes -----------------------------------------------------------------

  cat("\nNote: Full MCMC samples : $samples\n")
  cat("      Long-format samples: $tidy_samples\n")
  if (!is.null(x$activity_curve)) {
    cat("      Activity curve     : $activity_curve\n")
  }
  cat("      Trace plots        : MCMCvis::MCMCtrace(x$samples)\n")
  cat("\n")

  invisible(x)
}
