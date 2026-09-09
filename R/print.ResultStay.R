#' Print method for ResultStay objects
#'
#' Displays a formatted summary of the output from \code{bayes_stay_selection},
#' including model comparison (WAIC), posterior estimates of mean staying time,
#' convergence diagnostics, and a Bayesian p-value for the best model.
#'
#' @param x An object of class \code{"ResultStay"}, as returned by
#'   \code{bayes_stay_selection}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#' @export
print.ResultStay <- function(x, ...) {
  species_label <- if (!is.null(x$target_species)) x$target_species else "Unknown"
  family_label  <- if (!is.null(x$stay_family))   x$stay_family   else "Unknown"

  cat("\n=== ctrest: Staying time model selection ===\n")
  cat(sprintf("Species    : %s\n", species_label))
  cat(sprintf("Stay family: %s\n", family_label))

  # --- Model comparison (WAIC) -----------------------------------------------

  .section("Model comparison (WAIC)")
  waic_df   <- x$WAIC
  best_row  <- which.min(waic_df$WAIC)
  mark      <- rep("", nrow(waic_df))
  mark[best_row] <- "<- best"
  print_df  <- cbind(waic_df, Note = mark)
  print(print_df, row.names = FALSE, right = FALSE)

  best_model <- as.character(waic_df$Model[best_row])
  re_label   <- if ("Random_effect" %in% colnames(waic_df)) {
    as.character(waic_df$Random_effect[best_row])
  } else "NULL"
  cat(sprintf(
    "\nBest model: %s  (Random effect: %s)\n",
    best_model, re_label
  ))

  # --- Mean staying time estimates -------------------------------------------

  .section("Mean staying time (best model)")
  sr <- x$summary_result
  if (!is.null(sr) && nrow(sr) > 0) {
    print(sr, row.names = FALSE)
  }

  # --- Convergence -----------------------------------------------------------

  .section("Convergence")
  .check_convergence(x$summary_result)

  # --- Bayesian p-value ------------------------------------------------------

  .section("Bayesian p-value (best model)")
  p <- x$Bayesian_p_value
  if (!is.null(p) && !is.na(p)) {
    cat(sprintf("  %.3f  %s\n", p, .p_value_note(p)))
  } else {
    cat("  Not available.\n")
  }

  # --- Notes -----------------------------------------------------------------

  cat("\nNote: Full MCMC samples : $samples\n")
  cat("      Long-format samples: $tidy_samples\n")
  cat("      Trace plots        : MCMCvis::MCMCtrace(x$samples)\n")
  cat("      Full summary table : print(x$summary_result, n = Inf)\n")
  cat("\n")

  invisible(x)
}
