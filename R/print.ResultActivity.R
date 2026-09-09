#' Display a summary and optional activity curve for a ResultActivity object
#'
#' Prints a formatted summary of the output from \code{bayes_activity}, including
#' the estimated activity proportion, convergence diagnostics, and a Bayesian
#' p-value. When \code{plot = TRUE}, draws the posterior activity density curve
#' together with a kernel density estimate from the \code{activity} package.
#'
#' @param x An object of class \code{"ResultActivity"}, as returned by
#'   \code{bayes_activity}.
#' @param plot Logical; if \code{TRUE} (default), the activity density curve is
#'   plotted.
#' @param bw_adj A positive numeric bandwidth adjustment multiplier passed to
#'   the kernel density estimator. Default is \code{1.0}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#' @export
printResultActivity <- function(x, plot = TRUE, bw_adj = 1.0, ...) {
  if (!inherits(x, "ResultActivity"))
    stop("'x' must be a ResultActivity object returned by bayes_activity().", call. = FALSE)
  if (!is.logical(plot) || length(plot) != 1)
    stop("'plot' must be a single logical value (TRUE or FALSE).", call. = FALSE)
  if (!is.numeric(bw_adj) || length(bw_adj) != 1 || bw_adj <= 0)
    stop("'bw_adj' must be a single positive number.", call. = FALSE)

  species_label <- if (!is.null(x$target_species)) x$target_species else "Unknown"

  cat("\n=== ctrest: Activity estimation ===\n")
  cat(sprintf("Species: %s\n", species_label))

  # --- Activity proportion ----------------------------------------------------

  .section("Activity proportion")
  sr <- x$summary_result
  if (!is.null(sr) && nrow(sr) > 0) {
    print(sr, row.names = FALSE)
  }

  # --- Convergence -----------------------------------------------------------

  .section("Convergence")
  .check_convergence(x$summary_result)

  # --- Bayesian p-value ------------------------------------------------------

  .section("Bayesian p-value")
  p <- x$Bayesian_p_value
  if (!is.null(p) && !is.na(p)) {
    cat(sprintf("  %.3f  %s\n", p, .p_value_note(p)))
  } else {
    cat("  Not available.\n")
  }

  # --- Notes -----------------------------------------------------------------

  cat("\nNote: Full MCMC samples : $samples\n")
  cat("      Long-format samples: $tidy_samples\n")
  cat("      Activity curve data: $activity_curve\n")
  cat("      Trace plots        : MCMCvis::MCMCtrace(x$samples)\n")
  cat("\n")

  # --- Plot ------------------------------------------------------------------

  if (plot) {
    curve_df <- x$activity_curve
    act_vec  <- x$act_data

    # Kernel density via activity package for comparison
    bw_base <- tryCatch(
      activity::bwcalc(act_vec),
      error = function(e) NULL
    )
    if (!is.null(bw_base)) {
      kd    <- activity::fitact(act_vec, bw = bw_base * bw_adj, adj = 1)
      kd_df <- data.frame(x = kd@xvals, pdf = kd@pdf)
    } else {
      kd_df <- NULL
    }

    p_plot <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(
        data    = curve_df,
        mapping = ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper),
        fill    = "steelblue",
        alpha   = 0.3
      ) +
      ggplot2::geom_line(
        data    = curve_df,
        mapping = ggplot2::aes(x = .data$x, y = .data$mean),
        colour  = "steelblue",
        linewidth = 0.8
      ) +
      ggplot2::geom_rug(
        data    = data.frame(time = act_vec),
        mapping = ggplot2::aes(x = .data$time),
        sides   = "b",
        alpha   = 0.4
      ) +
      ggplot2::scale_x_continuous(
        name   = "Time of day (radians)",
        breaks = c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
        labels = c("0", "π/2", "π", "3π/2", "2π")
      ) +
      ggplot2::ylab("Activity density") +
      ggplot2::ggtitle(
        sprintf("Activity curve: %s", species_label),
        subtitle = "Posterior mean (blue) ± 95% CI (shaded); kernel density (dashed)"
      ) +
      ggplot2::theme_bw()

    if (!is.null(kd_df)) {
      p_plot <- p_plot +
        ggplot2::geom_line(
          data    = kd_df,
          mapping = ggplot2::aes(x = .data$x, y = .data$pdf),
          colour   = "grey40",
          linetype = "dashed",
          linewidth = 0.6
        )
    }

    print(p_plot)
  }

  invisible(x)
}

#' Print method for ResultActivity objects
#'
#' Calls \code{\link{printResultActivity}} with \code{plot = FALSE}. To display
#' the activity curve plot, use \code{printResultActivity(x, plot = TRUE)}.
#'
#' @param x An object of class \code{"ResultActivity"}.
#' @param ... Additional arguments passed to \code{printResultActivity}.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.ResultActivity <- function(x, ...) {
  printResultActivity(x, plot = FALSE, ...)
}
