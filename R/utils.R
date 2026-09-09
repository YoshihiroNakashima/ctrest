# -------------------------------------------------------------------------
# Internal Helper Functions (Do not export)
# -------------------------------------------------------------------------

#' Create all combinations of model terms
#' @noRd
.full_terms <- function(x) {
  if (length(x) == 0) return(list(stats::as.formula("~ 1")))
  combos <- unlist(lapply(0:length(x), function(i) utils::combn(x, i, simplify = FALSE)), recursive = FALSE)
  lapply(combos, function(vars) {
    if (length(vars) == 0) stats::as.formula("~ 1")
    else stats::as.formula(paste("~", paste(vars, collapse = " + ")))
  })
}

# Print section header
#' @noRd
.section <- function(title) {
  cat(sprintf("\n--- %s ---\n", title))
}

# Check MCMC convergence and print result
#' @noRd
.check_convergence <- function(df, threshold = 1.1) {
  if (is.null(df) || !"Rhat" %in% colnames(df)) {
    cat("  Rhat information not available.\n")
    return(invisible(NULL))
  }
  rhat    <- df$Rhat
  var_col <- if ("Variable" %in% colnames(df)) df$Variable else
             if ("variable" %in% colnames(df)) df$variable else
             as.character(seq_len(nrow(df)))
  bad     <- which(!is.na(rhat) & rhat > threshold)
  if (length(bad) > 0) {
    cat(sprintf("  Warning: %d parameter(s) with Rhat > %.1f:\n", length(bad), threshold))
    for (i in bad) {
      neff <- if ("n.eff" %in% colnames(df)) as.integer(df$n.eff[i]) else NA_integer_
      cat(sprintf("    %-35s Rhat = %.3f", var_col[i], rhat[i]))
      if (!is.na(neff)) cat(sprintf("  n.eff = %d", neff))
      cat("\n")
    }
    cat("  -> Consider increasing 'iter' or 'chains'.\n")
  } else {
    cat(sprintf("  All %d monitored parameter(s): Rhat <= %.1f.\n",
                sum(!is.na(rhat)), threshold))
  }
  invisible(NULL)
}

# Interpret a Bayesian p-value as a short string
#' @noRd
.p_value_note <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return("(not available)")
  if (p < 0.05 || p > 0.95) return("(model fit appears poor)")
  if (p < 0.10 || p > 0.90) return("(model fit is marginal)")
  "(model fit is adequate)"
}
