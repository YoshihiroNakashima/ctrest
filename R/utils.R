# -------------------------------------------------------------------------
# Internal Helper Functions (Do not export)
# -------------------------------------------------------------------------

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

# Design matrix of the staying-time model for the stations (one row per station)
#
# The staying time is observed for each video, but the density formula needs the mean staying
# time of each station. X_obs is the design matrix of the staying time records; the matrix of the
# stations has the same columns (same factor levels and contrasts), computed from the covariates
# of the station data.
#' @noRd
.station_design_matrix <- function(formula, model_frame_obs, X_obs, station_data) {
  terms_rhs <- stats::delete.response(stats::terms(stats::as.formula(formula)))
  vars <- all.vars(terms_rhs)
  miss <- setdiff(vars, names(station_data))
  if (length(miss) > 0) {
    stop("The covariate(s) of `formula_stay` must also be columns of `station_effort_data`, ",
         "because the mean staying time of each station is used for the density: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  isF  <- vapply(model_frame_obs, function(z) is.factor(z) || is.character(z), logical(1))
  xlev <- lapply(model_frame_obs[isF], function(z) levels(factor(z)))
  # na.pass: rows with missing values must not be dropped (the matrix has one row per station)
  mf <- tryCatch(
    stats::model.frame(terms_rhs, data = station_data,
                       xlev = if (length(xlev) > 0) xlev else NULL, na.action = stats::na.pass),
    error = function(e) {
      stop("The covariates of `formula_stay` in `station_effort_data` do not match those of `stay_data`: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  X <- stats::model.matrix(terms_rhs, mf, contrasts.arg = attr(X_obs, "contrasts"))
  if (!identical(colnames(X), colnames(X_obs))) {
    stop("The design matrix of `formula_stay` for the stations (", paste(colnames(X), collapse = ", "),
         ") does not match that of `stay_data` (", paste(colnames(X_obs), collapse = ", "), ").",
         call. = FALSE)
  }
  if (anyNA(X)) {
    stop("The covariates of `formula_stay` have missing values in `station_effort_data` (stations: ",
         paste(utils::head(station_data$Station[!stats::complete.cases(X)], 5), collapse = ", "),
         "). The mean staying time of every station is needed for the density.", call. = FALSE)
  }
  rownames(X) <- NULL
  X
}

# Interpret a Bayesian p-value as a short string
#' @noRd
.p_value_note <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return("(not available)")
  if (p < 0.05 || p > 0.95) return("(model fit appears poor)")
  if (p < 0.10 || p > 0.90) return("(model fit is marginal)")
  "(model fit is adequate)"
}
