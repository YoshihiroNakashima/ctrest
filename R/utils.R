# -------------------------------------------------------------------------
# Internal Helper Functions (Do not export)
# -------------------------------------------------------------------------

#' Create all combinations of model terms
#' @param x A character vector of predictor names
#' @return A list of formula objects
#' @noRd
.full_terms <- function(x) {
  if (length(x) == 0) return(list(stats::as.formula("~ 1")))

  combos <- unlist(lapply(0:length(x), function(i) utils::combn(x, i, simplify = FALSE)), recursive = FALSE)
  lapply(combos, function(vars) {
    if (length(vars) == 0) {
      stats::as.formula("~ 1")
    } else {
      stats::as.formula(paste("~", paste(vars, collapse = " + ")))
    }
  })
}
