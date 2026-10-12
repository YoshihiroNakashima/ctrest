#' Predict density from the covariates of a fitted REST model
#'
#' Computes posterior summaries of the density predicted by the density
#' sub-model (\code{formula_density}) of \code{\link{bayes_rest}} or
#' \code{\link{bayes_rest_multi}} at arbitrary covariate values. Covariates are
#' given on their original scale; standardization used during estimation is
#' applied internally.
#'
#' The density model is \eqn{\log D = \eta + a_{camera} + a_{transect}}, where
#' \eqn{\eta} is the linear predictor (for \code{bayes_rest_multi}, including the
#' species-level random effects, i.e. the species-specific coefficients),
#' \eqn{a_{camera} \sim N(0, \sigma_{camera}^2)} is the camera-level random
#' effect and \eqn{a_{transect} \sim N(0, \sigma_{transect}^2)} is the
#' transect-level random effect (if \code{random_effect_density} was given).
#' Two kinds of density are returned, and the column \code{Type} tells which:
#' \describe{
#'   \item{\code{"typical"}}{The density of a \emph{typical} camera, i.e. with
#'     the camera-level random effect set to zero: \eqn{\exp(\eta)}. The
#'     transect-level effect is also zero unless \code{newdata} contains the
#'     transect column (see below). This is the definition of \code{density}
#'     in \code{$summary_result}. Because the random effects act on the
#'     multiplicative scale, a typical camera is below the average camera when
#'     cameras differ much (by a factor of \eqn{\exp(-\sigma^2/2)}).}
#'   \item{\code{"average"}}{The \emph{mean} density over cameras (and over
#'     transects when the transect is not given): \eqn{\exp(\eta +
#'     \sigma_{camera}^2/2 [+ \sigma_{transect}^2/2])}. Use this to obtain the
#'     mean density of an area (e.g., to multiply by the area to get the number
#'     of animals) or to compare densities with those of other studies. It
#'     assumes that the camera-level variation is a real variation in density
#'     among locations that are represented by the cameras.}
#' }
#' For comparing covariate values, years or transects, the two types give
#' similar ratios, provided the variation among cameras is similar; the
#' absolute level differs by \eqn{\exp(\sigma^2/2)}.
#'
#' The mean density of an area depends on the covariates of the whole area:
#' set \code{average = TRUE} to average the predicted density over the rows of
#' \code{newdata} (for each posterior draw, so that the uncertainty is carried
#' over), e.g., over the covariate values of all locations of the study area,
#' or of the cameras (the default when \code{newdata} is \code{NULL}). The mean
#' of the density predicted at the mean covariate value is not the mean density
#' when the effect of a covariate is nonlinear on the density scale. If the
#' cameras do not represent the area (e.g., they are placed in a small part of
#' it), the result is the mean density of the places they represent.
#'
#' Prediction uses the posterior draws stored in \code{object$samples}, so no
#' re-estimation is needed. Formulas with data-dependent transformations such as
#' \code{poly()} or \code{scale()} are not supported; use \code{I()} or
#' \code{log()} instead.
#'
#' @param object An object of class \code{"ResultDensity"} returned by
#'   \code{\link{bayes_rest}} or \code{\link{bayes_rest_multi}}. Objects created
#'   by older versions of ctrest do not contain the information required for
#'   prediction and must be re-estimated.
#' @param newdata A data frame containing every covariate used in
#'   \code{formula_density}, on the original scale. If the density model has
#'   random effects for transects, adding the column named in
#'   \code{random_effect_density} gives the prediction for each transect (the
#'   effect of the transect is added). If \code{NULL} (default) and
#'   \code{average = FALSE}, a profile for each covariate is returned: the
#'   covariate is varied over its observed values (or a 25-point grid between
#'   its minimum and maximum when it is continuous) while the other covariates
#'   are held at their most frequent value (discrete covariates) or mean
#'   (continuous covariates). The column \code{Profile} gives the covariate
#'   that was varied. If the density model has no covariates and \code{newdata}
#'   is \code{NULL}, one row per species is returned. If \code{NULL} and
#'   \code{average = TRUE}, the covariates and transects of the cameras used in
#'   the estimation are used.
#' @param species A character vector of species names to predict (only for
#'   \code{bayes_rest_multi}). Default is all species in the model.
#' @param prob Probability mass of the credible interval. Default is 0.95.
#' @param type \code{"both"} (default), \code{"typical"} or \code{"average"}:
#'   the kind of density to return. See Details.
#' @param average Logical. If \code{TRUE}, the predicted densities are averaged
#'   over the rows of \code{newdata} (separately for each posterior draw), and
#'   one row per species (and per group of \code{by}) is returned, with the
#'   number of rows averaged in the column \code{n}. All rows have the same
#'   weight; repeat rows to give a larger weight (e.g., proportional to area).
#' @param by A character vector of column names of \code{newdata}. With
#'   \code{average = TRUE}, the average is taken within each combination of
#'   these columns (e.g., \code{by = "Year"} or a habitat type).
#' @param draws Logical. If \code{TRUE}, the posterior draws of the predicted
#'   density are returned instead of their summaries (see \code{Value}). The
#'   draws are joint across species, rows and the two types: the same posterior
#'   draw (column \code{draw}) is used for all of them, so that quantities
#'   combining several species or transects (e.g., biomass, community
#'   composition, the mean density of a habitat type) can be computed draw by
#'   draw and summarized afterwards.
#' @param n_draws Number of posterior draws to return when \code{draws = TRUE}
#'   (a random subset, the same for all rows; use \code{set.seed()} for
#'   reproducibility). Default is \code{NULL}: all draws.
#' @param ... Currently unused.
#'
#' @return A tibble with the column \code{Species}, the columns of
#'   \code{newdata} (or the columns of \code{by} and \code{n} when
#'   \code{average = TRUE}), the column \code{Type} (\code{"typical"} or
#'   \code{"average"}) and the posterior \code{mean}, \code{sd}, \code{lower},
#'   \code{median} and \code{upper} of the predicted density (individuals per
#'   km^2). The meaning of \code{Type} is printed below the table. With
#'   \code{draws = TRUE}, the posterior summaries are replaced by the columns
#'   \code{draw} (the index of the posterior draw) and \code{density} (the
#'   density in that draw), with one row per species, row of \code{newdata}
#'   (or group), type and draw.
#'
#' @examples
#' \dontrun{
#' fit <- bayes_rest_multi(
#'   formula_stay = Stay ~ 1,
#'   formula_density = ~ 1 + Year,
#'   formula_enter = ~ 1,
#'   random_effect_density = "Transect",
#'   station_effort_data = station_effort_data,
#'   stay_data = stay_data,
#'   activity_data = activity_data,
#'   focal_area = 2.67,
#'   target_species = c("SP01", "SP02", "SP03")
#' )
#'
#' # Density of each species in each year: typical camera and average
#' predict(fit, newdata = data.frame(Year = c(0, 1)))
#'
#' # Density of each transect in each year
#' predict(fit, newdata = expand.grid(Transect = unique(station_effort_data$Transect),
#'                                    Year = c(0, 1)))
#'
#' # Mean density over the cameras (of each year)
#' predict(fit, average = TRUE, by = "Year", type = "average")
#'
#' # Posterior draws of the density of each transect, e.g., to compute the
#' # biomass of each habitat type draw by draw
#' set.seed(1)
#' d <- predict(fit, newdata = expand.grid(Transect = unique(station_effort_data$Transect),
#'                                         Year = c(0, 1)),
#'              type = "average", draws = TRUE, n_draws = 500)
#' }
#'
#' @importFrom stats predict
#' @export
predict.ResultDensity <- function(object, newdata = NULL, species = NULL, prob = 0.95,
                                  type = c("both", "typical", "average"),
                                  average = FALSE, by = NULL,
                                  draws = FALSE, n_draws = NULL, ...) {

  dm <- object$density_model
  if (is.null(dm)) {
    stop("`object` does not contain the information required for prediction ",
         "(`density_model`). Please re-estimate the model with the current version of ctrest.",
         call. = FALSE)
  }
  if (!is.numeric(prob) || length(prob) != 1 || is.na(prob) || prob <= 0 || prob >= 1) {
    stop("`prob` must be a single number between 0 and 1.", call. = FALSE)
  }
  type  <- match.arg(type)
  types <- if (type == "both") c("typical", "average") else type
  if (!is.logical(average) || length(average) != 1 || is.na(average)) {
    stop("`average` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(by)) {
    if (!average) stop("`by` can only be used with `average = TRUE`.", call. = FALSE)
    if (!is.character(by) || length(by) == 0) {
      stop("`by` must be a character vector of column names of `newdata`.", call. = FALSE)
    }
  }
  if (!is.logical(draws) || length(draws) != 1 || is.na(draws)) {
    stop("`draws` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(n_draws)) {
    if (!draws) stop("`n_draws` can only be used with `draws = TRUE`.", call. = FALSE)
    if (!is.numeric(n_draws) || length(n_draws) != 1 || is.na(n_draws) || n_draws < 1) {
      stop("`n_draws` must be a positive integer.", call. = FALSE)
    }
  }

  # Species -------------------------------------------------------------------
  sp_all <- dm$species
  if (is.null(species)) {
    species <- sp_all
  } else {
    if (!dm$multi && !identical(as.character(species), sp_all)) {
      stop("`species` must be NULL or '", sp_all, "' for this single-species model.", call. = FALSE)
    }
    unknown <- setdiff(species, sp_all)
    if (length(unknown) > 0) {
      stop("Species not in the model: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
  }

  # Covariate values ----------------------------------------------------------
  if (is.null(newdata)) {
    if (average) {
      newdata <- dm$camera_data
      if (is.null(newdata)) {
        stop("`newdata` is required: this result does not contain the covariates of the cameras. ",
             "Please re-estimate the model with the current version of ctrest.", call. = FALSE)
      }
    } else {
      newdata <- .default_density_newdata(dm)
    }
  } else if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data frame or NULL.", call. = FALSE)
  }
  reserved <- intersect(c("Species", "Type"), names(newdata))
  if (length(reserved) > 0) {
    stop("`newdata` must not contain the column(s): ", paste(reserved, collapse = ", "),
         " (use the `species` and `type` arguments).", call. = FALSE)
  }
  if (!is.null(by) && !all(by %in% names(newdata))) {
    stop("`by` must be column names of `newdata`: ",
         paste(setdiff(by, names(newdata)), collapse = ", "), " not found.", call. = FALSE)
  }

  X      <- .density_design_new(dm, newdata)   # standardized design matrix (n x p)
  tr_idx <- .transect_index(dm, newdata)       # NULL or integer index (length n)
  n_new  <- nrow(X)
  p      <- ncol(X)

  # Groups of rows to average -------------------------------------------------
  if (average) {
    nd_all <- as.data.frame(newdata)
    key <- if (is.null(by)) rep("all", n_new) else do.call(paste, c(unname(as.list(nd_all[by])), sep = "\r"))
    grp_key  <- unique(key)
    grp_rows <- lapply(grp_key, function(k) which(key == k))
    grp_df   <- if (is.null(by)) NULL else nd_all[match(grp_key, key), by, drop = FALSE]
  }

  # Posterior draws -----------------------------------------------------------
  S <- .density_draws(object$samples)
  cn <- colnames(S)

  beta_cols <- .find_columns(cn, if (!dm$multi && p == 1L) c("beta_density", "beta_density[1]")
                                 else paste0("beta_density[", seq_len(p), "]"), p)
  if (is.null(beta_cols)) {
    stop("The posterior draws of `beta_density` are not stored in `object$samples`. ",
         "Please re-estimate the model with the current version of ctrest.", call. = FALSE)
  }
  beta <- S[, beta_cols, drop = FALSE]

  # Random-effect variances: the "average" density is typical * exp(sigma^2 / 2)
  log_adj <- 0
  if ("average" %in% types) {
    if (!"sigma_camera_density" %in% cn) {
      stop("The posterior draws of `sigma_camera_density` are not stored in `object$samples`. ",
           "Please re-estimate the model with the current version of ctrest.", call. = FALSE)
    }
    log_adj <- S[, "sigma_camera_density"]^2 / 2
    if (!is.null(dm$random_effect_density) && is.null(tr_idx)) {
      # the transect is not given: average over transects as well
      if (!"sigma_transect_density" %in% cn) {
        stop("The posterior draws of `sigma_transect_density` are not stored in `object$samples`. ",
             "Please re-estimate the model with the current version of ctrest.", call. = FALSE)
      }
      log_adj <- log_adj + S[, "sigma_transect_density"]^2 / 2
    }
  }

  # Posterior draws to return (the same subset for every row, species and type)
  n_total  <- nrow(S)
  draw_idx <- if (draws && !is.null(n_draws) && n_draws < n_total) {
    sort(sample.int(n_total, as.integer(n_draws)))
  } else {
    seq_len(n_total)
  }

  a  <- (1 - prob) / 2
  qs <- c(a, 0.5, 1 - a)

  summarize <- function(v) {                                  # v: draws x k
    q <- apply(v, 2, stats::quantile, probs = qs, names = FALSE)
    if (is.null(dim(q))) q <- matrix(q, nrow = length(qs))
    data.frame(mean = colMeans(v), sd = apply(v, 2, stats::sd),
               lower = q[1, ], median = q[2, ], upper = q[3, ], row.names = NULL)
  }

  out <- lapply(species, function(sp) {
    m <- match(sp, sp_all)

    if (dm$multi) {
      eff_cols <- match(paste0("species_effect_density[", m, ",", seq_len(p), "]"), cn)
      if (anyNA(eff_cols)) {
        stop("The posterior draws of `species_effect_density` are not stored in `object$samples`.",
             call. = FALSE)
      }
      coef_m <- beta + S[, eff_cols, drop = FALSE]
    } else {
      coef_m <- beta
    }

    eta <- coef_m %*% t(X)                                   # draws x n_new

    if (!is.null(tr_idx)) {
      nT <- length(dm$transect_levels)
      tr_names <- if (dm$multi) paste0("transect_effect_density[", seq_len(nT), ",", m, "]")
                  else paste0("transect_effect_density[", seq_len(nT), "]")
      tr_cols <- match(tr_names, cn)
      if (anyNA(tr_cols)) {
        stop("The posterior draws of `transect_effect_density` are not stored in `object$samples`.",
             call. = FALSE)
      }
      eta <- eta + S[, tr_cols, drop = FALSE][, tr_idx, drop = FALSE]
    }

    pieces <- lapply(types, function(ty) {
      dens <- if (ty == "typical") exp(eta) else exp(eta + log_adj)   # log_adj: one value per draw
      if (average) {
        v <- vapply(grp_rows, function(i) rowMeans(dens[, i, drop = FALSE]), numeric(nrow(dens)))
        if (is.null(dim(v))) v <- matrix(v, ncol = length(grp_rows))
        left <- data.frame(Species = rep(sp, length(grp_rows)), stringsAsFactors = FALSE)
        if (!is.null(grp_df)) left <- cbind(left, grp_df)
        left$n <- lengths(grp_rows)
      } else {
        v <- dens
        left <- data.frame(Species = rep(sp, n_new), stringsAsFactors = FALSE)
        nd <- as.data.frame(newdata)
        if (ncol(nd) > 0L) left <- cbind(left, nd)
      }
      left$Type <- ty
      if (draws) {
        v <- v[draw_idx, , drop = FALSE]                       # the same draws for every row, species and type
        long <- left[rep(seq_len(ncol(v)), each = nrow(v)), , drop = FALSE]
        long$draw    <- rep(draw_idx, times = ncol(v))
        long$density <- as.vector(v)
        rownames(long) <- NULL
        long
      } else {
        cbind(left, summarize(v))
      }
    })
    res_sp <- do.call(rbind, pieces)
    if (draws) return(res_sp)
    # keep the two types of the same row (or group) next to each other
    n_unit <- nrow(pieces[[1]])
    res_sp[order(rep(seq_len(n_unit), times = length(types)), rep(seq_along(types), each = n_unit)), ]
  })

  res <- tibble::as_tibble(do.call(rbind, out))
  if (draws) {
    attr(res, "n_draws") <- length(draw_idx)
    return(res)
  }
  class(res) <- c("ctrest_prediction", class(res))
  attr(res, "prob") <- prob
  res
}


#' Print method for predicted density (explains the column `Type`)
#' @param x A prediction returned by \code{predict.ResultDensity}.
#' @param ... Passed to the tibble print method.
#' @noRd
#' @export
print.ctrest_prediction <- function(x, ...) {
  NextMethod()
  if ("Type" %in% names(x)) {
    cat("# Type: typical = density of a typical camera (camera/transect random effects = 0);\n")
    cat("#       average = mean density over cameras (typical * exp(sigma^2 / 2)). See ?predict.ResultDensity\n")
  }
  invisible(x)
}


# -------------------------------------------------------------------------
# Internal helpers (not exported)
# -------------------------------------------------------------------------

# Maximum number of distinct values for a numeric covariate to be treated as discrete
.density_max_discrete <- 5L
# Number of grid points for continuous covariates in the default profiles
.density_n_grid <- 25L

# Summarize the observed values of each covariate (used for default profiles)
#' @noRd
.predictor_info <- function(data, vars) {
  stats::setNames(lapply(vars, function(v) {
    x <- data[[v]]
    if (is.null(x)) {
      stop("Covariate '", v, "' used in `formula_density` is not found in the data.", call. = FALSE)
    }
    x <- x[!is.na(x)]
    mode_of <- function(z) {
      u <- sort(unique(z))
      u[which.max(tabulate(match(z, u)))]
    }
    if (is.factor(x)) {
      x <- droplevels(x)
      list(type = "factor", levels = levels(x), values = levels(x), ref = mode_of(as.character(x)))
    } else if (is.character(x)) {
      list(type = "character", values = sort(unique(x)), ref = mode_of(x))
    } else if (is.logical(x)) {
      list(type = "logical", values = sort(unique(x)), ref = mode_of(x))
    } else if (is.numeric(x)) {
      u <- sort(unique(x))
      if (length(u) <= .density_max_discrete) {
        list(type = "numeric", values = u, ref = mode_of(x))
      } else {
        list(type = "numeric", values = seq(min(x), max(x), length.out = .density_n_grid), ref = mean(x))
      }
    } else {
      stop("Covariate '", v, "' has an unsupported class (", class(x)[1], ").", call. = FALSE)
    }
  }), vars)
}

# Collect everything required to predict density later (stored in the result)
# formula is kept as a character string: a formula created inside a function
# carries its environment, which would be saved together with the result.
# camera_data: the covariates (and transect) of the cameras used in the estimation,
# used by predict(average = TRUE) without newdata.
#' @noRd
.make_density_model <- function(formula, data, center = NULL, scale = NULL,
                                species, multi,
                                random_effect_density = NULL,
                                transect_levels = NULL,
                                camera_data = NULL) {
  f  <- stats::as.formula(formula)
  tr <- stats::delete.response(stats::terms(f))
  f_chr <- paste(deparse(f, width.cutoff = 500L), collapse = " ")
  if (grepl("(^|[^A-Za-z0-9_.])(poly|scale|ns|bs)\\(", f_chr)) {
    warning("`formula_density` contains a data-dependent transformation (poly, scale, ns or bs); ",
            "predicted density for new covariate values is not supported. ",
            "Use I() or log() instead.", call. = FALSE)
  }

  mf  <- stats::model.frame(tr, data = data)
  X   <- stats::model.matrix(tr, mf)
  isF <- vapply(mf, function(z) is.factor(z) || is.character(z), logical(1))
  xlevels <- lapply(mf[isF], function(z) levels(factor(z)))

  if (is.null(center)) center <- rep(0, ncol(X))
  if (is.null(scale))  scale  <- rep(1, ncol(X))
  if (length(center) != ncol(X) || length(scale) != ncol(X)) {
    stop("Scaling parameters do not match the design matrix of `formula_density`.", call. = FALSE)
  }
  names(center) <- names(scale) <- colnames(X)

  vars <- all.vars(tr)

  cam <- NULL
  if (!is.null(camera_data)) {
    cols <- intersect(c("Station", vars, random_effect_density), names(camera_data))
    if (all(c(vars, random_effect_density) %in% cols)) {
      cam <- as.data.frame(camera_data)[, cols, drop = FALSE]
      check_cols <- setdiff(cols, "Station")
      if (length(check_cols) > 0) {
        cam <- cam[stats::complete.cases(cam[, check_cols, drop = FALSE]), , drop = FALSE]
      }
      rownames(cam) <- NULL
    }
  }

  list(
    formula               = f_chr,
    predictors            = vars,
    predictor_info        = .predictor_info(data, vars),
    col_names             = colnames(X),
    xlevels               = xlevels,
    contrasts             = attr(X, "contrasts"),
    center                = center,
    scale                 = scale,
    species               = as.character(species),
    multi                 = multi,
    random_effect_density = random_effect_density,
    transect_levels       = transect_levels,
    camera_data           = cam
  )
}

# Used after estimation: a failure here must not discard the fitted model
#' @noRd
.make_density_model_safe <- function(...) {
  tryCatch(
    .make_density_model(...),
    error = function(e) {
      warning("The information required to predict density from covariates could not be stored: ",
              conditionMessage(e), call. = FALSE)
      NULL
    }
  )
}

# Default covariate values: one profile per covariate (others held at reference)
#' @noRd
.default_density_newdata <- function(dm) {
  info <- dm$predictor_info
  if (length(info) == 0L) return(NULL)
  refs <- lapply(info, function(z) {
    if (z$type == "factor") factor(z$ref, levels = z$levels) else z$ref
  })
  out <- lapply(names(info), function(v) {
    z <- info[[v]]
    vals <- if (z$type == "factor") factor(z$values, levels = z$levels) else z$values
    nd <- as.data.frame(lapply(refs, function(r) rep(r, length(vals))), stringsAsFactors = FALSE)
    names(nd) <- names(refs)
    nd[[v]] <- vals
    if (length(info) > 1L) nd$Profile <- v
    nd
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

# Standardized design matrix for new covariate values
#' @noRd
.density_design_new <- function(dm, newdata) {
  if (length(dm$predictors) == 0L) {
    n <- if (is.null(newdata)) 1L else nrow(newdata)
    X <- matrix(1, nrow = n, ncol = 1L, dimnames = list(NULL, dm$col_names))
  } else {
    if (is.null(newdata)) stop("`newdata` is required.", call. = FALSE)
    miss <- setdiff(dm$predictors, names(newdata))
    if (length(miss) > 0) {
      stop("`newdata` must contain the covariate(s): ", paste(miss, collapse = ", "), call. = FALSE)
    }
    nd <- as.data.frame(newdata)[, dm$predictors, drop = FALSE]
    bad <- vapply(nd, anyNA, logical(1))
    if (any(bad)) {
      stop("`newdata` contains missing values in: ", paste(names(nd)[bad], collapse = ", "), call. = FALSE)
    }
    if (nrow(nd) == 0L) stop("`newdata` has no rows.", call. = FALSE)
    tr <- stats::delete.response(stats::terms(stats::as.formula(dm$formula, env = globalenv())))
    mf <- stats::model.frame(tr, data = nd,
                             xlev = if (length(dm$xlevels) > 0) dm$xlevels else NULL,
                             na.action = stats::na.fail)
    X <- stats::model.matrix(tr, mf, contrasts.arg = dm$contrasts)
  }
  if (!identical(colnames(X), dm$col_names)) {
    stop("The design matrix of `newdata` does not match the fitted model (columns: ",
         paste(colnames(X), collapse = ", "), "; expected: ",
         paste(dm$col_names, collapse = ", "), ").", call. = FALSE)
  }
  X <- sweep(X, 2, dm$center, "-")
  X <- sweep(X, 2, dm$scale, "/")
  X
}

# Index of the transect of each row of newdata (NULL if the column is absent)
#' @noRd
.transect_index <- function(dm, newdata) {
  col <- dm$random_effect_density
  if (is.null(col) || is.null(newdata) || !col %in% names(newdata)) return(NULL)
  idx <- match(as.character(newdata[[col]]), dm$transect_levels)
  if (anyNA(idx)) {
    stop("Level(s) of '", col, "' not found in the fitted model: ",
         paste(unique(as.character(newdata[[col]])[is.na(idx)]), collapse = ", "), call. = FALSE)
  }
  idx
}

# Posterior draws (rows) of the parameters used for prediction
#' @noRd
.density_draws <- function(samples) {
  keep_pattern <- paste0("^(beta_density|species_effect_density|transect_effect_density|",
                         "sigma_camera_density|sigma_transect_density)(\\[|$)")
  chains <- if (inherits(samples, "mcmc.list")) unclass(samples) else list(samples)
  mats <- lapply(chains, function(ch) {
    cn   <- gsub("\\s+", "", colnames(ch))
    keep <- grepl(keep_pattern, cn)
    m <- as.matrix(ch[, keep, drop = FALSE])
    colnames(m) <- cn[keep]
    m
  })
  do.call(rbind, mats)
}

# Positions of the columns named in `candidates` (all of them), or NULL
#' @noRd
.find_columns <- function(cn, candidates, p) {
  if (p == 1L) {
    hit <- match(candidates, cn)
    hit <- hit[!is.na(hit)]
    return(if (length(hit) >= 1L) hit[1L] else NULL)
  }
  hit <- match(candidates, cn)
  if (anyNA(hit)) NULL else hit
}

# Early check before a long MCMC run: the information needed for prediction can be
# built from the data, and `newdata_density` (if given) is usable
#' @noRd
.check_newdata_density <- function(newdata, formula, data, random_effect_density = NULL) {
  if (!is.null(newdata)) {
    if (!is.data.frame(newdata)) {
      stop("`newdata_density` must be a data frame or NULL.", call. = FALSE)
    }
    reserved <- intersect(c("Species", "Type"), names(newdata))
    if (length(reserved) > 0) {
      stop("`newdata_density` must not contain the column(s): ", paste(reserved, collapse = ", "),
           call. = FALSE)
    }
  }
  if (!is.null(random_effect_density) && !random_effect_density %in% names(data)) {
    stop("`random_effect_density` ('", random_effect_density,
         "') is not a column of `station_effort_data`.", call. = FALSE)
  }
  tr_levels <- if (!is.null(random_effect_density)) levels(factor(data[[random_effect_density]])) else NULL
  dm <- .make_density_model(formula, data, species = "x", multi = FALSE,
                            random_effect_density = random_effect_density,
                            transect_levels = tr_levels)
  if (!is.null(newdata)) {
    tryCatch({
      .density_design_new(dm, newdata)
      .transect_index(dm, newdata)
    }, error = function(e) {
      stop(gsub("`newdata`", "`newdata_density`", conditionMessage(e), fixed = TRUE), call. = FALSE)
    })
  }
  invisible(NULL)
}

# Prediction table stored in the result (NULL if the model has no covariates and
# no newdata was requested; a failure never discards the fitted model)
#' @noRd
.summary_prediction <- function(object, newdata = NULL) {
  dm <- object$density_model
  if (is.null(newdata) && length(dm$predictors) == 0L) return(NULL)
  tryCatch(
    predict.ResultDensity(object, newdata = newdata),
    error = function(e) {
      warning("Predicted density could not be computed: ", conditionMessage(e),
              " (the fitted model is returned; call predict() on it with a valid `newdata`).",
              call. = FALSE)
      NULL
    }
  )
}
