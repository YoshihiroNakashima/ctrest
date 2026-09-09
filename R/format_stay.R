#' Prepare data for analyzing animal staying time within the focal area
#'
#' Joins detection records with station metadata and validates staying time data
#' for use with \code{bayes_rest} and \code{bayes_stay_selection}.
#'
#' @param detection_data A data frame containing individual detection records with
#'   at least a station ID column, a species name column, a staying time column, and
#'   a censoring indicator column.
#' @param station_data A data frame with one row per camera station, containing
#'   station-level covariates to be joined to the detection records.
#' @param col_name_station String. Column name for station ID (must be present in
#'   both \code{detection_data} and \code{station_data}). Default is \code{"Station"}.
#' @param col_name_species String. Column name for species names in
#'   \code{detection_data}. Default is \code{"Species"}.
#' @param col_name_stay String. Column name for staying time (numeric, must be > 0).
#'   Default is \code{"Stay"}.
#' @param col_name_cens String. Column name for the censoring indicator (0 = observed,
#'   1 = right-censored). Default is \code{"Cens"}.
#'
#' @return A joined data frame (tibble) with standardised column names
#'   (\code{Station}, \code{Species}, \code{Stay}, \code{Cens}) plus any additional
#'   columns from \code{station_data}, sorted by \code{Species} and \code{Station}.
#'
#' @import dplyr
#' @importFrom rlang sym !! .data
#' @export
format_stay <- function(detection_data,
                        station_data,
                        col_name_station = "Station",
                        col_name_species = "Species",
                        col_name_stay    = "Stay",
                        col_name_cens    = "Cens") {

  # --- Type checks ------------------------------------------------------------

  if (!is.data.frame(detection_data))
    stop("'detection_data' must be a data frame.", call. = FALSE)
  if (!is.data.frame(station_data))
    stop("'station_data' must be a data frame.", call. = FALSE)

  for (nm in c("col_name_station", "col_name_species", "col_name_stay", "col_name_cens")) {
    if (!is.character(get(nm)) || length(get(nm)) != 1)
      stop(sprintf("'%s' must be a single character string.", nm), call. = FALSE)
  }

  # --- Column presence --------------------------------------------------------

  req_det  <- c(col_name_station, col_name_species, col_name_stay, col_name_cens)
  miss_det <- setdiff(req_det, colnames(detection_data))
  if (length(miss_det) > 0)
    stop(
      sprintf(
        "Column(s) not found in 'detection_data': %s\n  Available columns: %s",
        paste(miss_det, collapse = ", "),
        paste(colnames(detection_data), collapse = ", ")
      ),
      call. = FALSE
    )

  if (!col_name_station %in% colnames(station_data))
    stop(
      sprintf(
        "Column '%s' not found in 'station_data'.\n  Available columns: %s",
        col_name_station,
        paste(colnames(station_data), collapse = ", ")
      ),
      call. = FALSE
    )

  # --- station_data integrity -------------------------------------------------

  if (any(duplicated(station_data[[col_name_station]])))
    stop(
      sprintf(
        "'station_data' must have exactly one row per station. Duplicate IDs found in column '%s'.\n  Did you accidentally pass an already-merged data frame?",
        col_name_station
      ),
      call. = FALSE
    )

  if (nrow(detection_data) == 0)
    stop("'detection_data' is empty (0 rows).", call. = FALSE)

  # --- Select and coerce columns ----------------------------------------------

  res <- detection_data %>%
    dplyr::select(
      Station = !!rlang::sym(col_name_station),
      Species = !!rlang::sym(col_name_species),
      Stay    = !!rlang::sym(col_name_stay),
      Cens    = !!rlang::sym(col_name_cens)
    ) %>%
    dplyr::mutate(
      Station = as.character(.data$Station),
      Stay    = suppressWarnings(as.numeric(.data$Stay)),
      Cens    = suppressWarnings(as.integer(.data$Cens))
    )

  # --- Content validation -----------------------------------------------------

  # Cens must be 0 or 1
  cens_obs <- res$Cens[!is.na(res$Cens)]
  bad_cens  <- sort(unique(cens_obs[!cens_obs %in% c(0L, 1L)]))
  if (length(bad_cens) > 0)
    stop(
      sprintf(
        "Column '%s' must contain only 0 (observed) and 1 (right-censored). Unexpected value(s) found: %s",
        col_name_cens,
        paste(bad_cens, collapse = ", ")
      ),
      call. = FALSE
    )

  # Stay must be positive
  n_nonpos <- sum(res$Stay <= 0, na.rm = TRUE)
  if (n_nonpos > 0) {
    warning(
      sprintf(
        "%d record(s) with Stay <= 0 in column '%s' will be excluded. Staying time must be positive (> 0).",
        n_nonpos, col_name_stay
      ),
      call. = FALSE
    )
    res <- res %>% dplyr::filter(is.na(.data$Stay) | .data$Stay > 0)
  }

  # Stay must be finite
  n_inf <- sum(!is.finite(res$Stay) & !is.na(res$Stay))
  if (n_inf > 0) {
    warning(
      sprintf(
        "%d record(s) with non-finite Stay values (Inf or -Inf) in column '%s' will be excluded.",
        n_inf, col_name_stay
      ),
      call. = FALSE
    )
    res <- res %>% dplyr::filter(is.na(.data$Stay) | is.finite(.data$Stay))
  }

  # Remove rows with NA Stay or NA Cens
  n_before <- nrow(res)
  res <- res %>% dplyr::filter(!is.na(.data$Stay), !is.na(.data$Cens))
  n_removed <- n_before - nrow(res)
  if (n_removed > 0)
    message(
      sprintf(
        "%d record(s) with missing values in '%s' or '%s' were excluded.",
        n_removed, col_name_stay, col_name_cens
      )
    )

  if (nrow(res) == 0)
    stop(
      "No valid records remain in 'detection_data' after input validation. Check Stay values and censoring indicators.",
      call. = FALSE
    )

  # --- Station ID mismatch warnings -------------------------------------------

  det_stations <- unique(res$Station)
  sta_stations <- unique(as.character(station_data[[col_name_station]]))

  extra_in_det <- setdiff(det_stations, sta_stations)
  if (length(extra_in_det) > 0)
    warning(
      sprintf(
        "Station(s) present in 'detection_data' but not in 'station_data' — covariate columns will be NA: %s",
        paste(extra_in_det, collapse = ", ")
      ),
      call. = FALSE
    )

  not_in_det <- setdiff(sta_stations, det_stations)
  if (length(not_in_det) > 0)
    message(
      sprintf(
        "Note: %d station(s) in 'station_data' have no stay records in 'detection_data': %s",
        length(not_in_det),
        paste(not_in_det, collapse = ", ")
      )
    )

  # --- Prepare station_data for join ------------------------------------------

  clean_station <- station_data
  if (col_name_station != "Station" && "Station" %in% colnames(clean_station))
    clean_station <- clean_station %>% dplyr::select(-.data$Station)

  clean_station <- clean_station %>%
    dplyr::rename(Station = !!rlang::sym(col_name_station)) %>%
    dplyr::mutate(Station = as.character(.data$Station))

  reserved_names <- c("Species", "Stay", "Cens")
  dup_cols <- intersect(colnames(clean_station), reserved_names)
  if (length(dup_cols) > 0)
    clean_station <- clean_station %>% dplyr::select(-dplyr::all_of(dup_cols))

  # --- Join and return --------------------------------------------------------

  res %>%
    dplyr::left_join(clean_station, by = "Station") %>%
    dplyr::arrange(.data$Species, .data$Station)
}

utils::globalVariables(c("Station", "Species", "Stay", "Cens"))
