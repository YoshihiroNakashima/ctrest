#' Aggregate the number of animal passes within focal areas for each station
#'
#' Summarises per-video pass counts from detection records into station-level
#' totals (REST) or pass-category distributions (RAD-REST) and joins with
#' station metadata.
#'
#' @param detection_data A data frame containing individual detection records.
#' @param station_data A data frame with one row per camera station containing
#'   station-level covariates. Must have a station ID column matching
#'   \code{col_name_station}.
#' @param col_name_station A string specifying the column name for station IDs
#'   (must be present in both data frames).
#' @param col_name_species A string specifying the column name for species names
#'   in \code{detection_data}.
#' @param col_name_y A string specifying the column name for the number of animal
#'   passes per video in \code{detection_data}. Values must be non-negative. For
#'   RAD-REST, values must be non-negative integers.
#' @param model A string specifying the model type: \code{"REST"} or
#'   \code{"RAD-REST"}.
#'
#' @return A data frame with aggregated detection counts joined with station
#'   metadata, sorted by \code{Species} and \code{Station}. For \code{"REST"},
#'   contains a column \code{Y} (total passes per station-species pair). For
#'   \code{"RAD-REST"}, contains \code{N} (total detections) and columns
#'   \code{y_0}, \code{y_1}, ... (counts per pass category).
#'
#' @export
#' @import dplyr tidyr stringr rlang
format_station_data <- function(detection_data,
                                station_data,
                                col_name_station,
                                col_name_species,
                                col_name_y,
                                model) {

  # --- Type checks ------------------------------------------------------------

  if (!is.data.frame(detection_data))
    stop("'detection_data' must be a data frame.", call. = FALSE)
  if (!is.data.frame(station_data))
    stop("'station_data' must be a data frame.", call. = FALSE)

  for (nm in c("col_name_station", "col_name_species", "col_name_y")) {
    if (!is.character(get(nm)) || length(get(nm)) != 1)
      stop(sprintf("'%s' must be a single character string.", nm), call. = FALSE)
  }

  if (!is.character(model) || length(model) != 1 || !model %in% c("REST", "RAD-REST"))
    stop("'model' must be either \"REST\" or \"RAD-REST\".", call. = FALSE)

  # --- Column presence --------------------------------------------------------

  req_cols  <- c(col_name_station, col_name_species, col_name_y)
  miss_cols <- setdiff(req_cols, colnames(detection_data))
  if (length(miss_cols) > 0)
    stop(
      sprintf(
        "Column(s) not found in 'detection_data': %s\n  Available columns: %s",
        paste(miss_cols, collapse = ", "),
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

  # --- Pass count validation --------------------------------------------------

  y_raw <- detection_data[[col_name_y]]
  y_num <- suppressWarnings(as.numeric(y_raw))

  if (any(!is.na(y_raw) & is.na(y_num)))
    stop(
      sprintf(
        "Column '%s' contains values that cannot be coerced to numeric. Pass counts must be numbers.",
        col_name_y
      ),
      call. = FALSE
    )

  n_neg <- sum(y_num < 0, na.rm = TRUE)
  if (n_neg > 0)
    stop(
      sprintf(
        "Column '%s' contains %d negative value(s). Pass counts must be non-negative (>= 0).",
        col_name_y, n_neg
      ),
      call. = FALSE
    )

  if (model == "RAD-REST") {
    non_int <- y_num[!is.na(y_num) & y_num != floor(y_num)]
    if (length(non_int) > 0)
      warning(
        sprintf(
          "Column '%s' contains non-integer value(s) (e.g., %g). For RAD-REST, pass counts should be whole numbers.",
          col_name_y, non_int[1]
        ),
        call. = FALSE
      )
  }

  # --- Station ID mismatch warning --------------------------------------------

  det_stations <- unique(as.character(detection_data[[col_name_station]]))
  sta_stations <- unique(as.character(station_data[[col_name_station]]))

  extra_in_det <- setdiff(det_stations, sta_stations)
  if (length(extra_in_det) > 0)
    warning(
      sprintf(
        "Station(s) in 'detection_data' not found in 'station_data' — these records will be excluded from the output: %s",
        paste(extra_in_det, collapse = ", ")
      ),
      call. = FALSE
    )

  # --- Select and coerce columns ----------------------------------------------

  det_clean <- detection_data %>%
    dplyr::select(
      Station = !!rlang::sym(col_name_station),
      Species = !!rlang::sym(col_name_species),
      y       = !!rlang::sym(col_name_y)
    ) %>%
    dplyr::mutate(
      Station = as.character(.data$Station),
      Species = as.character(.data$Species),
      y       = as.numeric(.data$y)
    )

  # --- Prepare station_data for join ------------------------------------------

  clean_station <- station_data
  if (col_name_station != "Station" && "Station" %in% colnames(clean_station))
    clean_station <- clean_station %>% dplyr::select(-.data$Station)
  clean_station <- clean_station %>%
    dplyr::rename(Station = !!rlang::sym(col_name_station)) %>%
    dplyr::mutate(Station = as.character(.data$Station))

  reserved_regex <- "^(Species|Y|N|y|y_\\d+)$"
  dup_cols <- setdiff(grep(reserved_regex, colnames(clean_station), value = TRUE), "Station")
  if (length(dup_cols) > 0)
    clean_station <- clean_station %>% dplyr::select(-dplyr::all_of(dup_cols))

  # --- Build full station x species grid --------------------------------------

  all_stations <- unique(clean_station$Station)
  all_species  <- unique(stats::na.omit(det_clean$Species))

  detection_grid <- tidyr::crossing(Station = all_stations, Species = all_species)

  # --- Model-specific aggregation ---------------------------------------------

  if (model == "REST") {

    n_na_y <- sum(is.na(det_clean$y) & !is.na(det_clean$Species))
    if (n_na_y > 0)
      warning(
        sprintf(
          "%d record(s) with NA in '%s' were assumed to be 1 pass (minimum detectable). Set explicit counts to suppress this warning.",
          n_na_y, col_name_y
        ),
        call. = FALSE
      )

    det_agg <- det_clean %>%
      dplyr::filter(!is.na(.data$Species)) %>%
      dplyr::mutate(y = ifelse(is.na(.data$y), 1, .data$y)) %>%
      dplyr::group_by(.data$Station, .data$Species) %>%
      dplyr::summarize(Y = sum(.data$y), .groups = "drop")

    detection_final <- detection_grid %>%
      dplyr::left_join(det_agg, by = c("Station", "Species")) %>%
      dplyr::mutate(Y = tidyr::replace_na(.data$Y, 0))

  } else {  # RAD-REST

    det_N <- det_clean %>%
      dplyr::filter(!is.na(.data$Species)) %>%
      dplyr::group_by(.data$Station, .data$Species) %>%
      dplyr::summarise(N = dplyr::n(), .groups = "drop")

    det_y <- det_clean %>%
      dplyr::filter(!is.na(.data$Species), !is.na(.data$y)) %>%
      dplyr::group_by(.data$Station, .data$Species, .data$y) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      tidyr::pivot_wider(
        names_from   = .data$y,
        values_from  = .data$n,
        names_prefix = "y_",
        values_fill  = 0
      )

    detection_final <- detection_grid %>%
      dplyr::left_join(det_N, by = c("Station", "Species")) %>%
      dplyr::left_join(det_y, by = c("Station", "Species"))

    y_cols       <- grep("^y_\\d+$", colnames(detection_final), value = TRUE)
    vars_to_fill <- c("N", y_cols)
    detection_final <- detection_final %>%
      dplyr::mutate(dplyr::across(dplyr::any_of(vars_to_fill), ~ tidyr::replace_na(.x, 0)))

    indices         <- if (length(y_cols) > 0) as.integer(stringr::str_extract(y_cols, "\\d+")) else integer(0)
    max_idx         <- max(c(0L, indices), na.rm = TRUE)
    missing_indices <- setdiff(0:max_idx, indices)
    for (col in paste0("y_", missing_indices)) detection_final[[col]] <- 0L

    sorted_y_cols   <- stringr::str_sort(grep("^y_\\d+$", colnames(detection_final), value = TRUE), numeric = TRUE)
    detection_final <- detection_final %>%
      dplyr::select(.data$Station, .data$Species, .data$N, dplyr::all_of(sorted_y_cols))
  }

  # --- Final join and sort ----------------------------------------------------

  final_joined <- detection_final %>%
    dplyr::left_join(clean_station, by = "Station") %>%
    dplyr::arrange(.data$Species, .data$Station)

  final_joined
}

utils::globalVariables(c("Station", "Species", "y", "Y", "N", "n"))
