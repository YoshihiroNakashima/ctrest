#' Randomly select videos with even distribution across camera stations
#'
#' Selects \code{N_sampled} videos from a detection dataset. Detections are
#' first filtered for independence (one per station per species per
#' \code{Indep_criteria}-minute window), and the remaining videos are then
#' sampled as evenly as possible across all stations.
#'
#' @param detection_data A data frame containing video detection records.
#' @param col_name_species A single character string. Column name for species.
#'   Default is \code{"Species"}.
#' @param col_name_station A single character string. Column name for camera
#'   station IDs. Default is \code{"Station"}.
#' @param col_name_datetime A single character string. Column name for
#'   detection datetimes. Accepts character or \code{POSIXct} values; common
#'   date-time string formats are supported automatically. Default is
#'   \code{"DateTime"}.
#' @param N_sampled A positive integer. Total number of videos to select.
#' @param Indep_criteria A non-negative number. Minimum time interval (in
#'   minutes) between two detections of the same species at the same station
#'   for them to be considered independent. Default is \code{30}.
#' @param seed An integer for reproducibility, or \code{NULL} (default).
#' @param target_species A character string specifying a single species to
#'   select. If \code{NULL} (default), all species are processed together.
#'
#' @return A data frame of selected video records. If fewer independent
#'   detections than \code{N_sampled} are available, all independent detections
#'   are returned with a warning. Returns \code{NULL} if no valid data remain
#'   after filtering.
#'
#' @details
#' The sampling strategy prevents bias towards highly active stations by
#' interleaving videos across stations: all "1st-choice" videos from every
#' station are collected first, then "2nd-choice" videos, and so on. The
#' first \code{N_sampled} rows of this interleaved list are returned.
#'
#' @export
#' @import dplyr
#' @importFrom lubridate parse_date_time
#' @importFrom rlang sym
#'
#' @examples
#' \dontrun{
#' selected <- select_videos(
#'   detection_data  = detection_data,
#'   col_name_species  = "Species",
#'   col_name_station  = "Station",
#'   col_name_datetime = "DateTime",
#'   N_sampled       = 100,
#'   Indep_criteria  = 30,
#'   target_species  = "SP01",
#'   seed            = 42
#' )
#' }
select_videos <- function(detection_data,
                          col_name_species  = "Species",
                          col_name_station  = "Station",
                          col_name_datetime = "DateTime",
                          N_sampled,
                          Indep_criteria    = 30,
                          seed              = NULL,
                          target_species    = NULL) {

  # --- Input validation -------------------------------------------------------

  if (!is.data.frame(detection_data))
    stop("'detection_data' must be a data frame.", call. = FALSE)
  if (nrow(detection_data) == 0)
    stop("'detection_data' is empty (0 rows).", call. = FALSE)

  for (nm in c("col_name_species", "col_name_station", "col_name_datetime")) {
    val <- get(nm)
    if (!is.character(val) || length(val) != 1)
      stop(sprintf("'%s' must be a single character string.", nm), call. = FALSE)
  }

  req_cols  <- c(col_name_species, col_name_station, col_name_datetime)
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

  if (missing(N_sampled) || is.null(N_sampled))
    stop("'N_sampled' must be specified.", call. = FALSE)
  if (!is.numeric(N_sampled) || length(N_sampled) != 1 ||
      N_sampled <= 0 || N_sampled != floor(N_sampled))
    stop("'N_sampled' must be a positive integer.", call. = FALSE)

  if (!is.numeric(Indep_criteria) || length(Indep_criteria) != 1 ||
      !is.finite(Indep_criteria)  || Indep_criteria < 0)
    stop(
      "'Indep_criteria' must be a single non-negative number (minutes).",
      call. = FALSE
    )

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || seed != floor(seed)))
    stop("'seed' must be a single integer or NULL.", call. = FALSE)

  if (!is.null(target_species)) {
    if (!is.character(target_species) || length(target_species) != 1)
      stop("'target_species' must be a single character string or NULL.", call. = FALSE)
    avail_sp <- unique(detection_data[[col_name_species]])
    if (!target_species %in% avail_sp)
      stop(
        sprintf(
          "Species '%s' not found in column '%s'.\n  Available species: %s",
          target_species, col_name_species,
          paste(sort(avail_sp), collapse = ", ")
        ),
        call. = FALSE
      )
  }

  # --- Setup ------------------------------------------------------------------

  if (!is.null(seed)) set.seed(seed)

  species_sym  <- rlang::sym(col_name_species)
  station_sym  <- rlang::sym(col_name_station)
  datetime_sym <- rlang::sym(col_name_datetime)

  # --- Parse datetimes --------------------------------------------------------

  df_processed <- detection_data %>%
    dplyr::filter(!is.na(!!datetime_sym)) %>%
    dplyr::mutate(
      .parsed_dt = lubridate::parse_date_time(
        !!datetime_sym,
        orders = c("ymd HMS", "ymd HM", "ymd", "dmy HMS", "dmy HM", "dmy",
                   "ymd IMS p"),
        quiet  = TRUE
      )
    ) %>%
    dplyr::filter(!is.na(.data$.parsed_dt))

  if (nrow(df_processed) == 0)
    stop(
      sprintf(
        "Could not parse any datetime value in column '%s'.\n  Example value seen: \"%s\"\n  Supported formats: YYYY-MM-DD HH:MM:SS, YYYY/MM/DD HH:MM, YYYY-MM-DD, DD-MM-YYYY HH:MM:SS",
        col_name_datetime,
        as.character(detection_data[[col_name_datetime]][1])
      ),
      call. = FALSE
    )

  # --- Filter species ---------------------------------------------------------

  if (!is.null(target_species)) {
    df_processed <- df_processed %>%
      dplyr::filter(!!species_sym == target_species)
  }

  if (nrow(df_processed) == 0) {
    message(sprintf("No records found for species '%s'. Returning NULL.", target_species))
    return(NULL)
  }

  # --- Identify independent detections ----------------------------------------

  independent_data <- df_processed %>%
    dplyr::arrange(!!species_sym, !!station_sym, .data$.parsed_dt) %>%
    dplyr::group_by(!!species_sym, !!station_sym) %>%
    dplyr::mutate(
      .diff_min = as.numeric(
        difftime(.data$.parsed_dt, dplyr::lag(.data$.parsed_dt), units = "mins")
      ),
      .indep = is.na(.data$.diff_min) | .data$.diff_min >= Indep_criteria
    ) %>%
    dplyr::filter(.data$.indep) %>%
    dplyr::ungroup()

  n_avail <- nrow(independent_data)

  if (n_avail == 0) {
    message("No independent detections remain after applying Indep_criteria. Returning NULL.")
    return(NULL)
  }

  # --- Stratified random sample -----------------------------------------------

  n_select <- min(N_sampled, n_avail)
  if (n_avail < N_sampled) {
    warning(
      sprintf(
        "Only %d independent detection(s) available; requested N_sampled = %d. Returning all %d.",
        n_avail, N_sampled, n_avail
      ),
      call. = FALSE
    )
  }

  sampled_result <- independent_data %>%
    dplyr::group_by(!!species_sym, !!station_sym) %>%
    dplyr::mutate(.pick_order = sample(dplyr::row_number())) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$.pick_order, sample(dplyr::row_number())) %>%
    dplyr::slice_head(n = n_select) %>%
    dplyr::select(-dplyr::any_of(c(".parsed_dt", ".diff_min", ".indep", ".pick_order")))

  message(sprintf(
    "Selected %d independent detection(s) from %d available (across %d station(s)).",
    nrow(sampled_result),
    n_avail,
    length(unique(sampled_result[[col_name_station]]))
  ))

  sampled_result
}
