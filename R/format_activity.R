#' Prepare data for estimating the proportion of time an animal is active
#'
#' Filters detection data to retain only independent detections and transforms
#' detection time into a circular format (radians) for activity modeling.
#'
#' @param detection_data A data frame containing individual detection records with at
#'   least a station ID column, a species name column, and a datetime column.
#' @param col_name_station A string specifying the column name for camera station IDs.
#' @param col_name_species A string specifying the column name for species names.
#' @param col_name_datetime A string specifying the column name for detection datetimes.
#'   Accepted formats: \code{"YYYY-MM-DD HH:MM:SS"}, \code{"YYYY/MM/DD HH:MM"}, or
#'   \code{"YYYY-MM-DD"} (time assumed 00:00). A \code{POSIXct} column is also accepted.
#' @param indep_time A non-negative number specifying the minimum time interval (in
#'   minutes) required between two detections of the same species at the same station
#'   to consider them independent. Default is 30.
#'
#' @return A data frame with columns \code{Species}, \code{Station}, and \code{time}
#'   (detection time in radians, 0–2π), containing only independent detections.
#'
#' @export
#' @import dplyr lubridate stringr
#' @importFrom rlang sym !!
#' @examples
#' \dontrun{
#' format_activity(
#'   detection_data    = detection_data,
#'   col_name_station  = "Station",
#'   col_name_species  = "Species",
#'   col_name_datetime = "DateTime",
#'   indep_time        = 30
#' )
#' }
format_activity <- function(detection_data,
                            col_name_station,
                            col_name_species,
                            col_name_datetime,
                            indep_time = 30) {

  # --- Input validation -------------------------------------------------------

  if (!is.data.frame(detection_data))
    stop("'detection_data' must be a data frame.", call. = FALSE)
  if (nrow(detection_data) == 0)
    stop("'detection_data' is empty (0 rows).", call. = FALSE)

  for (nm in c("col_name_station", "col_name_species", "col_name_datetime")) {
    if (!is.character(get(nm)) || length(get(nm)) != 1)
      stop(sprintf("'%s' must be a single character string.", nm), call. = FALSE)
  }

  if (!is.numeric(indep_time) || length(indep_time) != 1 ||
      !is.finite(indep_time)  || indep_time < 0)
    stop(
      "'indep_time' must be a single non-negative number (minimum independence interval in minutes).",
      call. = FALSE
    )

  required_cols <- c(col_name_station, col_name_species, col_name_datetime)
  missing_cols  <- setdiff(required_cols, colnames(detection_data))
  if (length(missing_cols) > 0)
    stop(
      sprintf(
        "Column(s) not found in 'detection_data': %s\n  Available columns: %s",
        paste(missing_cols, collapse = ", "),
        paste(colnames(detection_data), collapse = ", ")
      ),
      call. = FALSE
    )

  # --- Rename for internal use ------------------------------------------------

  detection_data <- detection_data %>%
    rename(
      Station  = !!sym(col_name_station),
      Species  = !!sym(col_name_species),
      DateTime = !!sym(col_name_datetime)
    )

  # --- Datetime parsing -------------------------------------------------------

  if (!inherits(detection_data$DateTime, "POSIXt")) {
    dt_raw    <- as.character(detection_data$DateTime)
    dt_filled <- ifelse(str_detect(dt_raw, ":"), dt_raw, paste0(dt_raw, " 00:00"))
    parsed    <- suppressWarnings(
      lubridate::parse_date_time(
        dt_filled,
        orders = c("Ymd HMS", "Ymd HM", "Ymd H", "Ymd"),
        quiet  = TRUE
      )
    )
    n_fail <- sum(is.na(parsed))
    if (n_fail == length(parsed)) {
      ex <- dt_raw[!is.na(dt_raw)][1]
      stop(
        sprintf(
          paste0(
            "Could not parse any datetime value in column '%s'.\n",
            "  Example value seen: \"%s\"\n",
            "  Supported formats: 'YYYY-MM-DD HH:MM:SS', 'YYYY/MM/DD HH:MM', 'YYYY-MM-DD'.\n",
            "  Alternatively, convert the column to POSIXct before calling this function."
          ),
          col_name_datetime, ex
        ),
        call. = FALSE
      )
    }
    if (n_fail > 0)
      warning(
        sprintf(
          "%d of %d datetime value(s) in column '%s' could not be parsed and will be excluded.",
          n_fail, length(parsed), col_name_datetime
        ),
        call. = FALSE
      )
    detection_data$DateTime <- parsed
    detection_data <- detection_data %>% filter(!is.na(DateTime))
  }

  # --- Independence filter ----------------------------------------------------

  activity_data <- detection_data %>%
    arrange(Station, DateTime) %>%
    group_by(Species, Station) %>%
    mutate(Indep = case_when(
      is.na(lag(DateTime))                                            ~ TRUE,
      difftime(DateTime, lag(DateTime), units = "mins") > indep_time ~ TRUE,
      TRUE                                                            ~ FALSE
    )) %>%
    filter(Indep) %>%
    mutate(
      time = 2 * pi *
        (hour(DateTime) * 3600 + minute(DateTime) * 60 + second(DateTime)) /
        (24 * 3600)
    ) %>%
    select(Species, Station, time) %>%
    ungroup()

  if (nrow(activity_data) == 0)
    stop(
      sprintf(
        "No independent detections remain after applying indep_time = %g min. Consider reducing 'indep_time'.",
        indep_time
      ),
      call. = FALSE
    )

  activity_data
}

time <- NULL
