#' Prepare data for estimating the proportion of time an animal is active
#'
#' This function filters detection data to retain only independent detections
#' of the specified species and transforms detection time into a circular format
#' for use in activity modeling.
#'
#' @param detection_data A data frame containing information for individual videos, including:
#'   - Camera station ID (character)
#'   - Capture datetime (character)
#'   - Species name (character)
#'   - Other columns (e.g., number of passes, staying time, censoring flag) are ignored by this function.
#' @param col_name_station A string specifying the column name that contains camera station IDs.
#' @param col_name_species A string specifying the column name that contains species names.
#' @param col_name_datetime A string specifying the column name that contains datetime information (character).
#'   If time is missing (e.g., only "2022/03/01"), "00:00" is assumed.
#' @param target_species A character vector specifying the species name(s) for which the activity proportion is to be estimated.
#' @param indep_time A numeric value specifying the minimum time interval (in minutes) to consider two detections as independent. Default is 30.
#'
#' @return A data frame with three columns:
#'   - `Species` (character): Species name.
#'   - `Station` (character): Camera station ID.
#'   - `time` (numeric): Time of detection, scaled from 0 to 2π (radians) for circular activity analysis.
#'   Only independent detections are included.
#'
#' @export
#' @import dplyr lubridate stringr
#' @examples
#' format_activity(
#'   detection_data = detection_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_datetime = "DateTime",
#'   target_species = "SP01",
#'   indep_time = 30
#' )
format_activity <- function(detection_data,
                            col_name_station,
                            col_name_species,
                            col_name_datetime,
                            target_species,
                            indep_time = 30) {

  # Check if required columns exist
  required_cols <- c(col_name_station, col_name_species, col_name_datetime)
  missing_cols <- required_cols[!required_cols %in% colnames(detection_data)]

  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing from 'detection_data':", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Check if target species exist in the dataset
  if (!any(detection_data[[col_name_species]] %in% target_species)) {
    warning("None of the target species are present in 'detection_data'. Returning an empty data frame.", call. = FALSE)
    return(data.frame())
  }

  # Rename columns for internal consistency
  detection_data <- detection_data %>%
    rename(
      Station = !!sym(col_name_station),
      Species = !!sym(col_name_species),
      DateTime = !!sym(col_name_datetime)
    )

  # Convert datetime column to POSIXct format
  if(is.character(detection_data$DateTime)){
    detection_data <- detection_data %>%
      mutate(
        DateTime_clean = if_else(
          str_detect(DateTime, ":"),
          DateTime,
          paste0(DateTime, " 00:00")
        ),
        DateTime = ymd_hm(DateTime_clean)
      )
  }

  # Filter and process independent detections
  activity_data <- detection_data %>%
    filter(Species %in% target_species) %>%
    arrange(Station, DateTime) %>%
    group_by(Species, Station) %>%
    mutate(Indep = case_when(
      is.na(lag(DateTime)) ~ TRUE,
      difftime(DateTime, lag(DateTime), units = "mins") > indep_time ~ TRUE,
      TRUE ~ FALSE
    )) %>%
    filter(Indep) %>%
    mutate(
      time = 2 * pi * (hour(DateTime) * 3600 + minute(DateTime) * 60 + second(DateTime)) / (24 * 3600)
    ) %>%
    select(Species, Station, time) %>%
    ungroup()

  return(activity_data)
}
