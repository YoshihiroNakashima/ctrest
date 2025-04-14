#' Calculate Camera Trapping Effort (Days) from Detection Data
#'
#' @description Computes the total camera trapping effort (in days) for each station based on detection data.
#' If monitoring periods (terms) are specified, effort is calculated separately for each term and then aggregated.
#' @param detection_data A data frame containing information on individual video records, with six columns:
#'   - Camera station ID (character)
#'   - Capture datetime (character)
#'   - Species name (character; not used in this function)
#'   - Number of passes through the focal area (numeric; not used in this function)
#'   - Staying time (in seconds) within the focal area per pass in each video (numeric; not used in this function)
#'   - Whether the observation of the pass was censored (1 = censored, 0 = observed; not used in this function)
#' @param station_data_formatted A data frame returned by the `format_station_data()` function. This will be updated with camera effort information.
#' @param col_name_station A string specifying the column name containing camera station IDs.
#' @param col_name_datetime A string specifying the column name containing datetime information. Ensure all timestamps are correct.
#' @param col_name_term (Optional) A string specifying the column name indicating different survey periods if cameras were checked during the study period. Set to NULL if only deployment and retrieval dates are available.
#' @param plot Logical. If TRUE, a plot showing the camera operation periods will be displayed. This is useful for identifying timestamp inconsistencies in `detection_data`.
#' @param font_size Font size for the station labels in the plot when `plot = TRUE`. Default is 8.
#' @return A data frame with the same structure as `station_data_formatted`, with an additional column `Effort` representing the camera trapping effort (in days) for each station. Stations with `Effort == 0` are removed from the final output.
#' @export
#' @import dplyr ggplot2 lubridate
#' @examples
#' station_data_rest <- format_station_data(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_y = "y",
#'   model = "REST",
#'   target_species = "SP01"
#' )
#' station_effort_rest <- add_effort(
#'   detection_data = detection_data,
#'   station_data_formatted = station_data_rest,
#'   col_name_station = "Station",
#'   col_name_term = "Term",
#'   col_name_datetime = "DateTime",
#'   plot = TRUE
#' )
add_effort <- function(detection_data,
                       station_data_formatted,
                       col_name_station,
                       col_name_datetime,
                       col_name_term = NULL,
                       plot = FALSE,
                       font_size = 8) {
  # Check if required columns exist
  required_cols <- c(col_name_station, col_name_datetime)
  if (!all(required_cols %in% colnames(detection_data))) {
    stop("Error: One or more specified column names do not exist in detection_data.")
  }
  detection_data <- detection_data %>%
    rename(Station = !!sym(col_name_station), DateTime = !!sym(col_name_datetime))
  if(!is.null(col_name_term)) {
    detection_data <- detection_data %>%
      rename(Term = !!sym(col_name_term))
  }

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

  if (any(is.na(detection_data[[col_name_datetime]]))) {
    stop("Error: Failed to convert datetime column. Ensure the format is YYYY-MM-DD HH:MM:SS. NA is not allowed.")
  }

  # Calculate effort per station (or per term if applicable)
  if (is.null(col_name_term)) {
    effort_temp <- detection_data %>%
      group_by(Station) %>%
      summarize(Start = min(DateTime, na.rm = TRUE), End = max(DateTime, na.rm = TRUE), .groups = 'drop') %>%
      mutate(effort = as.numeric(difftime(End, Start, units = "days")))
  } else {
    effort_temp <- detection_data %>%
      group_by(Station, Term) %>%
      summarize(Start = min(DateTime, na.rm = TRUE), End = max(DateTime, na.rm = TRUE), .groups = 'drop') %>%
      mutate(effort = as.numeric(difftime(End, Start, units = "days"))) %>%
      ungroup()
  }

  # Aggregate efforts per station and merge with station data
  effort_temp2 <- effort_temp %>%
    group_by(Station) %>%
    summarize(Effort = sum(effort), .groups = 'drop') %>%
    left_join(station_data_formatted, by = "Station") %>%
    filter(Effort != 0)

  # Warning if any station had zero effort
  if (any(effort_temp2$Effort == 0)) {
    warning("Warning: Sampling effort for some stations was calculated as 0, possibly due to a single detection. These stations were removed from the output.")
  }

  # Generate plot if requested
  if (plot) {
    g <- ggplot(data = effort_temp) +
      geom_segment(aes(x = Start, xend = End, y = Station, yend = Station),
                   alpha = 0.3,
                   linewidth = 2.5) +
      geom_point(aes(x = Start, y = Station), shape = 4, col = "red", alpha = 0.8) +
      theme(axis.text.y = element_text(size = nrow(effort_temp2) /25)) +
      xlab("Survey Period")
    print(g)
  }

  return(effort_temp2)
}

Start <- End <- effort <- Effort <- Station <- DateTime <- Term <- Species <- NULL
