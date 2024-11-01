#' Prepare data for estimating the proportion of time animal being active
#'
#' @param detection_data A data frame containing information for individual videos, with six columns: camera station name (station), capture datetime (datetime), species name (species), number of passes through the focal area (nfocal), duration of stay within the focal area for each pass in a video (stay), and whether the observation of the pass was censored (cens).
#' @param col_name_station Column name containing station id info
#' @param col_name_species Column name containing species name detected
#' @param col_name_datetime Column name containing datetime info
#' @param target_species The species name for which you want to estimate his density
#' @param indep_time The interval (in seconds) to consider detection as independent
#' @return A vector containing information about the times of an animal detection
#' @export
#' @import dplyr lubridate ggplot2
#' @examples format_activity(
#' detection_data = detection_data,
#' col_name_station = "Station",
#' col_name_species = "Species",
#' col_name_datetime = "DateTime",
#' target_species = "SP01",
#' indep_time = 30)
format_activity <- function(detection_data,
                            col_name_station,
                            col_name_species,
                            col_name_datetime,
                            target_species,
                            indep_time = 30){
  activity_data <- detection_data %>%
    rename(Station = !!sym(col_name_station), Species = !!sym(col_name_species), DateTime = !!sym(col_name_datetime)) %>%
    mutate(Indep = case_when(
      Station != lag(Station) ~ TRUE,
      difftime(DateTime, lag(DateTime), units = "mins") > indep_time ~ TRUE,
      TRUE ~ FALSE
    )) %>%
    filter(Species == target_species) %>%
    mutate(time = 2 * pi * (hour(DateTime) * 60 * 60 + minute(DateTime) * 60 + second(DateTime))/(24 * 60 * 60)) %>%
    filter(!is.na(time)) %>%
    pull(time)
  return(activity_data)
}

time <- Species <- DateTime <- NULL
