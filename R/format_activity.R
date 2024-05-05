#' Prepare data for estimating the proportion of time animal being active
#'
#' @param detection_data A data frame containing information for individual videos, with six columns: camera station name (station), capture datetime (datetime), species name (species), number of passes through the focal area (nfocal), duration of stay within the focal area for each pass in a video (stay), and whether the observation of the pass was censored (cens).
#' @param col_name_species XXX
#' @param col_name_datetime Column name containing datetime info
#' @param target_species The species name for which you want to estimate his density
#' @return A vector containing information about the times of an animal detection
#' @export
#' @import dplyr lubridate ggplot2
#' @examples format_activity(detection_data = detection_data,
#' col_name_species = "Species",
#' col_name_datetime = "DateTime",
#' target_species = "A")
format_activity <- function(detection_data,
                            col_name_species,
                            col_name_datetime,
                            target_species){
  activity_data <- detection_data %>%
    filter(!!sym(col_name_species) == target_species) %>%
    mutate(time = 2 * pi * (hour(!!sym(col_name_datetime)) * 60 * 60 + minute(!!sym(col_name_datetime)) * 60 + second(!!sym(col_name_datetime)))/(24 * 60 * 60)) %>%
    filter(!is.na(time)) %>%
    pull(time)
}

time <- NULL
