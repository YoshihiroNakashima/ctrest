#' Title
#'
#' @param detection_data A data frame containing information for individual videos, with six columns: camera station name (station), capture datetime (datetime), species name (species), number of passes through the focal area (nfocal), duration of stay within the focal area for each pass in a video (stay), and whether the observation of the pass was censored (cens).
#' @param target_species The species name for which you want to estimate his density
#' @param col_name_datetime Column name containing datetime info
#'
#' @return A vector containing information about the times of an animal detection.
#' @export
#' @import lubridate
#' @examples format_activity(detection_data = detection_data,
#' col_name_datetime = "datetime")
format_activity <-
  function(detection_data = detection_data,
           col_name_datetime = "datetime",
           target_species = "A"){
    act_temp <- detection_data %>%
      filter(species == target_species) %>%
      mutate(time = 2 * pi * (hour(!!sym(col_name_datetime)) * 60 * 60 + minute(!!sym(col_name_datetime)) * 60 + second(!!sym(col_name_datetime)))/(24 * 60 * 60)) %>%
      filter(!is.na(time)) %>%
      pull(time)
    print(act_temp)
  }

time <- species <- NULL
