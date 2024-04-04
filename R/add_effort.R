#' Caluculating camera trapping efforts from the detection data
#'
#' @param detection_data A data frame containing information for individual videos, with six columns: camera station name (station), capture datetime (datetime), species name (species), number of passes through the focal area (nfocal), duration of stay within the focal area for each pass in a video (stay), and whether the observation of the pass was censored (cens).
#' @param npass_data ""
#' @param col_name_station Column name containing station id info
#' @param col_name_datetime Column name containing datetime info
#'
#' @return data frame containing camera trapping effort (days) for each station
#' @details 'start' and 'end' represent the dates when the camera trap first and last captured an image, respectively.
#' @export
#' @import dplyr
#' @examples
#' add_effort(detection_data = detection_data,
#' npass_data = station_data,
#' col_name_station = "station",
#' col_name_datetime = "datetime")
add_effort <-
  function(detection_data,
           npass_data,
           col_name_station = "station",
           col_name_datetime = "datetime"){
    effort_temp <- detection_data %>%
      group_by(!!sym(col_name_station)) %>%
      summarize(start = min(!!sym(col_name_datetime), na.rm = TRUE), end = max(!!sym(col_name_datetime), na.rm = TRUE)) %>%
      mutate(effort = as.numeric(difftime(end, start))) %>%
      right_join(npass_data, by = col_name_station) %>%
      select(- start, - end)
    print(effort_temp)
  }

start <- end <- NULL
