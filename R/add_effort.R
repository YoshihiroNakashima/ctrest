#' Caluculating camera trapping efforts from the detection data
#'
#' @param detection_data A data frame containing information for individual videos, with six columns: camera station name (station), capture datetime (datetime), species name (species), number of passes through the focal area (nfocal), duration of stay within the focal area for each pass in a video (stay), and whether the observation of the pass was censored (cens).
#' @param station_data_formatted ""
#' @param col_name_station Column name containing station id info
#' @param col_name_datetime Column name containing datetime info
#' @param col_name_term Column name containing datetime info
#' @param plot Column name containing datetime info
#' @param font_size Column name containing datetime info
#' @return data frame containing camera trapping effort (days) for each station
#' @details 'Start' and 'End' represent the dates when the camera trap first and last captured an image, respectively.
#' @export
#' @import dplyr ggplot2
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
#'
#' station_data_complete <- add_effort(
#'   detection_data = detection_data,
#'   station_data_formatted = station_data_rest,
#'   col_name_station = "Station",
#'   col_name_term = "Term",
#'   col_name_datetime = "DateTime",
#'   plot = TRUE,
#'   font_size = 5
#' )
add_effort <-
  function(detection_data,
           station_data_formatted,
           col_name_station,
           col_name_datetime,
           col_name_term = NULL,
           plot = FALSE,
           font_size = 8){
    if(is.null(col_name_term)){
      effort_temp <- detection_data %>%
        rename(Station = !!sym(col_name_station), DateTime = !!sym(col_name_datetime), Term = !!sym(col_name_term)) %>%
        group_by(Station) %>%
        summarize(Start = min(DateTime, na.rm = TRUE), End = max(DateTime, na.rm = TRUE)) %>%
        mutate(effort = as.numeric(difftime(End, Start, units = "days")))
    } else {
      effort_temp <- detection_data %>%
        rename(Station = !!sym(col_name_station), DateTime = !!sym(col_name_datetime), Term = !!sym(col_name_term)) %>%
        group_by(Station, Term) %>%
        summarize(Start = min(DateTime, na.rm = TRUE), End = max(DateTime, na.rm = TRUE)) %>%
        mutate(effort = as.numeric(difftime(End, Start, units = "days"))) %>%
        ungroup()
    }

    effort_temp2 <- effort_temp %>%
      group_by(Station) %>%
      summarize(Effort = sum(effort)) %>%
      left_join(station_data_formatted, by = "Station") %>%
      arrange(Species) %>%
      filter(Effort != 0)

    if(sum(effort_temp2$Effort == 0)){
      print("Sampling efforts of some stations were calculated as 0 possibly due to only one detection. The stations' records were removed" , call. = FALSE)
    }

    if(plot == TRUE){
      g <- ggplot(data = effort_temp) +
        geom_segment(aes(x = Start, xend = End, y = Station, yend = Station),
                     alpha = 0.3,
                     linewidth = 2.5) +
        geom_point(aes(x = Start, y = Station), shape = 4, col = "red", alpha = 0.8) +
        theme(axis.text.y = element_text(size = font_size)) +
        xlab("Survey Period")
      print(g)
      return(effort_temp2)
    } else{
      return(effort_temp2)
    }

  }
Start <- End <- effort <- Effort <- Station <- DateTime <- Term <- Species <- NULL
