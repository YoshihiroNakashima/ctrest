#' Aggregate the number of passes within focal areas for each station from the detection data
#'
#' @param detection_data A data frame containing information for individual videos, with six columns: camera station name (station), capture datetime (datetime), species name (species), number of passes through the focal area (nfocal), duration of stay within the focal area for each pass in a video (stay), and whether the observation of the pass was censored (cens).
#' @param station_data A data frame containing information for each camera station, with one row per station. It includes the station name (station), location coordinates, covariates, and so forth. Unlike the detection_data, it lists all the camera stations that have been set up.
#' @param col_name_station Column name containing station id info
#' @param col_name_species Column name containing species name detected
#' @param col_name_y Column name containing number of animal passage within a focal area
#' @param target_species The species name for which you want to estimate his density
#' @param model Model name used ("REST" or "RAD-REST")
#' @return data frame containing at least 3 columns (station, species and npass)
#' @export
#' @import dplyr
#' @examples
#' format_station_data(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_y = "y",
#'   model = "REST",
#'   target_species = "A"
#' )

format_station_data <- function(detection_data,
                                station_data,
                                col_name_station,
                                col_name_species,
                                col_name_y,
                                target_species,
                                model) {
  if(!(model == "REST" | model == "RAD-REST")){
    stop("check model name! 'model' must be 'REST' or 'RAD-REST'" , call. = FALSE)
  }
  if(model == "REST"){
    detection_temp <- detection_data %>%
      filter(!!sym(col_name_species) == target_species) %>%
      group_by(!!sym(col_name_station), !!sym(col_name_species)) %>%
      summarize(Y = sum(!!sym(col_name_y))) %>%
      left_join(station_data, by = col_name_station) %>%
      rename(Species = !!sym(col_name_species),
             Station = !!sym(col_name_station)) %>%
      arrange(!!sym(col_name_species))
  } else if(model == "RAD-REST") {
    detection_temp_1 <- detection_data %>%
      filter(!is.na(col_name_y)) %>%
      filter(!!sym(col_name_species) == target_species) %>%
      group_by(Station, !!sym(col_name_y)) %>%
      summarise(n = n(), .groups = 'drop') %>%
      left_join(station_data, by = col_name_station) %>%
      tidyr::complete(!!sym(col_name_y), fill = list(n = 0)) %>%
      tidyr::pivot_wider(
        names_from = !!sym(col_name_y),
        values_from = n,
        names_prefix = "y_",
        values_fill = list(n = 0)
      )
    detection_temp <- detection_data %>%
      filter(!!sym(col_name_species) == target_species) %>%
      group_by(Station) %>%
      summarise(N = n(), .groups = 'drop') %>%
      left_join(detection_temp_1, by = "Station")
  }
}
Station <- NULL
