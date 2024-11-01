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
#' @import dplyr tidyr
#' @examples
#' format_station_data(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_y = "y",
#'   model = "REST",
#'   target_species = "SP01"
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
      rename(Station = !!sym(col_name_station), Species = !!sym(col_name_species), y = !!sym(col_name_y)) %>%
      filter(Species %in% target_species) %>%
      group_by(Station, Species) %>%
      summarize(Y = sum(y)) %>%
      right_join(station_data, by = "Station") %>%
      replace_na(list(Y = 0)) %>%
      arrange(Station) %>%
      arrange(Species)

  } else if (model == "RAD-REST") {
    detection_temp_1 <- detection_data %>%
      rename(Station = !!sym(col_name_station), Species = !!sym(col_name_species), y = !!sym(col_name_y)) %>%
      filter(!is.na(y)) %>%
      filter(Species == target_species) %>%
      group_by(Station, y) %>%
      summarise(n = n(), .groups = 'drop') %>%
      tidyr::complete(y, fill = list(n = 0)) %>%
      drop_na("y") %>%
      tidyr::pivot_wider(
        names_from = y,
        values_from = n,
        names_prefix = "y_",
        values_fill = list(n = 0)
      ) %>%
      right_join(station_data, by = "Station") %>%
      mutate(Species = target_species)

    detection_temp <- detection_data %>%
      rename(Station = !!sym(col_name_station), Species = !!sym(col_name_species), y = !!sym(col_name_y)) %>%
      filter(Species == target_species) %>%
      group_by(Station) %>%
      summarise(N = n(), .groups = 'drop') %>%
      right_join(detection_temp_1, by = "Station")  %>%
      dplyr::mutate(across(everything(), \(x) tidyr::replace_na(x, 0)))
  }
}
Station <- Species <- y <- NULL
