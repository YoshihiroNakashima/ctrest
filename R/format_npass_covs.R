#' Aggregate the number of passes within focal areas for each station from the detection data
#'
#' @param detection_data A data frame containing information for individual videos, with six columns: camera station name (station), capture datetime (datetime), species name (species), number of passes through the focal area (nfocal), duration of stay within the focal area for each pass in a video (stay), and whether the observation of the pass was censored (cens).
#' @param station_data A data frame containing information for each camera station, with one row per station. It includes the station name (station), location coordinates, covariates, and so forth. Unlike the detection_data, it lists all the camera stations that have been set up.
#' @param col_name_station Column name containing station id info
#' @param col_name_species Column name containing species name detected
#' @param target_species The species name for which you want to estimate his density
#'
#' @return data frame containing at least 3 columns (station, species and npass)
#' @export
#' @import dplyr
#' @examples
#' format_npass_covs(
#' detection_data = detection_data,
#' station_data = station_data,
#' col_name_station = "station",
#' col_name_species = "species",
#' target_species = "A"
#' )

format_npass_covs <- function(detection_data,
                              station_data,
                              col_name_station = "station",
                              col_name_species = "species",
                              target_species = "A") {
  detection_temp <- detection_data %>%
    group_by(!!sym(col_name_station), !!sym(col_name_species)) %>%
    summarize(npass = sum(nfocal)) %>%
    left_join(station_data, by = col_name_station) %>%
    filter(!!sym(col_name_species) == target_species) %>%
    arrange(!!sym(col_name_species))

  print(detection_temp)
}

