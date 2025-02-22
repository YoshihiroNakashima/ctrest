#' Pre
#'
#' @param detection_data A data frame containing information for individual videos, with six columns: camera station name (station), capture datetime (datetime), species name (species), number of passes through the focal area (nfocal), duration of stay within the focal area for each pass in a video (stay), and whether the observation of the pass was censored (cens).
#' @param detection_data A data frame containing information for each camera station, with one row per station. It includes the station name (station), location coordinates, covariates, and so forth. Unlike the detection_data, it lists all the camera stations that have been set up.
#' @param station_data ""
#' @param col_name_station Column name containing station id info
#' @param col_name_species Column name containing species name detected
#' @param col_name_stay Column name containing staying time
#' @param col_name_cens Column name containing infomation on censored (1) or not (0)
#' @param target_species The species name for which you want to estimate his density
#' @return Dataframe contaning staying time
#' @import dplyr
#' @importFrom stats reformulate
#' @importFrom magrittr %>%
#' @export
#' @examples format_stay(
#'  detection_data = detection_data,
#'  station_data = station_data,
#'  col_name_station = "Station",
#'  col_name_species = "Species",
#'  col_name_stay = "Stay",
#'  col_name_cens = "Cens",
#'  target_species = "SP01")
format_stay <-
  function(detection_data,
           station_data,
           col_name_station = "Station",
           col_name_species = "Species",
           col_name_stay = "Stay",
           col_name_cens = "Cens",
           target_species = "A") {
    if(sum(detection_data %>% pull(col_name_cens), na.rm = TRUE) / nrow(detection_data) > 0.5){
      warning("too many censored data! Check 1 indicates censored (unobserved) data!", call. = FALSE)
    }
    stay_data <- detection_data %>%
      rename(Stay = !!sym(col_name_stay), Cens = !!sym(col_name_cens), Species = !!sym(col_name_species)) %>%
      filter(!is.na(Stay) & !is.na(Cens) & Species %in% target_species) %>%
      dplyr::select(Station, Species, Stay, Cens) %>%
      arrange(Species) %>%
      left_join(station_data, by = "Station")
    return(stay_data)
  }
Species <- NULL
