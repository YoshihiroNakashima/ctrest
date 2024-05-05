#' Prepare staying
#'
#' @param detection_data A data frame containing information for individual videos, with six columns: camera station name (station), capture datetime (datetime), species name (species), number of passes through the focal area (nfocal), duration of stay within the focal area for each pass in a video (stay), and whether the observation of the pass was censored (cens).
#' @param col_name_station Column name containing station id info
#' @param col_name_species Column name containing species name detected
#' @param col_name_stay Column name containing staying time
#' @param col_name_cens Column name containing infomation on censored (1) or not (0)
#' @param target_species The species name for which you want to estimate his density
#'
#' @return Dataframe contaning staying time
#' @import dplyr
#' @importFrom stats reformulate
#' @export
#' @examples format_stay(detection_data = detection_data,
#'  col_name_station = "Station",
#'  col_name_species = "Species",
#'  col_name_stay = "Stay",
#'  col_name_cens = "Cens",
#'  target_species = "A")
format_stay <-
  function(detection_data,
           col_name_station = "Station",
           col_name_species = "Species",
           col_name_stay = "Stay",
           col_name_cens = "Cens",
           target_species = "A") {
    if(sum(detection_data %>% pull(col_name_cens), na.rm = TRUE) / nrow(detection_data) > 0.5){
      warning("too many censored data! Check 1 indicates censored (unobserved) data!", call. = FALSE)
    }
    stay_data <- detection_data %>%
      filter(!is.na(!!sym(col_name_stay)) &
               !is.na(!!sym(col_name_cens)) &
               !!sym(col_name_species) == target_species) %>%
      select(
        Station = !!sym(col_name_station),
        Species = !!sym(col_name_species),
        Stay = !!sym(col_name_stay),
        Cens = !!sym(col_name_cens)
      )
  }
