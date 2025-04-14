#' Prepares data for analyzing animal staying time within the focal area
#'
#' @param detection_data A data frame containing information for individual videos, with six columns:
#'   - Camera station ID (character)
#'   - Capture datetime (character, not used in this function)
#'   - Species name (character)
#'   - Number of passes through the focal area (numeric, not used in this function)
#'   - Staying time (seconds) within the focal area for each pass in a video (numeric)
#'   - Whether the observation of the pass was censored (1 = censored, 0 = observed)
#' @param station_data A data frame containing information for each camera station, with one row per station.
#'   It must include at least the camera station IDs (character) and may contain additional columns such as
#'   location coordinates, covariates, and other relevant information. Unlike `detection_data`, this data frame
#'   lists all the camera stations that have been set up.
#' @param col_name_station A string specifying the column name containing station ID information.
#' @param col_name_species A string specifying the column name containing detected species names.
#' @param col_name_stay A string specifying the column name containing staying time.
#' @param col_name_cens A string specifying the column name indicating whether the observation is censored (1) or not (0).
#' @param target_species A character vector specifying the species name(s) for which density estimation is desired.
#'
#' @return A data frame with staying time information for the specified species.
#'   The returned data frame contains the following columns:
#'   - `Station` (character): Camera station ID where the detection occurred.
#'   - `Species` (character): Species name.
#'   - `Stay` (numeric): Staying time (in seconds) within the focal area for each pass.
#'   - `Cens` (numeric): Indicator of whether the observation was censored (1 = censored, 0 = observed).
#'   Additionally, all columns from `station_data` are included, matched by `Station`.
#'   If no matching `Station` is found in `station_data`, those values will be `NA`.
#'   If none of the `target_species` are found in `detection_data`, an empty data frame with
#'   the same column structure will be returned.
#'
#' @import dplyr
#' @importFrom stats reformulate na.omit
#' @importFrom magrittr %>%
#' @export
#' @examples
#' format_stay(
#'  detection_data = detection_data,
#'  station_data = station_data,
#'  col_name_station = "Station",
#'  col_name_species = "Species",
#'  col_name_stay = "Stay",
#'  col_name_cens = "Cens",
#'  target_species = "SP01")
format_stay <- function(detection_data,
                        station_data,
                        col_name_station = "Station",
                        col_name_species = "Species",
                        col_name_stay = "Stay",
                        col_name_cens = "Cens",
                        target_species = "A") {

  if (!is.data.frame(detection_data)) {
    stop("detection_data must be a data frame.", call. = FALSE)
  }
  if (!is.data.frame(station_data)) {
    stop("station_data must be a data frame.", call. = FALSE)
  }

  required_cols <- c(col_name_station, col_name_species, col_name_stay, col_name_cens)
  missing_cols <- setdiff(required_cols, colnames(detection_data))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing from detection_data:", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  if (!(col_name_station %in% colnames(station_data))) {
    stop(paste("Column", col_name_station, "is missing from station_data."), call. = FALSE)
  }

  cens_values <- unique(na.omit(detection_data[[col_name_cens]]))
  if (!is.numeric(detection_data[[col_name_cens]]) || !all(cens_values %in% c(0, 1))) {
    stop(paste("Column", col_name_cens, "must be numeric and contain only 0 (observed) and 1 (censored)."), call. = FALSE)
  }

  if (!is.character(target_species)) {
    stop("target_species must be a character vector.", call. = FALSE)
  }

  if (!any(detection_data[[col_name_species]] %in% target_species)) {
    warning("None of the specified target_species were found in detection_data.", call. = FALSE)
  }

  censored_ratio <- sum(detection_data[[col_name_cens]], na.rm = TRUE) / nrow(detection_data)
  if (censored_ratio > 0.5) {
    warning("Too many censored data! Check that 1 indicates censored (unobserved) data!", call. = FALSE)
  }

  stay_data <- detection_data %>%
    rename(Stay = !!sym(col_name_stay),
           Cens = !!sym(col_name_cens),
           Species = !!sym(col_name_species)) %>%
    filter(!is.na(Stay) & !is.na(Cens) & Species %in% target_species) %>%
    dplyr::select(Station = !!sym(col_name_station), Species, Stay, Cens) %>%
    arrange(Species) %>%
    left_join(station_data, by = c("Station" = col_name_station))

  return(stay_data)
}

Species <- NULL
