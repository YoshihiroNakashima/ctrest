#' Aggregate the number of animal passes within focal areas for each station from detection data.
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
#' @param col_name_y A string specifying the column name containing the number of animal passes within a focal area.
#'   - If `model = "REST"`, this column must have values for all detections. Any `NA` values will be replaced with `1`.
#'   - If `model = "RAD-REST"`, only measured detections have count values, and unmeasured ones should be `NA`.
#' @param target_species A character vector specifying the species name(s) for which density estimation is desired.
#' @param model A string specifying the model type, either `"REST"` or `"RAD-REST"`.
#' @return A data frame with aggregated detection counts. The structure varies depending on the model:
#'   - If `model = "REST"`, the data frame contains:
#'     - `Station` (character): Camera station ID.
#'     - `Species` (character): Species name.
#'     - `Y` (numeric): The total number of passes through the focal area for each station-species pair.
#'     - All columns from `station_data` merged by `Station`.
#'   - If `model = "RAD-REST"`, the data frame contains:
#'     - `Station` (character): Camera station ID.
#'     - `Species` (character): Species name.
#'     - `N` (numeric): The total number of detected videos for each station-species pair.
#'     - `y_X` columns (`y_0`, `y_1`, ..., `y_max`): Number of detections classified by pass count.
#'     - All columns from `station_data` merged by `Station`.
#'   If no matching `Station` is found in `station_data`, the corresponding values will be `NA`.
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

  # Check model name
  if (!(model %in% c("REST", "RAD-REST"))) {
    stop("Invalid model name! 'model' must be either 'REST' or 'RAD-REST'.", call. = FALSE)
  }

  # Check if required columns exist
  required_cols <- c(col_name_station, col_name_species, col_name_y)
  missing_cols <- required_cols[!required_cols %in% colnames(detection_data)]
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing from 'detection_data':", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Check if target species exist in detection_data
  if (!any(detection_data[[col_name_species]] %in% target_species)) {
    warning("None of the target species are present in 'detection_data'. Returning an empty data frame.", call. = FALSE)
    return(data.frame())
  }

  # Check NA proportion in col_name_y for "REST"
  if (model == "REST") {
    na_ratio <- sum(is.na(detection_data[[col_name_y]])) / nrow(detection_data)
    if (na_ratio > 0.5) {
      warning("More than 50% of 'col_name_y' values are NA in 'REST' model. These values will be replaced with 1.", call. = FALSE)
    }
  }

  if (model == "REST") {
    detection_temp_0 <- crossing(
      Species = target_species,
      Station = station_data[[col_name_station]]
    )

    detection_temp <- detection_data %>%
      rename(Station = !!sym(col_name_station), Species = !!sym(col_name_species), y = !!sym(col_name_y)) %>%
      mutate(y = ifelse(is.na(y), 1, y)) %>%
      filter(Species %in% target_species) %>%
      group_by(Station, Species) %>%
      summarize(Y = sum(y), .groups = 'drop') %>%
      select(- Species) %>%
      right_join(detection_temp_0, by = "Station") %>%
      replace_na(list(Y = 0)) %>%
      arrange(Station) %>%
      left_join(station_data, by = c("Station"))

  } else if (model == "RAD-REST") {
    detection_temp_0 <- crossing(
      Species = target_species,
      Station = station_data[[col_name_station]]
    )

    detection_temp_1 <- detection_data %>%
      rename(Station = !!sym(col_name_station), Species = !!sym(col_name_species), y = !!sym(col_name_y)) %>%
      filter(!is.na(y)) %>%
      filter(Species %in% target_species) %>%
      group_by(Station, y, Species) %>%
      summarise(n = n(), .groups = 'drop') %>%
      tidyr::complete(y, fill = list(n = 0)) %>%
      drop_na("y") %>%
      tidyr::pivot_wider(
        names_from = y,
        values_from = n,
        names_prefix = "y_",
        values_fill = list(n = 0)
      ) %>%
      right_join(detection_temp_0, by = c("Station", "Species"))

    # Dynamically adjust missing columns
    index <- detection_temp_1 %>%
      select(starts_with("y_")) %>%
      names() %>%
      str_extract("\\d+") %>%
      as.integer()

    if (length(index) != max(index) + 1) {
      add_col <- paste0("y_", (0:max(index))[!(0:max(index)) %in% index])
      detection_temp_1 <- detection_temp_1 %>%
        mutate(!!!rlang::set_names(lapply(add_col, function(col) rep(0, nrow(detection_temp_1))), add_col)) %>%
        select(Station, Species,
               gtools::mixedsort(names(select(., starts_with("y_")))),
               everything())
    }

    detection_temp <- detection_data %>%
      rename(Station = !!sym(col_name_station), Species = !!sym(col_name_species), y = !!sym(col_name_y)) %>%
      filter(Species %in% target_species) %>%
      group_by(Station, Species) %>%
      summarise(N = n(), .groups = 'drop') %>%
      left_join(detection_temp_1, by = c("Station", "Species")) %>%
      right_join(detection_temp_0, by = c("Station", "Species")) %>%
      mutate(across(everything(), ~ replace_na(.x, 0))) %>%
      arrange(Species, Station) %>%
      left_join(station_data, by = "Station")
  }

  return(detection_temp)
}
Station <- Species <- y <- NULL
