#' Aggregate the number of animal passes within focal areas for each station from detection data.
#'
#' @param detection_data A data frame containing information for individual videos.
#' @param station_data A data frame containing information for each camera station.
#' @param col_name_station A string specifying the column name containing station ID information.
#' @param col_name_species A string specifying the column name containing detected species names.
#' @param col_name_y A string specifying the column name containing the number of animal passes.
#' @param model A string specifying the model type, either "REST" or "RAD-REST".
#' @return A data frame with aggregated detection counts joined with station data.
#' @export
#' @import dplyr tidyr stringr rlang

format_station_data <- function(detection_data,
                                station_data,
                                col_name_station,
                                col_name_species,
                                col_name_y,
                                model) {

  # Check model name
  if (!(model %in% c("REST", "RAD-REST"))) {
    stop("Invalid model name! 'model' must be either 'REST' or 'RAD-REST'.", call. = FALSE)
  }

  # Check if required columns exist in detection_data
  required_cols <- c(col_name_station, col_name_species, col_name_y)
  missing_cols <- required_cols[!required_cols %in% colnames(detection_data)]
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing from 'detection_data':", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Check if required column exists in station_data
  if (!col_name_station %in% colnames(station_data)) {
    stop(paste("The station ID column", col_name_station, "is missing from 'station_data'."), call. = FALSE)
  }

  # Prepare a full grid of Station x Species
  # Ensure we use the exact station list from station_data (to include stations with 0 detections)
  all_stations <- unique(station_data[[col_name_station]])
  all_species  <- unique(detection_data[[col_name_species]])

  detection_temp_0 <- tidyr::crossing(
    Station = all_stations,
    Species = all_species
  )

  # --- REST Model Logic ---
  if (model == "REST") {

    # Check NA proportion
    na_ratio <- sum(is.na(detection_data[[col_name_y]])) / nrow(detection_data)
    if (na_ratio > 0.5) {
      warning("More than 50% of 'col_name_y' values are NA in 'REST' model. These values will be replaced with 1.", call. = FALSE)
    }

    detection_aggregated <- detection_data %>%
      rename(Station = !!sym(col_name_station),
             Species = !!sym(col_name_species),
             y = !!sym(col_name_y)) %>%
      mutate(y = ifelse(is.na(y), 1, y)) %>%
      group_by(Station, Species) %>%
      summarize(Y = sum(y), .groups = 'drop') # Calculate Sum per Station-Species

    # Join with the full grid
    detection_final <- detection_temp_0 %>%
      left_join(detection_aggregated, by = c("Station", "Species")) %>%
      mutate(Y = replace_na(Y, 0)) %>% # Fill 0 only for the count column
      arrange(Station, Species)

    # --- RAD-REST Model Logic ---
  } else if (model == "RAD-REST") {

    # 1. Calculate N (Total detections per Station-Species)
    detection_counts_N <- detection_data %>%
      rename(Station = !!sym(col_name_station),
             Species = !!sym(col_name_species)) %>%
      group_by(Station, Species) %>%
      summarise(N = n(), .groups = 'drop')

    # 2. Calculate y_X (Frequency of each pass count)
    detection_counts_y <- detection_data %>%
      rename(Station = !!sym(col_name_station),
             Species = !!sym(col_name_species),
             y = !!sym(col_name_y)) %>%
      filter(!is.na(y)) %>%
      group_by(Station, Species, y) %>%
      summarise(n = n(), .groups = 'drop') %>%
      # Use pivot_wider to create columns y_0, y_1, ...
      pivot_wider(
        names_from = y,
        values_from = n,
        names_prefix = "y_",
        values_fill = list(n = 0)
      )

    # Join parts together
    detection_final <- detection_temp_0 %>%
      left_join(detection_counts_N, by = c("Station", "Species")) %>%
      left_join(detection_counts_y, by = c("Station", "Species"))

    # Fill NAs with 0 only for N and y_ columns
    # Identify y_ columns dynamically
    y_cols <- colnames(detection_final)[str_detect(colnames(detection_final), "^y_\\d+$")]
    vars_to_fill <- c("N", y_cols)

    detection_final <- detection_final %>%
      mutate(across(all_of(vars_to_fill), ~ replace_na(.x, 0)))

    # Ensure all intermediate y columns exist (e.g., if we have y_0 and y_2, we need y_1)
    # Extract numeric suffixes
    if (length(y_cols) > 0) {
      indices <- as.integer(str_extract(y_cols, "\\d+"))
      max_idx <- max(indices, na.rm = TRUE)
      full_indices <- 0:max_idx

      # Identify missing columns
      missing_y_cols <- paste0("y_", full_indices[!full_indices %in% indices])

      # Add missing columns filled with 0
      if (length(missing_y_cols) > 0) {
        new_cols <- rep(0, length(missing_y_cols))
        names(new_cols) <- missing_y_cols
        detection_final <- detection_final %>%
          tibble::add_column(!!!new_cols)
      }

      # Sort columns: Station, Species, N, then y_0, y_1 ...
      # Use stringr::str_sort for natural sorting (y_2 before y_10)
      sorted_y_cols <- str_sort(colnames(detection_final)[str_detect(colnames(detection_final), "^y_\\d+$")], numeric = TRUE)

      detection_final <- detection_final %>%
        select(Station, Species, N, all_of(sorted_y_cols))
    }
  }

  # --- Final Merge with Station Data ---
  # We renamed the ID column in detection_final to "Station".
  # We join with station_data using the user-provided column name mapped to "Station".

  join_by_clause <- setNames(col_name_station, "Station")

  detection_final <- detection_final %>%
    left_join(station_data, by = join_by_clause) %>%
    arrange(Station, Species)

  return(detection_final)
}
