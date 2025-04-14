#' A function to randomly select N_sampled videos such that the number of selected videos per camera station is as evenly distributed as possible.
#'
#' @param detection_data A data frame containing information for individual videos, with six columns:
#'   - Camera station ID (character)
#'   - Capture datetime (character, not used in this function)
#'   - Species name (character)
#'   - Number of passes through the focal area (numeric, not used in this function)
#'   - Staying time (seconds) within the focal area for each pass in a video (numeric)
#'   - Whether the observation of the pass was censored (1 = censored, 0 = observed)
#' @param col_name_species A string specifying the column name in `stay_data` that indicates the species name.
#' @param col_name_station A string specifying the column name containing camera station IDs.
#' @param N_sampled The total number of videos to be selected.
#' @param Indep_criteria A numeric value representing the interval (in minutes) used to determine whether detections are independent. The default is 30 minutes.
#' @param seed An integer value for setting the random seed.
#' @param target_species A character string specifying the species of interest. Only a single species can be specified.
#' @return A data frame containing the selected videos.
#' @export
#' @import dplyr
#' @importFrom lubridate ymd_hms
#' @importFrom tibble tibble
#' @importFrom rlang sym
#' @examples
#' select_videos(detection_data = detection_data,
#'               col_name_species = "Species",
#'               col_name_station = "Station",
#'               col_name_datetime = "DateTime",
#'               N_sampled = 100,
#'               Indep_criteria = 30,
#'               target_species = "SP01")


select_videos <- function(detection_data,
                          col_name_species = "Species",
                          col_name_station = "Station",
                          col_name_datetime = "DateTimeCorrected",
                          N_sampled,
                          Indep_criteria,
                          seed = NULL,
                          target_species = NULL) {

  # Set random seed
  if (!is.null(seed)) set.seed(seed)

  # Convert column names to symbols
  species_sym <- sym(col_name_species)
  station_sym <- sym(col_name_station)
  datetime_sym <- sym(col_name_datetime)

  # Filter by target_species if specified
  if (!is.null(target_species)) {
    detection_data <- detection_data %>%
      filter(!!species_sym %in% target_species)
  }

  # Process per species
  detection_data %>%
    group_by(!!species_sym) %>%
    group_split() %>%
    lapply(function(species_data) {

      # Process per station
      species_data %>%
        arrange(!!datetime_sym) %>%
        group_by(!!station_sym) %>%
        group_split() %>%
        lapply(function(station_data) {

          # Select videos that meet independence criteria
          selected_indices <- c()
          last_datetime <- ymd_hms("1900-01-01 00:00:00") # Initialize with an old datetime

          for (i in 1:nrow(station_data)) {
            current_datetime <- station_data[[col_name_datetime]][i]
            if (difftime(current_datetime, last_datetime, units = "mins") >= Indep_criteria) {
              selected_indices <- c(selected_indices, i)
              last_datetime <- current_datetime
            }
          }
          if (length(selected_indices) > 0){
            station_data[selected_indices, ]
          } else {
            NULL
          }
        }) %>%
        bind_rows() -> independent_data

      if (nrow(independent_data) == 0) return(NULL) # If no independent data available

      unique_stations <- unique(independent_data[[col_name_station]])
      num_unique_stations <- length(unique_stations)

      if (num_unique_stations >= N_sampled) {
        # If the number of stations is greater than or equal to N_sampled
        sampled_stations <- sample(unique_stations, N_sampled)
        sampled_data <- independent_data %>%
          filter(!!station_sym %in% sampled_stations) %>%
          group_by(!!station_sym) %>%
          slice_sample(n = 1) %>%
          ungroup()
        return(sampled_data)

      } else {
        # If the number of stations is less than N_sampled
        sampled_data <- tibble()
        available_data <- independent_data # Store available data for selection
        while(nrow(sampled_data) < N_sampled && nrow(available_data) > 0){ # Add condition to stop when available_data is empty
          temp_sampled <- available_data %>%
            group_by(!!station_sym) %>%
            slice_sample(n=1) %>%
            ungroup()
          sampled_data <- bind_rows(sampled_data, temp_sampled)
          # Remove selected data from available_data
          available_data <- anti_join(available_data, temp_sampled, by = c(col_name_station, col_name_datetime))
        }
        sampled_data <- sampled_data[1:N_sampled,] # Trim excess rows beyond N_sampled
        return(sampled_data)
      }
    }) %>%
    bind_rows()
}
