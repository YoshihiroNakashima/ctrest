#' Randomly select videos with even distribution across camera stations
#'
#' This function selects `N_sampled` videos from a detection dataset. It filters for
#' independent events based on a time interval and ensures that the selection is
#' as evenly distributed across camera stations as possible.
#'
#' @param detection_data A data frame containing video information. Must include columns for:
#'   \itemize{
#'     \item Camera station ID
#'     \item Capture datetime (Character or POSIXct; handles Excel-formatted strings)
#'     \item Species name
#'   }
#' @param col_name_species String. The column name for species. Default is "Species".
#' @param col_name_station String. The column name for station IDs. Default is "Station".
#' @param col_name_datetime String. The column name for datetime. Default is "DateTimeCorrected".
#' @param N_sampled Integer. The total number of videos to be selected.
#' @param Indep_criteria Numeric. The time interval (in minutes) to define independent detections. Default is 30 minutes.
#' @param seed Integer. Optional random seed for reproducibility. Default is NULL.
#' @param target_species Character. A specific species to filter for. If NULL, all species are processed.
#'
#' @return A data frame containing the selected videos.
#'         If the available independent data is less than `N_sampled`, all available data is returned.
#'         Returns NULL if no valid data is found.
#'
#' @details
#' \strong{1. Date Parsing:}
#' Uses \code{lubridate::parse_date_time} to handle various date formats, including cases
#' where Excel strips the time component from "00:00:00" entries (leaving only "YYYY-MM-DD").
#'
#' \strong{2. Independence Criteria:}
#' Detections at the same station are considered independent if they occur at least
#' `Indep_criteria` minutes after the previous detection.
#'
#' \strong{3. Sampling Strategy:}
#' To prevent bias towards stations with high activity:
#' \enumerate{
#'   \item Assigns a random priority rank (1st, 2nd, 3rd...) to videos \emph{within} each station.
#'   \item Collects all "1st choice" videos from all stations.
#'   \item Then collects all "2nd choice" videos, and so on.
#'   \item Slices the top `N_sampled` records from this interleaved list.
#' }
#'
#' @export
#' @import dplyr
#' @importFrom lubridate parse_date_time
#' @importFrom rlang sym
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' selected_df <- select_videos(
#'   detection_data = my_camera_data,
#'   N_sampled = 100,
#'   Indep_criteria = 30,
#'   target_species = "Cervus nippon",
#'   seed = 123
#' )
#' }
select_videos <- function(detection_data,
                          col_name_species = "Species",
                          col_name_station = "Station",
                          col_name_datetime = "DateTimeCorrected",
                          N_sampled,
                          Indep_criteria = 30,
                          seed = NULL,
                          target_species = NULL) {

  # 1. Setup: Set seed and convert column names to symbols
  if (!is.null(seed)) set.seed(seed)

  species_sym <- rlang::sym(col_name_species)
  station_sym <- rlang::sym(col_name_station)
  datetime_sym <- rlang::sym(col_name_datetime)

  # 2. Preprocessing: Parse datetime and filter species
  # parse_date_time is used to handle mixed formats (e.g., "ymd HMS" vs "ymd" caused by Excel)
  df_processed <- detection_data %>%
    filter(!is.na(!!datetime_sym)) %>%
    mutate(
      parsed_dt = lubridate::parse_date_time(
        !!datetime_sym,
        orders = c("ymd HMS", "ymd HM", "ymd", "dmy HMS", "dmy", "ymd IMS p"),
        quiet = TRUE
      )
    ) %>%
    filter(!is.na(parsed_dt)) # Remove rows where date parsing failed

  # Filter by target species if specified
  if (!is.null(target_species)) {
    df_processed <- df_processed %>%
      filter(!!species_sym %in% target_species)
  }

  # Check if data exists after filtering
  if (nrow(df_processed) == 0) {
    message("No valid data found matching the criteria or date format.")
    return(NULL)
  }

  # 3. Identify Independent Events
  # Calculate time difference from the previous row per station
  independent_data <- df_processed %>%
    arrange(!!species_sym, !!station_sym, parsed_dt) %>%
    group_by(!!species_sym, !!station_sym) %>%
    mutate(
      # Calculate difference in minutes
      diff_time = as.numeric(difftime(parsed_dt, lag(parsed_dt), units = "mins")),
      # Mark as independent if it's the first record (NA) or exceeds the criteria
      is_independent = if_else(is.na(diff_time) | diff_time >= Indep_criteria, TRUE, FALSE)
    ) %>%
    filter(is_independent) %>%
    ungroup()

  # Return NULL if no independent events remain
  if (nrow(independent_data) == 0) {
    return(NULL)
  }

  # 4. Stratified Random Sampling
  # Select videos evenly across stations
  sampled_result <- independent_data %>%
    group_by(!!species_sym, !!station_sym) %>%
    mutate(
      # Assign a random rank to each video within its station
      # (e.g., if a station has 5 videos, they are randomly ranked 1 to 5)
      pick_order = sample(row_number())
    ) %>%
    ungroup() %>%
    # Sort the entire dataset:
    # 1. By 'pick_order': Ensures we take the "1st" video from every station before taking any "2nd" videos.
    # 2. By 'sample(row_number())': Randomly shuffles stations within the same rank level.
    arrange(pick_order, sample(row_number())) %>%
    # Select the top N videos
    slice_head(n = N_sampled) %>%
    # Remove temporary calculation columns
    select(-parsed_dt, -diff_time, -is_independent, -pick_order)

  return(sampled_result)
}
