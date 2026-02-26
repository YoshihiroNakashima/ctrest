#' Calculate Camera Trapping Effort (Days) from Detection Data
#'
#' @description Computes the total camera trapping effort (in days) for each station based on detection data.
#' @param detection_data A data frame containing information on individual video records.
#' @param station_data_formatted A data frame returned by the `format_station_data()` function.
#' @param col_name_station A string specifying the column name containing camera station IDs.
#' @param col_name_datetime A string specifying the column name containing datetime information.
#' @param col_name_term (Optional) A string specifying the column name indicating different survey periods.
#' @param plot Logical. If TRUE, a plot showing the camera operation periods will be displayed.
#' @return A data frame with an additional column `Effort`. Stations with `Effort == 0` or no detections are removed.
#' @export
#' @import dplyr ggplot2 lubridate stringr rlang
add_effort <- function(detection_data,
                       station_data_formatted,
                       col_name_station,
                       col_name_datetime,
                       col_name_term = NULL,
                       plot = FALSE) {

  if (!is.data.frame(detection_data)) stop("'detection_data' must be a data frame.", call. = FALSE)
  if (!is.data.frame(station_data_formatted)) stop("'station_data_formatted' must be a data frame.", call. = FALSE)

  req_cols <- c(col_name_station, col_name_datetime)
  if (!is.null(col_name_term)) req_cols <- c(req_cols, col_name_term)

  missing_cols <- setdiff(req_cols, colnames(detection_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in detection_data:", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  if (!"Station" %in% colnames(station_data_formatted)) {
    stop("'station_data_formatted' must contain a 'Station' column. Did you run format_station_data() first?", call. = FALSE)
  }

  det_clean <- detection_data %>%
    dplyr::select(
      Station  = !!rlang::sym(col_name_station),
      DateTime = !!rlang::sym(col_name_datetime),
      dplyr::any_of(if (!is.null(col_name_term)) c(Term = col_name_term) else character())
    ) %>%
    dplyr::mutate(Station = as.character(.data$Station))

  if (!inherits(det_clean$DateTime, "POSIXt")) {
    det_clean <- det_clean %>%
      dplyr::mutate(
        DateTime = lubridate::parse_date_time(as.character(.data$DateTime), orders = c("Ymd HMS", "Ymd HM", "Ymd", "Ymd H"))
      )
  }

  if (any(is.na(det_clean$DateTime))) {
    stop("Error: Failed to parse some datetimes. Ensure the format is YYYY-MM-DD HH:MM:SS or similar.", call. = FALSE)
  }

  if (is.null(col_name_term)) {
    effort_temp <- det_clean %>% dplyr::group_by(.data$Station)
  } else {
    effort_temp <- det_clean %>% dplyr::group_by(.data$Station, .data$Term)
  }

  effort_temp <- effort_temp %>%
    dplyr::summarize(
      Start = min(.data$DateTime, na.rm = TRUE),
      End   = max(.data$DateTime, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(effort = as.numeric(difftime(.data$End, .data$Start, units = "days")))

  effort_agg <- effort_temp %>%
    dplyr::group_by(.data$Station) %>%
    dplyr::summarize(Effort = sum(.data$effort, na.rm = TRUE), .groups = 'drop')

  clean_station <- station_data_formatted
  if ("Effort" %in% colnames(clean_station)) {
    clean_station <- clean_station %>% dplyr::select(-.data$Effort)
  }

  final_data <- clean_station %>%
    dplyr::left_join(effort_agg, by = "Station")

  invalid_stations <- final_data %>%
    dplyr::filter(is.na(.data$Effort) | .data$Effort == 0) %>%
    dplyr::pull(.data$Station) %>%
    unique()

  if (length(invalid_stations) > 0) {
    warning(sprintf(
      "Sampling effort for the following stations was calculated as 0 (or had no detections). They will be removed from the output: %s",
      paste(invalid_stations, collapse = ", ")
    ), call. = FALSE)
  }

  final_data <- final_data %>%
    dplyr::filter(!is.na(.data$Effort), .data$Effort > 0)

  if (plot) {
    g <- ggplot2::ggplot(data = effort_temp) +
      ggplot2::geom_segment(
        ggplot2::aes(x = .data$Start, xend = .data$End, y = .data$Station, yend = .data$Station),
        alpha = 0.3, linewidth = 2.5
      ) +
      ggplot2::geom_point(ggplot2::aes(x = .data$Start, y = .data$Station), shape = 4, col = "red", alpha = 0.8) +
      ggplot2::geom_point(ggplot2::aes(x = .data$End, y = .data$Station), shape = 4, col = "blue", alpha = 0.8) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = max(5, 12 - nrow(effort_temp)/10))) +
      ggplot2::xlab("Survey Period") +
      ggplot2::ylab("Station") +
      ggplot2::ggtitle("Camera Trapping Operation Periods")
    print(g)
  }

  return(final_data)
}

utils::globalVariables(c("Station", "DateTime", "Term", "Start", "End", "effort", "Effort"))
