#' Calculate Camera Trapping Effort (Days) from Detection Data
#'
#' @description
#' Computes the total camera trapping effort (in days) for each station based
#' on detection records. When \code{col_name_term} is specified, effort is
#' calculated per term as the difference between the last and first detection
#' within that term; values across terms are then summed to give the total
#' effort per station. When \code{col_name_term} is \code{NULL}, effort is
#' calculated as the difference between the last and first detection across
#' the entire survey period for each station.
#'
#' @param detection_data A data frame containing individual detection records.
#' @param station_data_formatted A data frame returned by \code{format_station_data()}.
#'   Must contain a \code{Station} column.
#' @param col_name_station A string specifying the column name containing camera
#'   station IDs.
#' @param col_name_datetime A string specifying the column name containing
#'   datetime information. The format \code{"YYYY-MM-DD HH:MM:SS"} or similar
#'   is recommended.
#' @param col_name_term (Optional) A string specifying the column name indicating
#'   different survey terms. When provided, effort is computed separately for
#'   each term and summed per station.
#' @param plot Logical. If \code{TRUE}, a Gantt-style plot showing the camera
#'   operation periods per station is displayed. Default is \code{FALSE}.
#'
#' @return A data frame based on \code{station_data_formatted} with an additional
#'   \code{Effort} column (in days). Stations with \code{Effort == 0} or no
#'   detections are removed from the output with a warning.
#'
#' @details
#' Effort is approximated as \emph{last detection - first detection} (in days)
#' within each grouping. This approach is intended for use when actual camera
#' operation logs are unavailable. When \code{col_name_term} is specified,
#' inactive periods between terms are not counted toward total effort.
#'
#' @export
#'
#' @importFrom dplyr select mutate group_by summarize left_join filter pull
#'   any_of all_of across
#' @importFrom ggplot2 ggplot aes geom_segment geom_point theme element_text
#'   xlab ylab ggtitle
#' @importFrom lubridate parse_date_time
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' station_with_effort <- add_effort(
#'   detection_data         = my_detections,
#'   station_data_formatted = my_stations,
#'   col_name_station       = "Station",
#'   col_name_datetime      = "DateTime",
#'   col_name_term          = "Term",
#'   plot                   = TRUE
#' )
#' }
add_effort <- function(detection_data,
                       station_data_formatted,
                       col_name_station,
                       col_name_datetime,
                       col_name_term = NULL,
                       plot = FALSE) {

  # --- Input validation -------------------------------------------------------

  if (!is.data.frame(detection_data)) {
    stop("'detection_data' must be a data frame.", call. = FALSE)
  }
  if (!is.data.frame(station_data_formatted)) {
    stop("'station_data_formatted' must be a data frame.", call. = FALSE)
  }

  req_cols <- c(col_name_station, col_name_datetime)
  if (!is.null(col_name_term)) req_cols <- c(req_cols, col_name_term)

  missing_cols <- setdiff(req_cols, colnames(detection_data))
  if (length(missing_cols) > 0) {
    stop(
      paste("Missing columns in detection_data:", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!"Station" %in% colnames(station_data_formatted)) {
    stop(
      "'station_data_formatted' must contain a 'Station' column. Did you run format_station_data() first?",
      call. = FALSE
    )
  }

  # --- Column selection and renaming ------------------------------------------

  # Build a named vector for dplyr::all_of() renaming;
  # include Term only when col_name_term is provided
  select_cols <- c(
    Station  = col_name_station,
    DateTime = col_name_datetime
  )
  if (!is.null(col_name_term)) {
    select_cols <- c(select_cols, Term = col_name_term)
  }

  det_clean <- detection_data |>
    dplyr::select(dplyr::all_of(select_cols)) |>
    dplyr::mutate(Station = as.character(.data$Station))

  # --- Datetime parsing -------------------------------------------------------

  if (!inherits(det_clean$DateTime, "POSIXt")) {
    det_clean <- det_clean |>
      dplyr::mutate(
        DateTime = lubridate::parse_date_time(
          as.character(.data$DateTime),
          orders = c("Ymd HMS", "Ymd HM", "Ymd", "Ymd H")
        )
      )
  }

  if (anyNA(det_clean$DateTime)) {
    stop(
      "Failed to parse some datetime values. Ensure the format is 'YYYY-MM-DD HH:MM:SS' or similar.",
      call. = FALSE
    )
  }

  # --- Effort calculation -----------------------------------------------------
  # With col_name_term : group by Station x Term, compute (last - first) per term,
  #                       then sum across terms per station
  # Without col_name_term: group by Station only, compute (last - first) overall

  group_vars <- if (!is.null(col_name_term)) c("Station", "Term") else "Station"

  effort_by_group <- det_clean |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarize(
      Start  = min(.data$DateTime, na.rm = TRUE),
      End    = max(.data$DateTime, na.rm = TRUE),
      effort = as.numeric(difftime(.data$End, .data$Start, units = "days")),
      .groups = "drop"
    )

  # Aggregate per station by summing effort across terms
  effort_agg <- effort_by_group |>
    dplyr::group_by(.data$Station) |>
    dplyr::summarize(Effort = sum(.data$effort, na.rm = TRUE), .groups = "drop")

  # --- Join with station_data_formatted ---------------------------------------

  # Drop any pre-existing Effort column to avoid duplication
  clean_station <- station_data_formatted
  if ("Effort" %in% colnames(clean_station)) {
    clean_station <- clean_station |> dplyr::select(-"Effort")
  }

  final_data <- clean_station |>
    dplyr::left_join(effort_agg, by = "Station")

  # --- Warn and remove stations with zero or missing effort -------------------

  invalid_stations <- final_data |>
    dplyr::filter(is.na(.data$Effort) | .data$Effort == 0) |>
    dplyr::pull(.data$Station) |>
    unique()

  if (length(invalid_stations) > 0) {
    warning(
      sprintf(
        "Sampling effort for the following stations was 0 or had no detections and will be removed: %s",
        paste(invalid_stations, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  final_data <- final_data |>
    dplyr::filter(!is.na(.data$Effort), .data$Effort > 0)

  # --- Optional plot ----------------------------------------------------------

  if (plot) {
    # Scale y-axis label size by the number of station-term groups
    label_size <- max(5, 12 - nrow(effort_by_group) / 10)

    # Colour segments by Term when col_name_term is provided
    if (!is.null(col_name_term)) {
      g <- ggplot2::ggplot(data = effort_by_group) +
        ggplot2::geom_segment(
          ggplot2::aes(
            x      = .data$Start,
            xend   = .data$End,
            y      = .data$Station,
            yend   = .data$Station,
            colour = .data$Term
          ),
          alpha = 0.6, linewidth = 2.5
        )
    } else {
      g <- ggplot2::ggplot(data = effort_by_group) +
        ggplot2::geom_segment(
          ggplot2::aes(
            x    = .data$Start,
            xend = .data$End,
            y    = .data$Station,
            yend = .data$Station
          ),
          alpha = 0.3, linewidth = 2.5
        )
    }

    g <- g +
      ggplot2::geom_point(
        ggplot2::aes(x = .data$Start, y = .data$Station),
        shape = 4, colour = "red", alpha = 0.8
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = .data$End, y = .data$Station),
        shape = 4, colour = "blue", alpha = 0.8
      ) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = label_size)
      ) +
      ggplot2::xlab("Survey Period") +
      ggplot2::ylab("Station") +
      ggplot2::ggtitle("Camera Trapping Operation Periods")

    print(g)
  }

  return(final_data)
}

utils::globalVariables(c("Station", "DateTime", "Term", "Start", "End", "effort", "Effort"))
