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

  # 1. 基本的な入力チェック
  if (!is.data.frame(detection_data)) stop("'detection_data' must be a data frame.", call. = FALSE)
  if (!is.data.frame(station_data))   stop("'station_data' must be a data frame.", call. = FALSE)
  if (!(model %in% c("REST", "RAD-REST"))) {
    stop("Invalid model name! 'model' must be either 'REST' or 'RAD-REST'.", call. = FALSE)
  }

  req_cols <- c(col_name_station, col_name_species, col_name_y)
  missing_cols <- setdiff(req_cols, colnames(detection_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in 'detection_data':", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  if (!col_name_station %in% colnames(station_data)) {
    stop(paste("The station ID column", col_name_station, "is missing from 'station_data'."), call. = FALSE)
  }

  # 【重要】 station_data の重複チェック（デカルト積の防止）
  if (any(duplicated(station_data[[col_name_station]]))) {
    stop(sprintf("Error: 'station_data' must have exactly one row per station. Duplicate IDs detected in column '%s'. Did you accidentally pass an already merged dataset?", col_name_station), call. = FALSE)
  }

  # 2. データの分離と型定義
  det_clean <- detection_data %>%
    dplyr::select(
      Station = !!rlang::sym(col_name_station),
      Species = !!rlang::sym(col_name_species),
      y       = !!rlang::sym(col_name_y)
    ) %>%
    dplyr::mutate(
      Station = as.character(.data$Station),
      Species = as.character(.data$Species),
      y       = as.numeric(.data$y)
    )

  # 3. station_data の整形と衝突回避
  clean_station <- station_data
  if (col_name_station != "Station" && "Station" %in% colnames(clean_station)) {
    clean_station <- clean_station %>% dplyr::select(-.data$Station)
  }
  clean_station <- clean_station %>%
    dplyr::rename(Station = !!rlang::sym(col_name_station)) %>%
    dplyr::mutate(Station = as.character(.data$Station))

  reserved_regex <- "^(Species|Y|N|y|y_\\d+)$"
  dup_cols <- grep(reserved_regex, colnames(clean_station), value = TRUE)
  dup_cols <- setdiff(dup_cols, "Station")
  if (length(dup_cols) > 0) {
    clean_station <- clean_station %>% dplyr::select(-dplyr::all_of(dup_cols))
  }

  # 4. 全ステーション×種のグリッド作成
  all_stations <- unique(clean_station$Station)
  all_species  <- unique(stats::na.omit(det_clean$Species))

  detection_grid <- tidyr::crossing(
    Station = all_stations,
    Species = all_species
  )

  # 5. モデル別の集計ロジック
  if (model == "REST") {
    det_agg <- det_clean %>%
      dplyr::filter(!is.na(.data$Species)) %>%
      dplyr::mutate(y = ifelse(is.na(.data$y), 1, .data$y)) %>%
      dplyr::group_by(.data$Station, .data$Species) %>%
      dplyr::summarize(Y = sum(.data$y), .groups = 'drop')

    detection_final <- detection_grid %>%
      dplyr::left_join(det_agg, by = c("Station", "Species")) %>%
      dplyr::mutate(Y = tidyr::replace_na(.data$Y, 0))

  } else if (model == "RAD-REST") {
    det_N <- det_clean %>%
      dplyr::filter(!is.na(.data$Species)) %>%
      dplyr::group_by(.data$Station, .data$Species) %>%
      dplyr::summarise(N = dplyr::n(), .groups = 'drop')

    det_y <- det_clean %>%
      dplyr::filter(!is.na(.data$Species), !is.na(.data$y)) %>%
      dplyr::group_by(.data$Station, .data$Species, .data$y) %>%
      dplyr::summarise(n = dplyr::n(), .groups = 'drop') %>%
      tidyr::pivot_wider(
        names_from = .data$y,
        values_from = .data$n,
        names_prefix = "y_",
        values_fill = 0
      )

    detection_final <- detection_grid %>%
      dplyr::left_join(det_N, by = c("Station", "Species")) %>%
      dplyr::left_join(det_y, by = c("Station", "Species"))

    y_cols <- grep("^y_\\d+$", colnames(detection_final), value = TRUE)
    vars_to_fill <- c("N", y_cols)
    detection_final <- detection_final %>%
      dplyr::mutate(dplyr::across(dplyr::any_of(vars_to_fill), ~ tidyr::replace_na(.x, 0)))

    indices <- if (length(y_cols) > 0) as.integer(stringr::str_extract(y_cols, "\\d+")) else numeric(0)
    max_idx <- max(c(0, indices), na.rm = TRUE)
    full_indices <- 0:max_idx
    missing_indices <- setdiff(full_indices, indices)

    if (length(missing_indices) > 0) {
      missing_cols <- paste0("y_", missing_indices)
      for (col in missing_cols) {
        detection_final[[col]] <- 0
      }
    }

    sorted_y_cols <- stringr::str_sort(grep("^y_\\d+$", colnames(detection_final), value = TRUE), numeric = TRUE)
    detection_final <- detection_final %>%
      dplyr::select(.data$Station, .data$Species, .data$N, dplyr::all_of(sorted_y_cols))
  }

  # 6. station_data との最終結合
  final_joined <- detection_final %>%
    dplyr::left_join(clean_station, by = "Station") %>%
    dplyr::arrange(.data$Species, .data$Station)

  return(final_joined)
}

utils::globalVariables(c("Station", "Species", "y", "Y", "N", "n"))
