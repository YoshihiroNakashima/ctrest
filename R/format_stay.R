#' Prepares data for analyzing animal staying time within the focal area
#'
#' @param detection_data A data frame containing information for individual videos.
#' @param station_data A data frame containing information for each camera station.
#' @param col_name_station String. Column name for station ID.
#' @param col_name_species String. Column name for species names.
#' @param col_name_stay String. Column name for staying time (numeric).
#' @param col_name_cens String. Column name for censoring indicator (0 or 1).
#'
#' @return A joined data frame (tibble).
#'
#' @import dplyr
#' @importFrom rlang sym !! .data
#' @export
format_stay <- function(detection_data,
                        station_data,
                        col_name_station = "Station",
                        col_name_species = "Species",
                        col_name_stay = "Stay",
                        col_name_cens = "Cens") {

  # 1. 基本的な入力チェック
  if (!is.data.frame(detection_data)) stop("'detection_data' must be a data frame.", call. = FALSE)
  if (!is.data.frame(station_data))   stop("'station_data' must be a data frame.", call. = FALSE)

  if (!(col_name_station %in% colnames(station_data))) {
    stop(sprintf("Column '%s' missing from 'station_data'.", col_name_station), call. = FALSE)
  }

  # 【重要】 station_data の重複チェック（デカルト積の防止）
  if (any(duplicated(station_data[[col_name_station]]))) {
    stop(sprintf("Error: 'station_data' must have exactly one row per station. Duplicate IDs detected in column '%s'.", col_name_station), call. = FALSE)
  }

  # 2. 列の存在確認
  req_det <- c(col_name_station, col_name_species, col_name_stay, col_name_cens)
  missing_det <- setdiff(req_det, colnames(detection_data))
  if (length(missing_det) > 0) {
    stop(paste("Missing columns in 'detection_data':", paste(missing_det, collapse = ", ")), call. = FALSE)
  }

  # 3. データ処理（メイン）
  res <- detection_data %>%
    dplyr::select(
      Station = !!rlang::sym(col_name_station),
      Species = !!rlang::sym(col_name_species),
      Stay    = !!rlang::sym(col_name_stay),
      Cens    = !!rlang::sym(col_name_cens)
    ) %>%
    dplyr::mutate(
      Station = as.character(.data$Station),
      Stay    = as.numeric(.data$Stay),
      Cens    = as.integer(.data$Cens)
    ) %>%
    dplyr::filter(!is.na(.data$Stay), !is.na(.data$Cens))

  # 4. Cens列のバリデーション
  if (nrow(res) > 0) {
    cens_values <- unique(res$Cens)
    if (!all(cens_values %in% c(0, 1))) {
      stop("Censoring column must contain only 0 and 1.", call. = FALSE)
    }
  }

  # 5. station_data の整形
  clean_station_data <- station_data
  if (col_name_station != "Station" && "Station" %in% colnames(clean_station_data)) {
    clean_station_data <- clean_station_data %>% dplyr::select(-.data$Station)
  }

  clean_station_data <- clean_station_data %>%
    dplyr::rename(Station = !!rlang::sym(col_name_station)) %>%
    dplyr::mutate(Station = as.character(.data$Station))

  # 予約名がstation_dataにあれば除去して重複を防ぐ
  reserved_names <- c("Species", "Stay", "Cens")
  dup_cols <- intersect(colnames(clean_station_data), reserved_names)
  if (length(dup_cols) > 0) {
    clean_station_data <- clean_station_data %>% dplyr::select(-dplyr::all_of(dup_cols))
  }

  # 6. 最終的な結合
  final_data <- res %>%
    dplyr::left_join(clean_station_data, by = "Station") %>%
    dplyr::arrange(.data$Species, .data$Station)

  return(final_data)
}

utils::globalVariables(c("Station", "Species", "Stay", "Cens"))
