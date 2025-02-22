#' Aggregate the number of passes within focal areas for each station from the detection data
#'
#' @param detection_data A data frame containing information for individual videos, with six columns: camera station name (Station), capture datetime (Datetime), species name (Species), number of passes through the focal area (y), duration of stay within the focal area for each pass in a video (Stay), and whether the observation of the pass was censored (Cens).
#' @param station_data A data frame containing information for each camera station, with one row per station. It includes the station name (Station), location coordinates, covariates, and so forth. Unlike the detection_data, it lists all the camera stations that have been set up.
#' @param col_name_station Column name containing station ID info
#' @param col_name_species Column name containing species name detected
#' @param col_name_y Column name containing number of animal passage within a focal area
#' @param target_species The species names for which you want to estimate thier density
#' @param model Model name used ("REST" or "RAD-REST")
#' @return data frame containing at least 3 columns (station, species and npass)
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

  if(!(model == "REST" | model == "RAD-REST")){
    stop("check model name! 'model' must be 'REST' or 'RAD-REST.'" , call. = FALSE)
  }
  if(model == "REST"){
    detection_temp <- detection_data %>%
      rename(Station = !!sym(col_name_station), Species = !!sym(col_name_species), y = !!sym(col_name_y)) %>%
      filter(Species %in% target_species) %>%
      group_by(Station, Species) %>%
      summarize(Y = sum(y)) %>%
      right_join(station_data, by = c("Station")) %>%
      replace_na(list(Y = 0)) %>%
      arrange(Station)

  } else if (model == "RAD-REST") {
    detection_temp_0 <- crossing(
      Species = target_species,
      Station = station_data$Station
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

    # # 最大値を動的に取得し、欠損している列を補完
    index <- detection_temp_1 %>%
      select(starts_with("y_")) %>%            # y_ で始まる列を選択
      names() %>%                              # 列名を取得
      str_extract("\\d+") %>%                  # 数値部分を抽出
      as.integer()

    if(length(index) != max(index) + 1) {
      add_col <- paste0("y_", (0:max(index))[!(0:max(index)) %in% index])
      detection_temp_1 <- detection_temp_1 %>%
        mutate(!!!rlang::set_names(lapply(add_col, function(col) rep(as.integer(0), nrow(detection_temp_1))), add_col)) %>%
        select(Station, Species,
               gtools::mixedsort(names(select(., starts_with("y_")))),
               everything())
    }

    detection_temp <- detection_data %>%
      rename(Station = !!sym(col_name_station), Species = !!sym(col_name_species), y = !!sym(col_name_y)) %>%
      filter(Species %in% target_species) %>%
      group_by(Station, Species) %>%
      summarise(N = n(), .groups = 'drop') %>%
      left_join(detection_temp_1, by = c("Station", "Species"))  %>%
      right_join(detection_temp_0, by = c("Station", "Species")) %>%
      dplyr::mutate(across(everything(), \(x) tidyr::replace_na(x, 0))) %>%
      arrange(Species, Station) %>%
      left_join(station_data, by = c("Station")) %>%
      arrange(Species) %>%
      arrange(Station)
  }
  return(detection_temp)
}
Station <- Species <- y <- NULL
