select_videos <- function(detection_data,
                          col_name_species = "Species",
                          col_name_station = "Station",
                          col_name_datetime = "DateTimeCorrected",
                          N_sampled,
                          Indep_criteria,
                          seed = NULL,
                          target_species = NULL) {

  # seed値の設定
  if (!is.null(seed)) set.seed(seed)

  # 列名をシンボルに変換
  species_sym <- sym(col_name_species)
  station_sym <- sym(col_name_station)
  datetime_sym <- sym(col_name_datetime)

  # target_speciesが指定されている場合はフィルタリング
  if (!is.null(target_species)) {
    detection_data <- detection_data %>%
      filter(!!species_sym %in% target_species)
  }

  # 種ごとの処理
  detection_data %>%
    group_by(!!species_sym) %>%
    group_split() %>%
    lapply(function(species_data) {

      # ステーションごとの処理
      species_data %>%
        arrange(!!datetime_sym) %>%
        group_by(!!station_sym) %>%
        group_split() %>%
        lapply(function(station_data) {

          # 独立性基準を満たす動画を選択
          selected_indices <- c()
          last_datetime <- ymd_hms("1900-01-01 00:00:00") # 非常に過去の日時で初期化

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

      if (nrow(independent_data) == 0) return(NULL) # 独立データがない場合

      unique_stations <- unique(independent_data[[col_name_station]])
      num_unique_stations <- length(unique_stations)

      if (num_unique_stations >= N_sampled) {
        # ステーション数がN_sampled以上の場合
        sampled_stations <- sample(unique_stations, N_sampled)
        sampled_data <- independent_data %>%
          filter(!!station_sym %in% sampled_stations) %>%
          group_by(!!station_sym) %>%
          slice_sample(n = 1) %>%
          ungroup()
        return(sampled_data)

      } else {
        # ステーション数がN_sampled未満の場合
        sampled_data <- tibble()
        available_data <- independent_data # 選択可能なデータを保持
        while(nrow(sampled_data) < N_sampled && nrow(available_data) > 0){ # available_dataが空になったら終了条件を追加
          temp_sampled <- available_data %>%
            group_by(!!station_sym) %>%
            slice_sample(n=1) %>%
            ungroup()
          sampled_data <- bind_rows(sampled_data, temp_sampled)
          # 選択されたデータをavailable_dataから除外
          available_data <- anti_join(available_data, temp_sampled, by = c(col_name_station, col_name_datetime))
        }
        sampled_data <- sampled_data[1:N_sampled,] #N_sampledを超えた分を削除
        return(sampled_data)
      }
    }) %>%
    bind_rows()
}
