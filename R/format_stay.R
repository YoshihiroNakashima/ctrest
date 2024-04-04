format_stay <-
  function(detection_data,
           col_name_station = "station",
           col_name_species = "species",
           col_name_stay = "stay",
           col_name_cens = "cens",
           target_species = "A") {
    stay_temp <- detection_data %>%
      filter(!is.na(!!sym(col_name_stay)),
             !is.na(!!sym(col_name_cens)),
             !!sym(col_name_species) == target_species) %>%
      select(
        !!sym(col_name_station),
        !!sym(col_name_species),
        !!sym(col_name_stay),
        !!sym(col_name_cens)
      )
    print(stay_temp)
  }
