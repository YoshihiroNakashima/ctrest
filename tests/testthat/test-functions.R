test_that("fromat_stay functions work", {
  expect_type(format_stay(
    detection_data = detection_data, # 元データ
    station_data = station_data,
    col_name_station = "Station",    # カメラ設置点IDを含む列名
    col_name_species = "Species",    # 種名を含む列名
    col_name_stay = "Stay",          # 滞在時間を含む列名
    col_name_cens = "Cens",          # 打ち切りの有無を含む列名
    target_species = "SP01"             # 対象種名
  ), "list")
})
