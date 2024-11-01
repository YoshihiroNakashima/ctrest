library(tidyverse)
# 動物の密度は対数正規分布に従う
rlnorm2 <- function(n, mean=1, sd=1) {
  sdlog <- sqrt(log((sd/mean)^2 + 1))
  meanlog <- log(mean) - (sdlog^2) / 2
  rlnorm(n, meanlog = meanlog, sdlog = sdlog)
}
# 先に大量の撮影タイミングデータを作成
library("circular")
rv1 <- rvonmises(n = 10000, circular(pi /2), 8, control.circular=list(units="radian"))
rv2 <- rvonmises(n = 8000, circular(3 * pi) / 2, 8,  control.circular=list(units="radian"))
rv3 <- rvonmises(n = 6000, circular(pi), 4,  control.circular=list(units="radian"))
t.temp <- as.numeric(sample(c(rv1, rv2, rv3)))

N_station <- 100
station <- paste0("ST", formatC(1:N_station, flag = "0", width = 3))

# Term 1 ---------------------------------------------------------------
start_date_1 <- ymd_hms("2024-05-07 09:00:00") + minutes(round(runif(N_station, 0, 180))) # 1回目の調査開始
start_date_2 <- start_date_1 + minutes(24 * 60 * 30) # 2回目の調査開始 (1回目の1か月後)

# 1回目の調査終了
smp1 <- sample(0:1, N_station, replace = T) # 調査終了まで電池が切れたかどうか(1 = 切れた)
end_date_1 <- start_date_2 - minutes(round(runif(N_station, 0, 24 * 60 * 30))) * smp1 # 切れた場合は乱数を引く

# 2回目の調査終了
smp2 <- sample(0:1, N_station, replace = T) # 調査終了まで電池が切れたかどうか(1 = 切れた)
end_date_2 <- start_date_2 + minutes(24 * 60 * 30) - minutes(round(runif(N_station, 0, 24 * 60 * 30))) * smp2 # 切れた場合は乱数を引く

d.check <- tibble(Station = rep(station, 2), Start = c(start_date_1, start_date_2),
                  End = c(end_date_1, end_date_2))

# ggplot(data = d.check) +
#   geom_segment(aes(x = Start, xend = End, y = Station, yend = Station),
#                alpha = 0.3,
#                linewidth = 2.5) +
#   geom_point(aes(x = Start, y = Station), shape = 4, col = 2, alpha = 0.8)

effort <- d.check %>%
  mutate(Effort = difftime(End, Start, unit = "days")) %>%
  mutate(Term = rep(c("term1", "term2"), each = N_station))

total_effort <- effort %>%
  group_by(Station) %>%
  summarize(Effort = sum(Effort))

N_species <- 12
species <- paste0("SP", formatC(1:N_species, flag = "0", width = 2))
density <- rlnorm2(N_species, mean = 8, sd = 15)
density[1] <- 10
sp <- list(0)
for(k in 1:N_species) {
  prob <- c(0.4, 0.4, 0.1, 0.1)
  e_y <- sum(prob * (1:length(prob) - 1))
  mean_stay <- rlnorm2(N_species, mean = 2, sd = 0.25)


  library(tidyverse)
  N_station <- 100
  station <- paste0("ST", formatC(1:N_station, flag = "0", width = 3))
  focal <- 2.0
  censored <- 20

  # Assuming that cameras were checked once before survey terminated --------
  # Species A parameters
  mean <- 4
  sd <- 3

  # ラジアンから時分秒に変換
  time_in_seconds <- t.temp * (24 * 3600 / (2 * pi))
  time_as_hms <- seconds_to_period(time_in_seconds)

  library(activity)
  model<-fitact(t.temp, bw = 1.5*bwcalc(t.temp, K = 3),reps=1)

  activity <- model@act[1]

  # 撮影枚数の期待値
  p <- c(0.5, 0.4, 0.1)
  exp.n <- sum(p * 0:2)

  mu <- density[[k]] * focal * 24 * 60 * 60 * effort %>% pull(Effort) / mean / 1000000 * activity / exp.n
  library(MASS)
  N_detection  <- rnegbin(n = N_station * 2, mu, theta = 1.0) # 動画枚数

  # E_y
  prob <- MCMCpack::rdirichlet(n = N_station, alpha = p * 20)
  prob <- rbind(prob, prob)
  y_temp <- matrix(0, N_station * 2, length(p))
  y_list <- list(0)
  for(i in 1:(N_station * 2)){
    y_temp[i, ] <- t(stats::rmultinom(n = 1, size = N_detection[i], prob = prob[i, ]))
  }
  y <- rep(rep(0:2, N_station * 2), as.vector(t(y_temp)))

  # 検出時間の作成
  temp.time <- list(0)
  for(i in 1:nrow(effort)){
    temp.time[[i]] <- sort(as.POSIXct(runif(N_detection[i], effort$Start[i], effort$End[i]), origin = "1970-01-01"))
  }
  time <- as.Date((as.POSIXct(unlist(temp.time)))) + time_as_hms[1:sum(N_detection)]

  time_new <- 2 * pi * (hour(time) * 60 * 60 + minute(time) * 60 + second(time))/(24 * 60 * 60)

  model<-fitact(time_new, bw = 1.5*bwcalc(t.temp, K = 3),reps=1)

  # 滞在時間の生成
  stay <- round(rlnorm2(sum(N_detection), mean, sd), 1)
  hist(stay)
  stay[stay == 0] <- 0.1
  cens <- rep(0, sum(N_detection))
  cens[stay > 20] <- 1
  stay[stay > 20] <- 20

  # 動物の観測データ
  temp1 <- tibble(
    Station = rep(effort$Station, N_detection),
    DateTime = time,
    Term = rep(effort$Term, N_detection),
    Species = species[k],
    y = y
  ) %>%
    mutate(Stay = stay, Cens = cens) %>%
    group_by(Station, Term) %>%                # StationとTermごとにグループ化
    mutate(is_max = if_else(DateTime == max(DateTime), 1, 0)) %>% # 最大DateTimeに1、それ以外に0
    ungroup()  # グループ化解除

  sp[[k]] <- temp1
}

temp2 <- do.call("rbind", sp)


# 調査員のデータ(以下が、調査員が撮影された動画)
start_data <- effort %>%
  rename(DateTime = Start) %>%
  mutate(Species = "Surveyor") %>%
  dplyr::select(-End, -Effort) %>%
  mutate(y = NA, Stay = NA, Cens= NA, is_max = NA) %>%
  rbind(temp2) %>%
  arrange(Station, DateTime)

# カメラが途中で死ななかったものに関しては、普通に回収の日時を付ける
mal <- tibble(Station = rep(station, 2), Term = rep(c("term1", "term2"), each = N_station), Malfunction = c(smp1, smp2))

end_data <- effort %>%
  rename(DateTime = Start) %>%
  mutate(Species = "Surveyor") %>%
  mutate(Malfunction = c(smp1, smp2)) %>%
  filter(Malfunction == 0) %>%
  dplyr::select(-End, -Effort) %>%
  mutate(y = NA, Stay = NA, Cens= NA, is_max = NA) %>%
  dplyr::select(- Malfunction) %>%
  rbind(start_data) %>%
  arrange(Station, DateTime) %>%
  left_join(mal, by = c("Station", "Term")) %>%
  left_join(effort, by = c("Station", "Term")) %>%
  mutate(DateTime = if_else(!is.na(is_max) & Malfunction == 1 & is_max == 1, End, DateTime)) # 最後の撮影をカメラが死んだ日時に変換


detection_data <- end_data %>%
  dplyr::select(-Malfunction, -Start, -End, -Effort, -is_max)


station_data <- tibble(Station = station) %>%
  mutate(x1 = rnorm(n(), 0, 1)) %>%
  mutate(x2 = sample(c("A", "B", "C"), n(), replace = TRUE))


#usethis::use_mit_license()
effort_temp <- detection_data %>%
  group_by(Station, Term) %>%
  summarize(start = min(DateTime, na.rm = TRUE), end = max(DateTime, na.rm = TRUE)) %>%
  mutate(effort = as.numeric(difftime(end, start, units = "days"))) %>%
  ungroup() %>%
  group_by(Station) %>%
  summarize(effort = sum(effort)) %>%
  print(n = Inf)

library(devtools)
usethis::use_data(detection_data, overwrite = TRUE)
usethis::use_data(station_data, overwrite = TRUE)


devtools::document()
load_all()

