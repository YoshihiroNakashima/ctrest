library(tidyverse)
# 動物の密度は対数正規分布に従う
rlnorm2 <- function(n, mean=1, sd=1) {
  sdlog <- sqrt(log((sd/mean)^2 + 1))
  meanlog <- log(mean) - (sdlog^2) / 2
  rlnorm(n, meanlog = meanlog, sdlog = sdlog)
}
# ラジアンから ymd_hms に変換する関数（ベクトル対応）
radian_to_time <- function(rad_vec, date = "1970-01-01") {
  # ラジアンを0〜1の割合に変換（1日は2π）
  frac_day <- (rad_vec %% (2 * pi)) / (2 * pi)

  # 秒に変換（1日 = 86400秒）
  seconds <- frac_day * 86400

  # 基準日をPOSIXct形式で作成
  base_time <- ymd_hms(paste(date, "00:00:00"), tz = "UTC")

  # 各ラジアンに対応する時刻を生成
  result_times <- base_time + seconds

  return(result_times)
}

# 先に大量の撮影タイミングデータを作成
library("circular")
rv1 <- rvonmises(n = 10000, circular(pi /2), 8, control.circular=list(units="radian"))
rv2 <- rvonmises(n = 8000, circular(3 * pi) / 2, 8,  control.circular=list(units="radian"))
rv3 <- rvonmises(n = 6000, circular(pi), 4,  control.circular=list(units="radian"))
t.temp <- as.numeric(sample(c(rv1, rv2, rv3)))

library(activity)
model <- fitact(t.temp, bw = 1.5 * bwcalc(t.temp, K = 3),reps = 1)
plot(model)
activity <- model@act[1]
activity

N_station <- 100
station <- paste0("ST", formatC(1:N_station, flag = "0", width = 3))

start_date_1 <- ymd_hms("2024-05-07 09:00:00") + minutes(round(runif(N_station, 0, 180))) # 1回目の調査開始
end_date_1 <- start_date_1 +  minutes(24 * 60 * 30) - sample(0:1, N_station, prob = c(3/4, 1/4), rep = T) * minutes(round(runif(N_station, 24 * 60 * 5, 24 * 60 * 20)))

start_date_2 <- start_date_1 + minutes(24 * 60 * 30) # 2回目の調査開始 (1回目の1か月後)
end_date_2 <- start_date_2 + minutes(24 * 60 * 30) - sample(0:1, N_station, prob = c(3/4, 1/4), rep = T) * minutes(round(runif(N_station, 24 * 60 * 5, 24 * 60 * 20)))

d.surveyor <- tibble(Station = rep(station, 4),
                     DateTime = c(start_date_1, end_date_1, start_date_2, end_date_2),
                     Species = "Surveyor",
                     y = NA, Stay = NA, Cens= NA,
                     Term = rep(c("term1", "term2"), each = 2 * N_station))


d.check <- tibble(Station = rep(station, 2), Start = c(start_date_1, start_date_2),
                  End = c(end_date_1, end_date_2))

effort <- d.check %>%
  mutate(Effort = difftime(End, Start, unit = "days")) %>%
  mutate(Term = rep(c("term1", "term2"), each = N_station))

total_effort <- effort %>%
  group_by(Station) %>%
  summarize(Effort = sum(Effort))

N_species <- 12
species <- paste0("SP", formatC(1:N_species, flag = "0", width = 2))
density <- rlnorm2(N_species, mean = 10, sd = 10)
density[1] <- 10

probability <- c(0.4, 0.4, 0.1, 0.1)
e_y <- sum(probability * (1:length(probability) - 1))

library(tidyverse)
focal <- 1.96
censored <- 20

# Assuming that cameras were checked once before survey terminated --------
# Species A parameters
mean <- 4
sd <- 3


sp <- list(0)
for(k in 1:N_species) {

  # 撮影枚数の期待値
  mu <- density[k] * focal * 24 * 60 * 60 * effort %>% pull(Effort) / mean / 1000000 * activity / e_y
  library(MASS)
  N_detection  <- rnegbin(n = N_station * 2, mu, theta = 1.0) # 動画枚数

  # E_y
  prob <- MCMCpack::rdirichlet(n = N_station, alpha = probability * 20)
  prob <- rbind(prob, prob)
  y_temp <- matrix(0, nrow = N_station * 2, ncol = ncol(prob))
  y_list <- list(0)
  for(i in 1:(N_station * 2)){
    y_temp[i, ] <- t(stats::rmultinom(n = 1, size = N_detection[i], prob = prob[i, ]))
  }
  y <- rep(rep(0:(ncol(prob) - 1), N_station * 2), as.vector(t(y_temp)))

  # 検出時刻の作成
  time <- hms::as_hms(radian_to_time(t.temp))[1:sum(N_detection)]

  start_time <- effort %>% pull(Start)
  end_time <- effort %>% pull(End)

  # 乱数で日時を生成（数値として秒単位で間の値をランダムに取得）
  random_times_list <- list(0)
  for(i in 1:nrow(effort)) {
    random_times_list[[i]] <- as_datetime(runif(N_detection[i], as.numeric(start_time[i]) + 24 * 60 * 60, as.numeric(end_time[i]) - 24 * 60 * 60))
  }

  # 時刻情報を切り落とす（日付だけを取得）
  random_dates <- as_date(as_datetime(unlist(random_times_list)))

  detection_time <- as_datetime(random_dates) + time

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
    DateTime = detection_time,
    Term = rep(effort$Term, N_detection),
    Species = species[k],
    y = y
  ) %>%
    mutate(Stay = stay, Cens = cens) %>%
    group_by(Station, Term) %>%                # StationとTermごとにグループ化
    ungroup()  # グループ化解除

  sp[[k]] <- temp1
}

detection_data <- do.call("rbind", sp) %>%
  bind_rows(d.surveyor) %>%
  arrange(Station, DateTime)



station_data <- tibble(Station = station) %>%
  mutate(x1 = rnorm(n(), 0, 1)) %>%
  mutate(x2 = sample(c("A", "B", "C"), n(), replace = TRUE))


library(devtools)
usethis::use_data(detection_data, overwrite = TRUE)
usethis::use_data(station_data, overwrite = TRUE)


devtools::document()
load_all()

