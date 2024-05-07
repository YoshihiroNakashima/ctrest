rlnorm2 <- function(n, mean, sd) {
  sdlog <- sqrt(log((sd/mean)^2 + 1))
  meanlog <- log(mean) - (sdlog^2) / 2
  rlnorm(n, meanlog = meanlog, sdlog = sdlog)
}

library(tidyverse)
set.seed(1)
# Number of camera-traps
nstation <- 100
# Station Name
station <- paste0("ST", formatC(1:nstation, flag = "0", width = 3))

# setting basic variables

focal <- 2.0
censored <- 20


# Assuming that cameras were checked once before survey terminated --------
# Species A parameters
mean <- 4
sd <- 5
density <- 10

# activity
library("Rfast")
rv1 <- rvonmises(n = 10000, pi /2, 8, rads = TRUE)
rv2 <- rvonmises(n = 8000, 3 * pi / 2, 8, rads = TRUE)
rv3 <- rvonmises(n = 6000, pi, 4, rads = TRUE)

t.temp <- sample(c(rv1, rv2, rv3))

# t.temp <- rnorm(1000, pi, 0.6)
# t.temp <- t.temp[t.temp < 2 * pi]
# ラジアンから時分秒に変換
time_in_seconds <- t.temp * (24 * 3600 / (2 * pi))
time_as_hms <- seconds_to_period(time_in_seconds)

library(activity)

model<-fitact(t.temp, bw = 1.5*bwcalc(t.temp, K = 3),reps=1)
plot(model)

activity <- model@act[1]

# Term 1 ---------------------------------------------------------------

# Generate date and time
start_date_1 <- ymd_hms("2024-05-07 09:00:00") + minutes(round(runif(nstation, 0, 180)))
end_date_1_0 <- ymd_hms("2024-06-10 10:00:00")  + minutes(round(runif(nstation, 0, 180)))
smp1 <- sample(0:1, nstation, replace = T)
temp <- minutes(round(runif(nstation, 0, 24 * 60 * 30))) * smp1
end_date_1 <- end_date_1_0 - temp

start_date_2 <-  end_date_1_0 + minutes(round(runif(nstation, 0, 10)))
end_date_2_0 <- ymd_hms("2024-07-12 10:00:00")  + minutes(round(runif(nstation, 0, 180)))
smp2 <- sample(0:1, nstation, replace = T)
temp <- minutes(round(runif(nstation, 0, 24 * 60 * 30)))* smp2
end_date_2 <- end_date_2_0 - temp

d.check <- tibble(Station = rep(station, 2), Start = c(start_date_1, start_date_2),
       End = c(end_date_1, end_date_2))

ggplot(data = d.check) +
  geom_segment(aes(x = Start, xend = End, y = Station, yend = Station),
               alpha = 0.3,
               linewidth = 2.5) +
  geom_point(aes(x = Start, y = Station), shape = 4, col = 2, alpha = 0.8)

effort <- d.check %>%
  mutate(Effort = difftime(End, Start, unit = "days")) %>%
  mutate(Term = rep(c("term1", "term2"), each = nstation))

total_effort <- effort %>%
  group_by(Station) %>%
  summarize(Effort = sum(Effort))


# 撮影枚数の期待値

p <- c(0.5, 0.4, 0.1)
exp.n <- sum(p * 0:2)

mu <- density * focal * 24 * 60 * 60 * effort %>% pull(Effort) / mean / 1000000 * activity / exp.n
library(MASS)
ndetection  <- rnegbin(n = nstation * 2, mu, theta = 1.0) # 動画枚数

# E_y
prob <- MCMCpack::rdirichlet(n = nstation, alpha = p * 20)
prob <- rbind(prob, prob)
y_temp <- matrix(0, nstation * 2, length(p))
y_list <- list(0)
for(i in 1:(nstation * 2)){
  y_temp[i, ] <- t(stats::rmultinom(n = 1, size = ndetection[i], prob = prob[i, ]))
}
y <- rep(rep(0:2, nstation * 2), as.vector(t(y_temp)))


temp.time <- list(0)
for(i in 1:nrow(effort)){
  temp.time[[i]] <- sort(as.POSIXct(runif(ndetection[i], effort$Start[i], effort$End[i]), origin = "1970-01-01"))
}

time <- as.Date((as.POSIXct(unlist(temp.time)))) + time_as_hms[1:sum(ndetection)]

# 滞在時間の生成
stay <- round(rlnorm2(sum(ndetection), mean, sd), 1)
hist(stay)
stay[stay == 0] <- 0.1
cens <- rep(0, sum(ndetection))
cens[stay > 20] <- 1
stay[stay > 20] <- 20

# 動物の観測データ
temp1 <- tibble(Station = rep(effort$Station, ndetection),
       DateTime = time,
       Term = rep(effort$Term, ndetection),
       Species = "A",
       y = y) %>%
       mutate(Stay = stay, Cens = cens)

# 調査員のデータ
temp2 <- effort %>%
  mutate(End = c(end_date_1_0, end_date_2_0), smp = c(smp1, smp2)) %>%
  pivot_longer(cols = Start:End, names_to = "Start_End", values_to = "DateTime") %>%
  mutate(Species = "Surveyor") %>%
  filter(Start_End == "Start" | (Start_End == "Start" & smp == 1)) %>%
  dplyr::select(Station, Term, DateTime, Species) %>%
  mutate(y = NA, Stay = NA, Cens = NA)


detection_data <- rbind(temp1, temp2) %>%
  arrange(Station, DateTime)
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
  summarize(effort = sum(effort))


usethis::use_data(detection_data, overwrite = TRUE)
usethis::use_data(station_data, overwrite = TRUE)


devtools::document()
load_all()

