library(tidyverse)
set.seed(123)
# Number of camera-traps
nstation <- 100
# Station Name
station <- paste0("ST", formatC(1:nstation, flag = "0", width = 3))
# a covariate (e.g. habitat quality) for each station
covariate <- rnorm(nstation, 0, 1)

# Survey parameters
survey <- 50
effort <- tibble(station = station, effort = rep(24 * 60 * 60 * survey, nstation))
focal <- 3.0
censored <- 50

rlnorm2 <- function(n, mean, sd) {
  sdlog <- sqrt(log((sd/mean)^2 + 1))
  meanlog <- log(mean) - (sdlog^2) / 2
  rlnorm(n, meanlog = meanlog, sdlog = sdlog)
}

# Species A ---------------------------------------------------------------

# Species A parameters
rate <- 0.2
density <- 10
activity <- 0.5

# Expected number of videos per camera-trap
mu <- density * focal * 24 * 60 * 60 * survey * rate / 1000000 * activity * exp(0.5 * covariate) / mean(exp(0.5 * covariate))
ndetection <- rpois(nstation, mu)

# Random pass within focal area
nfocal <- rbinom(sum(ndetection), size = 1, p = 0.5)

# Generate date and time
start_date <- as.Date("2024-03-02")
end_date <- as.Date("2024-04-29")
date_sequence <- seq(start_date, end_date, by = "day")
date <- sample(date_sequence, sum(ndetection), replace = TRUE)

radian_time <- rnorm(sum(ndetection), mean = pi, sd = 0.5)
degree_time <- radian_time * (180 / pi)
hour <- as.integer(degree_time / 15)
minute <- (degree_time %% 15) * 4
time <- sprintf("%02d:%02.0f", hour, minute)
datetime <- paste(date, time)

# Create detection data
detection_data_a <- tibble(station = rep(station, ndetection), datetime, species = "A", nfocal) %>%
  mutate(stay = rlnorm2(n(), mean = 1 / rate, sd = 5)) %>%
  mutate(stay = if_else(nfocal >= 1, round(stay, 3), NA_real_)) %>%
  group_by(station) %>%
  mutate(rank = ifelse(!is.na(stay), 1:n(), NA_integer_)) %>%
  mutate(rank = ifelse(!is.na(stay), rank(rank), NA_integer_)) %>%
  mutate(stay = ifelse(rank <= 5, stay, NA_real_)) %>%
  mutate(rank = ifelse(rank <= 5, rank, NA_real_)) %>%
  select(-rank) %>%
  ungroup() %>%
  mutate(cens = ifelse(stay > censored, 1, 0), stay = ifelse(stay > censored, censored, stay)) %>%
  arrange(station, datetime)

print(detection_data_a, n = 100)

detection_data_a %>% filter(!is.na(stay) & cens == 0) %>% pull(stay) %>% mean()

# Species B ---------------------------------------------------------------

# Species B parameters
rate <- 0.1
density <- 5
activity <- 0.5

# Expected number of videos per camera-trap
mu <- density * focal * 24 * 60 * 60 * survey * rate / 10 ^ 6 * activity * exp(0.5 * covariate) / mean(exp(0.5 * covariate))
ndetection <- rpois(nstation, mu)

# Random pass within focal area
nfocal <- rbinom(sum(ndetection), size = 1, p = 0.3)
nfocal <- rbinom(sum(ndetection), size = 2, p = 0.1)

# Generate date and time
date <- sample(date_sequence, sum(ndetection), replace = TRUE)

radian_time <- rnorm(sum(ndetection), mean = 0, sd = 0.2)
degree_time <- radian_time * (180 / pi)
hour <- as.integer(degree_time / 15)
minute <- (degree_time %% 15) * 4
time <- sprintf("%02d:%02.0f", hour, minute)
datetime <- paste(date, time)

# Create detection data
detection_data_b <- tibble(station = rep(station, ndetection), datetime, species = "B", nfocal) %>%
  mutate(stay = rexp(n(), rate)) %>%
  mutate(stay = if_else(nfocal >= 1, round(stay, 3), NA_real_)) %>%
  group_by(station) %>%
  mutate(rank = ifelse(!is.na(stay), 1:n(), NA_integer_)) %>%
  mutate(rank = ifelse(!is.na(stay), rank(rank), NA_integer_)) %>%
  mutate(stay = ifelse(rank <= 5, stay, NA_real_)) %>%
  mutate(rank = ifelse(rank <= 5, rank, NA_real_)) %>%
  select(-rank) %>%
  ungroup() %>%
  mutate(cens = ifelse(stay > censored, 1, 0), stay = ifelse(stay > censored, censored, stay)) %>%
  arrange(station, datetime)

print(detection_data_a, n = 100)

detection_data_a %>% filter(!is.na(stay) & cens == 0) %>% pull(stay) %>% mean()

print(detection_data_b, n = 100)


# Integrate dataset -------------------------------------------------------

detection_data_temp <- rbind(detection_data_a, detection_data_b) %>% mutate(datetime = ymd_hm(datetime))
station_data <- tibble(station, covariate)

# add surveyor data -------------------------------------------------------
surveyor_start <- tibble(
  station = station,
  datetime = ymd_hms("2024-03-01 09:00:00") + minutes(round(runif(nstation, 0, 180))),
  species = "surveyor",
  nfocal = NA,
  stay = NA,
  cens = NA
)
surveyor_end <- tibble(
  station = station,
  datetime = ymd_hms("2024-04-30 09:00:00") + minutes(round(runif(nstation, 0, 180))),
  species = "surveyor",
  nfocal = NA,
  stay = NA,
  cens = NA
)
malfunc <- c("ST003", sample(station, 5, replace = FALSE))
surveyor_end[surveyor_end$station %in% malfunc, "datetime"] <- NA
surveyor_end <- surveyor_end %>% filter(!is.na(datetime))

detection_data <- rbind(detection_data_temp, surveyor_start, surveyor_end) %>%
  arrange(station, datetime)

#usethis::use_mit_license()

usethis::use_data(detection_data, overwrite = TRUE)
usethis::use_data(station_data, overwrite = TRUE)
