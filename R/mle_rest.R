#' Maximum likelihood parameter estimation for the REST/RAD-REST model based on the implementation of Nelder and Mead (1965).
#'
#' @param station_effort_data A data frame containing information for each camera station. Typically, this is the output of the `add_effort` function. Alternatively, a manually prepared data frame may be provided, with the required columns as follows:
#'   - If `model = "REST"`, the data frame must contain:
#'     - `Station` (character): Camera station ID.
#'     - `Effort` (numeric): Camera trapping effort (days) for each station.
#'     - `Species` (character): Species name.
#'     - `Y` (numeric): Total number of passes through the focal area for each station-species pair.
#'   - If `model = "RAD-REST"`, the data frame must contain:
#'     - `Station` (character): Camera station ID.
#'     - `Effort` (numeric): Camera trapping effort (days) for each station.
#'     - `Species` (character): Species name.
#'     - `N` (numeric): Total number of detected videos for each station-species pair.
#'     - `y_X` columns (`y_0`, `y_1`, ..., `y_max`): Number of videos classified by the number of passes observed.
#' @param stay_data A data frame returned by the `format_stay` function, containing:
#'   - `Station` (character): Camera station ID.
#'   - `Species` (character): Species name.
#'   - `Stay` (numeric): Staying time (in seconds) within the focal area for each detected pass.
#'   - `Censored` (binary): Whether the staying time was censored (1 = censored, 0 = observed).
#' @param stay_data A data frame returned by the `format_stay` function, containing:
#'   - `Station` (character): Camera station ID.
#'   - `Species` (character): Species name.
#'   - `Stay` (numeric): Staying time (in seconds) within the focal area for each detected pass.
#'   - `Censored` (binary): Whether the staying time was censored (1 = censored, 0 = observed).
#' @param activity_data A numeric vector of detection times (in radians), returned by the `format_activity` function.
#' @param focal_area The area of the focal area, in square meters.
#' @param model Specify the model to use. Choose either "REST" or "RAD-REST".
#' @param target_species A character string specifying the species to be analyzed. Only a single species can be specified.
#' @return A data frame containing the estimated density, its standard deviation (SD), the 95% lower and upper confidence limits, AIC, and delta AIC value for each model.
#' @export
#' @import activity
#' @importFrom stats dexp dgamma dlnorm dnbinom dpois dweibull optim pexp pgamma plnorm pweibull na.omit rexp rnorm step time
#' @examples
#' station_data_rest <- format_station_data(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_y = "y",
#'   model = "REST",
#'   target_species = "SP01"
#' )
#' station_effort_rest <- add_effort(
#'   detection_data = detection_data,
#'   station_data_formatted = station_data_rest,
#'   col_name_station = "Station",
#'   col_name_term = "Term",
#'   col_name_datetime = "DateTime",
#'   plot = TRUE,
#'   font_size = 5
#' )
#'
#' stay_data <- format_stay(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_stay = "Stay",
#'   col_name_cens = "Cens",
#'   target_species = "SP01"
#' )
#'
#' activity_data <- format_activity(
#'   detection_data = detection_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_datetime = "DateTime",
#'   target_species = "SP01",
#'   indep_time = 30
#' )
#'
#' model <- mle_rest(station_effort_data = station_effort_rest,
#'                   stay_data = stay_data,
#'                   activity_data = activity_data,
#'                   focal_area = 1.56,
#'                   model = "REST",
#'                   target_species = "SP01")
#'
mle_rest <- function(station_effort_data,
                     stay_data,
                     activity_data,
                     focal_area,
                     model,
                     target_species) {

  vector_time <- activity_data %>% pull(time)
  model_act <- activity::fitact(
    dat = vector_time,
    bw = activity::bwcalc(vector_time, K = 3),
    adj = 1.5,
    sample  = "none"
  )
  activity_proportion <- model_act@act
  stay <- stay_data %>% pull(Stay)
  cens <- stay_data %>% pull(Cens)
  effort <- station_effort_data %>%  pull(Effort)

  # RESTモデル -----------------------------------------------------------------
  if(model == "REST") {
    Y <- station_effort_data %>%  pull(Y)

    res <- vector("list", length = 2)
    names(res) <- c("Density", "Stay")
    res$Actvity <- model_act@act

    res$Density <- data.frame(ModelName = NA, Density = NA, SD = NA, Lower = NA, Upper = NA, AIC = NA, delta_AIC = NA)
    res$Stay <- vector("list", length = 4)
    names(res$Stay) <- c("Exponential", "Gammma", "Weibull", "Lognormal")
    res$Density[1:8, 1] <- c(paste("Poisson", c("Exponential", "Gammma", "Weibull", "Lognormal"), sep = "-"),
                             paste("NegBinom", c("Exponential", "Gammma", "Weibull", "Lognormal"), sep = "-"))

    # options(warn = -1)
    # Poisson-exponential model
    minus.log.lik.exp <- function(params) {
      density <- params[1]
      rate <- params[2]
      - (sum(dexp(stay[cens == 0], rate = rate, log = TRUE)) +
           sum(1 - pexp(stay[cens == 1], rate = rate, log.p = TRUE)) +
           sum(dpois(Y, lambda = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion * rate, log = TRUE)
           ))
    }
    fit_exp <- optim(par = c(density = 1, rate = 0.1), fn = minus.log.lik.exp, hessian = TRUE)

    fit_exp$AIC <- 2 * fit_exp$value + length(fit_exp$par) * 2
    sd_exp <- sqrt(diag(solve(fit_exp$hessian)))
    conf_exp <- 1.96 * sqrt(diag(solve(fit_exp$hessian)))

    res$Density[1, 2:6] <- c(fit_exp$par[1], sd_exp[1], fit_exp$par[1] - conf_exp[1], fit_exp$par[1] + conf_exp[1], fit_exp$AIC)
    res$Stay$Exponential <- fit_exp$par[2]

    # Poisson-Gamma model
    minus.log.lik.gam <- function(params) {
      density <- params[1]
      rate <- params[2]
      shape <- params[3]
      - (sum(dgamma(stay[cens == 0], shape = shape, rate = rate, log = TRUE)) +
           sum(1 - pgamma(stay[cens == 1], shape = shape, rate = rate, log.p = TRUE)) +
           sum(dpois(Y, lambda = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion * rate /shape, log = TRUE)
           ))
    }

    fit_gam <- optim(par = c(density = 5, shape = 0.5, rate = 1.0), fn = minus.log.lik.gam, hessian = TRUE)

    fit_gam$AIC <- 2 * fit_gam$value + length(fit_gam$par) * 2
    sd_gam <- sqrt(diag(solve(fit_gam$hessian)))
    conf_gam <- 1.96 * sqrt(diag(solve(fit_gam$hessian)))

    res$Density[2, 2:6] <- c(fit_gam$par[1], sd_gam[1], fit_gam$par[1] - conf_gam[1], fit_gam$par[1] + conf_gam[1], fit_gam$AIC)
    res$Stay$Gammma <- fit_gam$par[2:3]

    minus.log.lik.wei <- function(params) {
      density <- params[1]
      shape <- params[2]
      scale <- params[3]
      - (sum(dweibull(stay[cens == 0], shape = shape, scale = scale, log = TRUE)) +
           sum(1 - pweibull(stay[cens == 1], shape = shape, scale = scale, log.p = TRUE)) +
           sum(dpois(Y, lambda = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion / (scale * gamma(1 + 1 / shape)), log = TRUE)
           ))
    }

    fit_wei <- optim(par = c(density = 5, shape = 0.5, scale = 10), fn = minus.log.lik.wei, hessian = TRUE)

    fit_wei$AIC <- 2 * fit_wei$value + length(fit_wei$par) * 2
    sd_wei <- sqrt(diag(solve(fit_gam$hessian)))
    conf_wei <- 1.96 * sqrt(diag(solve(fit_wei$hessian)))

    res$Density[3, 2:6] <- c(fit_wei$par[1], sd_wei[1], fit_wei$par[1] - conf_wei[1], fit_wei$par[1] + conf_wei[1], fit_wei$AIC)
    res$Stay$Weibull <- fit_wei$par[2:3]

    # Poisson-lognormal model
    minus.log.lik.log <- function(params) {
      density <- params[1]
      meanlog <- params[2]
      sdlog <- params[3]
      - (sum(dlnorm(stay[cens == 0], meanlog = meanlog, sdlog = sdlog, log = TRUE)) +
           sum(1 - plnorm(stay[cens == 1], meanlog = meanlog, sdlog = sdlog, log.p = TRUE)) +
           sum(dpois(Y, lambda = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion / exp(meanlog + sdlog ^ 2 / 2), log = TRUE)
           ))
    }

    fit_log <- optim(par = c(density = 5, meanlog = 2, sdlog = 1), fn = minus.log.lik.log, hessian = TRUE)

    fit_log$AIC <- 2 * fit_log$value + length(fit_log$par) * 2
    sd_log <- sqrt(diag(solve(fit_log$hessian)))
    conf_log <- 1.96 * sqrt(diag(solve(fit_log$hessian)))

    res$Density[4, 2:6] <- c(fit_log$par[1], sd_log[1], fit_log$par[1] - conf_log[1], fit_log$par[1] + conf_log[1], fit_log$AIC)
    res$Stay$Lognormal <- fit_log$par[2:3]

    res$Density <- res$Density[order(res$Density$AIC), ]
    res$Density$delta_AIC <- res$Density$AIC - min(res$Density$AIC)


    # Negbin-exponential model
    minus.log.lik.exp <- function(params) {
      density <- params[1]
      rate <- params[2]
      size <- params[3]
      - (sum(dexp(stay[cens == 0], rate = rate, log = TRUE)) +
           sum(1 - pexp(stay[cens == 1], rate = rate, log.p = TRUE)) +
           sum(dnbinom(Y, size = size, mu = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion * rate, log = TRUE)
           ))
    }
    fit_exp <- optim(par = c(5, 1, 10), fn = minus.log.lik.exp, hessian = TRUE)

    fit_exp$AIC <- 2 * fit_exp$value + length(fit_exp$par) * 2
    sd_exp <- sqrt(diag(solve(fit_exp$hessian)))
    conf_exp <- 1.96 * sqrt(diag(solve(fit_exp$hessian)))

    res$Density[5, 2:6] <- c(fit_exp$par[1], sd_exp[1], fit_exp$par[1] - conf_exp[1], fit_exp$par[1] + conf_exp[1], fit_exp$AIC)
    res$Stay$Exponential <- fit_exp$par[2]


    # Negbin-Gamma model
    minus.log.lik.gam.neg <- function(params) {
      density <- params[1]
      rate <- params[2]
      shape <- params[3]
      size <- params[4]
      - (sum(dgamma(stay[cens == 0], shape = shape, rate = rate, log = TRUE)) +
           sum(1 - pgamma(stay[cens == 1], shape = shape, rate = rate, log.p = TRUE)) +
           sum(dnbinom(Y, size = size, mu = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion * rate / shape, log = TRUE)
           ))
    }

    fit_gam <- optim(par = c(5, 1, 1, 10), fn = minus.log.lik.gam.neg, hessian = TRUE)

    fit_gam$AIC <- 2 * fit_gam$value + length(fit_gam$par) * 2
    conf_gam <- 1.96 * sqrt(diag(solve(fit_gam$hessian)))
    sd_gam <- sqrt(diag(solve(fit_gam$hessian)))

    res$Density[6, 2:6] <- c(fit_gam$par[1], sd_gam[1], fit_gam$par[1] - conf_gam[1], fit_gam$par[1] + conf_gam[1], fit_gam$AIC)

    # Negbin-Weibull model
    minus.log.lik.wei.neg <- function(params) {
      density <- params[1]
      shape <- params[2]
      scale <- params[3]
      size <- params[4]
      - (sum(dweibull(stay[cens == 0], shape = shape, scale = scale, log = TRUE)) +
           sum(1 - pweibull(stay[cens == 1], shape = shape, scale = scale, log.p = TRUE)) +
           sum(dnbinom(Y, size = size, mu = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion / (scale * gamma(1 + 1 / shape)), log = TRUE)
           ))
    }

    fit_wei <- optim(par = c(5, 0.5, 10, 5), fn = minus.log.lik.wei.neg, hessian = TRUE)

    fit_wei$AIC <- 2 * fit_wei$value + length(fit_wei$par) * 2
    conf_wei <- 1.96 * sqrt(diag(solve(fit_wei$hessian)))
    sd_wei <- sqrt(diag(solve(fit_wei$hessian)))
    res$Density[7, 2:6] <- c(fit_wei$par[1], sd_wei[1], fit_wei$par[1] - conf_wei[1], fit_wei$par[1] + conf_wei[1], fit_wei$AIC)

    # Negbin-lognormal model
    minus.log.lik.log.neg <- function(params) {
      density <- params[1]
      meanlog <- params[2]
      sdlog <- params[3]
      size <- params[4]
      - (sum(dlnorm(stay[cens == 0], meanlog = meanlog, sdlog = sdlog, log = TRUE)) +
           sum(1 - plnorm(stay[cens == 1], meanlog = meanlog, sdlog = sdlog, log.p = TRUE)) +
           sum(dnbinom(Y, size = size, mu = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion / exp(meanlog + sdlog ^ 2 / 2), log = TRUE)
           ))
    }

    fit_log <- optim(par = c(5, 1, 1, 5), fn = minus.log.lik.log.neg, hessian = TRUE)

    fit_log$AIC <- 2 * fit_log$value + length(fit_log$par) * 2
    conf_log <- 1.96 * sqrt(diag(solve(fit_log$hessian)))
    sd_log <- sqrt(diag(solve(fit_log$hessian)))

    res$Density[8, 2:6] <- c(fit_log$par[1], sd_log[1], fit_log$par[1] - conf_log[1], fit_log$par[1] + conf_log[1], fit_log$AIC)
    res$Density <- res$Density[order(res$Density$AIC), ]
    res$Density$delta_AIC <- res$Density$AIC - min(res$Density$AIC)

    return(res)
  }

  # RAD-RESTモデル -------------------------------------------------------------

  if(model == "RAD-REST") {

    res <- vector("list", length = 2)
    names(res) <- c("Density", "Stay")
    res$Actvity <- model_act@act

    Y <- data.frame(station_effort_data)[, grep("y_", colnames(station_effort_data))]
    N <- station_effort_data %>% pull(N)
    trials <- apply(Y, 1, sum)
    nparams <- ncol(Y)


    res$Density <- data.frame(ModelName = NA, Density = NA, SD = NA, Lower = NA, Upper = NA, AIC = NA, delta_AIC = NA)
    res$Stay <- vector("list", length = 4)
    names(res$Stay) <- c("Exponential", "Gammma", "Weibull", "Lognormal")
    res$Density[1:8, 1] <- c(paste("Poisson", c("Exponential", "Gammma", "Weibull", "Lognormal"), sep = "-"),
                             paste("NegBinom", c("Exponential", "Gammma", "Weibull", "Lognormal"), sep = "-"))


    # Poisson-exponential model
    minus.log.lik.exp <- function(params) {
      density <- params[1]
      rate <- params[2]
      alpha <- params[3:(3 + (nparams) - 1)]
      - (sum(dexp(stay[cens == 0], rate = rate, log = TRUE)) +
           sum(1 - pexp(stay[cens == 1], rate = rate, log.p = TRUE)) +
           sum(ddirmnom(Y, size = trials, alpha = alpha, log = TRUE)) +
           sum(dpois(N, lambda = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion * rate
                     / (sum(alpha / sum(alpha) * 0:(length(alpha) - 1)))
                     , log = TRUE)
           ))
    }
    fit_exp <- optim(par = c(density = 5, rate = 1, alpha = rep(10, nparams)), fn = minus.log.lik.exp, hessian = TRUE,
                     method = "L-BFGS-B",
                     lower = c(0.01, 0.01, rep(0.01, nparams)))

    sum(fit_exp$par[3:5] / sum(fit_exp$par[3:5]) * 0:(length(fit_exp$par[3:5]) - 1))
    sum(apply(Y, 2, sum) / sum(apply(Y, 2, sum)) * 0:2)

    fit_exp$AIC <- 2 * fit_exp$value + length(fit_exp$par) * 2
    sd_exp <- sqrt(diag(solve(fit_exp$hessian)))
    conf_exp <- 1.96 * sqrt(diag(solve(fit_exp$hessian)))

    res$Density[1, 2:6] <- c(fit_exp$par[1], sd_exp[1], fit_exp$par[1] - conf_exp[1], fit_exp$par[1] + conf_exp[1], fit_exp$AIC)
    res$Stay$Exponential <- fit_exp$par[2]

    # Poisson-Gamma model
    minus.log.lik.gam <- function(params) {
      density <- params[1]
      rate <- params[2]
      shape <- params[3]
      alpha <- params[4:(4 + (nparams) - 1)]
      - (sum(dgamma(stay[cens == 0], shape = shape, rate = rate, log = TRUE)) +
           sum(1 - pgamma(stay[cens == 1], shape = shape, rate = rate, log.p = TRUE)) +
           sum(ddirmnom(Y, size = trials, alpha = alpha, log = TRUE)) +
           sum(dpois(N, lambda = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion * rate /shape
                     / (sum(alpha / sum(alpha) * 0:(length(alpha) - 1))), log = TRUE)
           ))
    }

    fit_gam <- optim(par = c(density = 5, shape = 1, rate = 1, alpha = rep(10, nparams)), fn = minus.log.lik.gam, hessian = TRUE,
                     method = "L-BFGS-B",
                     lower = c(0.01, 0.01, 0.01, rep(0.01, nparams)))

    fit_gam$AIC <- 2 * fit_gam$value + length(fit_gam$par) * 2
    sd_gam <- sqrt(diag(solve(fit_gam$hessian)))
    conf_gam <- 1.96 * sqrt(diag(solve(fit_gam$hessian)))

    res$Density[2, 2:6] <- c(fit_gam$par[1], sd_gam[1], fit_gam$par[1] - conf_gam[1], fit_gam$par[1] + conf_gam[1], fit_gam$AIC)
    res$Stay$Gammma <- fit_gam$par[2:3]

    # Poisson-Weibull model
    minus.log.lik.wei <- function(params) {
      density <- params[1]
      shape <- params[2]
      scale <- params[3]
      alpha <- params[4:(4 + (nparams) - 1)]
      - (sum(dweibull(stay[cens == 0], shape = shape, scale = scale, log = TRUE)) +
           sum(1 - pweibull(stay[cens == 1], shape = shape, scale = scale, log.p = TRUE)) +
           sum(ddirmnom(Y, size = trials, alpha = alpha, log = TRUE)) +
           sum(dpois(N, lambda = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion / (scale * gamma(1 + 1 / shape))
                     / (sum(alpha / sum(alpha) * 0:(length(alpha) - 1))), log = TRUE)
           ))
    }

    fit_wei <- optim(par = c(density = 5, shape = 0.5, scale = 10, alpha = rep(10, nparams)), fn = minus.log.lik.wei, hessian = TRUE,
                     method = "L-BFGS-B",
                     lower = c(0.01, 0.01, 0.01, rep(0.01, nparams)))

    fit_wei$AIC <- 2 * fit_wei$value + length(fit_wei$par) * 2
    sd_wei <- sqrt(diag(solve(fit_gam$hessian)))
    conf_wei <- 1.96 * sqrt(diag(solve(fit_wei$hessian)))

    res$Density[3, 2:6] <- c(fit_wei$par[1], sd_wei[1], fit_wei$par[1] - conf_wei[1], fit_wei$par[1] + conf_wei[1], fit_wei$AIC)
    res$Stay$Weibull <- fit_wei$par[2:3]

    # Poisson-lognormal model
    minus.log.lik.log <- function(params) {
      density <- params[1]
      meanlog <- params[2]
      sdlog <- params[3]
      alpha <- params[4:(4 + (nparams) - 1)]
      - (sum(dlnorm(stay[cens == 0], meanlog = meanlog, sdlog = sdlog, log = TRUE)) +
           sum(1 - plnorm(stay[cens == 1], meanlog = meanlog, sdlog = sdlog, log.p = TRUE)) +
           sum(ddirmnom(Y, size = trials, alpha = alpha, log = TRUE)) +
           sum(dpois(N, lambda = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion / exp(meanlog + sdlog ^ 2 / 2)
                     / (sum(alpha / sum(alpha) * 0:(length(alpha) - 1))), log = TRUE)
           ))
    }

    fit_log <- optim(par = c(density = 5, meanlog = 0.1, sdlog = 0.1, alpha = rep(10, nparams)), fn = minus.log.lik.log, hessian = TRUE,
                     method = "L-BFGS-B",
                     lower = c(0.01, 0.01, 0.01, rep(0.01, nparams)))

    fit_log$AIC <- 2 * fit_log$value + length(fit_log$par) * 2
    sd_log <- sqrt(diag(solve(fit_log$hessian)))
    conf_log <- 1.96 * sqrt(diag(solve(fit_log$hessian)))

    res$Density[4, 2:6] <- c(fit_log$par[1], sd_log[1], fit_log$par[1] - conf_log[1], fit_log$par[1] + conf_log[1], fit_log$AIC)
    res$Stay$Lognormal <- fit_log$par[2:3]

    res$Density <- res$Density[order(res$Density$AIC), ]
    res$Density$delta_AIC <- res$Density$AIC - min(res$Density$AIC)


    # Negbin-exponential model
    minus.log.lik.exp <- function(params) {
      density <- params[1]
      rate <- params[2]
      size <- params[3]
      alpha <- params[4:(4 + (nparams) - 1)]
      - (sum(dexp(stay[cens == 0], rate = rate, log = TRUE)) +
           sum(1 - pexp(stay[cens == 1], rate = rate, log.p = TRUE)) +
           sum(ddirmnom(Y, size = trials, alpha = alpha, log = TRUE)) +
           sum(dnbinom(N, size = size, mu = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion * rate
                       / (sum(alpha / sum(alpha) * 0:(length(alpha) - 1))), log = TRUE)
           ))
    }

    fit_exp <- optim(par = c(5, 0.1, 10, alpha = rep(5, nparams)), fn = minus.log.lik.exp, hessian = TRUE, method = "L-BFGS-B",
                     lower = c(0.01, 0.01, 0.01, rep(0.01, nparams)))

    fit_exp$AIC <- 2 * fit_exp$value + length(fit_exp$par) * 2
    sd_exp <- sqrt(diag(solve(fit_exp$hessian)))
    conf_exp <- 1.96 * sqrt(diag(solve(fit_exp$hessian)))

    res$Density[5, 2:6] <- c(fit_exp$par[1], sd_exp[1], fit_exp$par[1] - conf_exp[1], fit_exp$par[1] + conf_exp[1], fit_exp$AIC)
    res$Stay$Exponential <- fit_exp$par[2]


    # Negbin-Gamma model
    minus.log.lik.gam.neg <- function(params) {
      density <- params[1]
      rate <- params[2]
      shape <- params[3]
      size <- params[4]
      alpha <- params[5:(5 + (nparams) - 1)]
      - (sum(dgamma(stay[cens == 0], shape = shape, rate = rate, log = TRUE)) +
           sum(1 - pgamma(stay[cens == 1], shape = shape, rate = rate, log.p = TRUE)) +
           sum(ddirmnom(Y, size = trials, alpha = alpha, log = TRUE)) +
           sum(dnbinom(N, size = size, mu = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion * rate / shape
                       / (sum(alpha / sum(alpha) * 0:(length(alpha) - 1))), log = TRUE)
           ))
    }

    fit_gam <- optim(par = c(5, 0.5, 0.1, 10, alpha = rep(10, nparams)), fn = minus.log.lik.gam.neg, hessian = TRUE,
                     method = "L-BFGS-B",
                     lower = c(0.01, 0.01, 0.01, 0.01, rep(0.01, nparams)))

    fit_gam$AIC <- 2 * fit_gam$value + length(fit_gam$par) * 2
    conf_gam <- 1.96 * sqrt(diag(solve(fit_gam$hessian)))
    sd_gam <- sqrt(diag(solve(fit_gam$hessian)))

    res$Density[6, 2:6] <- c(fit_gam$par[1], sd_gam[1], fit_gam$par[1] - conf_gam[1], fit_gam$par[1] + conf_gam[1], fit_gam$AIC)

    # Negbin-Weibull model
    minus.log.lik.wei.neg <- function(params) {
      density <- params[1]
      shape <- params[2]
      scale <- params[3]
      size <- params[4]
      alpha <- params[5:(5 + (nparams) - 1)]
      - (sum(dweibull(stay[cens == 0], shape = shape, scale = scale, log = TRUE)) +
           sum(1 - pweibull(stay[cens == 1], shape = shape, scale = scale, log.p = TRUE)) +
           sum(ddirmnom(Y, size = trials, alpha = alpha, log = TRUE)) +
           sum(dnbinom(N, size = size, mu = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion / (scale * gamma(1 + 1 / shape))
                       / (sum(alpha / sum(alpha) * 0:(length(alpha) - 1))), log = TRUE)
           ))
    }

    fit_wei <- optim(par = c(5, 1, 5, 10, alpha = rep(10, nparams)), fn = minus.log.lik.wei.neg, hessian = TRUE,
                     method = "L-BFGS-B",
                     lower = c(0.01, 0.01, 0.01, 0.01, rep(0.01, nparams)))

    fit_wei$AIC <- 2 * fit_wei$value + length(fit_wei$par) * 2
    conf_wei <- 1.96 * sqrt(diag(solve(fit_wei$hessian)))
    sd_wei <- sqrt(diag(solve(fit_wei$hessian)))
    res$Density[7, 2:6] <- c(fit_wei$par[1], sd_wei[1], fit_wei$par[1] - conf_wei[1], fit_wei$par[1] + conf_wei[1], fit_wei$AIC)

    # Negbin-lognormal model
    minus.log.lik.log.neg <- function(params) {
      density <- params[1]
      meanlog <- params[2]
      sdlog <- params[3]
      size <- params[4]
      alpha <- params[5:(5 + (nparams) - 1)]
      - (sum(dlnorm(stay[cens == 0], meanlog = meanlog, sdlog = sdlog, log = TRUE)) +
           sum(1 - plnorm(stay[cens == 1], meanlog = meanlog, sdlog = sdlog, log.p = TRUE)) +
           sum(ddirmnom(Y, size = trials, alpha = alpha, log = TRUE)) +
           sum(dnbinom(N, size = size, mu = density * focal_area / 1000000 * effort * 60 * 60 * 24 * activity_proportion / exp(meanlog + sdlog ^ 2 / 2)
                       / (sum(alpha / sum(alpha) * 0:(length(alpha) - 1))), log = TRUE)
           ))
    }

    fit_log <- optim(par = c(5, 2, 1, 10, alpha = rep(10, nparams)), fn = minus.log.lik.log.neg, hessian = TRUE,
                     method = "L-BFGS-B",
                     lower = c(0.01, 0.01, 0.01, 0.01, rep(0.01, nparams)))

    fit_log$AIC <- 2 * fit_log$value + length(fit_log$par) * 2
    conf_log <- 1.96 * sqrt(diag(solve(fit_log$hessian)))
    sd_log <- sqrt(diag(solve(fit_log$hessian)))

    res$Density[8, 2:6] <- c(fit_log$par[1], sd_log[1], fit_log$par[1] - conf_log[1], fit_log$par[1] + conf_log[1], fit_log$AIC)
    res$Density <- res$Density[order(res$Density$AIC), ]
    res$Density$delta_AIC <- res$Density$AIC - min(res$Density$AIC)
    return(res)
  }
}
Effort <- Cens <- Stay <- NULL
