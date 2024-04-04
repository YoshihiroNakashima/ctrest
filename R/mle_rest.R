#' Maximu-likelihood estimation of animal density using the REST model
#'
#' @param npass_effort_data ##
#' @param stay_data ##
#' @param activity_data ##
#' @param col_name_npass ##
#' @param col_name_effort #E#
#' @param col_name_stay ##
#' @param col_name_cens ##
#' @param focal_area ##
#'
#' @return estimated density
#' @export
#' @import activity
#' @examples #' mle_rest(
#' npass_effort_data = npass_effort_data,
#' stay_data = stay_data,
#' activity_data = activity_data,
#' col_name_npass = "npass",
#' col_name_effort = "effort",
#' col_name_stay = "stay",
#' col_name_cens = "cens",
#' focal_area)
#' )
mle_rest <- function(npass_effort_data,
                     stay_data,
                     activity_data,
                     col_name_npass = "npass",
                     col_name_effort = "effort",
                     col_name_stay = "stay",
                     col_name_cens = "cens",
                     focal_area) {

  act <- activity::fitact(
      activity_data,
      bw = bwcalc(activity_data, K = 3),
      adj = 1.5,
      reps = 1
    )@act  # package activity

  stay <- stay_data %>% pull(!!sym(col_name_stay))
  cens <- stay_data %>% pull(!!sym(col_name_cens))
  npass <- npass_effort_data %>%  pull(!!sym(col_name_npass))
  effort <- npass_effort_data %>%  pull(!!sym(col_name_effort))

  res <- vector("list", length = 2)
  names(res) <- c("Density", "Stay")

  res$Density <- data.frame(ModelName = NA, Density = NA, Lower = NA, Upper = NA, AIC = NA, delta_AIC = NA)
  res$Stay <- vector("list", length = 4)
  names(res$Stay) <- c("Exponential", "Gammma", "Weibull", "Lognormal")
  res$Density[1:8, 1] <- c(paste("Poisson", c("Exponential", "Gammma", "Weibull", "Lognormal"), sep = "-"),
                           paste("NegBinom", c("Exponential", "Gammma", "Weibull", "Lognormal"), sep = "-"))

  options(warn = -1)
  # Poisson-exponential model
  minus.log.lik.exp <- function(params) {
    density <- params[1]
    rate <- params[2]
    - (sum(dexp(stay[cens == 0], rate = rate, log = TRUE)) +
      sum(1 - pexp(stay[cens == 1], rate = rate, log = TRUE)) +
      sum(dpois(npass, lambda = density * focal_area / 1000000 * effort * 60 * 60 * 24 * act * rate, log = TRUE)
      ))
  }
  fit_exp <- optim(par = c(density = 5, rate = 0.1), fn = minus.log.lik.exp, hessian = TRUE)

  fit_exp$AIC <- 2 * fit_exp$value + length(fit_exp$par) * 2
  conf_exp <- 1.96 * sqrt(diag(solve(fit_exp$hessian)))

  res$Density[1, 2:5] <- c(fit_exp$par[1], fit_exp$par[1] - conf_exp[1], fit_exp$par[1] + conf_exp[1], fit_exp$AIC)
  res$Stay$Exponential <- fit_exp$par[2]

  # Poisson-Gamma model
  minus.log.lik.gam <- function(params) {
    density <- params[1]
    rate <- params[2]
    shape <- params[3]
    - (sum(dgamma(stay[cens == 0], shape = shape, rate = rate, log = TRUE)) +
         sum(1 - pgamma(stay[cens == 1], shape = shape, rate = rate, log = TRUE)) +
         sum(dpois(npass, lambda = density * focal_area / 1000000 * effort * 60 * 60 * 24 * act * rate /shape, log = TRUE)
         ))
  }

  fit_gam <- optim(par = c(density = 5, shape = 0.5, rate = 0.1), fn = minus.log.lik.gam, hessian = TRUE)

  fit_gam$AIC <- 2 * fit_gam$value + length(fit_gam$par) * 2
  conf_gam <- 1.96 * sqrt(diag(solve(fit_gam$hessian)))

  res$Density[2, 2:5] <- c(fit_gam$par[1], fit_gam$par[1] - conf_gam[1], fit_gam$par[1] + conf_gam[1], fit_gam$AIC)
  res$Stay$Gammma <- fit_gam$par[2:3]

  # Poisson-Weibull model
  minus.log.lik.wei <- function(params) {
    density <- params[1]
    shape <- params[2]
    scale <- params[3]
    - (sum(dweibull(stay[cens == 0], shape = shape, scale = scale, log = TRUE)) +
         sum(1 - pweibull(stay[cens == 1], shape = shape, scale = scale, log = TRUE)) +
         sum(dpois(npass, lambda = density * focal_area / 1000000 * effort * 60 * 60 * 24 * act / (scale * gamma(1 + 1 / shape)), log = TRUE)
         ))
  }

  fit_wei <- optim(par = c(density = 5, shape = 0.5, scale = 10), fn = minus.log.lik.wei, hessian = TRUE)

  fit_wei$AIC <- 2 * fit_wei$value + length(fit_wei$par) * 2
  conf_wei <- 1.96 * sqrt(diag(solve(fit_wei$hessian)))

  res$Density[3, 2:5] <- c(fit_wei$par[1], fit_wei$par[1] - conf_wei[1], fit_wei$par[1] + conf_wei[1], fit_wei$AIC)
  res$Stay$Weibull <- fit_wei$par[2:3]

  # Poisson-lognormal model
  minus.log.lik.log <- function(params) {
    density <- params[1]
    meanlog <- params[2]
    sdlog <- params[3]
    - (sum(dlnorm(stay[cens == 0], meanlog = meanlog, sdlog = sdlog, log = TRUE)) +
         sum(1 - plnorm(stay[cens == 1], meanlog = meanlog, sdlog = sdlog, log = TRUE)) +
         sum(dpois(npass, lambda = density * focal_area / 1000000 * effort * 60 * 60 * 24 * act / exp(meanlog + sdlog ^ 2 / 2), log = TRUE)
         ))
  }

  fit_log <- optim(par = c(density = 5, meanlog = 2, sdlog = 1), fn = minus.log.lik.log, hessian = TRUE)

  fit_log$AIC <- 2 * fit_log$value + length(fit_log$par) * 2
  conf_log <- 1.96 * sqrt(diag(solve(fit_log$hessian)))

  res$Density[4, 2:5] <- c(fit_log$par[1], fit_log$par[1] - conf_log[1], fit_log$par[1] + conf_log[1], fit_log$AIC)
  res$Stay$Lognormal <- fit_log$par[2:3]

  res$Density <- res$Density[order(res$Density$AIC), ]
  res$Density$delta_AIC <- res$Density$AIC - min(res$Density$AIC)


  # Negbin-exponential model
  minus.log.lik.exp <- function(params) {
    density <- params[1]
    rate <- params[2]
    size <- params[3]
    - (sum(dexp(stay[cens == 0], rate = rate, log = TRUE)) +
         sum(1 - pexp(stay[cens == 1], rate = rate, log = TRUE)) +
         sum(dnbinom(npass, size = size, mu = density * focal_area / 1000000 * effort * 60 * 60 * 24 * act * rate, log = TRUE)
         ))
  }
  fit_exp <- optim(par = c(5, 0.1, 10), fn = minus.log.lik.exp, hessian = TRUE)

  fit_exp$AIC <- 2 * fit_exp$value + length(fit_exp$par) * 2
  conf_exp <- 1.96 * sqrt(diag(solve(fit_exp$hessian)))

  res$Density[5, 2:5] <- c(fit_exp$par[1], fit_exp$par[1] - conf_exp[1], fit_exp$par[1] + conf_exp[1], fit_exp$AIC)
  res$Stay$Exponential <- fit_exp$par[2]


  # Negbin-Gamma model
  minus.log.lik.gam.neg <- function(params) {
    density <- params[1]
    rate <- params[2]
    shape <- params[3]
    size <- params[4]
    - (sum(dgamma(stay[cens == 0], shape = shape, rate = rate, log = TRUE)) +
         sum(1 - pgamma(stay[cens == 1], shape = shape, rate = rate, log = TRUE)) +
         sum(dnbinom(npass, size = size, mu = density * focal_area / 1000000 * effort * 60 * 60 * 24 * act * rate / shape, log = TRUE)
         ))
  }

  fit_gam <- optim(par = c(5, 0.5, 0.1, 10), fn = minus.log.lik.gam.neg, hessian = TRUE)

  fit_gam$AIC <- 2 * fit_gam$value + length(fit_gam$par) * 2
  conf_gam <- 1.96 * sqrt(diag(solve(fit_gam$hessian)))

  res$Density[6, 2:5] <- c(fit_gam$par[1], fit_gam$par[1] - conf_gam[1], fit_gam$par[1] + conf_gam[1], fit_gam$AIC)

  # Negbin-Weibull model
  minus.log.lik.wei.neg <- function(params) {
    density <- params[1]
    shape <- params[2]
    scale <- params[3]
    size <- params[4]
    - (sum(dweibull(stay[cens == 0], shape = shape, scale = scale, log = TRUE)) +
         sum(1 - pweibull(stay[cens == 1], shape = shape, scale = scale, log = TRUE)) +
         sum(dnbinom(npass, size = size, mu = density * focal_area / 1000000 * effort * 60 * 60 * 24 * act / (scale * gamma(1 + 1 / shape)), log = TRUE)
         ))
  }

  fit_wei <- optim(par = c(5, 1, 5, 10), fn = minus.log.lik.wei.neg, hessian = TRUE)

  fit_wei$AIC <- 2 * fit_wei$value + length(fit_wei$par) * 2
  conf_wei <- 1.96 * sqrt(diag(solve(fit_wei$hessian)))

  res$Density[7, 2:5] <- c(fit_wei$par[1], fit_wei$par[1] - conf_wei[1], fit_wei$par[1] + conf_wei[1], fit_wei$AIC)

  # Negbin-lognormal model
  minus.log.lik.log.neg <- function(params) {
    density <- params[1]
    meanlog <- params[2]
    sdlog <- params[3]
    size <- params[4]
    - (sum(dlnorm(stay[cens == 0], meanlog = meanlog, sdlog = sdlog, log = TRUE)) +
         sum(1 - plnorm(stay[cens == 1], meanlog = meanlog, sdlog = sdlog, log = TRUE)) +
         sum(dnbinom(npass, size = size, mu = density * focal_area / 1000000 * effort * 60 * 60 * 24 * act / exp(meanlog + sdlog ^ 2 / 2), log = TRUE)
         ))
  }

  fit_log <- optim(par = c(5, 2, 1, 10), fn = minus.log.lik.log.neg, hessian = TRUE)

  fit_log$AIC <- 2 * fit_log$value + length(fit_log$par) * 2
  conf_log <- 1.96 * sqrt(diag(solve(fit_log$hessian)))

  res$Density[8, 2:5] <- c(fit_log$par[1], fit_log$par[1] - conf_log[1], fit_log$par[1] + conf_log[1], fit_log$AIC)
  res$Density <- res$Density[order(res$Density$AIC), ]
  res$Density$delta_AIC <- res$Density$AIC - min(res$Density$AIC)
  return(res)
}

