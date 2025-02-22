#' Model selection for the residence time model based on WAIC for REST/REST-RAD models
#'
#' @param formula_stay Model formula for the residence time within a focal area. e.g. Stay ~ 1 + x1
#' @param random_effect random_effect effect on the parameter of staying time model (as a character). e.g. random_effect = "Station"
#' @param station_data A data frame containing information for each camera station in each row
#' @param stay_data A data frame containing the time spent within a focal area
#' @param col_name_cens Column name of stay_data containing information on censored (1) or not (0)
#' @param family The probability distribution of the time spent. Specify the probability distribution selected by the bayes_stay_selection function.
#' @param plot If TRUE, plots the expected values of residence times.
#' @param iter The number of iterations. The default value is 2000.
#' @param cores The number of cores used for calculations.
#' @param warmup The number of warmup iterations. The default value is NULL (half of the iterations will be used for warmup).
#' @param chains The number of chains. The default value is 2.
#' @param thin The thinning interval. The default value is 1.
#' @param all_comb If TRUE, models with all combinations of covariates are compared. If FALSE, only the designated model in model_formula is run.
#' @param target_species Species name of interest.
#' @return WAIC values for each model
#' @import dplyr ggplot2 nimble
#' @importFrom tidyr unite extract
#' @importFrom purrr map
#' @importFrom stringr str_extract
#' @importFrom stats median model.response quantile
#' @export
#' @examples
#' stay_data <- format_stay(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_stay = "Stay",
#'   col_name_cens = "Cens",
#'   target_species = "SP01"
#' )
#' library(nimble)
#' bayes_stay_selection(formula_stay = Stay ~ 1 + x1,
#'                     station_data = station_data,
#'                     stay_data = stay_data,
#'                     col_name_cens = "Cens",
#'                     family = "lognormal",
#'                     plot = TRUE,
#'                     iter = 2000,
#'                     warmup = 1000,
#'                     chains = 2,
#'                     thin = 1,
#'                     all_comb = FALSE,
#'                     target_species = "SP01")

bayes_stay_selection <- function(
    formula_stay = Stay ~ 1,
    random_effect = NULL,
    station_data = station_data,
    stay_data = stay_data,
    col_name_cens = "Cens",
    family = "lognormal",
    plot = TRUE,
    cores = 1,
    iter = 3000,
    warmup = NULL,
    chains = 3,
    thin = 1,
    all_comb = TRUE,
    target_species = NULL
) {

  if(family=="lognormal" || family=="log-normal" || family=="gamma" || family=="weibull" || family=="exponential"){
  } else{
    stop(paste0("Input family type(", family,") is incorrect."))
  }


  # Define functions --------------------------------------------------------

  bit.test <- function(number, n) {
    (number %/% (2^n)) %% 2
  }

  full_terms <- function(x) {
    lapply(1:2^length(x),
           function(i) {
             r <- ""
             for (j in 1:length(x)) {
               r <- paste(r,
                          ifelse(bit.test((i - 1), (j - 1)),
                                 paste(x[j], " + ", sep=""), ""),
                          sep="")
             }
             r <- paste("Stay ~ 1 + ", r, "1", sep = "") #切片がベータの場合
             r <- strsplit(r, " \\+ 1$")[[1]][1]
           }
    )
  }

  # Define data for stay ----------------------------------------------------

  stay_data_join <- stay_data %>%
    left_join(station_data, by = intersect(names(stay_data), names(station_data)))

  if (family == "lognormal" ||
      family == "gamma" || family == "weibull" || family == "exponential") {

  } else {
    stop(paste0("Input family type(", family, ") is incorrect."))
  }

  target_species <- sort(target_species)

  # Extract variable names
  vars_stay <- all.vars(formula_stay)
  response_stay <- vars_stay[1]
  predictors_stay <- vars_stay[-1]

  model_frame_stay <- model.frame(formula_stay, stay_data_join)
  X_stay <- model.matrix(as.formula(formula_stay), model_frame_stay) # model.matrix(formula_stay, model_frame_stay)
  t <- model.response(model_frame_stay)
  censored <- stay_data_join[["Cens"]]

  # N of random effects
  if (!is.null(random_effect)) {
    levels <- unique(stay_data_join[[random_effect]])
    nLevels_stay <- length(levels)
  } else {
    nLevels_stay <- 0
  }

  # N of covariates
  N_station_species <- nrow(station_data)
  N_station <- length(unique(station_data$Station))
  nPreds_stay <- ncol(X_stay)

  names(t) <- NULL
  c_time <- t
  c_time[censored == 0] <- c_time[censored == 0] + 1
  t[censored == 1] <- NA
  N_stay <- length(t)

  formula_stay_all <- list(0)
  if (all_comb == TRUE) {
    formula_stay_all <- full_terms(c(predictors_stay))
  } else {
    formula_stay_all[[1]] <- formula_stay
  }

  raw_samples <- tidy_samples <- list(0)
  waic <- numeric(0)

  for(k in 1:length(formula_stay_all)) {

    formula_stay <- formula_stay_all[[k]]

    model_frame_stay <- model.frame(formula_stay, stay_data)
    X_stay <- model.matrix(as.formula(formula_stay_all[[k]]), model_frame_stay)
    t <- model.response(model_frame_stay)
    censored <- stay_data[[col_name_cens]]

    if (!is.null(random_effect)) {
      levels <- unique(stay_data[[random_effect]])
      nLevels_stay <- length(levels)
    } else {
      nLevels_stay <- 0
    }

    nPreds_stay <- ncol(X_stay)
    N_station <- nrow(station_data)

    names(t) <- NULL
    c_time <- t
    c_time[censored == 0] <- c_time[censored == 0] + 1
    t[censored == 1] <- NA
    N_stay <- length(t)

    data_stay <- list(t = t, censored = censored)
    cons_stay <- list(N_stay = N_stay, nPreds_stay = nPreds_stay, N_station = N_station, c_time = c_time, family = family)

    if (nPreds_stay > 1) {
      data_stay$X_stay <- X_stay
    }
    if (!is.null(random_effect)) {
      cons_stay$group <- as.numeric(factor(stay_data[[random_effect]]))
      cons_stay$nLevels_stay <- nLevels_stay
    } else {
      cons_stay$nLevels_stay <- 0
    }
    if (!is.null(random_effect)) {
      cons_stay$group <- as.numeric(factor(stay_data_join[[random_effect]]))
      cons_stay$nLevels_stay <- nLevels_stay
    } else {
      cons_stay$nLevels_stay <- 0
    }
    if (nPreds_stay > 1) {
      data_stay$X_stay <- X_stay
    }

    nSpecies <- length(target_species)
    species_stay <- stay_data_join %>% pull(Species) %>% factor() %>% as.numeric()

    cons_stay$nSpecies <- nSpecies
    cons_stay$species_stay <- species_stay

    ni <- iter
    nt <- thin
    nc <- chains
    nb <- warmup
    nmcmc <- (iter - warmup) * chains / thin

    # MCMC
    code <- nimbleCode({
      for(i in 1:N_stay) {
        censored[i] ~ dinterval(t[i], c_time[i])

        if(family == "exponential") {
          t[i] ~ dexp(rate = 1/exp(log_scale[i]))
        }
        if(family == "gamma") {
          t[i] ~ dgamma(shape = theta[species_stay[i]], rate = 1/exp(log_scale[i]))  # theta[species_stay[i]]を指数変換
        }
        if(family == "lognormal") {
          t[i] ~ dlnorm(meanlog = log_scale[i], sdlog = theta[species_stay[i]])  # theta[species_stay[i]]を指数変換
        }
        if(family == "weibull") {
          t[i] ~ dweibull(shape = theta[species_stay[i]], scale = exp(log_scale[i]))  # thetaを指数変換
        }

        if (nPreds_stay > 1) {
          if (nLevels_stay == 0) {
            log_scale[i] <- inprod(beta_stay[1:nPreds_stay] + v_stay[species_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay])
          } else {
            log_scale[i] <- inprod(beta_stay[1:nPreds_stay] + v_stay[species_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + u[group[i]]
          }
        } else {
          if (nLevels_stay == 0) {
            log_scale[i] <- beta_stay[1] + v_stay[species_stay[i], 1]
          } else {
            log_scale[i] <- beta_stay[1] + v_stay[species_stay[i], 1] + u[group[i]]
          }
        }
      }

      # 事前分布の修正
      for(j in 1:nPreds_stay) {
        beta_stay[j] ~ dnorm(0, sd = 100)
      }

      for(m in 1:nSpecies) {
        theta[m] ~ dgamma(shape_theta, rate_theta)  # dt()からdnorm()に変更
      }
      shape_theta ~ dgamma(0.01, 0.01)
      rate_theta ~ dgamma(0.01, 0.01)

      if (nLevels_stay > 0) {
        for(k in 1:nLevels_stay) {
          u[k] ~ dnorm(0, sd = sigma_u)
        }
        sigma_u ~ T(dnorm(0, sd = 2.5), 0, 5)  # 非負制約を追加
      }

      for(m in 1:nSpecies) {
        for(j in 1:nPreds_stay) {
          v_stay[m, j] ~ dnorm(0, sd = sigma_v[j])
        }
      }
      for(j in 1:nPreds_stay) {
        sigma_v[j] ~ T(dnorm(0, sd = 2.5), 0, 5)  # 非負制約を追加
      }


      # 期待値の計算を修正
      if (nPreds_stay == 1) {
        for(m in 1:nSpecies) {
          for(i in 1:N_station) {
            if(family == "exponential") mean_stay[i, m] <- exp(beta_stay[1] + v_stay[m, 1])

            if(family == "gamma") mean_stay[i, m] <- theta[m] * exp(beta_stay[1] + v_stay[m, 1])

            if(family == "lognormal") mean_stay[i, m] <- exp(beta_stay[1] + v_stay[m, 1] + (theta[m] ^ 2) / 2)

            if(family == "weibull") mean_stay[i, m] <- exp(lgamma(1 + 1 / theta[m]) )* exp(beta_stay[1] + v_stay[m, 1])
          }
        }
      }

      if (nPreds_stay > 1) {
        for(m in 1:nSpecies) {
          for(i in 1:N_station) {
            if(family == "exponential") {
              mean_stay[i, m] <- exp(inprod(beta_stay[1:nPreds_stay] + v_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
            }
            if(family == "gamma") {
              mean_stay[i, m] <- theta[m] * exp(inprod(beta_stay[1:nPreds_stay] + v_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
            }
            if(family == "lognormal") {
              mean_stay[i, m] <- exp(inprod(beta_stay[1:nPreds_stay] + v_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + theta[m] ^ 2 / 2)
            }
            if(family == "weibull") {
              mean_stay[i, m] <- exp(lgamma(1 + 1 / theta[m])) * exp(inprod(beta_stay[1:nPreds_stay] + v_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
            }
          }
        }
      }

      # Posterior predictive check
      for (i in 1:N_stay) {
        pred_t[i] ~ dlnorm(meanlog = log_scale[i], sdlog = theta[species_stay[i]])

        # 観測データの適合度統計量：対数尤度を使用
        loglike_obs[i] <- (1 - step(censored[i] - 0.5)) *
          dlnorm(t[i], meanlog = log_scale[i], sdlog = theta[species_stay[i]], log = TRUE) +
          step(censored[i] - 0.5) *
          plnorm(c_time[i], meanlog = log_scale[i], sdlog = theta[species_stay[i]], lower.tail = TRUE, log.p = TRUE)

        # 予測データの適合度統計量：対数尤度を使用
        loglike_pred[i] <- dlnorm(pred_t[i], meanlog = log_scale[i], sdlog = theta[species_stay[i]], log = TRUE)
      }

      sum_loglike_obs <- sum(loglike_obs[1:N_stay])
      sum_loglike_pred <- sum(loglike_pred[1:N_stay])

    })

    cat("Compiling the model. This may take a moment...\n")

    this_cluster <- makeCluster(nc)

    run_MCMC_pre <- function(info, data, constants, code, params, ni, nt, nb) {
      myModel <- nimbleModel(code = code,
                             data = data,
                             constants = constants,
                             inits = info$inits)
      CmyModel <- compileNimble(myModel)
      myMCMC <- buildMCMC(CmyModel, monitors = params)
      CmyMCMC <- compileNimble(myMCMC)
      results <- runMCMC(
        CmyMCMC,
        niter = ni,
        nburnin = nb,
        thin = nt,
        nchains = 1,
        setSeed = info$seed,
        samplesAsCodaMCMC = TRUE
      )
      results
    }

    inits_f <- function() {
      common_inits <- list(
        beta_stay = runif(nPreds_stay, -1, 1),
        t = ifelse(censored == 0, NA, c_time + runif(N_stay, 0.1, 1.0)),
        theta = runif(nSpecies, 0.5, 4.0),
        log_scale = runif(N_stay, 0.5, 4.0),
        v_stay = matrix(runif(nSpecies * nPreds_stay, -1, 1), nrow = nSpecies, ncol = nPreds_stay)
      )

      if (!is.null(random_effect)) {
        c(
          common_inits,
          list(
            u = runif(nLevels_stay, -1, 1),
            sigma_u = runif(1, 0.5, 2.5)
          )
        )
      } else {
        c(
          common_inits,
          list(sigma_v = runif(nPreds_stay, 0.5, 2.5))
        )
      }
    }

    per_chain_info <- lapply(1:nc, function(i) {
      list(
        seed = sample(1:9999, 1),
        inits = inits_f()
      )
    })
    params <- c("beta_stay", "log_scale", "mean_stay", "theta", "sum_loglike_obs", "sum_loglike_pred")

    cat("Running MCMC sampling. Please wait...\n")

    clusterEvalQ(this_cluster, {
      library(nimble)
    })

    clusterExport(this_cluster,
                  varlist = c("run_MCMC_pre"),
                  envir = environment())

    chain_output <- parLapply(cl = this_cluster, X = per_chain_info,
                              fun = run_MCMC_pre,
                              data = data_stay, code = code,
                              constants = cons_stay, params = params,
                              ni = ni, nt = nt, nb = nb
    )
    stopCluster(this_cluster)

    tidy_samples[[k]] <- MCMCvis::MCMCchains(chain_output) %>%
      as_tibble() %>%
      pivot_longer(cols = everything(),
                   names_to = "parameter",
                   values_to = "value") %>%
      group_by(parameter) %>%
      mutate(iteration = row_number()) %>%
      ungroup()

    pattern_scale <- paste0("log_scale\\[", 1:N_stay, "\\]")
    pattern_theta <- paste0("theta\\[", 1:nSpecies, "\\]")



    log_scale.samp <- tidy_samples[[k]] %>%
      filter(str_detect(parameter, "log_scale"))

    log_scale.samp <- matrix(log_scale.samp$value, nrow = nmcmc, byrow = TRUE)
    # library(bayesplot)
    #
    # # トレースプロット（全チェーン比較）
    # mcmc_trace(chain_output, pars = c("mean_stay[1, 1]"))

    # ベイズP値
    chi2_obs_total_samples <- tidy_samples[[k]] %>%
      filter(str_detect(parameter, "sum_loglike_obs"))

    chi2_pred_total_samples <- tidy_samples[[k]] %>%
      filter(str_detect(parameter, "sum_loglike_pred"))

    p_value <- mean(chi2_pred_total_samples > chi2_obs_total_samples, na.rm = TRUE)

    if(family != "exponential") tau.samp <- chain_output %>%
      map(~ .[, grep(paste(pattern_theta, collapse = "|"), colnames(.))]) %>%
      do.call(rbind, .)


    loglfstay <- matrix(NA, nrow(log_scale.samp), N_stay)

    if(family == "exponential") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- log(1 - pexp(c_time[j], rate = 1 /exp(log_scale.samp[ , j])))
        } else {
          loglfstay[, j] <- dexp(t[j], rate = 1 / exp(log_scale.samp[ , j]), log = T)
        }
      }
    }
    if(family == "gamma") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- log(1 - pgamma(c_time[j], shape = tau.samp[ , j],  scale = exp(log_scale.samp[ , j])))
        } else{
          loglfstay[, j] <- dgamma(t[j], shape = tau.samp[ , j], scale = exp(log_scale.samp[ , j]), log = T)
        }
      }
    }
    if(family == "lognormal" | family == "log-normal") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- log(1 - plnorm(c_time[j], log_scale.samp[ , j], tau.samp[ , j]))
        } else{
          loglfstay[, j] <- dlnorm(t[j], log_scale.samp[ , j], tau.samp[ , j], log = T)
        }
      }
    }
    if(family == "weibull") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- log(1 - pweibull(c_time[j], shape = tau.samp[ , j], scale = exp(log_scale.samp[ , j])))
        } else{
          loglfstay[, j] <- dweibull(t[j], shape = tau.samp[ ,j], scale = exp(log_scale.samp[ , j]), log = T)
        }
      }
    }

    lppd <- sum(log(colMeans(exp(loglfstay))))
    p.waic <- sum(apply(loglfstay, 2, var))
    # waic_mat <- loo::waic(loglfstay)waic_mat$estimates[3, 1] #
    waic[k] <- (- 2) * lppd + 2 * p.waic


    raw_samples[[k]] <- chain_output
  }

  WAIC <- data.frame(Model = as.character(unlist(formula_stay_all)), WAIC = waic) %>%
    arrange(WAIC)

  tidy_samples  <- tidy_samples[[which(WAIC[1, 1]  == unlist(formula_stay_all))]]

  best_stay_tidy_samples <- tidy_samples %>%
    filter(str_detect(parameter, "mean_stay")) %>%
    mutate(Species = rep(rep(target_species, each = N_station), nmcmc)) %>%
    mutate(Station = rep(station_data %>% pull(Station), nmcmc * nSpecies))


  if(nPreds_stay > 1) {
    mean_stay_data <- best_stay_tidy_samples %>%
      group_by(Station, Species) %>%
      summarise(
        median = median(value),
        mean = mean(value),
        sd = sd(value),
        cv = sd(value) / mean(value),
        lower = quantile(value, 0.025),
        upper = quantile(value, 0.975),
        .groups = "drop"
      ) %>%
      arrange(Species, Station)
  }
  if(nPreds_stay == 1) {
    mean_stay_data <- best_stay_tidy_samples %>%
      group_by(Species) %>%
      summarise(
        median = median(value),
        mean = mean(value),
        sd = sd(value),
        cv = sd(value) / mean(value),
        lower = quantile(value, 0.025),
        upper = quantile(value, 0.975),
        .groups = "drop"
      ) %>%
      arrange(Species, Station)
  }

  if (plot == TRUE) {
    if(nPreds_stay > 1) {
      g <- ggplot(mean_stay_data, aes(x = Station, y = median)) +
        geom_pointrange(aes(x = Station, y = median,
                            ymin = lower, ymax = upper)) +
        labs(
          title = paste("formula =", formula_stay, "(Median and 95% Credible Interval for mean staying time)"),
          x = "Station",
          y = "MeanStayingTime"
        ) +
        facet_wrap(~Species,  nrow = 4) +
        ylim(0, NA)
      print(g)
    }
    if(nPreds_stay == 1) {
      g <- ggplot(mean_stay_data, aes(x = Species, y = median)) +
        geom_pointrange(aes(x = Species, y = median,
                            ymin = lower, ymax = upper)) +
        labs(
          title = paste("formula =", formula_stay, "(Median and 95% Credible Interval for mean staying time)"),
          x = "Station",
          y = "MeanStayingTime"
        ) +
        ylim(0, NA)
      print(g)
    }
  }
  return(list(
    WAIC = WAIC,
    g = g,
    bayes_p_value = p_value,
    tidy_samples = tidy_samples))
}

