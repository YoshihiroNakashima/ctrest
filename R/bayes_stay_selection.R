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

  ni <- iter
  nt <- thin
  nc <- chains
  nb <- warmup

  vars_stay <- all.vars(formula_stay)
  response_stay <- vars_stay[1]
  predictors_stay <- vars_stay[-1]

  formula_stay_all <- list(0)
  if (all_comb == TRUE) {
    formula_stay_all <- full_terms(c(predictors_stay))
  } else {
    formula_stay_all[[1]] <- formula_stay
  }

  samples <- list(0)
  waic <- numeric(0)
  mean_stay_data <- list(0)
  for(k in 1:length(formula_stay_all)) {
    formula_stay <- formula_stay_all[[k]]

    model_frame_stay <- model.frame(formula_stay, stay_data)
    X_stay <- model.matrix(as.formula(formula_stay_all[[k]]), model_frame_stay)
    t <- model.response(model_frame_stay)
    censored <- stay_data[[col_name_cens]]

    if (!is.null(random_effect)) {
      levels <- unique(stay_data[[random_effect]])
      nLevels <- length(levels)
    } else {
      nLevels <- 0
    }

    N <- nrow(X_stay)
    nPreds <- ncol(X_stay)
    N_station <- nrow(station_data)

    names(t) <- NULL
    c_time <- t
    c_time[censored == 0] <- c_time[censored == 0] + 1  # 観測されたら観測値よりも小さい値
    t[censored == 1] <- NA
    N_stay <- length(t)


    data_pre <- list(t = t, censored = censored)
    cons_pre <- list(N_stay = N_stay, nPreds = nPreds, N_station = N_station, c_time = c_time, family = family)

    if (!is.null(random_effect)) {
      cons_pre$group <- as.numeric(factor(stay_data[[random_effect]]))
      cons_pre$nLevels <- nLevels
    } else {
      cons_pre$nLevels <- 0
    }
    if (nPreds > 1) {
      data_pre$X_stay <- X_stay
    }

    code <- nimbleCode({
      for(i in 1:N_stay) {
        censored[i] ~ dinterval(t[i], c_time[i])
        if(family == "exponential") {
          t[i] ~ dexp(scale = exp(log_scale[i]))
        }
        if(family == "gamma") {
          t[i] ~ dgamma(shape = theta, scale = exp(log_scale[i]))
        }
        if(family == "lognormal" | family == "log-normal") {
          t[i] ~ dlnorm(meanlog = log_scale[i], sdlog = theta)
        }
        if(family == "weibull") {
          t[i] ~ dweibull(shape = theta, scale = exp(log_scale[i]))
        }

        if (nPreds > 1) {
          if (nLevels == 0) {
            log_scale[i] <- inprod(beta[1:nPreds], X_stay[i, 1:nPreds])
          }
          if (nLevels > 0) {
            log_scale[i] <- inprod(beta[1:nPreds], X_stay[i, 1:nPreds]) + u[group[i]]
          }
        }
        if (nPreds == 1) {
          if (nLevels == 0) {
            log_scale[i] <- beta[1]
          }
          if (nLevels > 0) {
            log_scale[i] <- beta[1] + u[group[i]]
          }
        }
      }

      theta ~ dt(0, 1 / pow(2.5, -2), 1)

      for(j in 1:nPreds) {
        beta[j] ~ dnorm(0, sd = 100)
      }

      if (nLevels > 0) {
        for(k in 1:nLevels) {
          u[k] ~ dnorm(0, sd = sigma_u)
        }
        sigma_u ~ dt(0, 1 / pow(2.5, -2), 1)
      }

      if (nPreds == 1) {
        for(i in 1:N_station){
          if(family == "exponential") mean_stay[i] <- exp(beta[1])
          if(family == "gamma") mean_stay[i] <- theta * exp(beta[1])
          if(family == "lognormal") mean_stay[i] <- exp(beta[1] + theta ^2 / 2)
          if(family == "weibull") mean_stay[i] <-  lgamma(1 + 1 / theta) + exp(beta[1])
        }
      }

      if (nPreds > 1) {
        if(family == "exponential") {
          for(i in 1:N_station){
            mean_stay[i] <- exp(inprod(beta[1:nPreds], X_stay[i, 1:nPreds]))
          }
        }
        if(family == "gamma") {
          for(i in 1:N_station){
            mean_stay[i] <- theta * exp(inprod(beta[1:nPreds], X_stay[i, 1:nPreds]))
          }
        }
        if(family == "lognormal") {
          for(i in 1:N_station){
            mean_stay[i] <- exp(inprod(beta[1:nPreds], X_stay[i, 1:nPreds]) + theta ^2 / 2)
          }
        }
        if(family == "weibull") {
          for(i in 1:N_station){
            mean_stay[i] <- lgamma(1 + 1 / theta) + exp(inprod(beta[1:nPreds], X_stay[i, 1:nPreds]))
          }
        }
      }

    })

    init.t <- ifelse(censored == 0, NA, c_time + 1)

    inits <- list(beta = rep(0, nPreds), theta = 0.1, t = init.t)
    if (!is.null(random_effect)) {
      inits$u <- rep(0, nLevels)
      inits$sigma_u <- 1
    }
    cat("Compiling the model. This may take a moment...\n")

    model <-  suppressMessages(nimbleModel(code, constants = cons_pre, data = data_pre, inits = inits))

    params <- c("beta", "log_scale", "mean_stay", "theta")

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
      list(
        list(beta = rep(0, nPreds), theta = 1.0)
      )
    }

    per_chain_info <- lapply(1:nc, function(i) {
      list(
        seed = sample(1:9999, 1),
        inits = inits_f()
      )
    })
    params <- c("beta", "log_scale", "mean_stay", "theta")

    cat("Running MCMC sampling. Please wait...\n")

    clusterEvalQ(this_cluster, {
      library(nimble)
    })

    clusterExport(this_cluster,
                  varlist = c("run_MCMC_pre"),
                  envir = environment())

    chain_output <- parLapply(cl = this_cluster, X = per_chain_info,
                              fun = run_MCMC_pre,
                              data = data_pre, code = code,
                              constants = cons_pre, params = params,
                              ni = ni, nt = nt, nb = nb
    )
    stopCluster(this_cluster)


    pattern <- paste0("log_scale\\[", 1:N_stay, "\\]")

    log_scale.samp <- chain_output %>%
      map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
      do.call(rbind, .)
    if(family != "exponential") tau.samp <- do.call(rbind, chain_output)[, "theta"]

    loglfstay<-matrix(NA, nrow(log_scale.samp), N_stay)

    if(family == "exponential") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- 1 - pexp(c_time[j], rate = 1 /exp(log_scale.samp[ , j]), lower.tail = F, log.p = T)
        } else {
          loglfstay[, j] <- dexp(t[j], rate = 1 / exp(log_scale.samp[ , j]), log = T)
        }
      }
    }
    if(family == "gamma") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- 1 - pgamma(c_time[j], shape = tau.samp,lower.tail = F, scale = exp(log_scale.samp[ , j]), log.p = T)
        } else{
          loglfstay[, j] <- dgamma(t[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), log = T)
        }
      }
    }
    if(family == "lognormal" | family == "log-normal") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- 1 - plnorm(c_time[j], log_scale.samp[ , j], tau.samp,lower.tail = F, log.p = T)
        } else{
          loglfstay[, j] <- dlnorm(t[j], log_scale.samp[ , j], tau.samp, log = T)
        }
      }
    }
    if(family == "weibull") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- 1 - pweibull(c_time[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), lower.tail = F, log.p = T)
        } else{
          loglfstay[, j] <- dweibull(t[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), log = T)
        }
      }
    }

    lppd <- sum(log(colMeans(exp(loglfstay))))
    p.waic <- sum(apply(loglfstay, 2, var))
    waic[k] <- (- 2) * lppd + 2 * p.waic

    tidy_samples <- MCMCvis::MCMCchains(chain_output) %>%
      as_tibble() %>%
      pivot_longer(cols = everything(),
                   names_to = "parameter",
                   values_to = "value") %>%
      group_by(parameter) %>%
      mutate(iteration = row_number())


    mean_stay_data[[k]] <- tidy_samples %>%
      filter(grepl("^mean_stay\\[\\d+\\]$", parameter)) %>%
      mutate(index = as.numeric(str_extract(parameter, "\\d+"))) %>%
      filter(index <= 100) %>%
      group_by(parameter, index) %>%
      summarise(
        median = median(value),
        lower = quantile(value, 0.025),
        upper = quantile(value, 0.975),
        .groups = "drop"
      ) %>%
      arrange(index)

    samples[[k]] <- chain_output
  }

  model_summary <- list(0)
  model_summary$WAIC <- data.frame(Model = as.character(unlist(formula_stay_all)), WAIC = waic) %>%
    arrange(WAIC)
  model_summary$all_samples <- samples[[which(model_summary$WAIC[1, 1]  == unlist(formula_stay_all))]]
  model_summary$MeanStay_samples <- mean_stay_data[[which(model_summary$WAIC[1, 1]  == unlist(formula_stay_all))]]

  if (plot == TRUE) {
    g <- ggplot(model_summary$MeanStay_samples, aes(x = index, y = median)) +
      geom_point() +
      geom_segment(aes(x = index, xend = index, y = lower, yend = upper)) +
      labs(
        title = paste("formula =", formula_stay, "(Median and 95% Credible Interval for mean staying time)"),
        x = "Station",
        y = "MeanStayingTime"
      ) +
      ylim(0, NA)
    print(g)
  }
  return(model_summary$WAIC)
}

