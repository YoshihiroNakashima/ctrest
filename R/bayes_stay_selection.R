#' Model selection for staying time analysis using WAIC for REST/REST-RAD models based on MCMC samplings with nimble
#'
#' @param formula_stay A model formula for staying time within a focal area. For example, `Stay ~ 1 + x1`.
#' @param random_effect A character string specifying the random effect on the parameter of the staying time model. For example, `random_effect = "Station"`.
#' @param stay_data A data frame returned by the `format_stay` function, containing the processed staying time data.
#' @param col_name_cens A string specifying the column name in `stay_data` that indicates whether the observation is censored (1) or not (0).
#' @param family A character string specifying the probability distribution of staying time. Choose from `"exponential"`, `"gamma"`, `"weibull"`, or `"lognormal"`.
#' @param cores The number of CPU cores to use for parallel computation. Default is 3.
#' @param iter The total number of MCMC iterations per chain. Default is 5000
#' @param warmup The number of warm-up (burn-in) iterations per chain. Default is 1000.
#' @param thin The thinning interval for MCMC sampling. Default is 4 (no thinning).
#' @param chains The number of MCMC chains. Default is 3.
#' @param all_comb A logical value indicating whether to compare models with all possible combinations of covariates. If `FALSE`, only the specified model in `formula_stay` is evaluated. Default is FALSE.
#' @param target_species A character string specifying the species of interest. Only a single species can be specified.
#' @return A list of class \code{"ResultStay"}, which includes the following components:
#' \describe{
#'   \item{\code{WAIC}}{An object containing WAIC (Widely Applicable Information Criterion) results for model comparison.}
#'   \item{\code{Bayesian_p_value}}{Bayesian p-value for the best model, used to assess model fit.}
#'   \item{\code{summary_result}}{A data frame summarizing posterior estimates of the mean staying time.}
#'   \item{\code{samples}}{A \code{coda::mcmc.list} object containing MCMC samples for all parameters.}
#' }
#' The returned object has a custom print method that displays WAIC, Bayesian p-value, and summary of staying time estimates.
#' You can access full MCMC samples via \code{$samples}, and analyze convergence using the \code{MCMCvis} package.

#' @export
#' @import dplyr nimble
#' @importFrom tidyr unite extract
#' @importFrom purrr map
#' @importFrom stringr str_extract str_detect
#' @importFrom stats as.formula formula model.frame model.matrix sd var runif median quantile model.response rexp rnorm step dexp pexp dgamma pgamma dlnorm plnorm dweibull pweibull dnbinom
#' @export
#' @examples
#' stay_data <- format_stay(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_stay = "Stay",
#'   col_name_cens = "Cens"
#' )
#' bayes_stay_selection(formula_stay = Stay ~ 1 + x1,
#'                     stay_data = stay_data,
#'                     col_name_cens = "Cens",
#'                     family = "lognormal",
#'                     cores = 2,
#'                     iter = 5000,
#'                     warmup = 1000,
#'                     chains = 2,
#'                     thin = 4,
#'                     all_comb = FALSE,
#'                     target_species = "SP01")

bayes_stay_selection <- function(
    formula_stay = Stay ~ 1,
    random_effect = NULL,
    stay_data = stay_data,
    col_name_cens = "Cens",
    family = "lognormal",
    cores = 3,
    iter = 5000,
    warmup = 1000,
    chains = 3,
    thin = 4,
    all_comb = FALSE,
    target_species = NULL
) {

  if(family=="lognormal" || family=="log-normal" || family=="gamma" || family=="weibull" || family=="exponential"){
  } else{
    stop(paste0("Input family type(", family,") is incorrect."))
  }

  if (!inherits(formula_stay, "formula")) {
    stop("formula_stay must be a formula.")
  }

  if (!is.null(random_effect) && !(random_effect %in% colnames(stay_data))) {
    stop(paste0("random_effect column '", random_effect, "' not found in stay_data."))
  }

  if (!is.data.frame(stay_data)) {
    stop("stay_data must be a data frame.")
  }

  if (!(col_name_cens %in% colnames(stay_data))) {
    stop(paste0("Column '", col_name_cens, "' not found in stay_data."))
  }

  if (!is.numeric(cores) || cores < 1 || (cores %% 1 != 0)) {
    stop("cores must be a positive integer.")
  }

  if (!is.numeric(iter) || iter < 1 || (iter %% 1 != 0)) {
    stop("iter must be a positive integer.")
  }

  if (!is.null(warmup) && (!is.numeric(warmup) || warmup < 1 || (warmup %% 1 != 0))) {
    stop("warmup must be a positive integer or NULL.")
  }

  if (!is.numeric(chains) || chains < 1 || (chains %% 1 != 0)) {
    stop("chains must be a positive integer.")
  }

  if (!is.numeric(thin) || thin < 1 || (thin %% 1 != 0)) {
    stop("thin must be a positive integer.")
  }

  if (!is.null(target_species) && !all(target_species %in% stay_data$Species)) {
    stop("Some values in target_species are not found in stay_data$Species.")
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
  stay_data <- stay_data %>%
    filter(Species == target_species) %>%
    arrange(Station)

  station.id <- stay_data %>% pull(Station) %>% unique(.)


  if (family == "lognormal" ||
      family == "gamma" || family == "weibull" || family == "exponential") {

  } else {
    stop(paste0("Input family type(", family, ") is incorrect."))
  }

  # Extract variable names
  vars_stay <- all.vars(formula_stay)
  predictors_stay <- vars_stay[-1]

  formula_stay_all <- list(0)
  if (all_comb == TRUE) {
    formula_stay_all <- full_terms(c(predictors_stay))
  } else {
    formula_stay_all[[1]] <- formula_stay
  }

  mcmc_samples <- tidy_samples <- list(0)
  waic <- nPreds <- p_value <- numeric(0)

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
    N_station <- length(unique(stay_data$Station))

    names(t) <- NULL
    c_time <- t
    c_time[censored == 0] <- c_time[censored == 0] + 0.1
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
      cons_stay$group <- as.numeric(factor(stay_data[[random_effect]]))
      cons_stay$nLevels_stay <- nLevels_stay
    } else {
      cons_stay$nLevels_stay <- 0
    }
    if (nPreds_stay > 1) {
      data_stay$X_stay <- X_stay
    }

    nmcmc <- (iter - warmup) * chains / thin

    # MCMC2
    # MCMC
    code <- nimbleCode({
      for (i in 1:N_stay) {
        censored[i] ~ dinterval(t[i], c_time[i])

        if (family == "exponential") {
          t[i] ~ dexp(rate = 1/scale[i])
          pred_t[i] ~ dexp(rate = 1/scale[i])
          loglike_obs[i] <- (1 - step(censored[i] - 0.5)) * dexp(t[i], rate = 1/scale[i], log = 1) +
          step(censored[i] - 0.5) * log(1 - pexp(c_time[i], rate = 1/scale[i]))
          loglike_pred[i] <- dexp(pred_t[i], rate = 1/scale[i], log = 1)
        }

        if (family == "gamma") {
          t[i] ~ dgamma(shape = theta_stay, rate = exp(-log(scale[i])))
          pred_t[i] ~ dgamma(shape = theta_stay, rate = exp(-log(scale[i])))
          loglike_obs[i] <- (1 - step(censored[i] - 0.5)) * dgamma(t[i], shape = theta_stay, rate = exp(-log(scale[i])), log = 1) +
            step(censored[i] - 0.5) * log(1 - pgamma(c_time[i], shape = theta_stay, rate = exp(-log(scale[i]))))
          loglike_pred[i] <- dgamma(pred_t[i], shape = theta_stay, rate = exp(-log(scale[i])), log = 1)
        }

        if (family == "lognormal") {
          t[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay)
          pred_t[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay)
          loglike_obs[i] <- (1 - step(censored[i] - 0.5)) * dlnorm(t[i], meanlog = log(scale[i]), sdlog = theta_stay, log = 1) +
            step(censored[i] - 0.5) * log(1 - plnorm(c_time[i], meanlog = log(scale[i]), sdlog = theta_stay))
          loglike_pred[i] <- dlnorm(pred_t[i], meanlog = log(scale[i]), sdlog = theta_stay, log = 1)
          meanlog[i] <- log(scale[i])
        }

        if (family == "weibull") {
          t[i] ~ dweibull(shape = theta_stay, scale = scale[i])
          pred_t[i] ~ dweibull(shape = theta_stay, scale = scale[i])
          loglike_obs[i] <- (1 - step(censored[i] - 0.5)) * dweibull(t[i], shape = theta_stay, scale = scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pweibull(c_time[i], shape = theta_stay, scale = scale[i]))
          loglike_pred[i] <- dweibull(pred_t[i], shape = theta_stay, scale = scale[i], log = 1)
        }

        # Calculate log(scale[i])
        if (nPreds_stay == 1) {
          if(nLevels_stay > 0){
            log(scale[i]) <- beta_stay[1] +  random_effect[group[i]]
          } else{
            log(scale[i]) <- beta_stay[1]
          }

        } else {
          if(nLevels_stay > 0){
            log(scale[i]) <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect[group[i]]
          } else {
            log(scale[i]) <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay])
          }
        }
      }
      if (family == "lognormal") {
        sdlog <- theta_stay
      } else {
        shape <- theta_stay
      }

      # Priors
      for (j in 1:nPreds_stay) {
        beta_stay[j] ~ dnorm(0, sd = 100)
      }
      theta_stay ~ dgamma(0.01, 0.01)

      if (nLevels_stay > 0) {
        for (k in 1:nLevels_stay) {
          random_effect[k] ~ dnorm(0, sd = sigma_u)
        }
        sigma_u ~ T(dnorm(0, sd = 2.5), 0, 5)
      }
      # Expected value calculation
      if (nPreds_stay == 1) {

          if (nLevels_stay == 0) {
            if(family == "exponential") {
              mean_stay <- exp(beta_stay[1])
            }
            if(family == "gamma") {
              mean_stay <- theta_stay * exp(beta_stay[1])
            }
            if(family == "lognormal") {
              mean_stay <- exp(beta_stay[1] + theta_stay ^ 2 / 2)
            }
            if(family == "weibull") {
              mean_stay <-  lgamma(1 + 1 / theta_stay) + exp(beta_stay[1])
            }
          }
          if (nLevels_stay > 0) {
            if(family == "exponential") {
              mean_stay <- exp(beta_stay[1])
              # for(i in 1:N_station) {
              #   mean_stay[i] <- exp(beta_stay[1] + random_effect[group[i]])
              # }
            }
            if(family == "gamma") {
              mean_stay <- theta_stay * exp(beta_stay[1])
              # for(i in 1:N_station) {
              #   mean_stay[i] <- theta_stay * exp(beta_stay[1] + random_effect[group[i]])
              # }
            }

            if(family == "lognormal") {
              mean_stay <- exp(beta_stay[1] + theta_stay ^ 2 / 2)
              # for(i in 1:N_station) {
              #   mean_stay[i] <- exp(beta_stay[1] + random_effect[group[i]] + theta_stay ^ 2 / 2)
              # }
            }
            if(family == "weibull") {
              mean_stay <-  lgamma(1 + 1 / theta_stay) + exp(beta_stay[1])
              # for(i in 1:N_station) {
              #   mean_stay[i] <-  lgamma(1 + 1 / theta_stay) + exp(beta_stay[1] + random_effect[group[i]])
              # }
            }
          }
      }
      if (nPreds_stay > 1) {
        if (nLevels_stay == 0) {
          if(family == "exponential") {
            for(i in 1:N_station){
              mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
            }
          }
          if(family == "gamma") {
            for(i in 1:N_station){
              mean_stay[i] <- theta_stay * exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
            }
          }
          if(family == "lognormal") {
            for(i in 1:N_station){
              mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + theta_stay ^ 2 / 2)
            }
          }
          if(family == "weibull") {
            for(i in 1:N_station){
              mean_stay[i] <- lgamma(1 + 1 / theta_stay) + exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
            }
          }
        }
        if (nLevels_stay > 0) {
          if(family == "exponential") {
            for(i in 1:N_station){
              mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
            }
          }
          if(family == "gamma") {
            for(i in 1:N_station){
              # mean_stay[i] <- theta_stay * exp(inprod(beta_stay[1:nPreds_stay] + random_effect[group[i]], X_stay[i, 1:nPreds_stay]))
              mean_stay[i] <- theta_stay * exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
            }
          }
          if(family == "lognormal") {
            for(i in 1:N_station){
              mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + theta_stay ^ 2 / 2)
            }
          }
          if(family == "weibull") {
            for(i in 1:N_station){
              mean_stay[i] <- lgamma(1 + 1 / theta_stay) + exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
            }
          }
        }
      }

      sum_loglike_obs <- sum(loglike_obs[1:N_stay])
      sum_loglike_pred <- sum(loglike_pred[1:N_stay])

      # Deviance Calculation
      deviance_obs <- -2 * sum_loglike_obs
      deviance_pred <- -2 * sum_loglike_pred
    })


    cat("Compiling the model. This may take a moment...\n")

    this_cluster <- makeCluster(chains)
    # on.exit(stopCluster(this_cluster), add = TRUE)

    run_MCMC_pre <- function(info, data, constants, code, params, iter, thin, warmup) {
      myModel <- nimbleModel(code = code,
                             data = data,
                             constants = constants,
                             inits = info$inits)
      CmyModel <- compileNimble(myModel)
      myMCMC <- buildMCMC(CmyModel, monitors = params)
      CmyMCMC <- compileNimble(myMCMC)
      results <- runMCMC(
        CmyMCMC,
        niter = iter,
        nburnin = warmup,
        thin = thin,
        nchains = 1,
        setSeed = info$seed,
        samplesAsCodaMCMC = TRUE
      )
      results
    }

    inits_f <- function() {
      common_inits <- list(
        beta_stay = runif(nPreds_stay, -0.1, 0.1),
        t = ifelse(censored == 0, NA, c_time + runif(N_stay, 0.1, 1.0)),
        theta_stay = runif(1, 0.2, 2.0),
        scale = runif(N_stay, 0.2, 2.0)
      )

      if (!is.null(random_effect)) {
        c(
          common_inits,
          list(
            random_effect = runif(nLevels_stay, -1, 1),
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

    per_chain_info <- lapply(1:chains, function(i) {
      list(
        seed = sample(1:9999, 1),
        inits = inits_f()
      )
    })
    if(family == "exponential") prms <- c("scale", "mean_stay")
    if(family == "gamma" | family == "weibull") prms <- c("scale", "shape", "mean_stay")
    if(family == "lognormal") prms <- c("meanlog", "sdlog", "mean_stay")
    if(family == "exponential") prms <- c("scale", "shape", "mean_stay")

    params <- c(prms, "beta_stay", "sum_loglike_obs", "sum_loglike_pred", "deviance_obs", "deviance_pred", "loglike_obs")

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
                              iter = iter, thin = thin, warmup = warmup
    )
    stopCluster(this_cluster)
    cat("Estimation is finished!\n")

    # Summarize results -------------------------------------------------------

    mcmc_samples[[k]] <- MCMCvis::MCMCchains(chain_output,
                                             mcmc.list = TRUE,
                                             params = prms)

    samples <- MCMCvis::MCMCchains(chain_output)
    tidy_samples[[k]] <- samples %>%
      as_tibble() %>%
      pivot_longer(cols = everything(),
                   names_to = "parameter",
                   values_to = "value") %>%
      group_by(parameter) %>%
      mutate(iteration = row_number()) %>%
      ungroup()

    # Bayesian p-value
    chi2_obs_total_samples <- samples[ , grep("sum_loglike_obs", colnames(mcmc_samples))]
    chi2_pred_total_samples <- samples[ , grep("sum_loglike_pred", colnames(mcmc_samples))]
    deviance_obs_samples <- samples[, "deviance_obs"]
    deviance_pred_samples <- samples[, "deviance_pred"]

    p_value[k] <- mean(deviance_pred_samples > deviance_obs_samples)

    # WAIC --------------------------------------------------------------------

    loglfstay <- MCMCvis::MCMCchains(chain_output,
                                     mcmc.list = FALSE,
                                     params = c("loglike_obs"))

    lppd <- sum(log(colMeans(exp(loglfstay))))
    p.waic <- sum(apply(loglfstay, 2, var))
    waic[k] <- (- 2) * lppd + 2 * p.waic

    nPreds[k] <- nPreds_stay
  }

  if(is.null(random_effect)) random_effect <- "NULL"
  WAIC <- data.frame(Model = as.character(unlist(formula_stay_all)), Family = rep(family, k),
                     Random_effect = random_effect, WAIC = waic) %>%
    arrange(WAIC)

  best.model <- which(WAIC[1, 1]  == unlist(formula_stay_all))

  mcmc_samples <- mcmc_samples[[best.model]]
  tidy_samples_best  <- tidy_samples[[best.model]]
  p_value <- p_value[best.model]

  mcmc_samples_mean <- MCMCvis::MCMCchains(chain_output, mcmc.list = TRUE, params = c("mean_stay"))
  summary_mean_temp <- MCMCsummary(mcmc_samples_mean, round = 2) %>%
    tibble::rownames_to_column(., var = "Variable") %>%
    tibble(.) %>%
    rename(lower = `2.5%`, median = `50%`, upper = `97.5%`)

  if(nPreds_stay > 1) {
    summary_mean <- tibble(Station = station.id) %>%
      mutate(summary_mean_temp)
  } else {
    summary_mean <- mutate(summary_mean_temp)
  }
  stay_result <- list(
    WAIC = WAIC,
    Bayesian_p_value = p_value,
    summary_result = summary_mean,
    samples = mcmc_samples
  )
  class(stay_result) <- "ResultStay"

  stay_result
}

