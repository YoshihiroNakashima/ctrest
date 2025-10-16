#' Bayesian parameter estimation for the REST/RAD-REST model via MCMC sampling using `nimble`.
#'
#' @param formula_stay A model formula for staying times within the focal area. For example, `Stay ~ 1 + x1`.
#' @param formula_density A model formula for animal density. For example, `~ 1 + x1` (note that the left-hand side must be omitted).
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
#' @param random_effect A specification of the random effect structure on staying time (e.g., `~ (1 | Station)`).
#' @param activity_data A data frame containing a "time" column, representing detection times transformed into radians. Typically, this is the output of the format_activity function.
#' @param activity_estimation Method used to estimate activity patterns. Choose `"kernel"` for fixed kernel density estimation (Rowcliffe et al. 2014) or `"mixture"` for nonparametric Bayesian estimation using von Mises mixture models (Nakashima et al. 2025).
#' @param bw_adj Bandwidth adjustment parameter for kernel density estimation. The default is 1.0. See Rowcliffe et al. (2014) for details.
#' @param C The maximum number of von Mises components used for the mixture model. Required only if `activity_estimation = "mixture"`. The default is 10.
#' @param stay_family The probability distribution used for modeling staying times. This should be selected based on the output of the `bayes_stay_selection` function.
#' @param focal_area The area of the focal area, in square meters.
#' @param cores The number of CPU cores to use for parallel computation. Default is 3.
#' @param iter The total number of MCMC iterations per chain. Default is 5000
#' @param warmup The number of warm-up (burn-in) iterations per chain. Default is 1000.
#' @param chains The number of MCMC chains. Default is 2.
#' @param thin The thinning interval for MCMC sampling. Default is 1 (no thinning).
#' @param all_comb Logical. If TRUE, all combinations of covariates in the density model will be evaluated and compared. If FALSE, only the specified model in `formula_density` will be run.
#' @param model The model to be used. Choose either `"REST"` or `"RAD-REST"`.
#' @param target_species A character string specifying the species to be analyzed. Only a single species can be specified.
#' @return A list of class \code{"ResultDensity"}, which includes the following components:
#' \describe{
#'   \item{\code{WAIC}}{An object containing WAIC (Widely Applicable Information Criterion) results.}
#'   \item{\code{summary_result}}{A data frame summarizing the posterior estimates of the mean staying time.}
#'   \item{\code{samples}}{A \code{coda::mcmc.list} object containing MCMC samples for all model parameters.}
#' }
#' The returned object has a custom print method that displays the WAIC and summary of staying time estimates.
#' You can access full MCMC samples via \code{$samples}, and analyze convergence using the \code{MCMCvis} package.

#' @export
#' @import nimble activity parallel MCMCvis
#' @importFrom stats as.formula formula model.frame model.matrix sd var runif median quantile model.response rexp rnorm step dexp pexp dgamma pgamma dlnorm plnorm dweibull pweibull dnbinom
#' @importFrom dplyr select
#' @importFrom extraDistr ddirmnom
#' @examples
#' station_data_rest <- format_station_data(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_y = "y",
#'   model = "REST"
#' )
#' station_effort_REST <- add_effort(
#'   detection_data = detection_data,
#'   station_data_formatted = station_data_rest,
#'   col_name_station = "Station",
#'   col_name_term = "Term",
#'   col_name_datetime = "DateTime",
#'   plot = TRUE
#' )
#' stay_data <- format_stay(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_stay = "Stay",
#'   col_name_cens = "Cens"
#' )
#' activity_data <- format_activity(
#'   detection_data = detection_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_datetime = "DateTime",
#'   indep_time = 30
#' )
#' rest_model <- bayes_rest(
#'   formula_stay = Stay ~ 1 + x1,
#'   formula_density = ~ 1,
#'   station_effort_data = station_effort_REST,
#'   stay_data = stay_data,
#'   activity_data = activity_data,
#'   activity_estimation = "kernel",
#'   bw_adj = 1.5,
#'   stay_family = "lognormal",
#'   focal_area = 1.96,
#'   cores = 2,
#'   iter = 3000,
#'   warmup = 1000,
#'   chains = 2,
#'   thin = 1,
#'   model = "REST",
#'   all_comb = FALSE,
#'   target_species = "SP01"
#' )
#'
bayes_rest <- function(formula_stay,
                       formula_density,
                       station_effort_data,
                       stay_data,
                       random_effect = NULL,
                       activity_data,
                       activity_estimation = "kernel",
                       C = 10,
                       bw_adj = 1.0,
                       stay_family = "lognormal",
                       focal_area,
                       cores = 3,
                       iter = 5000,
                       warmup = 1000,
                       chains = 3,
                       thin = 2,
                       model = "REST",
                       all_comb = FALSE,
                       target_species
){
  ni <- iter
  nt <- thin
  nc <- chains
  nb <- warmup
  nmcmc <- (iter - warmup) * chains / thin

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
             r <- paste("~ 1 + ", r, "1", sep = "")
             r <- strsplit(r, " \\+ 1$")[[1]][1]
           }
    )
  }

  # Check -----------------------------------------------------------------

  # Check stay_family
  if (!stay_family %in% c("lognormal", "gamma", "weibull", "exponential")) {
    stop(paste0("Input stay_family type (", stay_family, ") is incorrect."))
  }

  # Check single species
  if (length(target_species) > 1) {
    stop("The 'bayes_rest' function does not support multiple species. Please run the model separately for each species or use 'bayes_rest_multi' for multispecies analysis.")
  }

  # Check activity_estimation
  if (!activity_estimation %in% c("kernel", "mixture")) {
    stop("`activity_estimation` must be either 'kernel' or 'mixture'.")
  }

  # Check model
  if (!model %in% c("REST", "RAD-REST")) {
    stop("`model` must be either 'REST' or 'RAD-REST'.")
  }

  # Check target_species
  if (!is.character(target_species) || length(target_species) != 1) {
    stop("`target_species` must be a single character string.")
  }

  # Check focal_area
  if (!is.numeric(focal_area) || length(focal_area) != 1 || focal_area <= 0) {
    stop("`focal_area` must be a positive numeric value representing area in square meters.")
  }

  # Check cores
  if (!is.numeric(cores) || cores < 1) {
    stop("`cores` must be a positive integer.")
  }

  # Check iter
  if (!is.numeric(iter) || iter <= 0 || length(iter) != 1) {
    stop("`iter` must be a positive integer.")
  }

  # Check warmup
  if (!is.null(warmup) && (!is.numeric(warmup) || warmup <= 0 || length(warmup) != 1)) {
    stop("`warmup` must be NULL or a positive integer.")
  }

  # Check chains
  if (!is.numeric(chains) || chains < 1 || length(chains) != 1) {
    stop("`chains` must be a positive integer.")
  }

  # Check thin
  if (!is.numeric(thin) || thin < 1 || length(thin) != 1) {
    stop("`thin` must be a positive integer.")
  }

  # Check all_comb
  if (!is.logical(all_comb) || length(all_comb) != 1) {
    stop("`all_comb` must be either TRUE or FALSE.")
  }

  # Check formula_stay
  if (!inherits(formula_stay, "formula")) {
    stop("`formula_stay` must be a valid model formula.")
  }

  # Check formula_density
  if (!inherits(formula_density, "formula")) {
    stop("`formula_density` must be a valid model formula (e.g., `~ 1 + x1`).")
  }

  # Check stay_data columns
  required_cols_stay <- c("Station", "Species", "Stay", "Cens")
  if (!all(required_cols_stay %in% colnames(stay_data))) {
    stop("`stay_data` must contain the following columns: 'Station', 'Species', 'Stay', and 'Cens'.")
  }

  # # Check activity_data
  # if (!is.numeric(activity_data)) {
  #   stop("`activity_data` must be a numeric vector of detection times in radians.")
  # }


  # Prepare dataset ---------------------------------------------------------------

  act_data <- activity_data %>%
    filter(Species == target_species) %>%
    pull(time)
  if(activity_estimation == "kernel") {
    model_act <- fitact(act_data, bw = bw_adj*bwcalc(act_data, K = 3), reps=1)
    activity_proportion <- model_act@act
  } else {
    dens.x <- seq(0, 2 * pi, 0.02)
  }
  N_act <- activity_data %>% pull(time) %>% length(.)


  station_effort_data <- station_effort_data %>%
    filter(Species == target_species)
  stay_data_join <- stay_data %>%
    filter(Species == target_species) %>%
    left_join(station_effort_data, by = intersect(names(stay_data), names(station_effort_data))) %>%
    filter(!is.na(Stay))
  N_station <- nrow(station_effort_data)
  station.id <- unique(station_effort_data$Station)
  stay <- stay_data_join %>% pull(Stay)
  censored <- stay_data_join %>% pull(Cens)

  c_time <- stay
  c_time[censored == 0] <- c_time[censored == 0] + 1
  stay[censored == 1] <- NA
  N_stay <- length(stay)

  S <- focal_area * 10 ^ -6
  N_period <- station_effort_data$Effort[1:N_station] * 60 * 60 * 24

  # Extract variable names
  vars_stay <- all.vars(formula_stay)
  response_stay <- vars_stay[1]
  predictors_stay <- vars_stay[-1]

  model_frame_stay <- model.frame(formula_stay, stay_data_join)
  X_stay <- model.matrix(as.formula(formula_stay), model_frame_stay)
  nPreds_stay <- ncol(X_stay)

  # N of random effects
  if (!is.null(random_effect)) {
    levels <- unique(stay_data_join[[random_effect]])
    nLevels_stay  <- length(levels)
  } else {
    nLevels_stay  <- 0
  }

  predictors_density <- all.vars(formula_density)

  formula_density_all <- list(0)
  if (all_comb == TRUE) {
    formula_density_all <- full_terms(c(predictors_density))
  } else {
    formula_density_all[[1]] <- formula_density
  }

  # Activity proportion estimation ("mixture") ------------------------------


  if(activity_estimation == "mixture") {
    dens.x <- seq(0, 2 * pi, 0.02)
    ndens <- length(dens.x)
    N <- length(act_data)
    constants <- list(N = N, C = C, dens.x = dens.x, ndens = ndens)
    data <- list(act_data = act_data)

    code <- nimbleCode({
      for(k in 1:(C-1)) {
        v[k] ~ dbeta(1, alpha)
      }
      alpha ~ dgamma(1, 1)
      w[1:C] <- stick_breaking(v[1:(C-1)])
      for(k in 1:C) {
        mu_mix[k] ~ dunif(0, 2 * 3.141592654)
        kappa_mix[k] ~ dgamma(1, 0.01)
      }
      for(n in 1:N) {
        group[n] ~ dcat(w[1:C])
        act_data[n] ~ dvonMises(mu_mix[group[n]], kappa_mix[group[n]])

        act_data_pred[n] ~ dvonMises(mu_mix[group[n]], kappa_mix[group[n]])
        loglike_obs_act[n] <- dvonMises(act_data[n], mu_mix[group[n]], kappa_mix[group[n]], log = 1)
        loglike_pred_act[n] <- dvonMises(act_data_pred[n], mu_mix[group[n]], kappa_mix[group[n]], log = 1)
      }
      for (j in 1:ndens) {
        for (i in 1:C) {
          dens.cpt[i, j] <- w[i] * dvonMises(dens.x[j] , mu_mix[i], kappa_mix[i], log = 0)
        }
        activity_density[j] <- sum(dens.cpt[1:C, j])
      }
      activity_proportion <- 1.0 / (2 * 3.141592654 * max(activity_density[1:ndens]));
    })

    inits_f <- function() {
      list(
        mu_mix = runif(constants$C, 0, 2 * pi),
        kappa_mix = rgamma(constants$C, 1, 0.01),
        group = sample(1:constants$C, size = constants$N, replace = TRUE),
        v = rbeta(constants$C-1, 1, 1),
        alpha = 1
      )
    }

    # Define the von Mises distribution function using nimbleFunction
    dvonMises <- nimbleFunction(
      run = function(x = double(0), kappa = double(0), mu = double(0), log = integer(0)) {
        returnType(double(0))
        ccrit <- 1E-6
        s <- 1
        i <- 1
        inc <- 1
        x_2i <- 0
        satisfied <- FALSE
        while(!satisfied) {
          x_2i <- kappa / (2 * i)
          inc <- inc * x_2i * x_2i
          s <- s + inc
          i <- i + 1
          satisfied <- inc < ccrit
        }
        prob <- exp(kappa * cos(x - mu)) / (2 * pi * s)
        if (log) {
          return(log(prob))
        } else {
          return(prob)
        }
      }
    )

    rvonMises <- nimbleFunction(
      run = function(n = integer(0), kappa = double(0), mu = double(0)) {
        returnType(double(0))
        return(0)
      }
    )

    suppressMessages(registerDistributions(list(
      dvonMises = list(
        BUGSdist = "dvonMises(kappa, mu)",
        types = c('value = double(0)', 'kappa = double(0)', 'mu = double(0)'),
        pqAvail = FALSE
      )
    )))

    run_MCMC_vonMises <- function(info, data, constants, code, params, ni, nt, nb) {
      myModel <- nimbleModel(code = code,
                             data = data,
                             constants = constants,
                             inits = info$inits)

      CmyModel <- compileNimble(myModel)

      configModel <- configureMCMC(myModel, monitors = params)
      myMCMC <- buildMCMC(configModel, monitors = params)
      CmyMCMC <- compileNimble(myMCMC)

      results <- runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
    }

    per_chain_info <- lapply(1:nc, function(i) {
      list(
        seed = sample(1:9999, 1),
        inits = inits_f()
      )
    })

    params <- c("activity_density", "activity_proportion", "mu_mix", "kappa_mix", "w", "loglike_obs_act")
    cat("Running MCMC sampling. Please wait...\n")

    this_cluster <- makeCluster(nc)
    clusterEvalQ(this_cluster, {
      library(nimble)
    })
    clusterExport(this_cluster,
                  c("dvonMises", "rvonMises", "registerDistributions", "run_MCMC_vonMises"),
                  envir = environment())

    actv_chain_output <- parLapply(cl = this_cluster, X = per_chain_info,
                              fun = run_MCMC_vonMises,
                              data = data, code = code,
                              constants = constants, params = params,
                              ni = ni, nt = nt, nb = nb)
    stopCluster(this_cluster)

    actv_out_trace <- actv_chain_output %>%
      map(~ .[, grep(paste("activity_proportion", collapse = "|"), colnames(.))])

    summary <- MCMCsummary(actv_chain_output, round = 3)

    activity_density_estimates <- summary %>%
      as.data.frame() %>%
      tibble::rownames_to_column("variable") %>%
      filter(stringr::str_starts(variable, "activity_density")) %>%
      mutate(x = dens.x)

    loglact <- MCMCchains(actv_chain_output, params = c("loglike_obs_act"))
  }

  # Original-REST model ----------------------------------------------------------

  if(model == "REST") {

    tidy_samples <- mcmc_samples <- list(0)
    waic <- numeric(0)

    for(k in 1:length(formula_density_all)) {

      formula_density <- formula_density_all[[k]]
      model_frame_density <- model.frame(formula_density, station_effort_data)
      X_density <- model.matrix(as.formula(formula_density), model_frame_density)
      nPreds_density <- ncol(X_density)

      # prepare dataset for the original REST

      y <- station_effort_data %>% pull(Y)

      data_REST <- list(y = y,
                        X_density = X_density,
                        stay = stay,
                        censored = censored)

      cons_REST <- list(N_station = N_station,
                        S = S,
                        N_period = N_period,
                        nPreds_stay = nPreds_stay,
                        N_stay = N_stay,
                        c_time = c_time,
                        nPreds_density = nPreds_density,
                        stay_family = stay_family,
                        activity_estimation = activity_estimation)

      # Random effects
      if (!is.null(random_effect)) {
        cons_REST$group_stay <- as.numeric(factor(stay_data_join[[random_effect]]))
        cons_REST$nLevels_stay  <- nLevels_stay
      } else {
        cons_REST$nLevels_stay  <- 0
      }
      if(activity_estimation == "kernel") {
        cons_REST$activity_proportion <- activity_proportion
      }
      if (nPreds_stay > 1) {
        data_REST$X_stay <- X_stay
      }

      Model_REST <- nimbleCode(
        {
          for(i in 1:N_stay) {
            censored[i] ~ dinterval(stay[i], c_time[i])

            if (stay_family == "exponential") {
              stay[i] ~ dexp(rate = 1/scale[i])
              pred_t[i] ~ dexp(rate = 1/scale[i])
              loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dexp(stay[i], rate = 1/scale[i], log = 1) +
                step(censored[i] - 0.5) * log(1 - pexp(c_time[i], rate = 1/scale[i]))
              loglike_pred_stay[i] <- dexp(pred_t[i], rate = 1/scale[i], log = 1)
            }

            if (stay_family == "gamma") {
              stay[i] ~ dgamma(shape = theta_stay, rate = exp(-log(scale[i])))
              pred_t[i] ~ dgamma(shape = theta_stay, rate = exp(-log(scale[i])))
              loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dgamma(stay[i], shape = theta_stay, rate = exp(-log(scale[i])), log = 1) +
                step(censored[i] - 0.5) * log(1 - pgamma(c_time[i], shape = theta_stay, rate = exp(-log(scale[i]))))
              loglike_pred_stay[i] <- dgamma(pred_t[i], shape = theta_stay, rate = exp(-log(scale[i])), log = 1)
            }

            if (stay_family == "lognormal") {
              stay[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay)
              pred_t[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay)
              loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dlnorm(stay[i], meanlog = log(scale[i]), sdlog = theta_stay, log = 1) +
                step(censored[i] - 0.5) * log(1 - plnorm(c_time[i], meanlog = log(scale[i]), sdlog = theta_stay))
              loglike_pred_stay[i] <- dlnorm(pred_t[i], meanlog = log(scale[i]), sdlog = theta_stay, log = 1)
              meanlog[i] <- log(scale[i])
            }

            if (stay_family == "weibull") {
              stay[i] ~ dweibull(shape = theta_stay, scale = scale[i])
              pred_t[i] ~ dweibull(shape = theta_stay, scale = scale[i])
              loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dweibull(stay[i], shape = theta_stay, scale = scale[i], log = 1) +
                step(censored[i] - 0.5) * log(1 - pweibull(c_time[i], shape = theta_stay, scale = scale[i]))
              loglike_pred_stay[i] <- dweibull(pred_t[i], shape = theta_stay, scale = scale[i], log = 1)
            }

            if (nPreds_stay > 1) {
              if (nLevels_stay  == 0) {
                log(scale[i]) <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay])
              } else {
                log(scale[i]) <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]]
              }
            } else {
              if (nLevels_stay  == 0) {
                log(scale[i]) <- beta_stay[1]
              } else {
                log(scale[i]) <- beta_stay[1] + random_effect_stay[group_stay[i]]
              }
            }
          }

          theta_stay ~ dgamma(1, 1)

          if (stay_family == "lognormal") {
            sdlog <- theta_stay
          } else {
            shape <- theta_stay
          }

          for(j in 1:nPreds_stay) {
            beta_stay[j] ~ T(dnorm(0, 100), -10, 10)
          }

          if (nLevels_stay  > 0) {
            for(k in 1:nLevels_stay ) {
              random_effect_stay[k] ~ dnorm(0, sd = sigma_stay)
            }
            sigma_stay ~ T(dnorm(0, 100), 0, 10)
          }

          # Expected value calculation
          if (nPreds_stay == 1) {
            if (nLevels_stay == 0) {
              if(stay_family == "exponential") {
                mean_stay <- exp(beta_stay[1])
              }
              if(stay_family == "gamma") {
                mean_stay <- theta_stay * exp(beta_stay[1])
              }
              if(stay_family == "lognormal") {
                mean_stay <- exp(beta_stay[1] + theta_stay ^ 2 / 2)
              }
              if(stay_family == "weibull") {
                mean_stay <-  lgamma(1 + 1 / theta_stay) + exp(beta_stay[1])
              }
            }
            if (nLevels_stay > 0) {
              if(stay_family == "exponential") {
                mean_stay <- exp(beta_stay[1])
                # for(i in 1:N_station) {
                #   mean_stay[i] <- exp(beta_stay[1] + random_effect[group[i]])
                # }
              }
              if(stay_family == "gamma") {
                mean_stay <- theta_stay * exp(beta_stay[1])
                # for(i in 1:N_station) {
                #   mean_stay[i] <- theta_stay * exp(beta_stay[1] + random_effect[group[i]])
                # }
              }

              if(stay_family == "lognormal") {
                mean_stay <- exp(beta_stay[1] + theta_stay ^ 2 / 2)
                # for(i in 1:N_station) {
                #   mean_stay[i] <- exp(beta_stay[1] + random_effect[group[i]] + theta_stay ^ 2 / 2)
                # }
              }
              if(stay_family == "weibull") {
                mean_stay <-  lgamma(1 + 1 / theta_stay) + exp(beta_stay[1])
                # for(i in 1:N_station) {
                #   mean_stay[i] <-  lgamma(1 + 1 / theta_stay) + exp(beta_stay[1] + random_effect[group[i]])
                # }
              }
            }
          }
          if (nPreds_stay > 1) {
            if (nLevels_stay == 0) {
              if(stay_family == "exponential") {
                for(i in 1:N_station){
                  mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
                }
              }
              if(stay_family == "gamma") {
                for(i in 1:N_station){
                  mean_stay[i] <- theta_stay * exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
                }
              }
              if(stay_family == "lognormal") {
                for(i in 1:N_station){
                  mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + theta_stay ^ 2 / 2)
                }
              }
              if(stay_family == "weibull") {
                for(i in 1:N_station){
                  mean_stay[i] <- lgamma(1 + 1 / theta_stay) + exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
                }
              }
            }
            if (nLevels_stay > 0) {
              if(stay_family == "exponential") {
                for(i in 1:N_station){
                  mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay] + random_effect[group[i]], X_stay[i, 1:nPreds_stay]))
                }
              }
              if(stay_family == "gamma") {
                for(i in 1:N_station){
                  mean_stay[i] <- theta_stay * exp(inprod(beta_stay[1:nPreds_stay] + random_effect[group[i]], X_stay[i, 1:nPreds_stay]))
                }
              }
              if(stay_family == "lognormal") {
                for(i in 1:N_station){
                  mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay] + random_effect[group[i]], X_stay[i, 1:nPreds_stay]) + theta_stay ^ 2 / 2)
                }
              }
              if(stay_family == "weibull") {
                for(i in 1:N_station){
                  mean_stay[i] <- lgamma(1 + 1 / theta_stay) + exp(inprod(beta_stay[1:nPreds_stay] + random_effect[group[i]], X_stay[i, 1:nPreds_stay]))
                }
              }
            }
          }
          shape <- theta_stay

          # Posterior predictive check
          for (i in 1:N_station) {
            y[i] ~ dnbinom(size = size, prob = p[i])
            p[i] <- size / (size + mu[i])

            y_rep[i] ~ dnbinom(size = size, prob = p[i])

            loglike_obs_y[i] <- dnbinom(y[i], size, p[i], log = 1)
            loglike_pred_y[i] <- dnbinom(y_rep[i], size, p[i], log = 1)
          }
          size ~ dgamma(1, 1)


          # REST formula
          if(nPreds_density == 1) {
            if(nPreds_stay == 1) {
              for(i in 1:N_station) {
                log(mu[i]) <- log(density) + log(S) + log(N_period[i]) - log(mean_stay) + log(activity_proportion)
              }
              log(density) <- beta_density
              beta_density ~ dnorm(0, sd = 100)
            }
            if(nPreds_stay > 1) {
              for(i in 1:N_station) {
                log(mu[i]) <- log(density) + log(S) + log(N_period[i]) - log(mean_stay[i]) + log(activity_proportion)
              }
              log(density) <- beta_density
              beta_density ~ dnorm(0, sd = 100)
            }
          }
          if(nPreds_density > 1) {
            if(nPreds_stay == 1) {
              for(i in 1:N_station) {
                log(mu[i]) <- log_local_density[i] + log(S) + log(N_period[i]) - log(mean_stay) + log(activity_proportion)
                log_local_density[i] <- inprod(beta_density[1:nPreds_density], X_density[i, 1:nPreds_density])
              }
              for(j in 1:nPreds_density) {
                beta_density[j] ~ dnorm(0, sd = 100)
              }
              for(i in 1:N_station) {
                density[i] <- exp(inprod(beta_density[1:nPreds_density], X_density[i, 1:nPreds_density]))
              }
            }
            if(nPreds_stay > 1) {
              for(i in 1:N_station) {
                log(mu[i]) <- log_local_density[i] + log(S) + log(N_period[i]) - log(mean_stay[i]) + log(activity_proportion)
                log_local_density[i] <- inprod(beta_density[1:nPreds_density], X_density[i, 1:nPreds_density])
              }
              for(j in 1:nPreds_density) {
                beta_density[j] ~ dnorm(0, sd = 100)
              }
              for(i in 1:N_station) {
                density[i] <- exp(inprod(beta_density[1:nPreds_density], X_density[i, 1:nPreds_density]))
              }
            }
          }
        }#end
      )


      inits_f <- function() {
        common_inits <- list(
          beta_stay = runif(nPreds_stay, -1, 1),
          stay = ifelse(censored == 0, NA, c_time + runif(N_stay, 0.1, 1.0)),
          theta_stay = runif(1, 0.5, 1.5),
          scale = runif(N_stay, 0, 2),
          beta_density = rnorm(nPreds_density, 0, 1)
        )

        if (!is.null(random_effect)) {
          return(c(
            common_inits,
            list(
              random_effect_stay = runif(nLevels_stay , -1, 1),
              sigma_stay = runif(1, 0.5, 2.5)
            )
          ))
        }
        return(common_inits)
      }

      dvonMises <- nimbleFunction(
        run = function(x = double(0), kappa = double(0), mu = double(0), log = integer(0)) {
          returnType(double(0))
          ccrit <- 1E-6
          s <- 1
          i <- 1
          inc <- 1
          x_2i <- 0
          satisfied <- FALSE
          while(!satisfied) {
            x_2i <- kappa / (2 * i)
            inc <- inc * x_2i * x_2i
            s <- s + inc
            i <- i + 1
            satisfied <- inc < ccrit
          }
          prob <- exp(kappa * cos(x - mu)) / (2 * pi * s)
          if (log) {
            return(log(prob))
          } else {
            return(prob)
          }
        }
      )

      rvonMises <- nimbleFunction(
        run = function(n = integer(0), kappa = double(0), mu = double(0)) {
          returnType(double(0))
          return(0)
        }
      )

      suppressMessages(registerDistributions(list(
        dvonMises = list(
          BUGSdist = "dvonMises(kappa, mu)",
          types = c('value = double(0)', 'kappa = double(0)', 'mu = double(0)'),
          pqAvail = FALSE
        )
      )))
      if(activity_estimation != "mixture") {
        run_MCMC_REST <- function(info, data, constants, code, params, ni, nt, nb) {
          myModel <- nimbleModel(code = code,
                                 data = data,
                                 constants = constants,
                                 inits = info$inits)

          CmyModel <- compileNimble(myModel)

          configModel <- configureMCMC(myModel, monitors = params)
          myMCMC <- buildMCMC(configModel, monitors = params)
          CmyMCMC <- compileNimble(myMCMC)

          results <- runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
        }
      }

      if(activity_estimation == "mixture") {
        run_MCMC_REST <- function(info, data, constants, code, params, ni, nt, nb) {
          myModel <- nimbleModel(code = code,
                                 data = data,
                                 constants = constants,
                                 inits = info$inits)
          CmyModel <- compileNimble(myModel)

          configModel <- configureMCMC(myModel, monitors = params)
          configModel$removeSampler(c("activity_proportion"))
          configModel$addSampler(
            target = c("activity_proportion"),
            type = 'prior_samples',
            samples = info$actv_samples
          )

          myMCMC <- buildMCMC(configModel, monitors = params)
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
          return(results)
        }

        per_chain_info <- lapply(1:length(actv_out_trace), function(i) {
          list(
            seed = sample(1:9999, 1),
            inits = inits_f(),
            actv_samples = as.matrix(actv_out_trace[[i]])
          )
        })
      }

      if(activity_estimation == "kernel") {
        per_chain_info <- lapply(1:nc, function(i) {
          list(
            seed = sample(1:9999, 1),
            inits = inits_f()
          )
        })
      }
      if(activity_estimation == "mixture") {
        per_chain_info <- lapply(1:length(actv_out_trace), function(i) {
          list(
            seed = sample(1:9999, 1),
            inits = inits_f(),
            actv_samples = as.matrix(actv_out_trace[[i]])
          )
        })
      }

      if(stay_family == "exponential") prms <- c("scale", "mean_stay")
      if(stay_family == "gamma" | stay_family == "weibull") prms <- c("scale", "shape", "mean_stay")
      if(stay_family == "lognormal") prms <- c("meanlog", "sdlog", "mean_stay")
      if(stay_family == "exponential") prms <- c("scale", "shape", "mean_stay")

      prms <- c(prms, "density", "p", "size")

      params <- c(prms, "loglike_obs_stay", "loglike_obs_y", "loglike_pred_stay", "loglike_pred_y")

      if(activity_estimation == "mixture") {
        params <- c(params, "activity_proportion")
      }

      cat("Running MCMC sampling. Please wait...\n")

      this_cluster <- makeCluster(nc)
      clusterEvalQ(this_cluster, {
        library(nimble)
      })

      clusterExport(this_cluster,
                    c("dvonMises", "rvonMises", "registerDistributions", "run_MCMC_REST"),
                    envir = environment())
      chain_output <- parLapply(
        cl = this_cluster,
        X = per_chain_info,
        fun = run_MCMC_REST,
        data = data_REST,
        code = Model_REST,
        constants = cons_REST,
        params = params,
        ni = ni,
        nt = nt,
        nb = nb
      )
      stopCluster(this_cluster)
      cat("Estimation is finished!\n")

      loglfstay <- MCMCchains(chain_output, params = c("loglike_obs_stay"))
      loglfy <- MCMCchains(chain_output, params = c("loglike_obs_y"))

      if(activity_estimation == "mixture") {
        loglfall<-cbind(loglfstay, loglfy, loglact)
      } else {
        loglfall<-cbind(loglfstay, loglfy)
      }
      lppd<-sum(log(colMeans(exp(loglfall))))
      p.waic<-sum(apply(loglfall, 2, var))
      waic[k] <- (-2) * lppd+2 * p.waic

      # Summarize the estimates -------------------------------------------------
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

    }
  }

  # RAD-REST model ----------------------------------------------------------

  if(model == "RAD-REST") {

    model_frame_density <- model.frame(formula_density, station_effort_data)
    X_density <- model.matrix(as.formula(formula_density), model_frame_density)
    nPreds_density <- ncol(X_density)

    # prepare dataset for the original REST

    y <- station_effort_data %>% select(starts_with("y_")) %>% as.matrix()
    N_judge <- apply(y, 1, sum)
    N_group <- ncol(y)
    N_detection <- station_effort_data %>% pull(N)

    waic <- numeric(0)
    tidy_samples <- mcmc_samples <- list(0)

    for(k in 1:length(formula_density_all)) {
      formula_density <- formula_density_all[[k]]

      model_frame_density <- model.frame(formula_density, station_effort_data)
      X_density <- model.matrix(as.formula(formula_density), model_frame_density) # model.matrix(formula_stay, model_frame_stay)

      if (!is.null(random_effect)) {
        levels <- unique(station_effort_data[[random_effect]])
        nLevels_stay  <- length(levels)
      } else {
        nLevels_stay  <- 0
      }

      nPreds_density <- ncol(X_density)


      data_REST <- list(y = y,
                        N_detection = N_detection,
                        X_density = X_density,
                        stay = stay,
                        censored = censored,
                        N_judge = N_judge)

      cons_REST <- list(N_station = N_station,
                        S = S,
                        N_period = N_period,
                        nPreds_stay = nPreds_stay,
                        N_stay = N_stay,
                        c_time = c_time,
                        nPreds_density = nPreds_density,
                        stay_family = stay_family,
                        N_group = N_group,
                        activity_estimation = activity_estimation)

      if(activity_estimation == "kernel") {
        cons_REST$activity_proportion <- activity_proportion
      }

      # Random effects
      if (!is.null(random_effect)) {
        cons_REST$group_stay <- as.numeric(factor(stay_data_join[[random_effect]]))
        cons_REST$nLevels_stay  <- nLevels_stay
      } else {
        cons_REST$nLevels_stay  <- 0
      }
      # Covariates
      if (nPreds_stay > 1) {
        data_REST$X_stay <- X_stay
      }


      # REST model
      Model_REST <- nimbleCode(
        {
          for(i in 1:N_stay) {
            censored[i] ~ dinterval(stay[i], c_time[i])

            if (stay_family == "exponential") {
              stay[i] ~ dexp(rate = 1/scale[i])
              pred_t[i] ~ dexp(rate = 1/scale[i])
              loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dexp(stay[i], rate = 1/scale[i], log = 1) +
                step(censored[i] - 0.5) * log(1 - pexp(c_time[i], rate = 1/scale[i]))
              loglike_pred_stay[i] <- dexp(pred_t[i], rate = 1/scale[i], log = 1)
            }

            if (stay_family == "gamma") {
              stay[i] ~ dgamma(shape = theta_stay, rate = exp(-log(scale[i])))
              pred_t[i] ~ dgamma(shape = theta_stay, rate = exp(-log(scale[i])))
              loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dgamma(stay[i], shape = theta_stay, rate = exp(-log(scale[i])), log = 1) +
                step(censored[i] - 0.5) * log(1 - pgamma(c_time[i], shape = theta_stay, rate = exp(-log(scale[i]))))
              loglike_pred_stay[i] <- dgamma(pred_t[i], shape = theta_stay, rate = exp(-log(scale[i])), log = 1)
            }

            if (stay_family == "lognormal") {
              stay[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay)
              pred_t[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay)
              loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dlnorm(stay[i], meanlog = log(scale[i]), sdlog = theta_stay, log = 1) +
                step(censored[i] - 0.5) * log(1 - plnorm(c_time[i], meanlog = log(scale[i]), sdlog = theta_stay))
              loglike_pred_stay[i] <- dlnorm(pred_t[i], meanlog = log(scale[i]), sdlog = theta_stay, log = 1)
              meanlog[i] <- log(scale[i])
            }

            if (stay_family == "weibull") {
              stay[i] ~ dweibull(shape = theta_stay, scale = scale[i])
              pred_t[i] ~ dweibull(shape = theta_stay, scale = scale[i])
              loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dweibull(stay[i], shape = theta_stay, scale = scale[i], log = 1) +
                step(censored[i] - 0.5) * log(1 - pweibull(c_time[i], shape = theta_stay, scale = scale[i]))
              loglike_pred_stay[i] <- dweibull(pred_t[i], shape = theta_stay, scale = scale[i], log = 1)
            }

            if (nPreds_stay > 1) {
              if (nLevels_stay  == 0) {
                log(scale[i]) <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay])
              } else {
                log(scale[i]) <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]]
              }
            } else {
              if (nLevels_stay  == 0) {
                log(scale[i]) <- beta_stay[1]
              } else {
                log(scale[i]) <- beta_stay[1] + random_effect_stay[group_stay[i]]
              }
            }
          }

          theta_stay ~ dgamma(1, 1)

          if (stay_family == "lognormal") {
            sdlog <- theta_stay
          } else {
            shape <- theta_stay
          }

          for(j in 1:nPreds_stay) {
            beta_stay[j] ~ dnorm(0, 100)
          }

          if (nLevels_stay  > 0) {
            for(k in 1:nLevels_stay ) {
              random_effect_stay[k] ~ dnorm(0, sd = sigma_stay)
            }
            sigma_stay ~ T(dnorm(0, sd = 100), 0, 5)
          }

          # Expected value calculation
          if (nPreds_stay == 1) {

            if(stay_family == "exponential") {
              mean_stay <- exp(beta_stay[1])
            }
            if(stay_family == "gamma") {
              mean_stay <- theta_stay * exp(beta_stay[1])
            }
            if(stay_family == "lognormal") {
              mean_stay <- exp(beta_stay[1] + theta_stay ^ 2 / 2)
            }
            if(stay_family == "weibull") {
              mean_stay <-  lgamma(1 + 1 / theta_stay) + exp(beta_stay[1])
            }
          }
          if (nPreds_stay > 1) {
            if(stay_family == "exponential") {
              for(i in 1:N_station){
                mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
              }
            }
            if(stay_family == "gamma") {
              for(i in 1:N_station){
                mean_stay[i] <- theta_stay * exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
              }
            }
            if(stay_family == "lognormal") {
              for(i in 1:N_station){
                mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + theta_stay ^ 2 / 2)
              }
            }
            if(stay_family == "weibull") {
              for(i in 1:N_station){
                mean_stay[i] <- lgamma(1 + 1 / theta_stay) + exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
              }
            }
          }

          # mean_pass
          for (g in 1:N_group) {
            alpha_Dirichlet[g] ~ dgamma(shape = 0.01, rate = 0.01)
          }

          # Expected values
          alpha_sum <- sum(alpha_Dirichlet[1:N_group])
          for (g in 1:N_group) {
            p_expected[g] <- alpha_Dirichlet[g] / alpha_sum
            c_expected[g] <- p_expected[g] * (g - 1)
          }
          for (j in 1:N_station) {
            y[j, 1:N_group] ~ ddirchmulti(alpha_Dirichlet[1:N_group], N_judge[j])
            pred_y[j, 1:N_group] ~ ddirchmulti(alpha_Dirichlet[1:N_group], N_judge[j])
            loglike_obs_y[j] <- ddirchmulti(y[j, 1:N_group], alpha_Dirichlet[1:N_group], N_judge[j], log = 1)
            loglike_pred_y[j] <- ddirchmulti(pred_y[j, 1:N_group], alpha_Dirichlet[1:N_group], N_judge[j], log = 1)
          }
          mean_pass <- sum(c_expected[1:N_group])


          # likelihood
          for(i in 1:N_station){
            N_detection[i] ~ dnbinom(size = size, prob = p[i])
            p[i] <- size / (size + mu[i])

            N_detection_rep[i] ~ dnbinom(size = size, prob = p[i])
            loglike_obs_detection[i] <- dnbinom(N_detection[i], size, p[i])
            loglike_pred_detection[i] <- dnbinom(N_detection_rep[i], size, p[i])
          }
          size ~ dgamma(1, 1)

          # REST formula
          if(nPreds_density == 1) {
            for(i in 1:N_station) {
              if(nPreds_stay == 1) {
                log(mu[i]) <- log(density) + log(S) + log(N_period[i]) - log(mean_stay) + log(activity_proportion) - log(mean_pass)
              } else {
                log(mu[i]) <- log(density) + log(S) + log(N_period[i]) - log(mean_stay[i]) + log(activity_proportion) - log(mean_pass)
              }
            }
            log(density) <- beta_density
            beta_density ~ dnorm(0, sd = 100)
          }
          if(nPreds_density > 1) {
            for(i in 1:N_station) {
              if(nPreds_stay == 1) {
                log(mu[i]) <- log(density[i]) + log(S) + log(N_period[i]) - log(mean_stay) + log(activity_proportion) - log(mean_pass)
              } else {
                log(mu[i]) <- log(density[i]) + log(S) + log(N_period[i]) - log(mean_stay[i]) + log(activity_proportion) - log(mean_pass)
              }
              log(density[i]) <- inprod(beta_density[1:nPreds_density], X_density[i, 1:nPreds_density])
            }
            for(j in 1:nPreds_density) {
              beta_density[j] ~ dnorm(0, sd = 100)
            }
          }
        }#end
      )

      inits_f <- function() {
        common_inits <- list(
          beta_stay = rnorm(nPreds_stay, 0, 1),
          stay = ifelse(censored == 0, c_time, c_time + rexp(N_stay, rate = 1 / mean(c_time))),
          theta_stay = runif(1, 0.1, 1.5),
          log_scale = runif(N_stay, -2.0, 2.0),
          alpha_Dirichlet = rexp(N_group, rate = 1),
          beta_density = rnorm(nPreds_density, 0, 1),
          size = runif(1, 1, 10),
          mu = rep(mean(N_detection, na.rm = TRUE), N_station)
        )

        if (!is.null(random_effect)) {
          random_effect_inits <- list(
            random_effect_stay = rnorm(nLevels_stay , 0, 1),
            sigma_stay = runif(1, 0.1, 0.5)
          )
          return(c(common_inits, random_effect_inits))
        }
        return(common_inits)
      }
      if(activity_estimation != "mixture") {
        run_MCMC_RAD <- function(info, data, constants, code, params, ni, nt, nb) {
          myModel <- nimbleModel(code = code,
                                 data = data,
                                 constants = constants,
                                 inits = info$inits)

          CmyModel <- compileNimble(myModel)

          configModel <- configureMCMC(myModel, monitors = params)
          myMCMC <- buildMCMC(configModel, monitors = params) # configModelを入れる
          CmyMCMC <- compileNimble(myMCMC)

          results <- runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
        }
        if(activity_estimation != "mixture") {
          per_chain_info <- lapply(1:nc, function(i) {
            list(
              seed = sample(1:9999, 1),
              inits = inits_f()
            )
          })
        }
      }
      if(activity_estimation == "mixture") {
        run_MCMC_RAD <- function(info, data, constants, code, params, ni, nt, nb) {
          myModel <- nimbleModel(code = code,
                                 data = data,
                                 constants = constants,
                                 inits = info$inits)
          CmyModel <- compileNimble(myModel)

          configModel <- configureMCMC(myModel, monitors = params)
          configModel$removeSampler(c("activity_proportion")) # (add_remove)
          configModel$addSampler(
            target = c("activity_proportion"), # add_remove,
            type = 'prior_samples',
            samples = info$actv_samples
          )

          myMCMC <- buildMCMC(configModel, monitors = params)
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
          return(results)
        }

        per_chain_info <- lapply(1:length(actv_out_trace), function(i) {
          list(
            seed = sample(1:9999, 1),
            inits = inits_f(),
            actv_samples = as.matrix(actv_out_trace[[i]])
          )
        })
      }

      ddirchmulti <- nimbleFunction(
        run = function(x = double(1),
                       alpha = double(1),
                       size = double(0),
                       log = integer(0)) {
          returnType(double(0))
          logProb <- lgamma(size + 1) - sum(lgamma(x + 1)) + lgamma(sum(alpha)) -
            sum(lgamma(alpha)) + sum(lgamma(alpha + x)) -
            lgamma(sum(alpha) + size)
          if (log)
            return(logProb)
          else
            return(exp(logProb))
        }
      )
      rdirchmulti <- nimbleFunction(
        run = function(n = integer(0),
                       alpha = double(1),
                       size = double(0)) {
          returnType(double(1))
          if (n != 1)
            print("rdirchmulti only allows n = 1; using n = 1.")
          p <- rdirch(1, alpha)
          return(rmulti(1, size = size, prob = p))
        }
      )
      dvonMises <- nimbleFunction(
        run = function(x = double(0), kappa = double(0), mu = double(0), log = integer(0)) {
          returnType(double(0))
          ccrit <- 1E-6
          s <- 1
          i <- 1
          inc <- 1
          x_2i <- 0
          satisfied <- FALSE
          while(!satisfied) {
            x_2i <- kappa / (2 * i)
            inc <- inc * x_2i * x_2i
            s <- s + inc
            i <- i + 1
            satisfied <- inc < ccrit
          }
          prob <- exp(kappa * cos(x - mu)) / (2 * pi * s)
          if (log) {
            return(log(prob))
          } else {
            return(prob)
          }
        }
      )
      rvonMises <- nimbleFunction(
        run = function(n = integer(0), kappa = double(0), mu = double(0)) {
          returnType(double(0))
          return(0)
        }
      )

      suppressMessages(registerDistributions(list(
        dvonMises = list(
          BUGSdist = "dvonMises(kappa, mu)",
          types = c('value = double(0)', 'kappa = double(0)', 'mu = double(0)'),
          pqAvail = FALSE
        ),
        ddirchmulti = list(
          BUGSdist = "ddirchmulti(alpha, size)",
          types = c('value = double(1)', 'alpha = double(1)', 'size = double(0)'),
          pqAvail = FALSE
        )
      )))

      if(stay_family == "exponential") prms <- c("scale", "mean_stay")
      if(stay_family == "gamma" | stay_family == "weibull") prms <- c("scale", "shape", "mean_stay")
      if(stay_family == "lognormal") prms <- c("meanlog", "sdlog", "mean_stay")
      if(stay_family == "exponential") prms <- c("scale", "shape", "mean_stay")

      prms <- c(prms, "density", "mean_stay", "mean_pass", "mu", "alpha_Dirichlet", "p", "size", "beta_stay", "beta_density")
      if(activity_estimation == "mixture") {
        prms <- c(prms, "activity_proportion")
      }
      params <- c(prms, "loglike_obs_stay", "loglike_obs_y", "loglike_obs_detection", "loglike_pred_detection", "loglike_pred_stay", "loglike_pred_y")

      cat("Running MCMC sampling. Please wait...\n")


      this_cluster <- makeCluster(nc)
      clusterEvalQ(this_cluster, {
        library(nimble)
      })

      clusterExport(this_cluster,
                    c("dvonMises", "rvonMises", "ddirchmulti", "rdirchmulti", "registerDistributions", "run_MCMC_RAD"),
                    envir = environment())

      chain_output <- parLapply(
        cl = this_cluster,
        X = per_chain_info,
        fun = run_MCMC_RAD,
        data = data_REST,
        code = Model_REST,
        constants = cons_REST,
        params = params,
        ni = ni,
        nt = nt,
        nb = nb
      )
      stopCluster(this_cluster)
      cat("Estimation is finished!\n")

      ## WAIC
      loglfy <- MCMCchains(chain_output, params = c("loglike_obs_y"))
      loglfstay <- MCMCchains(chain_output, params = c("loglike_obs_stay"))
      loglfN <- MCMCchains(chain_output, params = c("loglike_obs_detection"))

      if(activity_estimation == "mixture") loglfall <- cbind(loglfstay, loglfy, loglfN, loglact)
      if(activity_estimation == "kernel") loglfall <- cbind(loglfstay, loglfy, loglfN)

      lppd <- sum(log(colMeans(exp(loglfall))))
      p.waic <- sum(apply(loglfall, 2, var))
      waic[k] <- (-2) * lppd + 2 * p.waic

      # Summarize results -------------------------------------------------
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
    }
  }


  if(is.null(random_effect)) random_effect <- "NULL"
  WAIC <- data.frame(Model = as.character(unlist(formula_density_all)),
                     Random_effect = random_effect, WAIC = waic) %>%
    arrange(WAIC)

  best.model <- which(WAIC[1, 1]  == unlist(formula_density_all))

  mcmc_samples <- mcmc_samples[[best.model]]
  tidy_samples_best  <- tidy_samples[[best.model]]
  # p_value <- p_value[best.model]

  if(activity_estimation == "mixture") {
    sample_activity <- MCMCvis::MCMCchains(actv_chain_output,
                                           mcmc.list = TRUE,
                                           params = c("activity_proportion"))
    mcmc_samples <- lapply(seq_along(mcmc_samples), function(i) {
      cbind(mcmc_samples[[i]], sample_activity[[i]])
    })

    # mcmc_samples <- purrr::map2(mcmc_samples, sample_activity, ~ cbind(.x, .y))
  }
  prms <- c("density", "mean_stay")
  if(model == "RAD-REST") prms <- c(prms, "mean_pass")
  if(activity_estimation == "mixture") prms <- c(prms, "activity_proportion")

  mcmc_samples_mean <- MCMCvis::MCMCchains(chain_output, mcmc.list = TRUE, params = prms)
  summary_mean_temp <- MCMCsummary(mcmc_samples_mean, round = 2) %>%
    tibble::rownames_to_column(., var = "Variable") %>%
    tibble(.) %>%
    rename(lower = `2.5%`, median = `50%`, upper = `97.5%`) %>%
    mutate(cv = sd / mean)

  if(nPreds_density > 1) {
    summary_mean <- tibble(Species = target_species, Station = station.id) %>%
      mutate(summary_mean_temp)
  } else {
    summary_mean <- tibble(Species = target_species, summary_mean_temp)
  }

  density_result <- list(
    WAIC = WAIC,
    summary_result = summary_mean,
    samples = mcmc_samples
  )
  if(activity_estimation == "mixture") {
    density_result$activity_curve <- activity_density_estimates
  }
  class(density_result) <- "ResultDensity"

  density_result
}
time <- Stay <- Cens <- NULL
