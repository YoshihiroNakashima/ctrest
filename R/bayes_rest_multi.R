#' Bayesian Parameter Estimation of a Multispecies RAD-REST Model Based on MCMC Sampling with RStan
#' @param formula_stay A model formula for staying times within the focal area (e.g., Stay ~ 1 + x1). Random effects should be specified separately using the random_effect argument.
#' @param formula_density A model formula for animal density (e.g., ~ 1 + x1). Note that the formula must not have a left-hand side.
#' @param station_effort_data A data frame containing information for each camera station. Typically, this is the output of the add_effort function. Alternatively, a manually prepared data frame may be provided with the following required columns:
#' - For model = "REST":
#' - Station (character): Camera station ID.
#' - Effort (numeric): Camera trapping effort (in days) at each station.
#' - Species (character): Species name.
#' - Y (numeric): Total number of passes through the focal area for each station-species pair.
#' - For model = "RAD-REST":
#' - Station (character): Camera station ID.
#' - Effort (numeric): Camera trapping effort (in days) at each station.
#' - Species (character): Species name.
#' - N (numeric): Total number of detected videos for each station-species pair.
#' - y_X columns (y_0, y_1, ..., y_max): Number of videos categorized by the number of observed passes.
#' @param stay_data A data frame returned by the format_stay function, containing:
#' - Station (character): Camera station ID.
#' - Species (character): Species name.
#' - Stay (numeric): Staying time (in seconds) within the focal area for each detected pass.
#' - Censored (binary): Whether the staying time was censored (1 = censored, 0 = observed).
#' @param random_effect A character string specifying a random effect on mean staying time. (Note: Species-level random effects are automatically included and do not need to be specified.)
#' @param activity_data A data frame containing a "time" column, representing detection times transformed into radians. Typically, this is the output of the format_activity function.
#' @param activity_estimation The method used to estimate activity patterns. Choose "kernel" for fixed kernel density estimation (Rowcliffe et al. 2014), or "mixture" for nonparametric Bayesian estimation using von Mises mixture models (Nakashima et al. 2025).
#' @param bw_adj Bandwidth adjustment parameter for kernel density estimation. The default is 1.0. See Rowcliffe et al. (2014) for details.
#' @param C The maximum number of von Mises components to use in the mixture model. Required only if activity_estimation = "mixture".
#' @param stay_family The probability distribution used to model staying times. This should be selected based on the output of the bayes_stay_selection function.
#' @param focal_area The size of the focal area, in square meters.
#' @param cores The number of CPU cores to use for parallel computation. Default is 3.
#' @param iter The total number of MCMC iterations per chain. Default is 5,000.
#' @param warmup The number of warm-up (burn-in) iterations per chain. Default is 1,000.
#' @param chains The number of MCMC chains. Default is 3.
#' @param thin The thinning interval for MCMC sampling. Default is 4 (no thinning).
#' @param all_comb Logical. If TRUE, all combinations of covariates in the density model will be evaluated and compared. If FALSE, only the model specified in formula_density will be run.
#' @param target_species A vector specifying the species to be analyzed. Multiple species must be specified.
#' @return A list of class \code{"ResultDensity"}, which includes the following components:
#' \describe{
#'   \item{\code{WAIC}}{An object containing WAIC (Widely Applicable Information Criterion) results.}
#'   \item{\code{summary_result}}{A data frame summarizing the posterior estimates of the mean staying time.}
#'   \item{\code{samples}}{A \code{coda::mcmc.list} object containing MCMC samples for all model parameters.}
#' }
#' The returned object has a custom print method that displays the WAIC and summary of staying time estimates.
#' You can access full MCMC samples via \code{$samples}, and analyze convergence using the \code{MCMCvis} package.

#' @export
#' @import nimble activity parallel MCMCvis tibble
#' @importFrom stats as.formula formula model.frame model.matrix sd var runif median quantile model.response rexp rnorm step dexp pexp dgamma pgamma dlnorm plnorm dweibull pweibull dnbinom
#' @importFrom dplyr select
#' @importFrom extraDistr ddirmnom
#' @examples
#' station_data_RAD <- format_station_data(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_y = "y",
#'   model = "RAD-REST"
#' )
#' station_effort_RAD <- add_effort(
#'   detection_data = detection_data,
#'   station_data_formatted = station_data_RAD,
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
#' rest_model <- bayes_rest_multi(
#'   formula_stay = Stay ~ 1,
#'   formula_density = ~ 1,
#'   station_effort_data = station_effort_RAD,
#'   stay_data = stay_data,
#'   activity_data = activity_data,
#'   stay_family = "lognormal",
#'   focal_area = 1.96,
#'   cores = 2,
#'   iter = 5000,
#'   warmup = 1000,
#'   chains = 2,
#'   thin = 1,
#'   target_species = c("SP01", "SP02", "SP03", "SP04", "SP05"),
#'   all_comb = FALSE
#' )
#'
bayes_rest_multi <- function(formula_stay,
                             formula_density,
                             station_effort_data,
                             stay_data,
                             random_effect = NULL,
                             activity_data,
                             activity_estimation = "kernel",
                             bw_adj = 1.0,
                             C = 10,
                             stay_family = "lognormal",
                             focal_area,
                             cores = 3,
                             iter = 5000,
                             warmup = 1000,
                             chains = 3,
                             thin = 4,
                             target_species,
                             all_comb = FALSE) {

  # Define functions --------------------------------------------------------
  cat("Running MCMC sampling. Please wait...\n")
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


  # check options -----------------------------------------------------------

  if (!inherits(formula_stay, "formula")) {
    stop("`formula_stay` must be a valid formula (e.g., Stay ~ 1 + x1).")
  }

  if (!inherits(formula_density, "formula")) {
    stop("`formula_density` must be a valid formula (e.g., ~ 1 + x1).")
  }

  if (!is.data.frame(station_effort_data)) {
    stop("`station_effort_data` must be a data frame.")
  }

  if (!is.data.frame(stay_data)) {
    stop("`stay_data` must be a data frame.")
  }

  if (!is.character(random_effect) && !is.null(random_effect)) {
    stop("`random_effect` must be a character vector or NULL.")
  }

  if (!is.data.frame(activity_data)) {
    stop("`activity_data` must be a data frame.")
  }

  if (!(activity_estimation %in% c("kernel", "mixture"))) {
    stop("`activity_estimation` must be either 'kernel' or 'mixture'.")
  }

  if (!is.numeric(bw_adj) || length(bw_adj) != 1 || bw_adj <= 0) {
    stop("`bw_adj` must be a positive number.")
  }

  if (activity_estimation == "mixture" && (missing(C) || !is.numeric(C) || length(C) != 1 || C <= 0)) {
    stop("`C` must be a positive integer when `activity_estimation = 'mixture'`.")
  }

  if (!is.character(stay_family) || length(stay_family) != 1) {
    stop("`stay_family` must be a single string.")
  }

  if (!is.numeric(focal_area) || length(focal_area) != 1 || focal_area <= 0) {
    stop("`focal_area` must be a positive number.")
  }

  if (!is.numeric(cores) || length(cores) != 1 || cores < 1) {
    stop("`cores` must be a positive integer.")
  }

  if (!is.numeric(iter) || iter <= 0) {
    stop("`iter` must be a positive integer.")
  }

  if (!is.numeric(warmup) || warmup <= 0) {
    stop("`warmup` must be a positive integer.")
  }

  if (!is.numeric(chains) || chains <= 0) {
    stop("`chains` must be a positive integer.")
  }

  if (!is.numeric(thin) || thin <= 0) {
    stop("`thin` must be a positive integer.")
  }

  if (!is.logical(all_comb) || length(all_comb) != 1) {
    stop("`all_comb` must be TRUE or FALSE.")
  }

  if (!is.character(target_species) || length(target_species) < 2) {
    stop("`target_species` must be a character vector with at least two species.")
  }

  # Activity proportion estimation ------------------------------------------

  ni <- iter
  nt <- thin
  nc <- chains
  nb <- warmup
  nmcmc <- (iter - warmup) * chains / thin
  actv_out_trace <- list(0)

  if(activity_estimation == "mixture") {
    out_trace_0 <- list(0)
    for(m in 1:length(target_species)) {

      act_data <- activity_data %>%
        filter(Species == target_species[m]) %>%
        pull(time)
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

      out_trace_0[[m]] <- actv_chain_output %>%
        map(~ .[, grep(paste("activity_proportion", collapse = "|"), colnames(.))])

    }

    for (j in 1:length(actv_chain_output)) {
      actv_list <- lapply(out_trace_0, function(species_results) {
        as.vector(species_results[[j]])
      })

      actv_matrix <- do.call(cbind, actv_list)
      colnames(actv_matrix) <- paste0("activity_proportion[", 1:length(target_species), "]")
      actv_out_trace[[j]] <- coda::mcmc(actv_matrix)
    }
  }

  # Define data for stay ----------------------------------------------------

  target_species <- sort(target_species)

  station_effort_data <- station_effort_data %>%
    filter(Species %in% target_species) %>%
    arrange(Species, Station)

  stay_data_join <- stay_data %>%
    filter(Species %in% target_species) %>%
    left_join(station_effort_data, by = intersect(names(stay_data), names(station_effort_data))) %>%
    filter(!is.na(Stay)) %>%
    arrange(Species, Station)

  if (stay_family == "lognormal" ||
      stay_family == "gamma" || stay_family == "weibull" || stay_family == "exponential") {

  } else {
    stop(paste0("Input family type(", stay_family, ") is incorrect."))
  }


  # Extract variable names
  vars_stay <- all.vars(formula_stay)
  response_stay <- vars_stay[1]
  predictors_stay <- vars_stay[-1]

  model_frame_stay <- model.frame(formula_stay, stay_data_join)
  X_stay <- model.matrix(as.formula(formula_stay), model_frame_stay) # model.matrix(formula_stay, model_frame_stay)
  stay <- model.response(model_frame_stay)
  is.censored <- stay_data_join[["Cens"]]

  # N of random effects
  if (!is.null(random_effect)) {
    levels <- unique(stay_data_join[[random_effect]])
    nLevels_stay <- length(levels)
  } else {
    nLevels_stay <- 0
  }

  # N of covariates
  N_station_species <- nrow(station_effort_data)
  station.id <- station_effort_data %>%
    arrange(Species, Station) %>%
    pull(Station)

  N_station <- length(unique(station_effort_data$Station))
  nPreds_stay <- ncol(X_stay)

  names(stay) <- NULL
  c_time <- stay
  c_time[is.censored == 0] <- c_time[is.censored == 0] + 1
  stay[is.censored == 1] <- NA
  N_stay <- length(stay)

  y <- station_effort_data %>% dplyr::select(starts_with("y_")) %>% as.matrix()
  N_judge <- apply(y, 1, sum)
  N_group <- ncol(y)
  N_detection <- station_effort_data %>% pull(N)
  nSpecies <- length(target_species)

  if(activity_estimation == "kernel") {
    activity_proportion <- numeric(0)
    for(i in 1:nSpecies) {
      act_sp <- activity_data %>%
        filter(Species == target_species[i])
      model_act <- fitact(act_sp %>% pull(time), bw = NULL, adj = bw_adj, reps = 1)
      activity_proportion[i] <- model_act@act
    }
  }

  vars_density <- all.vars(formula_density)
  predictors_density <- vars_density


  formula_density_all <- list(0)
  if (all_comb == TRUE) {
    formula_density_all <- full_terms(c(predictors_density))
  } else {
    formula_density_all[[1]] <- formula_density
  }

  samples <- list(0)
  waic <- numeric(0)
  density_estimates <- list(0)


  # Estimate density --------------------------------------------------------
  waic <- numeric(0)
  tidy_samples <- mcmc_samples <- list(0)

  for(k in 1:length(formula_density_all)) {

    formula_density <- formula_density_all[[k]]

    model_frame_density <- model.frame(formula_density, station_effort_data)
    X_density <- model.matrix(as.formula(formula_density), model_frame_density) # model.matrix(formula_stay, model_frame_stay)

    if (!is.null(random_effect)) {
      levels <- unique(station_effort_data[[random_effect]])
      nLevels_stay <- length(levels)
    } else {
      nLevels_stay <- 0
    }

    nPreds_density <- ncol(X_density)

    S <- focal_area * 10^-6
    N_period <- station_effort_data$Effort[1:N_station] * 60 * 60 * 24
    N_detection_matrix <- matrix(N_detection, ncol = nSpecies, byrow = F)
    X_density <- X_density[1:N_station , ]

    species_id_stay <- stay_data_join %>% pull(Species) %>% factor() %>% as.numeric()
    species_id_ey <- station_effort_data %>% pull(Species) %>% factor() %>% as.numeric()

    data_density <- list(stay = stay, is.censored = is.censored, y = y, N_judge = N_judge, N_detection_matrix = N_detection_matrix)

    if(nPreds_density > 1) data_density$X_density <- X_density

    cons_density <- list(
      N_stay = N_stay,
      nPreds_stay = nPreds_stay,
      N_station_species = N_station_species,
      N_station = N_station,
      c_time = c_time,
      stay_family = stay_family,
      N_group = N_group,
      nSpecies = nSpecies,
      species_id_stay = species_id_stay,
      species_id_ey = species_id_ey,
      S = S,
      N_period = N_period,
      nPreds_density = nPreds_density
    )
    if(activity_estimation == "kernel") {
      cons_density$activity_proportion <- activity_proportion
    }

    if (!is.null(random_effect)) {
      cons_density$group <- as.numeric(factor(stay_data_join[[random_effect]]))
      cons_density$nLevels_stay <- nLevels_stay
    } else {
      cons_density$nLevels_stay <- 0
    }

    if (nPreds_stay > 1) {
      data_density$X_stay <- X_stay
    }



  # Nimble code -------------------------------------------------------------

    code <- nimbleCode({

      # model for stay
      for(i in 1:N_stay) {
        censored[i] ~ dinterval(stay[i], c_time[i])

        if (stay_family == "exponential") {
          stay[i] ~ dexp(rate = 1 / scale[i])
          pred_t[i] ~ dexp(rate = 1 / scale[i])
          loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dexp(stay[i], rate = 1/scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pexp(c_time[i], rate = 1/scale[i]))
          loglike_pred_stay[i] <- dexp(pred_t[i], rate = 1/scale[i], log = 1)
        }

        if (stay_family == "gamma") {
          stay[i] ~ dgamma(shape = theta_stay[species_id_stay[i]], rate = 1 / scale[i])
          pred_t[i] ~ dgamma(shape = theta_stay[species_id_stay[i]], rate = 1/ scale[i])
          loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dgamma(stay[i], shape = theta_stay[species_id_stay[i]], rate = 1/ scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pgamma(c_time[i], shape = theta_stay[species_id_stay[i]], rate = 1 / scale[i]))
          loglike_pred_stay[i] <- dgamma(pred_t[i], shape = theta_stay[species_id_stay[i]], rate = exp(-log(scale[i])), log = 1)
        }

        if (stay_family == "lognormal") {
          stay[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]])
          pred_t[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]])
          loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dlnorm(stay[i], meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]], log = 1) +
            step(censored[i] - 0.5) * log(1 - plnorm(c_time[i], meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]]))
          loglike_pred_stay[i] <- dlnorm(pred_t[i], meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]], log = 1)

          meanlog[i] <- log(scale[i])
        }

        if (stay_family == "weibull") {
          stay[i] ~ dweibull(shape = theta_stay[species_id_stay[i]], scale = scale[i])
          pred_t[i] ~ dweibull(shape = theta_stay[species_id_stay[i]], scale = scale[i])
          loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dweibull(stay[i], shape = theta_stay[species_id_stay[i]], scale = scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pweibull(c_time[i], shape = theta_stay[species_id_stay[i]], scale = scale[i]))
          loglike_pred_stay[i] <- dweibull(pred_t[i], shape = theta_stay[species_id_stay[i]], scale = scale[i], log = 1)
        }

        if (nPreds_stay > 1) {
          if (nLevels_stay == 0) {
            log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay])
          } else {
            log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]]
          }
        } else { # nPreds_stay == 1
          if (nLevels_stay == 0) {
            log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1]
          } else {
            log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] + random_effect_stay[group_stay[i]]
          }
        }
      }

      for(j in 1:nPreds_stay) {
        beta_stay[j] ~ dnorm(0, sd = 100)
      }

      if (nLevels_stay > 0) {
        for(k in 1:nLevels_stay) {
          random_effect_stay[k] ~ dnorm(0, sd = sigma_stay)
        }
        sigma_stay ~ T(dnorm(0, sd = 100), 0, 5)
      }


      for(m in 1:nSpecies) {
        theta_stay[m] ~ dgamma(shape_stay, rate_stay)
        if(stay_family == "lognormal") {
          sdlog[m] <- theta_stay[m]
        } else {
          shape[m] <- theta_stay[m]
        }
      }
      shape_stay ~ dgamma(0.1, 0.1)
      rate_stay ~ dgamma(0.1, 0.1)

      for(m in 1:nSpecies) {
        for(j in 1:nPreds_stay) {
          species_effect_stay[m, j] ~ dnorm(0, sd = sigma_species_stay)
        }
      }
      sigma_species_stay ~  T(dnorm(0, sd = 100), 0, 5)

      # Expected value calculation
      if (nPreds_stay == 1) {
        if(stay_family == "exponential") {
          for(m in 1:nSpecies) {
            mean_stay[m] <- exp(beta_stay[1] + species_effect_stay[m, 1])
          }
        }
        if(stay_family == "gamma") {
          for(m in 1:nSpecies) {
            mean_stay[m] <- theta_stay[m] * exp(beta_stay[1] + species_effect_stay[m, 1])
          }
        }

        if(stay_family == "lognormal") {
          for(m in 1:nSpecies) {
            mean_stay[m] <- exp(beta_stay[1] + species_effect_stay[m, 1] + theta_stay[m] ^ 2 / 2)
          }
        }
        if(stay_family == "weibull") {
          for(m in 1:nSpecies) {
            mean_stay[m] <-  lgamma(1 + 1 / theta_stay[m]) + exp(beta_stay[1] + species_effect_stay[m, 1])
          }
        }
      }
      if (nPreds_stay > 1) {
        if(stay_family == "exponential") {
          for(m in 1:nSpecies) {
            for(i in 1:N_station) {
              mean_stay[i, m] <- exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
            }
          }
        }
        if(stay_family == "gamma") {
          for(m in 1:nSpecies) {
            for(i in 1:N_station) {
              mean_stay[i, m] <- theta_stay[m] * exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
            }
          }
        }
        if(stay_family == "lognormal") {
          for(m in 1:nSpecies) {
            for(i in 1:N_station) {
              mean_stay[i, m] <- exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + theta_stay[m] ^ 2 / 2)
            }
          }

        }
        if(stay_family == "weibull") {
          for(m in 1:nSpecies) {
            for(i in 1:N_station) {
              mean_stay[i, m] <-  lgamma(1 + 1 / theta_stay[m]) * exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
            }
          }
        }
      }

      # model for y
      for (j in 1:N_station_species) {
        y[j, 1:N_group] ~ ddirchmulti(alpha_Dirichlet[species_id_ey[j], 1:N_group], N_judge[j])

        pred_y[j, 1:N_group] ~ ddirchmulti(alpha_Dirichlet[species_id_ey[j], 1:N_group], N_judge[j])
        loglike_obs_y[j] <- ddirchmulti(y[j, 1:N_group], alpha_Dirichlet[species_id_ey[j], 1:N_group], N_judge[j], log = 1)
        loglike_pred_y[j] <- ddirchmulti(pred_y[j, 1:N_group], alpha_Dirichlet[species_id_ey[j], 1:N_group], N_judge[j], log = 1)
      }

      for (m in 1:nSpecies) {
        for (g in 1:N_group) {
          alpha_Dirichlet[m, g] ~ dgamma(shape_alpha, rate_alpha)
        }
      }
      shape_alpha ~ dgamma(0.01, 0.01)
      rate_alpha ~ dgamma(0.01, 0.01)

      for(m in 1:nSpecies) {
        alpha_sum[m] <- sum(alpha_Dirichlet[m, 1:N_group])
        for(g in 1:N_group) {
          p_expected[m, g] <- alpha_Dirichlet[m, g] / alpha_sum[m]
          c_expected[m, g] <- p_expected[m, g] * (g - 1)
        }
        mean_pass[m] <- sum(c_expected[m, 1:N_group])
      }

      # model for N_detection
      for(m in 1:nSpecies) {
        for(i in 1:N_station){
          N_detection_matrix[i, m] ~ dnbinom(size = size[m], prob = p[i, m])
          p[i, m] <- size[m] / (size[m] + mu[i, m])

          N_detection_rep[i, m] ~ dnbinom(size = size[m], prob = p[i, m])
          loglike_obs_detection[i] <- dnbinom(N_detection_rep[i, m], size[m], p[i, m])
          loglike_pred_detection[i] <- dnbinom(N_detection_rep[i, m], size[m], p[i, m])
        }
        size[m] ~ dgamma(1, 1)
      }

      # REST formula
      if(nPreds_density == 1) {
        for (m in 1:nSpecies) {
          for(i in 1:N_station) {
            if(nPreds_stay == 1) {
              log(mu[i, m]) <- log(density[m]) + log(S) + log(N_period[i]) - log(mean_stay[m]) + log(activity_proportion[m]) - log(mean_pass[m])
            } else {
              log(mu[i, m]) <- log(density[m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m]) - log(mean_pass[m])
            }
          }
          log(density[m]) <- beta_dens + species_effect_density[m, 1]
        }
        beta_dens ~ dnorm(0, sd = 100)
      }


      if(nPreds_density > 1) {
        for(m in 1:nSpecies) {
          for(i in 1:N_station) {
            if(nPreds_stay == 1) {
              log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[m]) + log(activity_proportion[m]) - log(mean_pass[m])
            } else {
              log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m]) - log(mean_pass[m])
            }
            log(density[i, m]) <- inprod(beta_dens[1:nPreds_density] + species_effect_density[m, 1:nPreds_density], X_density[i, 1:nPreds_density])
          }
        }
        for(j in 1:nPreds_density) {
            beta_dens[j] ~ dnorm(0, sd = 100)
        }
      }

      # Random effects
      for(m in 1:nSpecies) {
        for(k in 1:nPreds_density) {
          species_effect_density[m, k] ~ dnorm(0, sd_species_density)
        }
      }
      sd_species_density ~ T(dnorm(0, sd = 100), 0, 5)
    })


    cat("Compiling the model. This may take a moment...\n")
    inits_f <- function() {
      common_inits <- list(
        beta_stay = runif(nPreds_stay, -1, 1),
        stay = ifelse(is.censored == 0, NA, c_time + runif(N_stay, 0.1, 1.0)),
        theta_stay = runif(nSpecies, 0.5, 4.0),
        # scale = runif(N_stay, 0.5, 4.0),
        species_effect_stay = matrix(runif(nSpecies * nPreds_stay, -1, 1), nrow = nSpecies, ncol = nPreds_stay),

        beta_dens = rnorm(nPreds_density, 0, 1),
        species_effect_density = matrix(rnorm(nSpecies * nPreds_density, 0, 0.5), nrow = nSpecies, ncol = nPreds_density),
        sd_density = runif(1, 0.01, 2),
        sd_species_density = runif(1, 0.01, 2),

        alpha_Dirichlet = matrix(runif(nSpecies * N_group, 0.01, 2), nrow = nSpecies, ncol = N_group),
        species_effect_ey = matrix(runif(nSpecies * N_group, 0.01, 2), nrow = nSpecies, ncol = N_group),
        sd_species_y = runif(1, 0.01, 2),
        shape_alpha = runif(1, 0.1, 0.2),
        rate_alpha = runif(1, 0.1, 0.2)
      )

      if (!is.null(random_effect)) {
        c(
          common_inits,
          list(
            random_effect_stay = runif(nLevels_stay, -1, 1),
            sigma_stay = runif(1, 0.5, 2.5)
          )
        )
      } else {
        c(
          common_inits,
          list(sigma_species_stay = runif(1, 0.5, 2.5))
        )
      }
    }

    this_cluster <- makeCluster(nc)

    clusterEvalQ(this_cluster, {
      library(nimble)
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
      registerDistributions(list(
        ddirchmulti = list(
          BUGSdist = "ddirchmulti(alpha, size)",
          types = c('value = double(1)', 'alpha = double(1)', 'size = double(0)'),
          pqAvail = FALSE
        )
      ))
    })

    if(activity_estimation == "kernel") {
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
    }
    if(activity_estimation != "mixture") {
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

    prms <- c(prms, "density", "mean_stay", "mean_pass") #, "mu", "alpha_Dirichlet", "p", "size", "beta_stay", "beta_dens")

    params <- c(prms, "loglike_obs_stay", "loglike_obs_y", "loglike_obs_detection", "loglike_pred_detection", "loglike_pred_stay", "loglike_pred_y")

    if(activity_estimation == "mixture") {
      params <- c(params, "activity_proportion")
    }
    clusterEvalQ(this_cluster, {
      library(nimble)
    })

    clusterExport(this_cluster,
                  c("ddirchmulti", "rdirchmulti", "registerDistributions", "run_MCMC_RAD"),
                  envir = environment())
    # cat("Running MCMC sampling. Please wait...\n")
    chain_output <- parLapply(
      cl = this_cluster,
      X = per_chain_info,
      fun = run_MCMC_RAD,
      data = data_density,
      code = code,
      constants = cons_density,
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

    if(activity_estimation == "mixture") loglfall <- cbind(loglfstay, loglfy, loglfN)
    if(activity_estimation == "kernel") loglfall <- cbind(loglfstay, loglfy, loglfN)

    lppd <- sum(log(colMeans(exp(loglfall))))
    p.waic <- sum(apply(loglfall, 2, var))
    waic[k] <- (-2) * lppd + 2 * p.waic

    # Summarize results -------------------------------------------------
    mcmc_samples[[k]] <- MCMCvis::MCMCchains(chain_output,
                                             mcmc.list = TRUE,
                                             params = prms)

    samples <- MCMCvis::MCMCchains(chain_output, params = prms)
    tidy_samples[[k]] <- samples %>%
      as_tibble() %>%
      pivot_longer(cols = everything(),
                   names_to = "parameter",
                   values_to = "value") %>%
      group_by(parameter) %>%
      mutate(iteration = row_number()) %>%
      ungroup()
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
                                           params = "activity_proportion")
    mcmc_samples <- lapply(seq_along(mcmc_samples), function(i) {
      cbind(mcmc_samples[[i]], sample_activity[[i]])
    })

    # mcmc_samples <- purrr::map2(mcmc_samples, sample_activity, ~ cbind(.x, .y))
  }

  mcmc_samples_mean <- MCMCvis::MCMCchains(chain_output, mcmc.list = TRUE, params = c("density"))
  summary_mean_temp <- MCMCsummary(mcmc_samples_mean, round = 2) %>%
    tibble::rownames_to_column(., var = "Variable") %>%
    tibble(.) %>%
    rename(lower = `2.5%`, median = `50%`, upper = `97.5%`)

  if(nPreds_density == 1) {
    summary_mean <- summary_mean_temp %>%
      mutate(Species = station_effort_data %>% pull(Species) %>%  unique(.)) %>%
      select(Species, Variable:n.eff)
  } else {
    summary_mean <- summary_mean_temp %>%
      mutate(Species = station_effort_data %>% pull(Species), Station = station_effort_data %>% pull(Station)) %>%
      select(Species, Station, Variable:n.eff)
  }

  if(nPreds_stay == 1) {
    prms <- c("mean_stay", "mean_pass")
    mcmc_samples_mean <- MCMCvis::MCMCchains(chain_output, mcmc.list = TRUE, params = prms)
    summary_mean_temp <- MCMCsummary(mcmc_samples_mean, round = 2) %>%
      tibble::rownames_to_column(., var = "Variable") %>%
      tibble(.) %>%
      rename(lower = `2.5%`, median = `50%`, upper = `97.5%`) %>%
      mutate(Species = rep(target_species, 2))
  } else {
    prms <- c("mean_stay")
    mcmc_samples_mean <- MCMCvis::MCMCchains(chain_output, mcmc.list = TRUE, params = prms)
    summary_mean_temp_1 <- MCMCsummary(mcmc_samples_mean, round = 2) %>%
      tibble::rownames_to_column(., var = "Variable") %>%
      tibble(.) %>%
      rename(lower = `2.5%`, median = `50%`, upper = `97.5%`) %>%
      mutate(Species = station_effort_data %>% pull(Species), Station = station_effort_data %>% pull(Station)) %>%
      select(Species, Station, Variable:n.eff)

    prms <- c("mean_pass")
    mcmc_samples_mean <- MCMCvis::MCMCchains(chain_output, mcmc.list = TRUE, params = prms)
    summary_mean_temp_2 <- MCMCsummary(mcmc_samples_mean, round = 2) %>%
      tibble::rownames_to_column(., var = "Variable") %>%
      tibble(.) %>%
      rename(lower = `2.5%`, median = `50%`, upper = `97.5%`) %>%
      mutate(Species = target_species, Station = NA) %>%
      select(Species, Station, Variable:n.eff)

    summary_mean_temp <- summary_mean_temp_1 %>%
      bind_rows(summary_mean_temp_2)
  }

  summary_mean <- summary_mean %>%
    bind_rows(summary_mean_temp) %>%
    mutate(cv = sd / mean)

  density_result <- list(
    WAIC = WAIC,
    summary_result = summary_mean,
    samples = mcmc_samples
  )
  class(density_result) <- "ResultDensity"

  density_result
}
time <- Stay <- NULL
