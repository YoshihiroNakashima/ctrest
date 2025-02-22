#' Bayesian parameter estimation of REST/RAD-REST model based on MCMC samplings with Rstan
#'
#' @param formula_stay Model formula for the staying time within a focal area. Grammar follows lme4::glmer function. e.g. Stay ~ 1 + (1|Group)
#' @param formula_density Model formula for animal density. e.g. ~ 1 + x1
#' @param station_effort_data A data frame containing information for each camera station (with sampling days) in each row
#' @param stay_data A data frame containing the time spent within a focal area
#' @param random_effect A random effect on mean staying time
#' @param activity_data A vector of detection times transformed into radians
#' @param activity_estimation "kernel" or "mixture"
#' @param add_prior "kernel" or "mixture"
#' @param bw_adj Bandwidth adjustment of kernel density estimation. The default value is 1.5. See Rowcliffe et al. () for details.
#' @param C ....
#' @param stay_family The probability distribution of staying times. Specify the probability distribution selected by the bayes_stay_selection function.
#' @param focal_area The area of a focal area
#' @param cores The number of cores to use for parallel computation. The default value is 1.
#' @param iter The length of iterations. The default value is 2000.
#' @param warmup The length of warmup. The default value is NULL (half of the iterations will be used for warmup).
#' @param chains The number of chains. The default value is 2.
#' @param thin The thinning interval. The default value is 1.
#' @param all_comb If TRUE, models with all combinations of covariates are compared. If FALSE, only the designated model in model_formula is run.
#' @param model Specify the model to use. Choose either "REST" or "RAD-REST".
#' @param target_species Specify the model to use. Choose either "REST" or "RAD-REST".
#'
#' @return XXXX
#' @export
#' @import nimble activity parallel MCMCvis
#' @importFrom stats as.formula formula model.frame model.matrix sd var runif median quantile model.response
#' @importFrom dplyr select
#' @importFrom extraDistr ddirmnom
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
#' station_effort_REST <- add_effort(
#'   detection_data = detection_data,
#'   station_data_formatted = station_data_rest,
#'   col_name_station = "Station",
#'   col_name_term = "Term",
#'   col_name_datetime = "DateTime",
#'   plot = TRUE,
#'   font_size = 5
#' )
#' stay_data <- format_stay(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_stay = "Stay",
#'   col_name_cens = "Cens",
#'   target_species = "SP01"
#' )
#' activity_data <- format_activity(
#'   detection_data = detection_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_datetime = "DateTime",
#'   target_species = "SP01",
#'   indep_time = 30
#' )
#' library(nimble)
#' rest_model <- bayes_rest(
#'   formula_stay = Stay ~ 1 + x1,
#'   formula_density = ~ 1,
#'   station_effort_data = station_effort_REST,
#'   stay_data = stay_data,
#'   activity_data = activity_data,
#'   activity_estimation = "kernel",
#'   bw_adj = 1.5,
#'   C = 10,
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
                       add_prior = actv_chain_output,
                       bw_adj = 1.5,
                       stay_family = "lognormal",
                       focal_area,
                       cores = 3,
                       iter = 3000,
                       warmup = 1000,
                       chains = 3,
                       thin = 5,
                       model = "REST",
                       all_comb = F,
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

  # 引数指定の確認 -----------------------------------------------------------------

  if (stay_family == "lognormal" ||
      stay_family == "gamma" || stay_family == "weibull" || stay_family == "exponential") {

  } else {
    stop(paste0("Input stay_family type(", stay_family, ") is incorrect."))
  }
  # bayes_rest関数は、複数種の指定には対応していません。種ごとに推定するか、bayes_rest_multi関数を使って下さい。

  # データセットの作成 ---------------------------------------------------------------

  act_data <- activity_data %>%
    filter(Species == target_species) %>%
    pull(time)
  if(activity_estimation == "kernel") {
    model_act <- fitact(act_data, bw = bw_adj*bwcalc(act_data, K = 3), reps=1)
    actv <- model_act@act
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

  stay <- stay_data_join %>% pull(Stay) # model.response(model_frame_stay)
  censored <- stay_data_join %>% pull(Cens)

  c_time <- stay
  c_time[censored == 0] <- c_time[censored == 0] + 1
  stay[censored == 1] <- NA
  N_stay <- length(stay)

  S <- focal_area * 10^-6
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
    nLevels <- length(levels)
  } else {
    nLevels <- 0
  }

  # 共変量
  predictors_density <- all.vars(formula_density)

  formula_density_all <- list(0)
  if (all_comb == TRUE) {
    formula_density_all <- full_terms(c(predictors_density))
  } else {
    formula_density_all[[1]] <- formula_density
  }


  # 活動時間割合を分離する場合 -----------------------------------------------------------
  if(activity_estimation == "prior_samples") {
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
        mu_mix[k] ~ dunif(0, 2 * 3.141592654)  # von Misesの平均
        kappa_mix[k] ~ dgamma(1, 0.01)  # von Misesの集中度パラメータ
      }
      for(n in 1:N) {
        group[n] ~ dcat(w[1:C])
        act_data[n] ~ dvonMises(mu_mix[group[n]], kappa_mix[group[n]])
      }
      # 確率密度の計算
      for (j in 1:ndens) {
        for (i in 1:C) {
          dens.cpt[i, j] <- w[i] * dvonMises(dens.x[j] , mu_mix[i], kappa_mix[i], log = 0)
        }
        activity_dens[j] <- sum(dens.cpt[1:C, j])
      }
      # 活動時間割合
      actv <- 1.0 / (2 * 3.141592654 * max(activity_dens[1:ndens]));
    })

    inits_f <- function() {
      list(
        mu_mix = runif(constants$C, 0, 2*pi),
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

    registerDistributions(list(
      dvonMises = list(
        BUGSdist = "dvonMises(kappa, mu)",
        types = c('value = double(0)', 'kappa = double(0)', 'mu = double(0)'),
        pqAvail = FALSE
      )
    ))

    run_MCMC_vonMises <- function(info, data, constants, code, params, ni, nt, nb) {
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

    per_chain_info <- lapply(1:nc, function(i) {
      list(
        seed = sample(1:9999, 1),
        inits = inits_f()
      )
    })

    params <- c("activity_dens", "actv", "mu_mix", "kappa_mix", "w")
    cat("Running MCMC sampling. Please wait...\n")

    this_cluster <- makeCluster(nc)
    clusterEvalQ(this_cluster, {
      library(nimble)
    })
    # クラスターに必要な関数と変数をエクスポート
    clusterExport(this_cluster,
                  c("dvonMises", "rvonMises", "registerDistributions", "run_MCMC_vonMises"),
                  envir = environment())
    # MCMCの設定
    actv_chain_output <- parLapply(cl = this_cluster, X = per_chain_info,
                              fun = run_MCMC_vonMises,
                              data = data, code = code,
                              constants = constants, params = params,
                              ni = ni, nt = nt, nb = nb)
    stopCluster(this_cluster)

    actv_out_trace <- actv_chain_output %>%
      map(~ .[, grep(paste("actv", collapse = "|"), colnames(.))])

    # activity_density_estimates
    out_vonMises <- coda::mcmc.list(actv_chain_output) %>% posterior::as_draws()
    summary <- posterior::summarize_draws(out_vonMises)
    activity_density_estimates <- summary %>%
      dplyr::filter(stringr::str_starts(variable, "activity_dens")) %>%
      dplyr::mutate(x = dens.x)
  }
  ######################################
  if(activity_estimation == "prior_add") {
    actv_chain_output <- add_prior
    actv_out_trace <- actv_chain_output %>%
      map(~ .[, grep(paste("actv", collapse = "|"), colnames(.))])

    activity_estimation <- "prior_samples"

    out_vonMises <- coda::mcmc.list(actv_chain_output) %>% posterior::as_draws()
    summary <- posterior::summarize_draws(out_vonMises)
    activity_density_estimates <- summary %>%
      dplyr::filter(stringr::str_starts(variable, "activity_dens")) %>%
      dplyr::mutate(x = dens.x)
  }
  ######################################


  # Original-REST model ----------------------------------------------------------

  if(model == "REST") {

    tidy_samples <- list(0)
    waic <- numeric(0)
    density_estimates <- list(0)

    for(k in 1:length(formula_density_all)) {

      formula_density <- formula_density_all[[k]]

      # 共変量
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
                        focal_area = S,
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
        cons_REST$nLevels <- nLevels
      } else {
        cons_REST$nLevels <- 0
      }
      if(activity_estimation == "kernel") {
        cons_REST$actv <- actv
      }
      if(activity_estimation == "mixture") {
        data_REST$act_data <- activity_data %>% pull(time)
        cons_REST$N_act <- N_act
        cons_REST$dens.x <- seq(0, 2 * pi, 0.02)
        cons_REST$ndens <- length(dens.x)
        cons_REST$C <- 10
      }

      # Covariates
      if (nPreds_stay > 1) {
        data_REST$X_stay <- X_stay
      }

      Model_REST <- nimbleCode(
        {
          if(activity_estimation == "mixture") {
            for(c in 1:(C-1)) {
              v[c] ~ dbeta(1, alpha)
            }
            alpha ~ dgamma(1, 1)
            w[1:C] <- stick_breaking(v[1:(C-1)])
            for(c in 1:C) {
              mu_mix[c] ~ dunif(0, 2 * 3.141592654)  # von Misesの平均
              kappa_mix[c] ~ dgamma(1, 0.01)  # von Misesの集中度パラメータ
            }
            for(n in 1:N_act) {
              group[n] ~ dcat(w[1:C])
              act_data[n] ~ dvonMises(mu_mix[group[n]], kappa_mix[group[n]])
            }
            # 確率密度の計算
            for (j in 1:ndens) {
              for (i in 1:C) {
                dens.cpt[i, j] <- w[i] * dvonMises(dens.x[j] , mu_mix[i], kappa_mix[i], log = 0)
              }
              activity_dens[j] <- sum(dens.cpt[1:C, j])
            }
            # 活動時間割合
            actv <- 1.0 / (2 * 3.141592654 * max(activity_dens[1:ndens]));
          }

          # mean_stay

          for(i in 1:N_stay) {
            censored[i] ~ dinterval(stay[i], c_time[i])

            if(stay_family == "exponential") {
              stay[i] ~ dexp(rate = 1/exp(log_scale[i]))
            }
            if(stay_family == "gamma") {
              stay[i] ~ dgamma(shape = theta_stay, rate = 1/exp(log_scale[i]))
            }
            if(stay_family == "lognormal") {
              stay[i] ~ dlnorm(meanlog = log_scale[i], sdlog = theta_stay)  # theta_stay
            }
            if(stay_family == "weibull") {
              stay[i] ~ dweibull(shape = theta_stay, scale = exp(log_scale[i]))  # theta_stayを指数変換
            }

            if (nPreds_stay > 1) {
              if (nLevels == 0) {
                log_scale[i] <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay])
              } else {
                log_scale[i] <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]]
              }
            } else {
              if (nLevels == 0) {
                log_scale[i] <- beta_stay[1]
              } else {
                log_scale[i] <- beta_stay[1] + random_effect_stay[group_stay[i]]
              }
            }
          }


          theta_stay ~ dgamma(1, 1)

          for(j in 1:nPreds_stay) {
            beta_stay[j] ~ T(dnorm(0, 100), -10, 10)
          }

          if (nLevels > 0) {
            for(k in 1:nLevels) {
              random_effect_stay[k] ~ dnorm(0, sd = sigma_stay)
            }
            sigma_stay ~ T(dnorm(0, 100), 0, 10)
          }

          if (nPreds_stay == 1) {
            for(i in 1:N_station){
              if(stay_family == "exponential") mean_stay[i] <- exp(beta_stay[1])
              if(stay_family == "gamma") mean_stay[i] <- theta_stay * exp(beta_stay[1])
              if(stay_family == "lognormal") mean_stay[i] <- exp(beta_stay[1] + theta_stay ^ 2 / 2)
              if(stay_family == "weibull") mean_stay[i] <-  lgamma(1 + 1 / theta_stay) + exp(beta_stay[1])
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

          # likelihood
          for(i in 1:N_station){
            y[i] ~ dnbinom(size = r, prob = p[i])
            p[i] <- r / (r + mu[i])
          }
          r ~ dgamma(1, 1)

          # prior
          if(nPreds_density == 1) {
            for(i in 1:N_station) {
              log(mu[i]) <- log(density[i]) + log(focal_area) + log(N_period[i]) - log(mean_stay[i]) + log(actv)
              log(density[i]) <- beta_density
            }
            beta_density ~ dnorm(0, sd = 100) # dnorm(0, sd = 100)
          }
          if(nPreds_density > 1) {
            for(i in 1:N_station) {
              log(mu[i]) <- log_local_density[i] + log(focal_area) + log(N_period[i]) - log(mean_stay[i]) + log(actv)
              log_local_density[i] <- inprod(beta_density[1:nPreds_density], X_density[i, 1:nPreds_density]) + eps[i]
              eps[i] ~ dnorm(0, sd_density)
            }
            for(j in 1:nPreds_density) {
              beta_density[j] ~ dnorm(0, sd = 100)
            }
            for(i in 1:N_station) {
              density[i] <- exp(inprod(beta_density[1:nPreds_density], X_density[i, 1:nPreds_density]))
            }
          }
          sd_density ~ T(dt(0, 1 / pow(2.5, -2), 1), 0, 10)
        }#end
      )


      inits_f <- function() {
        common_inits <- list(
          beta_stay = runif(nPreds_stay, -1, 1),
          stay = ifelse(censored == 0, NA, c_time + runif(N_stay, 0.1, 1.0)),
          theta_stay = runif(1, 0.5, 1.5),
          log_scale = runif(N_stay, 0.5, 4.0),

          beta_density = rnorm(nPreds_density, 0, 1),
          sd_density = runif(1, 0.01, 2)
        )

        if (!is.null(random_effect)) {
          return(c(
            common_inits,
            list(
              random_effect_stay = runif(nLevels, -1, 1),
              sigma_stay = runif(1, 0.5, 2.5)
            )
          ))
        }
        if(activity_estimation == "mixture") {
          return(c(
            common_inits,
            inits <- list(mu_mix = runif(cons_REST$C, 0, 2*pi),
                          kappa_mix = rgamma(cons_REST$C, 1, 0.01),
                          group = sample(1:cons_REST$C, size = cons_REST$N_act, replace = TRUE),
                          v = rbeta(cons_REST$C-1, 1, 1),
                          alpha = 1)
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

      registerDistributions(list(
        dvonMises = list(
          BUGSdist = "dvonMises(kappa, mu)",
          types = c('value = double(0)', 'kappa = double(0)', 'mu = double(0)'),
          pqAvail = FALSE
        )
      ))
      if(activity_estimation != "prior_samples") {
        run_MCMC_REST <- function(info, data, constants, code, params, ni, nt, nb) {
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

      if(activity_estimation == "prior_samples") {
        run_MCMC_REST <- function(info, data, constants, code, params, ni, nt, nb) {
          myModel <- nimbleModel(code = code,
                                 data = data,
                                 constants = constants,
                                 inits = info$inits)
          CmyModel <- compileNimble(myModel)
          # MCMCの設定
          configModel <- configureMCMC(myModel, monitors = params)
          configModel$removeSampler(c("actv")) # (add_remove)
          configModel$addSampler(
            target = c("actv"), # add_remove,
            type = 'prior_samples',
            samples = info$actv_samples
          )

          myMCMC <- buildMCMC(configModel, monitors = params)
          CmyMCMC <- compileNimble(myMCMC)

          # MCMCの実行
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
      if(activity_estimation != "prior_samples") {
        per_chain_info <- lapply(1:nc, function(i) {
          list(
            seed = sample(1:9999, 1),
            inits = inits_f()
          )
        })
      }
      if(activity_estimation == "prior_samples") {
        per_chain_info <- lapply(1:length(actv_out_trace), function(i) {
          list(
            seed = sample(1:9999, 1),
            inits = inits_f(),
            actv_samples = as.matrix(actv_out_trace[[i]])
          )
        })
      }

      params <- c("density", "mean_stay", "log_scale", "theta_stay", "p", "r")
      if(activity_estimation == "mixture") {
        params <- c(params, "actv", "activity_dens")
      }

      cat("Running MCMC sampling. Please wait...\n")

      this_cluster <- makeCluster(nc)
      clusterEvalQ(this_cluster, {
        library(nimble)
      })
      # クラスターに必要な関数と変数をエクスポート
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

      ### 結果のまとめ
      tidy_samples[[k]] <- MCMCvis::MCMCchains(chain_output) %>%
        as_tibble() %>%
        pivot_longer(cols = everything(),
                     names_to = "parameter",
                     values_to = "value") %>%
        group_by(parameter) %>%
        mutate(iteration = row_number())

      library(coda)
      library(posterior)
      library(bayesplot)

      # # チェーン結果を coda::mcmc.list に統合
      # combined_output <- mcmc.list(chain_output)
      #
      # # coda オブジェクトを posterior の配列形式に変換
      # mcmc_array <- as_draws_array(combined_output)
      #
      # # トレースプロット
      # mcmc_trace(mcmc_array, pars = "mean_stay[2]")



      ######WAIC
      # pattern <- paste0("log_scale\\[", 1:N_stay, "\\]")
      #
      # log_scale.samp <- chain_output %>%
      #   map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
      #   do.call(rbind, .)

      # log_scale の列名を生成
      log_scale_indices <- seq_len(N_stay)
      log_scale_cols <- paste0("log_scale[", log_scale_indices, "]")

      # 必要な列だけを抽出
      log_scale.samp <- chain_output %>%
        map(~ .[, colnames(.) %in% log_scale_cols, drop = FALSE]) %>%
        do.call(rbind, .)

      if(stay_family != "exponential") tau.samp <- do.call(rbind, chain_output)[, "theta_stay"]

      loglfstay<-matrix(0, nrow(log_scale.samp), N_stay)

      if(stay_family == "exponential") {
        for(j in 1:N_stay){
          if(is.na(stay[j])){
            loglfstay[, j] <- log(1 - pexp(c_time[j], rate = 1 /exp(log_scale.samp[ , j])))
          } else {
            loglfstay[, j] <- dexp(stay[j], rate = 1 / exp(log_scale.samp[ , j]), log = T)
          }
        }
      }
      if(stay_family == "gamma") {
        for(j in 1:N_stay){
          if(is.na(stay[j])){
            loglfstay[, j] <- log(1 - pgamma(c_time[j], shape = tau.samp,lower.tail = F, scale = exp(log_scale.samp[ , j])))
          } else{
            loglfstay[, j] <- dgamma(stay[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), log = T)
          }
        }
      }
      if(stay_family == "lognormal") {
        for(j in 1:N_stay){
          if(is.na(stay[j])){
            loglfstay[, j] <- log(1 - plnorm(c_time[j], log_scale.samp[ , j], tau.samp))
          } else{
            loglfstay[, j] <- dlnorm(stay[j], log_scale.samp[ , j], tau.samp, log = T)
          }
        }
      }
      if(stay_family == "weibull") {
        for(j in 1:N_stay){
          if(is.na(stay[j])){
            loglfstay[, j] <- log(1 - pweibull(c_time[j], shape = tau.samp, scale = exp(log_scale.samp[ , j])))
          } else{
            loglfstay[, j] <- dweibull(stay[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), log = T)
          }
        }
      }


      # p の列名を生成
      p_indices <- seq_len(N_station)
      p_cols <- paste0("p[", p_indices, "]")

      # 必要な列だけを抽出
      p.samp <- chain_output %>%
        map(~ .[, colnames(.) %in% p_cols, drop = FALSE]) %>%
        do.call(rbind, .)
      # pattern <- paste0("p\\[", 1:N_station, "\\]")
      # p.samp <- chain_output %>%
      #   map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
      #   do.call(rbind, .)

      loglfy <- matrix(NA, nrow(log_scale.samp), N_station)
      r.samp <- do.call(rbind, chain_output)[, "r"]
      for(i in 1:N_station) {
        loglfy[,i] <- dnbinom(y[i], prob = p.samp[,i], size = r.samp, log = T)
      }

      loglfall<-cbind(loglfstay, loglfy)
      lppd<-sum(log(colMeans(exp(loglfall))))
      p.waic<-sum(apply(loglfall,2,var))
      waic[k] <- (-2)*lppd+2*p.waic


      density_estimates[[k]] <- tidy_samples[[k]] %>%
        filter(grepl("^density\\[\\d+\\]$", parameter)) %>%
        mutate(index = as.numeric(str_extract(parameter, "\\d+"))) %>%
        filter(index <= 100) %>%
        group_by(parameter, index) %>%
        summarise(
          median = median(value),
          mean = mean(value),
          sd = sd(value),
          cv = sd / mean,
          lower = quantile(value, 0.025),
          upper = quantile(value, 0.975),
          .groups = "drop"
        ) %>%
        arrange(index)
    }
  }

  # RAD-REST model ----------------------------------------------------------

  if(model == "RAD-REST") {
    # 共変量
    model_frame_density <- model.frame(formula_density, station_effort_data)
    X_density <- model.matrix(as.formula(formula_density), model_frame_density)
    nPreds_density <- ncol(X_density)

    # prepare dataset for the original REST

    y <- station_effort_data %>% select(starts_with("y_")) %>% as.matrix()
    N_judge <- apply(y, 1, sum)
    N_group <- ncol(y)
    N_detection <- station_effort_data %>% pull(N)

    #model_mcmc <- suppressMessages(nimbleModel(code, constants = cons_pre, data = data_pre, inits = inits))

    samples <- list(0)
    waic <- numeric(0)
    density_estimates <- tidy_samples <- list(0)

    for(k in 1:length(formula_density_all)) {
      formula_density <- formula_density_all[[k]]

      model_frame_density <- model.frame(formula_density, station_effort_data)
      X_density <- model.matrix(as.formula(formula_density), model_frame_density) # model.matrix(formula_stay, model_frame_stay)

      if (!is.null(random_effect)) {
        levels <- unique(station_effort_data[[random_effect]])
        nLevels <- length(levels)
      } else {
        nLevels <- 0
      }

      nPreds_density <- ncol(X_density)


      data_REST <- list(y = y,
                        N_detection = N_detection,
                        X_density = X_density,
                        stay = stay,
                        censored = censored,
                        N_judge = N_judge)

      cons_REST <- list(N_station = N_station,
                        focal_area = S,
                        N_period = N_period,
                        nPreds_stay = nPreds_stay,
                        N_stay = N_stay,
                        c_time = c_time,
                        nPreds_density = nPreds_density,
                        stay_family = stay_family,
                        N_group = N_group,
                        activity_estimation = activity_estimation)

      if(activity_estimation == "kernel") {
        cons_REST$actv <- actv
      }
      if(activity_estimation == "mixture") {
        data_REST$act_data <- activity_data %>% pull(time)
        cons_REST$N_act <- N_act
        cons_REST$dens.x <- dens.x
        cons_REST$ndens <- length(dens.x)
        cons_REST$C <- C
      }

      # Random effects
      if (!is.null(random_effect)) {
        cons_REST$group_stay <- as.numeric(factor(stay_data_join[[random_effect]]))
        cons_REST$nLevels <- nLevels
      } else {
        cons_REST$nLevels <- 0
      }
      # Covariates
      if (nPreds_stay > 1) {
        data_REST$X_stay <- X_stay
      }


      # REST model
      Model_REST <- nimbleCode(
        {
          if(activity_estimation == "mixture") {
            for(c in 1:(C-1)) {
              v[c] ~ dbeta(1, alpha)
            }
            alpha ~ dgamma(1, 1)
            w[1:C] <- stick_breaking(v[1:(C-1)])
            for(c in 1:C) {
              mu_mix[c] ~ dunif(0, 2 * 3.141592654)  # von Misesの平均
              kappa_mix[c] ~ dgamma(1, 0.01)  # von Misesの集中度パラメータ
            }
            for(n in 1:N_act) {
              group[n] ~ dcat(w[1:C])
              act_data[n] ~ dvonMises(mu_mix[group[n]], kappa_mix[group[n]])
            }
            # 確率密度の計算
            for (j in 1:ndens) {
              for (i in 1:C) {
                dens.cpt[i, j] <- w[i] * dvonMises(dens.x[j], mu_mix[i], kappa_mix[i], log = 0)
              }
              activity_dens[j] <- sum(dens.cpt[1:C, j])
            }
            # 活動時間割合
            actv <- 1.0 / (2 * 3.141592654 * max(activity_dens[1:ndens]))
          }

          for(i in 1:N_stay) {
            censored[i] ~ dinterval(stay[i], c_time[i])
            if(stay_family == "exponential") {
              stay[i] ~ dexp(scale = exp(log_scale[i]))
            }
            if(stay_family == "gamma") {
              stay[i] ~ dgamma(shape = theta_stay, scale = exp(log_scale[i]))
            }
            if(stay_family == "lognormal") {
              stay[i] ~ dlnorm(meanlog = log_scale[i], sdlog = theta_stay)
            }
            if(stay_family == "weibull") {
              stay[i] ~ dweibull(shape = theta_stay, scale = exp(log_scale[i]))
            }

            if (nPreds_stay > 1) {
              if (nLevels == 0) {
                log_scale[i] <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay])
              }
              if (nLevels > 0) {
                log_scale[i] <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]]
              }
            }
            if (nPreds_stay == 1) {
              if (nLevels == 0) {
                log_scale[i] <- beta_stay[1]
              }
              if (nLevels > 0) {
                log_scale[i] <- beta_stay[1] + random_effect_stay[group_stay[i]]
              }
            }
          }

          theta_stay ~ dgamma(0.1, 0.1)

          for(j in 1:nPreds_stay) {
            beta_stay[j] ~ T(dnorm(0, sd = 20), -5, 5)
          }

          if (nLevels > 0) {
            for(k in 1:nLevels) {
              random_effect_stay[k] ~ dnorm(0, sd = sigma_stay)
            }
            sigma_stay ~ T(dnorm(0, sd = 100), 0, 5)
          }

          if (nPreds_stay == 1) {
            for(i in 1:N_station){
              if(stay_family == "exponential") mean_stay[i] <- exp(beta_stay[1])
              if(stay_family == "gamma") mean_stay[i] <- theta_stay * exp(beta_stay[1])
              if(stay_family == "lognormal") mean_stay[i] <- exp(beta_stay[1] + theta_stay ^ 2 / 2)
              if(stay_family == "weibull") mean_stay[i] <-  lgamma(1 + 1 / theta_stay) + exp(beta_stay[1])
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

          # e_y
          # 観測モデル
          for (j in 1:N_station) {
            y[j, 1:N_group] ~ ddirchmulti(alpha_dir[1:N_group], N_judge[j])
          }

          # alpha の生成
          for (g in 1:N_group) {
            alpha_dir[g] ~ dgamma(shape = 0.01, rate = 0.01)
          }

          # Expected values計算
          alpha_sum <- sum(alpha_dir[1:N_group])
          for (g in 1:N_group) {
            p_expected[g] <- alpha_dir[g] / alpha_sum
            c_expected[g] <- p_expected[g] * (g - 1)
          }
          e_y <- sum(c_expected[1:N_group])


          # likelihood
          for(i in 1:N_station){
            N_detection[i] ~ dnbinom(size = r, prob = p[i])
            p[i] <- r / (r + mu[i])
          }
          r ~ dgamma(1, 1)

          # prior
          if(nPreds_density == 1) {
            for(i in 1:N_station) {
              log(mu[i]) <- log(density[i]) + log(focal_area) + log(N_period[i]) - log(mean_stay[i]) + log(actv) - log(e_y)
              log(density[i]) <- beta_density
            }
            beta_density ~ dnorm(0, sd = 100)
          }
          if(nPreds_density > 1) {
            for(i in 1:N_station) {
              log(mu[i]) <- log_local_density[i] + log(focal_area) + log(N_period[i]) - log(mean_stay[i]) + log(actv) - log(e_y)
              log_local_density[i] <- inprod(beta_density[1:nPreds_density], X_density[i, 1:nPreds_density]) + eps[i]
              eps[i] ~ dnorm(0, sd_density)
            }

            for(j in 1:nPreds_density) {
              beta_density[j] ~ dnorm(0, sd = 100)
            }
            for(i in 1:N_station) {
              density[i] <- exp(inprod(beta_density[1:nPreds_density], X_density[i, 1:nPreds_density]))
            }
          }
          sd_density ~ T(dt(0, 1 / pow(2.5, -2), 1), 0, 10)


          # Posterior predictive check
          for (i in 1:N_station) {
            N_detection_rep[i] ~ dnbinom(size = r, prob = p[i])

            # 観測データの対数尤度
            loglike_obs[i] <- dnbinom(N_detection[i], r, p[i])

            # 複製データの対数尤度
            loglike_pred[i] <- dnbinom(N_detection_rep[i], r, p[i])
          }

          sum_loglike_obs <- sum(loglike_obs[1:N_station])
          sum_loglike_pred <- sum(loglike_pred[1:N_station])

        }#end
      )

      # 初期値関数
      inits_f <- function() {
        # 共通の初期値
        common_inits <- list(
          beta_stay = rnorm(nPreds_stay, 0, 1), # 小さな値で中心化
          stay = ifelse(censored == 0, c_time, c_time + rexp(N_stay, rate = 1 / mean(c_time))), # 分布に応じた初期値
          theta_stay = runif(1, 0.1, 1.5), # 汎用的な範囲
          log_scale = runif(N_stay, 0.5, 4.0), # X_stayから計算
          alpha_dir = rexp(N_group, rate = 1), # 正の値
          beta_density = rnorm(nPreds_density, 0, 1), # 小さめの正規分布
          sd_density = runif(1, 0.1, 0.5), # 小さめの範囲
          r = runif(1, 1, 10), # 一様分布で広く初期化
          mu = rep(mean(N_detection, na.rm = TRUE), N_station) # データから計算
        )

        # ランダム効果が存在する場合
        if (!is.null(random_effect)) {
          random_effect_inits <- list(
            random_effect_stay = rnorm(nLevels, 0, 1), # 標準正規分布
            sigma_stay = runif(1, 0.1, 0.5) # 小さめの初期値
          )
          return(c(common_inits, random_effect_inits))
        }
        if(activity_estimation == "mixture") {
          return(c(
            common_inits,
            inits <- list(mu_mix = runif(cons_REST$C, 0, 2*pi),
                          kappa_mix = rgamma(cons_REST$C, 1, 0.01),
                          group = sample(1:cons_REST$C, size = cons_REST$N_act, replace = TRUE),
                          v = rbeta(cons_REST$C-1, 1, 1),
                          alpha = 1)
          ))
        }
        return(common_inits)
      }

      if(activity_estimation != "prior_samples") {
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

      if(activity_estimation == "prior_samples") {
        run_MCMC_RAD <- function(info, data, constants, code, params, ni, nt, nb) {
          myModel <- nimbleModel(code = code,
                                 data = data,
                                 constants = constants,
                                 inits = info$inits)
          CmyModel <- compileNimble(myModel)
          # MCMCの設定
          configModel <- configureMCMC(myModel, monitors = params)
          configModel$removeSampler(c("actv")) # (add_remove)
          configModel$addSampler(
            target = c("actv"), # add_remove,
            type = 'prior_samples',
            samples = info$actv_samples
          )

          myMCMC <- buildMCMC(configModel, monitors = params)
          CmyMCMC <- compileNimble(myMCMC)

          # MCMCの実行
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
      if(activity_estimation != "prior_samples") {
        per_chain_info <- lapply(1:nc, function(i) {
          list(
            seed = sample(1:9999, 1),
            inits = inits_f()
          )
        })
      }
      if(activity_estimation == "prior_samples") {
        per_chain_info <- lapply(1:length(actv_out_trace), function(i) {
          list(
            seed = sample(1:9999, 1),
            inits = inits_f(),
            actv_samples = as.matrix(actv_out_trace[[i]])
          )
        })
      }

        # 必要な関数を定義・登録
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

        registerDistributions(list(
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
        ))
      # })


      # チェーンごとの初期情報設定
      # per_chain_info <- lapply(1:nc, function(i) {
      #   list(seed = sample(1:9999, 1), inits = inits_f())
      # })

      params <- c("mean_stay",
                  "e_y",
                  "density",
                  "log_scale",
                  "theta_stay",
                  "alpha_dir",
                  "density",
                  "beta_stay",
                  "beta_density",
                  "p",
                  "r",
                  "sum_loglike_obs",
                  "sum_loglike_pred",
                  "mu")

      if(activity_estimation == "mixture") {
        params <- c(params, "actv", "activity_dens")
      }
      if(activity_estimation == "prior_samples") {
        params <- c(params, "actv")
      }
      cat("Running MCMC sampling. Please wait...\n")


      this_cluster <- makeCluster(nc)
      clusterEvalQ(this_cluster, {
        library(nimble)
      })
      # クラスターに必要な関数と変数をエクスポート
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


      ### 結果のまとめ
      # out <- mcmc.list(chain_output) %>% as_draws()
      # summary <- summary(out)
      # summary

      tidy_samples[[k]] <- MCMCvis::MCMCchains(chain_output) %>%
        as_tibble() %>%
        pivot_longer(cols = everything(),
                     names_to = "parameter",
                     values_to = "value") %>%
        group_by(parameter) %>%
        mutate(iteration = row_number()) %>%
        ungroup()

      ## WAIC
      ###### 滞在時間
      # log_scale の列名を生成
      log_scale_indices <- seq_len(N_stay)
      log_scale_cols <- paste0("log_scale[", log_scale_indices, "]")

      # 必要な列だけを抽出
      log_scale.samp <- chain_output %>%
        map(~ .[, colnames(.) %in% log_scale_cols, drop = FALSE]) %>%
        do.call(rbind, .)

      if(stay_family != "exponential") tau.samp <- do.call(rbind, chain_output)[, "theta_stay"]

      loglfstay<-matrix(0, nrow(log_scale.samp), N_stay)

      if(stay_family == "exponential") {
        for(j in 1:N_stay){
          if(is.na(stay[j])){
            loglfstay[, j] <- log(1 - pexp(c_time[j], rate = 1 /exp(log_scale.samp[ , j])))
          } else {
            loglfstay[, j] <- dexp(stay[j], rate = 1 / exp(log_scale.samp[ , j]), log = T)
          }
        }
      }
      if(stay_family == "gamma") {
        for(j in 1:N_stay){
          if(is.na(stay[j])){
            loglfstay[, j] <- log(1 - pgamma(c_time[j], shape = tau.samp,lower.tail = F, scale = exp(log_scale.samp[ , j])))
          } else{
            loglfstay[, j] <- dgamma(stay[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), log = T)
          }
        }
      }
      if(stay_family == "lognormal") {
        for(j in 1:N_stay){
          if(is.na(stay[j])){
            loglfstay[, j] <- log(1 - plnorm(c_time[j], log_scale.samp[ , j], tau.samp))
          } else{
            loglfstay[, j] <- dlnorm(stay[j], log_scale.samp[ , j], tau.samp, log = T)
          }
        }
      }
      if(stay_family == "weibull") {
        for(j in 1:N_stay){
          if(is.na(stay[j])){
            loglfstay[, j] <- log(1 - pweibull(c_time[j], shape = tau.samp, scale = exp(log_scale.samp[ , j])))
          } else{
            loglfstay[, j] <- dweibull(stay[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), log = T)
          }
        }
      }

      ###### 侵入回数
      pattern <- paste0("alpha_dir\\[", 1:N_group, "\\]")

      alpha.samp <- chain_output %>%
        map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
        do.call(rbind, .)

      loglfN <- extraDistr::ddirmnom(y, size = N_detection, alpha = alpha.samp[1:N_group])

      ###### 観測回数
      loglfy <- matrix(NA, nrow(log_scale.samp), N_station)

      pattern <- paste0("p\\[", 1:N_station, "\\]")
      p.samp <- chain_output %>%
        map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
        do.call(rbind, .)

      r.samp <- do.call(rbind, chain_output)[, "theta_stay"]
      for(i in 1:N_station) {
        loglfy[,i] <- dnbinom(N_detection[i], prob = p.samp[,i], size = r.samp, log = T)
      }

      ###### 統合
      loglfall<-cbind(loglfstay, loglfy)
      lppd<-sum(log(colMeans(exp(loglfall))))
      p.waic<-sum(apply(loglfall,2,var))
      waic[k] <- (-2)*lppd+2*p.waic


      density_estimates[[k]] <- tidy_samples[[k]] %>%
        filter(grepl("^density\\[\\d+\\]$", parameter)) %>%
        mutate(index = as.numeric(str_extract(parameter, "\\d+"))) %>%
        group_by(parameter, index) %>%
        summarise(
          median = median(value),
          mean = mean(value),
          sd = sd(value),
          cv = sd / mean,
          lower = quantile(value, 0.025),
          upper = quantile(value, 0.975),
          .groups = "drop"
        ) %>%
        arrange(index)
    } # kのループ
  }　# モデルのループ


  WAIC <- data.frame(Model = as.character(unlist(formula_density_all)), WAIC = waic) %>%
    arrange(WAIC)
  density <- density_estimates[[which(WAIC[1, 1]  == unlist(formula_density_all))]]
  all_samples <- tidy_samples[[which(WAIC[1, 1]  == unlist(formula_density_all))]]

  station_id <- station_effort_data %>% arrange(Station) %>%  pull(Station)

  stay_estimates <- all_samples %>%
    filter(grepl("^mean_stay", parameter)) %>%
    rename(Station = parameter) %>%
    group_by(Station) %>% #
    summarise(
      median = median(value),
      mean = mean(value),
      sd = sd(value),
      cv = sd / mean,
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    mutate(parameter_index = as.numeric(gsub("mean_stay\\[|\\]", "", Station))) %>%
    arrange(parameter_index) %>%
    select(-parameter_index) %>%
    mutate(Station = station_id)

  if(nPreds_stay == 1) {
    stay_estimates <- stay_estimates %>% dplyr::slice(1) %>% mutate(Station = "Mean_stay")
  }
  ey_estimates <- NULL
  if(model == "RAD-REST") {
    ey_estimates <- all_samples %>%
      filter(grepl("^e_y", parameter)) %>%
      summarise(
        median = median(value),
        mean = mean(value),
        sd = sd(value),
        cv = sd / mean,
        lower = quantile(value, 0.025),
        upper = quantile(value, 0.975),
        .groups = "drop"
      ) %>%
      mutate(Station = "Meay_nEntry") %>% # 列を追加
      relocate(Station, .before = median) # 一番左に移動
  }

  density_estimates <- all_samples %>%
    filter(grepl("^density", parameter)) %>%
    rename(Station = parameter) %>%
    group_by(Station) %>% #
    summarise(
      median = median(value),
      mean = mean(value),
      sd = sd(value),
      cv = sd / mean,
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    mutate(parameter_index = as.numeric(gsub("density\\[|\\]", "", Station))) %>%
    arrange(parameter_index) %>%
    select(-parameter_index) %>%
    mutate(Station = station_id)



  dispersion_estimates <- all_samples %>%
    filter(parameter == "r") %>%
    rename(Station = parameter) %>%
    group_by(Station) %>% #
    summarise(
      median = median(value),
      mean = mean(value),
      sd = sd(value),
      cv = sd / mean,
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    )




  if(activity_estimation == "mixture" | activity_estimation == "prior_samples") {
    activity_estimates <- all_samples %>%
    filter(parameter == "actv") %>% #
    summarise(
      median = median(value),
      mean = mean(value),
      sd = sd(value),
      cv = sd / mean,
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    )
  }

  if(nPreds_density == 1) {
    density_estimates <- density_estimates %>% dplyr::slice(1) %>% mutate(Station = "Mean_density")
  }

  # 結果の処理
  if(activity_estimation == "mixture") {
    out_vonMises <- coda::mcmc.list(chain_output) %>% posterior::as_draws()
    summary <- posterior::summarize_draws(out_vonMises)
    activity_density_estimates <- summary %>%
      dplyr::filter(stringr::str_starts(variable, "activity_dens")) %>%
      dplyr::mutate(x = dens.x)


    g_act <- ggplot() +
      geom_histogram(data = data.frame(time = time),
                     mapping = aes(x = time, y = after_stat(density)),
                     bins = 30, fill = "gray") +
      geom_line(data = activity_density_estimates,
                mapping = aes(x = x, y = mean),
                colour = "red", linewidth = 1) +
      geom_ribbon(data = activity_density_estimates,
                  mapping = aes(x = x, ymin = q5, ymax = q95),
                  fill = "red", alpha = 0.3) +
      labs(x = "Activity level",
           y = "Probability density")

    print(g_act)
  }

  # ベイズP値
  chi2_obs_total_samples <- tidy_samples[[k]] %>%
    filter(str_detect(parameter, "sum_loglike_obs"))

  chi2_pred_total_samples <- tidy_samples[[k]] %>%
    filter(str_detect(parameter, "sum_loglike_pred"))

  p_value <- mean(chi2_pred_total_samples > chi2_obs_total_samples, na.rm = TRUE)

  result <- list(
    WAIC = waic,
    stay_estimates = stay_estimates,
    ey_estimates = ey_estimates,
    density_estimates = density_estimates,
    tidy_samples = tidy_samples[[1]],
    bayes_p_value = p_value)

  if(activity_estimation != "kernel") {
    pattern <- paste0("alpha_dir\\[", 1:N_group, "\\]")

    alpha.samp <- chain_output %>%
      map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
      do.call(rbind, .)

    result$activity_estimates <- activity_estimates
    result$activity_density_estimates <- activity_density_estimates

    result$dispersion <- dispersion_estimates
  }
  return(result)
}
