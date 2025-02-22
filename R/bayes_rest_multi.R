#' Bayesian parameter estimation of multispecies RAD-REST model based on MCMC samplings with Rstan
#'
#' @param formula_stay Model formula for the staying time within a focal area. Grammar follows lme4::glmer function. e.g. Stay ~ 1 + (1|Group)
#' @param formula_density Model formula for animal density. e.g. ~ 1 + x1
#' @param station_effort_data A data frame containing information for each camera station (with sampling days) in each row
#' @param stay_data A data frame containing the time spent within a focal area
#' @param random_effect A random effect on mean staying time#'
#' @param col_name_species Column name of stay_data containing information on species name
#' @param activity_data A vector of detection times transformed into radians
#' @param activity_estimation "kernel" or "mixture"
#' @param add_prior """
#' @param bw_adj Bandwidth adjustment of kernel density estimation. The default value is 1.5. See Rowcliffe et al. () for details.
#' @param stay_family The probability distribution of staying times. Specify the probability distribution selected by the bayes_stay_selection function.
#' @param focal_area The area of a focal area
#' @param cores The number of cores to use for parallel computation. The default value is 1.
#' @param iter The length of iterations. The default value is 2000.
#' @param warmup The length of warmup. The default value is NULL (half of the iterations will be used for warmup).
#' @param chains The number of chains. The default value is 2.
#' @param thin The thinning interval. The default value is 1.
#' @param all_comb If TRUE, models with all combinations of covariates are compared. If FALSE, only the designated model in model_formula is run.
#' @param model Specify the model to use. Choose either "REST" or "RAD-REST".
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
#'   bw_adj = 1.5,
#'   stay_family = "lognormal",
#'   focal_area = 1.96,
#'   cores = 2,
#'   iter = 3000,
#'   warmup = 1000,
#'   chains = 2,
#'   thin = 1,
#'   model = "REST",
#'   all_comb = FALSE
#' )
#'
bayes_rest_multi <- function(formula_stay,
                             formula_density,
                             station_effort_data,
                             stay_data,
                             col_name_cens = "Cens",
                             col_name_species = "Species",
                             random_effect = NULL,
                             activity_data,
                             activity_estimation = "kernel",
                             add_prior = out_trace,
                             bw_adj = 1.5,
                             C = 10,
                             stay_family = "lognormal",
                             focal_area,
                             cores = 2,
                             iter = 3000,
                             warmup = NULL,
                             chains = 2,
                             thin = 1,
                             target_species,
                             all_comb = F) {



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

  # 活動時間割合を分離する場合 -----------------------------------------------------------

  ni <- iter
  nt <- thin
  nc <- chains
  nb <- warmup
  nmcmc <- (iter - warmup) * chains / thin
  actv_out_trace <- list(0)
  if(activity_estimation == "prior_samples") {
    out_trace_0 <- list(0)
    for(s in 1:length(target_species)) {

      act_data <- activity_data %>%
        filter(Species == target_species[s]) %>%
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
        actv <- 1.0 / (2 * 3.141592654 * max(activity_dens[1:ndens]))

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
      chain_output <- parLapply(cl = this_cluster, X = per_chain_info,
                                fun = run_MCMC_vonMises,
                                data = data, code = code,
                                constants = constants, params = params,
                                ni = ni, nt = nt, nb = nb)
      stopCluster(this_cluster)

      out_trace_0[[s]] <- chain_output %>%
        map(~ .[, grep(paste("actv", collapse = "|"), colnames(.))])

    }


    for (j in 1:length(chain_output)) { # チェーンごとにループ
      actv_list <- lapply(out_trace_0, function(species_results) {
        as.vector(species_results[[j]]) # 特定のチェーンの結果を抽出
      })

      actv_matrix <- do.call(cbind, actv_list)
      colnames(actv_matrix) <- paste0("actv[", 1:length(target_species), "]")
      actv_out_trace[[j]] <- coda::mcmc(actv_matrix) # codaオブジェクトに変換して格納
    }
  }
  ######################################
  if(activity_estimation == "prior_add") {
    actv_chain_output <- add_prior
    actv_out_trace <- actv_chain_output %>%
      map(~ .[, grep(paste("actv", collapse = "|"), colnames(.))])

    activity_estimation <- "prior_samples"

    # out_vonMises <- coda::mcmc.list(actv_chain_output) %>% posterior::as_draws()
    # summary <- posterior::summarize_draws(out_vonMises)
    # activity_density_estimates <- summary %>%
    #   dplyr::filter(stringr::str_starts(variable, "activity_dens")) %>%
    #   dplyr::mutate(x = dens.x)

    # out_trace <- add_prior
    # activity_estimation <- "prior_samples"
  }
  ######################################

  # Define data for stay ----------------------------------------------------

  target_species <- sort(target_species)

  station_effort_data <- station_effort_data %>%
    filter(Species %in% target_species)

  stay_data_join <- stay_data %>%
    filter(Species %in% target_species) %>%
    left_join(station_effort_data, by = intersect(names(stay_data), names(station_effort_data))) %>%
    filter(!is.na(Stay))

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
    nLevels <- length(levels)
  } else {
    nLevels <- 0
  }

  # N of covariates
  N_station_species <- nrow(station_effort_data)
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
    actv <- numeric(0)
    for(i in 1:nSpecies) {
      act_sp <- activity_data %>%
        filter(Species == target_species[i])
      model_act <- fitact(act_sp %>% pull(time), bw = bw_adj * bwcalc(act_sp %>% pull(time), K = 3), reps = 1)
      actv[i] <- model_act@act
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

  # pattern_alpha <- paste0("e_y\\[", 1:nSpecies, "\\]")
  # pattern_stay <- c(outer(1:N_station, 1:nSpecies, function(s, g) {
  #   paste0("mean_stay\\[", s, ", ", g, "\\]")
  # }))
  # pattern <- c(pattern_alpha, pattern_stay)

  S <- focal_area * 10^-6
  N_period <- station_effort_data$Effort[1:N_station] * 60 * 60 * 24
  N_detection_matrix <- matrix(N_detection, ncol = nSpecies, byrow = F)
  X_density <- X_density[1:N_station , ]

  # 滞在時間の種IDを数値化
  species_id_stay <- stay_data_join %>% pull(Species) %>% factor() %>% as.numeric()

  # 侵入回数の種IDを数値化
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
    cons_density$actv <- actv
  }

  if (!is.null(random_effect)) {
    cons_density$group <- as.numeric(factor(stay_data_join[[random_effect]]))
    cons_density$nLevels <- nLevels
  } else {
    cons_density$nLevels <- 0
  }

  if (nPreds_stay > 1) {
    data_density$X_stay <- X_stay
  }


  # MCMC
  code <- nimbleCode({
    for(i in 1:N_stay) {
      is.censored[i] ~ dinterval(stay[i], c_time[i])

      if(stay_family == "exponential") {
        stay[i] ~ dexp(rate = 1/exp(log_scale[i]))
      }
      if(stay_family == "gamma") {
        stay[i] ~ dgamma(shape = theta_stay[species_id_stay[i]], rate = 1/exp(log_scale[i]))  # theta_stay[species_id_stay[i]]を指数変換
      }
      if(stay_family == "lognormal") {
        stay[i] ~ dlnorm(meanlog = log_scale[i], sdlog = theta_stay[species_id_stay[i]])  # theta_stay[species_id_stay[i]]を指数変換
      }
      if(stay_family == "weibull") {
        stay[i] ~ dweibull(shape = theta_stay[species_id_stay[i]], scale = exp(log_scale[i]))  # theta_stayを指数変換
      }

      if (nPreds_stay > 1) {
        if (nLevels == 0) {
          log_scale[i] <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay])
        } else {
          log_scale[i] <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group[i]]
        }
      } else {
        if (nLevels == 0) {
          log_scale[i] <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1]
        } else {
          log_scale[i] <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] + random_effect_stay[group[i]]
        }
      }
    }

    # 事前分布の修正
    for(j in 1:nPreds_stay) {
      beta_stay[j] ~ dnorm(0, sd = 100)
    }

    for(m in 1:nSpecies) {
      theta_stay[m] ~ dgamma(shape_stay, rate_stay)
    }
    shape_stay ~ dgamma(0.01, 0.01)
    rate_stay ~ dgamma(0.01, 0.01)

    if (nLevels > 0) {
      for(k in 1:nLevels) {
        random_effect_stay[k] ~ dnorm(0, sd = sigma_stay)
      }
      sigma_stay ~ T(dt(0, 1 / pow(2.5, -2), 1), 0, 10)  # 非負制約を追加
    }

    for(m in 1:nSpecies) {
      for(j in 1:nPreds_stay) {
        species_effect_stay[m, j] ~ dnorm(0, sd = sigma_species_stay) # [j]
      }
    }
    sigma_species_stay ~ T(dt(0, 1 / pow(2.5, -2), 1), 0, 10)  # 非負制約を追加


    # 期待値の計算を修正
    if (nPreds_stay == 1) {
      for(m in 1:nSpecies) {
        for(i in 1:N_station) {
          if(stay_family == "exponential") mean_stay[i, m] <- exp(beta_stay[1] + species_effect_stay[m, 1])

          if(stay_family == "gamma") mean_stay[i, m] <- theta_stay[m] * exp(beta_stay[1] + species_effect_stay[m, 1])

          if(stay_family == "lognormal") mean_stay[i, m] <- exp(beta_stay[1] + species_effect_stay[m, 1] + (theta_stay[m] ^ 2) / 2)

          if(stay_family == "weibull") mean_stay[i, m] <- lgamma(1 + 1 / theta_stay[m]) * exp(beta_stay[1] + species_effect_stay[m, 1])
        }
      }
    }

    if (nPreds_stay > 1) {
      for(m in 1:nSpecies) {
        for(i in 1:N_station) {
          if(stay_family == "exponential") {
            mean_stay[i, m] <- exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
          }
          if(stay_family == "gamma") {
            mean_stay[i, m] <- theta_stay[m] * exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
          }
          if(stay_family == "lognormal") {
            mean_stay[i, m] <- exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + theta_stay[m] ^ 2 / 2)
          }
          if(stay_family == "weibull") {
            mean_stay[i, m] <- lgamma(1 + 1 / theta_stay[m]) * exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
          }
        }
      }
    }

    # 観測モデル
    for (j in 1:N_station_species) {
      y[j, 1:N_group] ~ ddirchmulti(alpha_dir_species[species_id_ey[j], 1:N_group], N_judge[j])
    }

    # 種ランダム効果の事前分布
    for (s in 1:nSpecies) {
      for (g in 1:N_group) {
        # alpha_dir_species[s, g] <- exp(log_alpha_dir_species[s, g])
        # log_alpha_dir_species[s, g] ~ dnorm(mu_alpha[g], sd = sd_species_y)

        alpha_dir_species[s, g] ~ dgamma(shape_alpha, rate_alpha)
      }
    }
    shape_alpha ~ dgamma(0.01, 0.01)
    rate_alpha ~ dgamma(0.01, 0.01)
    # sd_species_y ~ T(dnorm(0, sd = 10), 0, 10)
    # ランダム効果の分散
    # for(g in 1:N_group) {
    #   mu_alpha[g] ~ dnorm(0, sd = 100)
    # }

    # Expected values計算
    for(s in 1:nSpecies) {
      alpha_sum[s] <- sum(alpha_dir_species[s, 1:N_group])
      for(g in 1:N_group) {
        p_expected[s, g] <- alpha_dir_species[s, g] / alpha_sum[s]
        c_expected[s, g] <- p_expected[s, g] * (g - 1)
      }
      e_y[s] <- sum(c_expected[s, 1:N_group])
    }


    # likelihood
    for(s in 1:nSpecies) {
      for(i in 1:N_station){
        N_detection_matrix[i, s] ~ dnbinom(size = r[s], prob = p[i, s])
        p[i, s] <- r[s] / (r[s] + mu[i, s])
        # N_detection_matrix[i, s] ~ dpois(mu[i, s])
      }
      r[s] ~ dgamma(1, 1) #dgamma(shape_r, scale_r)
    }
    #shape_r ~ dgamma(0.1, 0.1)
    #scale_r ~ dgamma(0.1, 0.1)

    # Prior with Species Random Effect
    if (nPreds_density == 1) {
      for (s in 1:nSpecies) {
        for (i in 1:N_station) {
          log(mu[i, s]) <- log(density[s]) + log(S) + log(N_period[i]) - log(mean_stay[i, s]) + log(actv[s]) - log(e_y[s])
          # log(local_density[i, s]) ~ dnorm(beta_density + species_effect_density[s, 1], sd_density) #<- beta_density + species_effect_density[s, 1] #
          #  eps[i, s] ~ dnorm(0, sd_density)
          # local_density[i, s] <- exp(beta_density + species_effect_density[s, 1])
          # local_density[i, s] ~ dgamma(sd_density, sd_density / density[s])
        }
        density[s] <- exp(beta_density + species_effect_density[s, 1])
      }
      beta_density ~ dnorm(0, 100)
    }

    if (nPreds_density > 1) {
      for(i in 1:N_station) {
        for(s in 1:nSpecies) {
          log(mu[i, s]) <- log_local_density[i, s] + log(focal_area) + log(N_period[i]) - log(mean_stay[i, s]) + log(actv[s]) - log(e_y[s])
          log_local_density[i, s] <- dnorm(inprod(beta_density[1:nPreds_density] + species_effect_density[s, 1:nPreds_density], X_density[i, 1:nPreds_density]), sd_density)
        }
      }
      for(j in 1:nPreds_density) {
        beta_density[j] ~ dnorm(0, 100) #dunif(-10, 10) #dnorm(0, 100) # T(dnorm(0, sd = 5), -5, 5)
      }
      # for(s in 1:nSpecies) {
      #   for(i in 1:N_station) {
      #     density[i, s] <- exp(inprod(beta_density[1:nPreds_density] + species_effect_density[s, 1:nPreds_density], X_density[i, 1:nPreds_density]))
      #   }
      # }
    }

    # Random effects
    for(m in 1:nSpecies) {
      for(k in 1:nPreds_density) {
        species_effect_density[m, k] ~ dnorm(0, sd_species_density)
      }
    }
    sd_species_density ~ T(dt(0, 1 / pow(2.5, -2), 1), 0, 10) # dunif(0, 5) #
  })


  cat("Compiling the model. This may take a moment...\n")
  inits_f <- function() {
    common_inits <- list(
      beta_stay = runif(nPreds_stay, -1, 1),
      stay = ifelse(is.censored == 0, NA, c_time + runif(N_stay, 0.1, 1.0)),
      theta_stay = runif(nSpecies, 0.5, 4.0),
      log_scale = runif(N_stay, 0.5, 4.0),
      species_effect_stay = matrix(runif(nSpecies * nPreds_stay, -1, 1), nrow = nSpecies, ncol = nPreds_stay),

      beta_density = rnorm(nPreds_density, 0, 1),
      species_effect_density = matrix(rnorm(nSpecies * nPreds_density, 0, 0.5), nrow = nSpecies, ncol = nPreds_density),
      sd_density = runif(1, 0.01, 2),
      sd_species_density = runif(1, 0.01, 2),

      alpha_dir_species = matrix(runif(nSpecies * N_group, 0.01, 2), nrow = nSpecies, ncol = N_group),
      species_effect_ey = matrix(runif(nSpecies * N_group, 0.01, 2), nrow = nSpecies, ncol = N_group),
      sd_species_y = runif(1, 0.01, 2),
      # mu_alpha = runif(N_group, 0.1, 0.2)
      shape_alpha = runif(1, 0.1, 0.2),
      rate_alpha = runif(1, 0.1, 0.2)
    )

    if (!is.null(random_effect)) {
      c(
        common_inits,
        list(
          random_effect_stay = runif(nLevels, -1, 1),
          sigma_stay = runif(1, 0.5, 2.5)
        )
      )
    } else {
      c(
        common_inits,
        list(sigma_species_stay = runif(nPreds_stay, 0.5, 2.5))
      )
    }
  }

  # クラスター初期化
  this_cluster <- makeCluster(nc)

  # カスタム分布の登録をクラスター全体で評価
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

  # # MCMC実行関数
  # run_MCMC_RAD <- function(info, data, constants, code, params, ni, nt, nb) {
  #   myModel <- nimbleModel(code = code, data = data, constants = constants, inits = info$inits)
  #   CmyModel <- compileNimble(myModel)
  #   myMCMC <- buildMCMC(CmyModel, monitors = params)
  #   CmyMCMC <- compileNimble(myMCMC)
  #   runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
  # }

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

  params <- c("mean_stay",
              "e_y",
              "density",
              "log_scale",
              "theta_stay",
              "alpha_dir_species",
              "beta_stay",
              "species_effect_density",
              "r")

  if(activity_estimation == "prior_samples") {
    params <- c(params, "actv")
  }
  clusterEvalQ(this_cluster, {
    library(nimble)
  })
  # クラスターに必要な関数と変数をエクスポート
  clusterExport(this_cluster,
                c("ddirchmulti", "rdirchmulti", "registerDistributions", "run_MCMC_RAD"),
                envir = environment())
  cat("Running MCMC sampling. Please wait...\n")
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


  # mcmcサンプルをtidy形式に変換
  samples <- MCMCvis::MCMCchains(chain_output)
  nparams <- ncol(samples)
  nmcmc <- nrow(samples)

  tidy_samples <- samples %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    mutate(chain = rep(1:chains, each = nparams * nmcmc / chains)) %>%
    group_by(parameter, chain) %>%
    mutate(iteration = row_number()) %>%
    ungroup()


  station_id <- station_effort_data %>% arrange(Station) %>%  pull(Station)
  species_id <- stay_data %>%
    filter(Species %in% target_species) %>%
    arrange(Species) %>%
    pull(Species) %>%
    unique(.)

  stay_estimates <- tidy_samples %>%
    filter(grepl("^mean_stay", parameter)) %>%
    mutate(Species = rep(rep(species_id, each = N_station), nmcmc)) %>%
    mutate(Station = rep(station_id, nmcmc)) %>%
    group_by(Species, Station) %>% #
    summarise(
      median = median(value),
      mean = mean(value),
      sd = sd(value),
      cv = sd / mean,
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    filter(!duplicated(Species))

  # true_stay <- mean_stay
  # true_stay_df <- data.frame(
  #   Species = paste0("SP", sprintf("%03d", 1:length(true_stay))),
  #   true_value = true_stay
  # )
  #
  g_stay <- ggplot(stay_estimates, aes(x = Station, y = median)) +
    # geom_hline(data = true_stay_df,
    #            aes(yintercept = true_value),
    #            color = "red", linetype = "dashed") +
    geom_point() +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    facet_wrap(~Species, nrow = 3, ncol = 4) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      y = "Mean_stay Estimate (Median)",
      title = "Mean_stay",
      caption = "Red dashed line: True density"
    ) +
    ylim(0, NA)

#
#   test <- vector("list", 12)
#   names(test) <- sp
#   for(i in 1:12) {
#     mtest <- station_effort_data %>% filter(Species == sp[i]) %>% select(y_0:y_16) %>% as.matrix()
#
#     test[[i]] <- sum(apply(mtest, 2, sum) / sum(apply(mtest, 2, sum)) * 0:14)
#   }
#   test

  ey_estimates <- tidy_samples %>%
    filter(grepl("^e_y", parameter)) %>%
    mutate(Species = rep(species_id, nmcmc)) %>%
    group_by(Species) %>% #
    summarise(
      median = median(value),
      mean = mean(value),
      sd = sd(value),
      cv = sd / mean,
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    )
  #
  # ey_estimates %>% mutate(Species = unique(stay_data %>% pull(Species)))
  # true_stay <- mean_stay
  # true_stay_df <- data.frame(
  #   Species = paste0("SP", sprintf("%03d", 1:length(true_stay))),
  #   true_value = true_stay
  # )
  #
  g_ey <- ggplot(ey_estimates, aes(x = Species, y = median)) +
    # geom_hline(data = true_stay_df,
    #            aes(yintercept = true_value),
    #            color = "red", linetype = "dashed") +
    geom_point() +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      y = "N of enterings (Median)",
      title = "N of enterings (E_y)",
      caption = "Red dashed line: True density"
    )+
    ylim(0, NA)
  ######WAIC
  #
  # pattern <- paste0("log_scale\\[", 1:N_stay, "\\]")
  #   if(chains > 1) {
  #     log_scale.samp <- chain_output %>%
  #       map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
  #       do.call(rbind, .)
  #     if(stay_family != "exponential") tau.samp <- do.call(rbind, chain_output)[, "theta_stay"]
  #   } else {
  #     log_scale.samp <- chain_output
  #     if(stay_family != "exponential") tau.samp <- log_scale.samp[, "theta_stay"]
  #   }
  #
  #   if(stay_family != "exponential") tau.samp <- do.call(rbind, chain_output)[, "theta_stay"]
  #
  #   loglfstay<-matrix(NA, nrow(log_scale.samp), N_stay)
  #
  #   if(stay_family == "exponential") {
  #     for(j in 1:N_stay){
  #       if(is.na(stay[j])){
  #         loglfstay[, j] <- 1 - pexp(c_time[j], rate = 1 /exp(log_scale.samp[ , j]), lower.tail = F, log.p = T)
  #       } else {
  #         loglfstay[, j] <- dexp(stay[j], rate = 1 / exp(log_scale.samp[ , j]), log = T)
  #       }
  #     }
  #   }
  #   if(stay_family == "gamma") {
  #     for(j in 1:N_stay){
  #       if(is.na(stay[j])){
  #         loglfstay[, j] <- 1 - pgamma(c_time[j], shape = tau.samp,lower.tail = F, scale = exp(log_scale.samp[ , j]), log.p = T)
  #       } else{
  #         loglfstay[, j] <- dgamma(stay[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), log = T)
  #       }
  #     }
  #   }
  #   if(stay_family == "lognormal") {
  #     for(j in 1:N_stay){
  #       if(is.na(stay[j])){
  #         loglfstay[, j] <- 1 - plnorm(c_time[j], log_scale.samp[ , j], tau.samp,lower.tail = F, log.p = T)
  #       } else{
  #         loglfstay[, j] <- dlnorm(stay[j], log_scale.samp[ , j], tau.samp, log = T)
  #       }
  #     }
  #   }
  #   if(stay_family == "weibull") {
  #     for(j in 1:N_stay){
  #       if(is.na(stay[j])){
  #         loglfstay[, j] <- 1 - pweibull(c_time[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), lower.tail = F, log.p = T)
  #       } else{
  #         loglfstay[, j] <- dweibull(stay[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), log = T)
  #       }
  #     }
  #   }
  # #
  # #   pattern <- paste0("alpha_dir\\[", 1:N_group, "\\]")
  #
  #   pattern <- c(outer(1:nSpecies, 1:N_group, function(s, g) {
  #     paste0("alpha_dir\\[", s, ", ", g, "\\]")
  #   }))
  #
  #   alpha.samp <- chain_output %>%
  #     map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
  #     do.call(rbind, .)
  #
  #   loglfN <- extraDistr::ddirmnom(y, size = N_detection, alpha = alpha.samp[1:N_group])

  # loglfy <- matrix(NA, nrow(log_scale.samp), N_station)
  #
  #
  # pattern <- paste0("p\\[", 1:N_station, "\\]")
  # p.samp <- chain_output %>%
  #   map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
  #   do.call(rbind, .)
  #
  # r.samp <- tau.samp <- do.call(rbind, chain_output)[, "theta_stay"]
  # for(i in 1:N_station) {
  #   loglfy[,i] <- dnbinom(N_detection[i], prob = p.samp[,i], size = r.samp, log = T)
  # }
  #
  # loglfall<-cbind(loglfstay, loglfN, loglfy)
  # lppd<-sum(log(colMeans(exp(loglfall))))
  # p.waic<-sum(apply(loglfall,2,var))
  # waic[k] <- (-2)*lppd+2*p.waic
  # samples[[k]] <- chain_output


  density_estimates <- tidy_samples %>%
    filter(grepl("^density", parameter)) %>%
    mutate(Species = rep(species_id, nmcmc)) %>%
    group_by(Species) %>% #
    summarise(
      median = median(value),
      mean = mean(value),
      sd = sd(value),
      cv = sd / mean,
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    )

  if(activity_estimation != "kernel") {
    activity_estimates <- tidy_samples %>%
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
  # true_density <- density
  # true_density_df <- data.frame(
  #   Species = paste0("SP", sprintf("%03d", 1:length(true_density))),
  #   true_value = true_density
  # )

  g_dens <- ggplot(density_estimates, aes(x = Species, y = median)) +
    # geom_hline(data = true_density_df,
    #            aes(yintercept = true_value),
    #            color = "red", linetype = "dashed") +
    geom_point() +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    # facet_wrap(~Species, nrow = 3, ncol = 4) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      y = "Density Estimate (Median)",
      title = "Species Density Estimates by Station",
      caption = "Red dashed line: True density"
    ) +
    scale_y_log10()

  }

  # density_estimates %>%
  #   mutate(true_dens = density) %>%
  #   mutate(include = ifelse(lower <= true_dens & upper >= true_dens, 0, 1)) %>%
  #   pull(include) %>%
  #   sum(.)

  tidy_samples_dens <- tidy_samples %>%
    filter(grepl("^density", parameter)) %>%
    mutate(Species = rep(species_id, nmcmc))


  g_mcmc_density <- ggplot(tidy_samples_dens, aes(x = iteration, y = value)) +
    geom_line(aes(col = as.character(chain))) +
    facet_wrap(~Species)

  print(g_mcmc_density)

  return(list(g_stay = g_stay,
              g_ey = g_ey,
              g_dens = g_dens,
              stay_estimates =stay_estimates,
              ey_estimates = ey_estimates,
              density_estimates = density_estimates,
              g_mcmc_density = g_mcmc_density,
              tidy_samples = tidy_samples
              ))
}

