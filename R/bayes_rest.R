#' Bayesian parameter estimation of REST/RAD-REST model based on MCMC samplings with Rstan
#'
#' @param formula_stay Model formula for the staying time within a focal area. Grammar follows lme4::glmer function. e.g. Stay ~ 1 + (1|Group)
#' @param formula_density Model formula for animal density. e.g. ~ 1 + x1
#' @param station_effort_data A data frame containing information for each camera station (with sampling days) in each row
#' @param stay_data A data frame containing the time spent within a focal area
#' @param col_name_cens Column name of stay_data containing information on censored (1) or not (0)
#' @param random_effect A random effect on mean staying time
#' @param activity_data A vector of detection times transformed into radians
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
#'   col_name_cens = "Cens",
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
bayes_rest <- function(formula_stay,
                       formula_density,
                       station_effort_data,
                       stay_data,
                       col_name_cens = "Cens",
                       random_effect = NULL,
                       activity_data,
                       bw_adj = 1.5,
                       stay_family = "lognormal",
                       focal_area,
                       cores = 2,
                       iter = 3000,
                       warmup = NULL,
                       chains = 2,
                       thin = 1,
                       model = "REST",
                       all_comb = F
){
  ni <- iter
  nt <- thin
  nc <- chains
  nb <- warmup


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


  model_act <- fitact(activity_data, bw = bw_adj*bwcalc(activity_data, K = 3), reps=1)
  actv <- model_act@act

  stay_data_join <- stay_data %>%
    left_join(station_effort_data, by = intersect(names(stay_data), names(station_effort_data)))

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
  t <- model.response(model_frame_stay)
  censored <- stay_data_join[["Cens"]]

  # N of random effects
  if (!is.null(random_effect)) {
    levels <- unique(stay_data_join[[random_effect]])
    nLevels <- length(levels)
  } else {
    nLevels <- 0
  }

  # N of covariates
  N_station <- nrow(station_effort_data)
  nPreds <- ncol(X_stay)

  names(t) <- NULL
  c_time <- t
  c_time[censored == 0] <- c_time[censored == 0] + 1
  t[censored == 1] <- NA
  N_stay <- length(t)


  # RAD-REST model ----------------------------------------------------------

  if(model == "RAD-REST") {
    y <- station_effort_data %>% dplyr::select(starts_with("y_")) %>% as.matrix()
    N_judge <- apply(y, 1, sum)
    N_group <- ncol(y)
    N_detection <- station_effort_data %>% pull(N)

    data_pre <- list(t = t, censored = censored, y = y, N_judge = N_judge)

    cons_pre <- list(
      N_stay = N_stay,
      nPreds = nPreds,
      N_station = N_station,
      c_time = c_time,
      stay_family = stay_family,
      N_group = N_group
    )

    if (!is.null(random_effect)) {
      cons_pre$group <- as.numeric(factor(stay_data_join[[random_effect]]))
      cons_pre$nLevels <- nLevels
    } else {
      cons_pre$nLevels <- 0
    }

    if (nPreds > 1) {
      data_pre$X_stay <- X_stay
    }

    # MCMC
    code <- nimbleCode(
      {

        for(i in 1:N_stay) {
          censored[i] ~ dinterval(t[i], c_time[i])
          if(stay_family == "exponential") {
            t[i] ~ dexp(scale = exp(log_scale[i]))
          }
          if(stay_family == "gamma") {
            t[i] ~ dgamma(shape = theta, scale = exp(log_scale[i]))
          }
          if(stay_family == "lognormal") {
            t[i] ~ dlnorm(meanlog = log_scale[i], sdlog = theta)
          }
          if(stay_family == "weibull") {
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
          beta[j] ~ T(dnorm(0, 100), -5, 5)
        }

        if (nLevels > 0) {
          for(k in 1:nLevels) {
            u[k] ~ dnorm(0, sd = sigma_u)
          }
          sigma_u ~ dt(0, 1 / pow(2.5, -2), 1)
        }

        if (nPreds == 1) {
          for(i in 1:N_station){
            if(stay_family == "exponential") mean_stay[i] <- exp(beta[1])
            if(stay_family == "gamma") mean_stay[i] <- theta * exp(beta[1])
            if(stay_family == "lognormal") mean_stay[i] <- exp(beta[1] + theta ^ 2 / 2)
            if(stay_family == "weibull") mean_stay[i] <-  lgamma(1 + 1 / theta) + exp(beta[1])
          }
        }

        if (nPreds > 1) {
          if(stay_family == "exponential") {
            for(i in 1:N_station){
              mean_stay[i] <- exp(inprod(beta[1:nPreds], X_stay[i, 1:nPreds]))
            }
          }
          if(stay_family == "gamma") {
            for(i in 1:N_station){
              mean_stay[i] <- theta * exp(inprod(beta[1:nPreds], X_stay[i, 1:nPreds]))
            }
          }
          if(stay_family == "lognormal") {
            for(i in 1:N_station){
              mean_stay[i] <- exp(inprod(beta[1:nPreds], X_stay[i, 1:nPreds]) + theta ^ 2 / 2)
            }
          }
          if(stay_family == "weibull") {
            for(i in 1:N_station){
              mean_stay[i] <- lgamma(1 + 1 / theta) + exp(inprod(beta[1:nPreds], X_stay[i, 1:nPreds]))
            }
          }
        }

        for(j in 1:N_station) {
          y[j, 1:N_group] ~ ddirchmulti(alpha_dir[1:N_group], N_judge[j])
        }
        for(i in 1:N_group) {
          alpha_dir[i] ~ dexp(1 / N_group)
        }
        alpha_sum <- sum(alpha_dir[1:N_group])
        for(i in 1:N_group) {
          p_expected[i] <- alpha_dir[i] / alpha_sum
          c_expected[i] <- p_expected[i] * (i - 1)
        }
        e_y <- sum(c_expected[1:N_group])
      } #end
    )

    cat("Compiling the model. This may take a moment...\n")

    # Define the von Mises distribution function using nimbleFunction
    # ddirchmulti <- nimbleFunction(
    #   run = function(x = double(1),
    #                  alpha = double(1),
    #                  size = double(0),
    #                  log = integer(0)) {
    #     returnType(double(0))
    #     logProb <- lgamma(size + 1) - sum(lgamma(x + 1)) + lgamma(sum(alpha)) -
    #       sum(lgamma(alpha)) + sum(lgamma(alpha + x)) -
    #       lgamma(sum(alpha) + size)
    #     if (log)
    #       return(logProb)
    #     else
    #       return(exp(logProb))
    #   }
    # )
    #
    # rdirchmulti <- nimbleFunction(
    #   run = function(n = integer(0),
    #                  alpha = double(1),
    #                  size = double(0)) {
    #     returnType(double(1))
    #     if (n != 1)
    #       print("rdirchmulti only allows n = 1; using n = 1.")
    #     p <- rdirch(1, alpha)
    #     return(rmulti(1, size = size, prob = p))
    #   }
    # )
    # assign("rdirchmulti", rdirchmulti, envir = .GlobalEnv)
    # assign("ddirchmulti", ddirchmulti, envir = .GlobalEnv)

    init.t <- ifelse(censored == 0, NA, c_time + 1)

    inits <- list(beta = rep(1.0, nPreds),
                  theta = 0.1, t = init.t,
                  alpha_dir = runif(N_group, 0, 1))
    if (!is.null(random_effect)) {
      inits$u <- rep(0, nLevels)
      inits$sigma_u <- 1
    }
    list(
      log_scale = runif(1, 0.5, 4.0),
      theta = runif(1, 0.5, 4.0),
      alpha_dir = runif(N_group, 0, 1)
    )


    model_mcmc <- suppressMessages(nimbleModel(code, constants = cons_pre, data = data_pre, inits = inits))

    # クラスター初期化
    this_cluster <- makeCluster(nc)

    # parallel::clusterEvalQ(this_cluster, {
    #   library(nimble)
    # })

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

    # MCMC実行関数
    run_MCMC_pre <- function(info, data, constants, code, params, ni, nt, nb) {
      myModel <- nimbleModel(code = code, data = data, constants = constants, inits = info$inits)
      CmyModel <- compileNimble(myModel)
      myMCMC <- buildMCMC(CmyModel, monitors = params)
      CmyMCMC <- compileNimble(myMCMC)
      runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
    }

    # クラスターに関数・変数をエクスポート
    clusterExport(this_cluster,
                  c("ddirchmulti", "rdirchmulti", "registerDistributions", "run_MCMC_pre"),
                  envir = environment())


    # 初期値関数
    inits_f <- function() {
      list(log_scale = runif(1, 0.5, 4.0), theta = runif(1, 0.5, 4.0), alpha_dir = runif(N_group, 0, 1))
    }

    # チェーンごとの初期情報設定
    per_chain_info <- lapply(1:nc, function(i) {
      list(seed = sample(1:9999, 1), inits = inits_f())
    })

    params <- c("mean_stay", "e_y", "log_scale", "theta", "alpha_dir")
    # MCMC実行
    cat("Running MCMC sampling. Please wait...\n")
    out_pre <- parLapply(
      cl = this_cluster,
      X = per_chain_info,
      fun = run_MCMC_pre,
      data = data_pre,
      code = code,
      constants = cons_pre,
      params = params,
      ni = ni,
      nt = nt,
      nb = nb
    )
    stopCluster(this_cluster)
    # this_cluster <- makeCluster(nc)
    #
    # suppressMessages(registerDistributions(list(
    #   ddirchmulti = list(
    #     BUGSdist = "ddirchmulti(alpha, size)",
    #     types = c('value = double(1)', 'alpha = double(1)', 'size = double(0)'),
    #     pqAvail = FALSE
    #   )
    # )))
    #
    # clusterEvalQ(this_cluster, {
    #   library(nimble)
    #   source("R/dirichlet_multinomial.R")
    #   registerDistributions(list(
    #     ddirchmulti = list(
    #       BUGSdist = "ddirchmulti(alpha, size)",
    #       types = c('value = double(1)', 'alpha = double(1)', 'size = double(0)'),
    #       pqAvail = FALSE
    #     )
    #   ))
    # })
    #
    # clusterExport(this_cluster, c("ddirchmulti", "rdirchmulti", "registerDistributions"),
    #               envir = environment())
    #
    #
    # run_MCMC_pre <- function(info, data, constants, code, params, ni, nt, nb) {
    #
    #   #library(nimble)
    #
    #   myModel <- nimbleModel(code = code,
    #                          data = data,
    #                          constants = constants,
    #                          inits = info$inits)
    #
    #   CmyModel <- compileNimble(myModel)
    #
    #   myMCMC <- buildMCMC(CmyModel, monitors = params)
    #   CmyMCMC <- compileNimble(myMCMC)
    #
    #   results <- runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
    #
    #   results
    # }
    #
    # # 初期値関数
    # inits_f <- function() {
    #   list(
    #     log_scale = runif(1, 0.5, 4.0),
    #     theta = runif(1, 0.5, 4.0),
    #     alpha_dir = runif(N_group, 0, 1)
    #   )
    # }
    #
    # # info
    # per_chain_info <- lapply(1:nc, function(i) {
    #   list(
    #     seed = sample(1:9999, 1),
    #     inits = inits_f()
    #   )
    # })
    # params <- c("mean_stay", "e_y", "log_scale", "theta", "alpha_dir")
    #
    # cat("Running MCMC sampling. Please wait...\n")
    #
    # out_pre <- parLapply(cl = this_cluster, X = per_chain_info,
    #                      fun = run_MCMC_pre,
    #                      data = data_pre, code = code,
    #                      constants = cons_pre, params = params,
    #                      ni = ni, nt = nt, nb = nb
    # )
    # stopCluster(this_cluster)

    pattern <- paste0("log_scale\\[", 1:N_stay, "\\]")

    ######WAIC

    if(chains > 1) {
      log_scale.samp <- out_pre %>%
        map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
        do.call(rbind, .)
      if(stay_family != "exponential") tau.samp <- do.call(rbind, out_pre)[, "theta"]
    } else {
      log_scale.samp <- out_pre
      if(stay_family != "exponential") tau.samp <- log_scale.samp[, "theta"]
    }

    if(stay_family != "exponential") tau.samp <- do.call(rbind, out_pre)[, "theta"]

    loglfstay<-matrix(NA, nrow(log_scale.samp), N_stay)

    if(stay_family == "exponential") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- 1 - pexp(c_time[j], rate = 1 /exp(log_scale.samp[ , j]), lower.tail = F, log.p = T)
        } else {
          loglfstay[, j] <- dexp(t[j], rate = 1 / exp(log_scale.samp[ , j]), log = T)
        }
      }
    }
    if(stay_family == "gamma") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- 1 - pgamma(c_time[j], shape = tau.samp,lower.tail = F, scale = exp(log_scale.samp[ , j]), log.p = T)
        } else{
          loglfstay[, j] <- dgamma(t[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), log = T)
        }
      }
    }
    if(stay_family == "lognormal") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- 1 - plnorm(c_time[j], log_scale.samp[ , j], tau.samp,lower.tail = F, log.p = T)
        } else{
          loglfstay[, j] <- dlnorm(t[j], log_scale.samp[ , j], tau.samp, log = T)
        }
      }
    }
    if(stay_family == "weibull") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- 1 - pweibull(c_time[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), lower.tail = F, log.p = T)
        } else{
          loglfstay[, j] <- dweibull(t[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), log = T)
        }
      }
    }

    pattern <- paste0("alpha_dir\\[", 1:N_group, "\\]")

    alpha.samp <- out_pre %>%
      map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
      do.call(rbind, .)

    loglfN <- extraDistr::ddirmnom(y, size = N_detection, alpha = alpha.samp[1:N_group])



    # Density model -----------------------------------------------------------


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

      nPreds <- ncol(X_density)
      N_station <- nrow(station_effort_data)


      # REST model
      Model_REST <- nimbleCode(
        {
          # likelihood
          for(i in 1:N_station){
            N_detection[i] ~ dnbinom(size = r, prob = p[i])
            p[i] <- r / (r + mu[i])
          }
          r ~ dgamma(0.01, 0.01)
          # prior
          if(nPreds == 1) {
            for(i in 1:N_station) {
              log(mu[i]) <- log(density[i]) + log(focal_area) + log(N_period[i]) - log(mean_stay[i]) + log(actv) - log(e_y)
              log(density[i]) <- beta
            }
            beta ~ dnorm(0, sd = 100)
          }
          if(nPreds > 1) {
            for(i in 1:N_station) {
              log(mu[i]) <- log_local_density[i] + log(focal_area) + log(N_period[i]) - log(mean_stay[i]) + log(actv) - log(e_y)
              log_local_density[i] <- inprod(beta[1:nPreds], X_density[i, 1:nPreds]) + eps[i]
              eps[i] ~ dnorm(0, sd_density)
            }

            for(j in 1:nPreds) {
              beta[j] ~ dnorm(0, sd = 100)
            }
            for(i in 1:N_station) {
              density[i] <- exp(inprod(beta[1:nPreds], X_density[i, 1:nPreds]))
            }
          }
          sd_density ~ dunif(0, 5)
        }#end
      )

      run_MCMC_REST <- function(info, data, constants, code, params, ni, nt, nb) {

        # library(nimble)
        myModel <- nimbleModel(code = code,
                               data = data,
                               constants = constants,
                               inits = info$inits)
        CmyModel <- compileNimble(myModel)
        configModel <- configureMCMC(myModel, monitors = params)
        configModel$removeSamplers(c("e_y", "mean_stay"))
        configModel$addSampler(
          target = c(c("e_y", "mean_stay")),
          type = 'prior_samples',
          samples = info$mean_stay_samples
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
      }

      inits_f <- function() {
        list(
          beta = rep(1, nPreds), r = 1
        )
      }
      if (nPreds > 1) {
        inits$sd_density <- 0.1
      }

      # pattern <- c(paste0("mean_stay\\[", 1:N_station, "\\]"), "e_y")
      #
      # out_trace <- out_pre %>%
      #   map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))])

      out_trace <- lapply(coda::mcmc.list(out_pre), function(chain) {
        chain[, grep("^mean_stay\\[", colnames(chain))]
      })


      S <- focal_area * 10^-6
      N_period <- station_effort_data$Effort[1:N_station] * 60 * 60 * 24

      data_REST <- list(N_detection = N_detection, X_density = X_density)
      cons_REST <- list(N_station = N_station, focal_area = S, N_period = N_period, actv = actv, nPreds = nPreds)


      params <- c("density", "beta", "p", "r")

      # info
      per_chain_info <- lapply(1:length(out_trace), function(i) {
        list(
          seed = sample(1:9999, 1),
          inits = inits_f(),
          mean_stay_samples = as.matrix(out_trace[[i]])
        )
      })
      }

      this_cluster <- makeCluster(nc)

      # カスタム分布の登録をクラスター全体で評価
      clusterEvalQ(this_cluster, {
        library(nimble)
      })
      clusterExport(this_cluster,
                    varlist = c("run_MCMC_REST"),
                    envir = environment())

      chain_output <- parLapply(cl = this_cluster, X = per_chain_info,
                                fun = run_MCMC_REST,
                                data = data_REST, code = Model_REST,
                                constants = cons_REST, params = params,
                                ni = ni, nt = nt, nb = nb
      )
      stopCluster(this_cluster)

      loglfy <- matrix(NA, nrow(log_scale.samp), N_station)


      pattern <- paste0("p\\[", 1:N_station, "\\]")
      p.samp <- chain_output %>%
        map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
        do.call(rbind, .)

      r.samp <- tau.samp <- do.call(rbind, out_pre)[, "theta"]
      for(i in 1:N_station) {
        loglfy[,i] <- dnbinom(N_detection[i], prob = p.samp[,i], size = r.samp, log = T)
      }

      loglfall<-cbind(loglfstay, loglfN, loglfy)
      lppd<-sum(log(colMeans(exp(loglfall))))
      p.waic<-sum(apply(loglfall,2,var))
      waic[k] <- (-2)*lppd+2*p.waic
      samples[[k]] <- chain_output

      tidy_samples <- MCMCvis::MCMCchains(chain_output) %>%
        as_tibble() %>%
        pivot_longer(cols = everything(),
                     names_to = "parameter",
                     values_to = "value") %>%
        group_by(parameter) %>%
        mutate(iteration = row_number())


      density_estimates[[k]] <- tidy_samples %>%
        filter(grepl("^density\\[\\d+\\]$", parameter)) %>%
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
    }

  # Original-REST model ----------------------------------------------------------

  if(model == "REST") {

    y <- station_effort_data %>% pull(Y)
    N_station <- length(y)

    data_pre <- list(t = t, censored = censored, y = y)

    cons_pre <- list(N_stay = N_stay, nPreds = nPreds, N_station = N_station, c_time = c_time, stay_family = stay_family,
                     N_station = N_station)

    # Random effects
    if (!is.null(random_effect)) {
      cons_pre$group <- as.numeric(factor(stay_data_join[[random_effect]]))
      cons_pre$nLevels <- nLevels
    } else {
      cons_pre$nLevels <- 0
    }
    # Covariates
    if (nPreds > 1) {
      data_pre$X_stay <- X_stay
    }

    code <- nimbleCode(
      {
        for(i in 1:N_stay) {
          censored[i] ~ dinterval(t[i], c_time[i])
          if(stay_family == "exponential") {
            t[i] ~ dexp(scale = exp(log_scale[i]))
          }
          if(stay_family == "gamma") {
            t[i] ~ dgamma(shape = theta, scale = exp(log_scale[i]))
          }
          if(stay_family == "lognormal") {
            t[i] ~ dlnorm(meanlog = log_scale[i], sdlog = theta)
          }
          if(stay_family == "weibull") {
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
          beta[j] ~ T(dnorm(0, 100), -5, 5)
        }

        if (nLevels > 0) {
          for(k in 1:nLevels) {
            u[k] ~ dnorm(0, sd = sigma_u)
          }
          sigma_u ~ dt(0, 1 / pow(2.5, -2), 1)
        }

        if (nPreds == 1) {
          for(i in 1:N_station){
            if(stay_family == "exponential") mean_stay[i] <- exp(beta[1])
            if(stay_family == "gamma") mean_stay[i] <- theta * exp(beta[1])
            if(stay_family == "lognormal") mean_stay[i] <- exp(beta[1] + theta ^ 2 / 2)
            if(stay_family == "weibull") mean_stay[i] <-  lgamma(1 + 1 / theta) + exp(beta[1])
          }
        }

        if (nPreds > 1) {
          if(stay_family == "exponential") {
            for(i in 1:N_station){
              mean_stay[i] <- exp(inprod(beta[1:nPreds], X_stay[i, 1:nPreds]))
            }
          }
          if(stay_family == "gamma") {
            for(i in 1:N_station){
              mean_stay[i] <- theta * exp(inprod(beta[1:nPreds], X_stay[i, 1:nPreds]))
            }
          }
          if(stay_family == "lognormal") {
            for(i in 1:N_station){
              mean_stay[i] <- exp(inprod(beta[1:nPreds], X_stay[i, 1:nPreds]) + theta ^ 2 / 2)
            }
          }
          if(stay_family == "weibull") {
            for(i in 1:N_station){
              mean_stay[i] <- lgamma(1 + 1 / theta) + exp(inprod(beta[1:nPreds], X_stay[i, 1:nPreds]))
            }
          }
        }
      } #end
    )

    init.t <- ifelse(censored == 0, NA, c_time + 1)
    inits <- list(beta = rep(1.0, nPreds), theta = 0.1, t = init.t)

    if (!is.null(random_effect)) {
      inits$u <- rep(0, nLevels)
      inits$sigma_u <- 1
    }

    cat("Compiling the model. This may take a moment...\n")

    model_mcmc <- suppressMessages(nimbleModel(code, constants = cons_pre, data = data_pre, inits = inits))

    this_cluster <- makeCluster(nc)
    run_MCMC_pre <- function(info, data, constants, code, params, ni, nt, nb) {
      # library(nimble)
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
        log_scale = runif(1, 0.5, 4.0),
        theta = runif(1, 0.5, 4.0)
      )
    }

    per_chain_info <- lapply(1:nc, function(i) {
      list(
        seed = sample(1:9999, 1),
        inits = inits_f()
      )
    })
    params <- c("mean_stay", "log_scale", "theta")

    cat("Running MCMC sampling. Please wait...\n")

    clusterEvalQ(this_cluster, {
      library(nimble)
    })

    clusterExport(this_cluster,
                  varlist = c("run_MCMC_pre"),
                  envir = environment())

    out_pre <- parLapply(
      cl = this_cluster,
      X = per_chain_info,
      fun = run_MCMC_pre,
      data = data_pre,
      code = code,
      constants = cons_pre,
      params = params,
      ni = ni,
      nt = nt,
      nb = nb
    )
    stopCluster(this_cluster)

    # tidy_samples <- MCMCvis::MCMCchains(out_pre) %>%
    #   as_tibble() %>%
    #   pivot_longer(cols = everything(),
    #                names_to = "parameter",
    #                values_to = "value") %>%
    #   group_by(parameter) %>%
    #   mutate(iteration = row_number()) %>%
    #   filter(grepl("^mean_stay\\[\\d+\\]$", parameter)) %>%
    #   mutate(index = as.numeric(str_extract(parameter, "\\d+"))) %>%
    #   filter(index <= 100) %>%
    #   group_by(parameter, index) %>%
    #   summarise(
    #     median = median(value),
    #     lower = quantile(value, 0.025),
    #     upper = quantile(value, 0.975),
    #     .groups = "drop"
    #   ) %>%
    #   arrange(index)
    #
    # g <- ggplot(tidy_samples, aes(x = index, y = median)) +
    #   geom_point() +
    #   geom_segment(aes(x = index, xend = index, y = lower, yend = upper)) +
    #   labs(
    #     title = paste("formula =", formula_stay, "(Median and 95% Credible Interval for mean staying time)"),
    #     x = "Station",
    #     y = "MeanStayingTime"
    #   ) +
    #   ylim(0, NA)
    # print(g)


    pattern <- paste0("log_scale\\[", 1:N_stay, "\\]")

    log_scale.samp <- out_pre %>%
      map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
      do.call(rbind, .)

    if(stay_family != "exponential") tau.samp <- do.call(rbind, out_pre)[, "theta"]

    loglfstay<-matrix(NA, nrow(log_scale.samp), N_stay)

    if(stay_family == "exponential") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- 1 - pexp(c_time[j], rate = 1 /exp(log_scale.samp[ , j]), lower.tail = F, log.p = T)
        } else {
          loglfstay[, j] <- dexp(t[j], rate = 1 / exp(log_scale.samp[ , j]), log = T)
        }
      }
    }
    if(stay_family == "gamma") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- 1 - pgamma(c_time[j], shape = tau.samp,lower.tail = F, scale = exp(log_scale.samp[ , j]), log.p = T)
        } else{
          loglfstay[, j] <- dgamma(t[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), log = T)
        }
      }
    }
    if(stay_family == "lognormal") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- 1 - plnorm(c_time[j], log_scale.samp[ , j], tau.samp,lower.tail = F, log.p = T)
        } else{
          loglfstay[, j] <- dlnorm(t[j], log_scale.samp[ , j], tau.samp, log = T)
        }
      }
    }
    if(stay_family == "weibull") {
      for(j in 1:N_stay){
        if(is.na(t[j])){
          loglfstay[, j] <- 1 - pweibull(c_time[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), lower.tail = F, log.p = T)
        } else{
          loglfstay[, j] <- dweibull(t[j], shape = tau.samp, scale = exp(log_scale.samp[ , j]), log = T)
        }
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

      N <- nrow(X_density)
      nPreds <- ncol(X_density)
      N_station <- nrow(station_effort_data)

      Model_REST <- nimbleCode(
        {
          # mean_stay

          # likelihood
          for(i in 1:N_station){
            y[i] ~ dnbinom(size = r, prob = p[i])
            p[i] <- r / (r + mu[i])
          }
          r ~ dgamma(0.01, 0.01)
          # prior
          if(nPreds == 1) {
            for(i in 1:N_station) {
              log(mu[i]) <- log(density[i]) + log(focal_area) + log(N_period[i]) - log(mean_stay[i]) + log(actv)
              log(density[i]) <- beta
            }
            beta ~ dnorm(0, sd = 100)
          }
          if(nPreds > 1) {
            for(i in 1:N_station) {
              log(mu[i]) <- log_local_density[i] + log(focal_area) + log(N_period[i]) - log(mean_stay[i]) + log(actv)
              log_local_density[i] <- inprod(beta[1:nPreds], X_density[i, 1:nPreds]) + eps[i]
              eps[i] ~ dnorm(0, sd_density)
            }
            for(j in 1:nPreds) {
              beta[j] ~ dnorm(0, sd = 100)
            }
            for(i in 1:N_station) {
              density[i] <- exp(inprod(beta[1:nPreds], X_density[i, 1:nPreds]))
            }
          }
          sd_density ~ dunif(0, 5)
        }#end
      )

      run_MCMC_REST <- function(info, data, constants, code, params, ni, nt, nb) {

        # library(nimble)
        myModel <- nimbleModel(code = code,
                               data = data,
                               constants = constants,
                               inits = info$inits)

        CmyModel <- compileNimble(myModel)

        configModel <- configureMCMC(myModel, monitors = params)
        configModel$removeSamplers(c("mean_stay"))
        configModel$addSampler(
          target = c(c("mean_stay")),
          type = 'prior_samples',
          samples = info$mean_stay_samples
        )

        myMCMC <- buildMCMC(configModel, monitors = params) # configModelを入れる
        CmyMCMC <- compileNimble(myMCMC)

        results <- runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt,nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
      }

      inits_f <- function() {
        list(
          beta = rep(1, nPreds), r = 1
        )
      }
      if (nPreds > 1) {
        inits$sd_density <- 0.1
      }


      #out_pre <- mcmc.list(chain_output)

      # pattern <- c(paste0("mean_stay\\[", 1:N_station, "\\]"))
      # out_trace <- out_pre %>%
      #   map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))])

      out_trace <- lapply(coda::mcmc.list(out_pre), function(chain) {
        chain[, grep("^mean_stay\\[", colnames(chain))]
      })

      per_chain_info <- lapply(1:length(out_trace), function(i) {
        list(
          seed = sample(1:9999, 1),
          inits = inits_f(),
          mean_stay_samples = as.matrix(out_trace[[i]])
        )
      })

      S <- focal_area * 10^-6
      N_period <- station_effort_data$Effort[1:N_station] * 60 * 60 * 24

      data_REST <- list(y = y, X_density = X_density)
      cons_REST <- list(N_station = N_station, focal_area = S, N_period = N_period, actv = actv, nPreds = nPreds)


      params <- c("density", "beta", "p", "r")

      cat("Running MCMC sampling. Please wait...\n")

      this_cluster <- makeCluster(nc)
      clusterEvalQ(this_cluster, {
        library(nimble)
      })

      chain_output <- parLapply(cl = this_cluster, X = per_chain_info,
                                fun = run_MCMC_REST,
                                data = data_REST, code = Model_REST,
                                constants = cons_REST, params = params,
                                ni = ni, nt = nt, nb = nb
      )
      stopCluster(this_cluster)


      ######WAIC
      loglfy <- matrix(NA, nrow(log_scale.samp), N_station)


      pattern <- paste0("p\\[", 1:N_station, "\\]")
      p.samp <- chain_output %>%
        map(~ .[, grep(paste(pattern, collapse = "|"), colnames(.))]) %>%
        do.call(rbind, .)

      r.samp <- tau.samp <- do.call(rbind, out_pre)[, "theta"]
      for(i in 1:N_station) {
        loglfy[,i] <- dnbinom(y[i], prob = p.samp[,i], size = r.samp, log = T)
      }

      loglfall<-cbind(loglfstay, loglfy)
      lppd<-sum(log(colMeans(exp(loglfall))))
      p.waic<-sum(apply(loglfall,2,var))
      waic[k] <- (-2)*lppd+2*p.waic
      samples[[k]] <- chain_output

      tidy_samples <- MCMCvis::MCMCchains(chain_output) %>%
        as_tibble() %>%
        pivot_longer(cols = everything(),
                     names_to = "parameter",
                     values_to = "value") %>%
        group_by(parameter) %>%
        mutate(iteration = row_number())


      density_estimates[[k]] <- tidy_samples %>%
        filter(grepl("^density\\[\\d+\\]$", parameter)) %>%
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
    }
  }

  model_summary <- list(0)
  model_summary$WAIC <- data.frame(Model = as.character(unlist(formula_density_all)), WAIC = waic) %>%
    arrange(WAIC)
  model_summary$density <- density_estimates[[which(model_summary$WAIC[1, 1]  == unlist(formula_density_all))]]
  model_summary$all_samples <- samples[[which(model_summary$WAIC[1, 1]  == unlist(formula_density_all))]]

  return(model_summary)
}


