utils::globalVariables(c("time", "Stay", "Cens", "Species", "Effort", "Station"))

#' Bayesian parameter estimation for the REST/RAD-REST model via MCMC sampling using `nimble`.
#'
#' @param formula_stay A model formula for staying times within the focal area. For example, \code{Stay ~ 1 + x1}. The left-hand side must specify the column name for staying time.
#' @param formula_density A model formula for animal density. For example, \code{~ 1 + x1}. Note that the left-hand side must be omitted as density is a latent parameter.
#' @param formula_enter A model formula for the total number of passes (\code{Y}) through the focal area for each station-species pair. For example, \code{~ 1 + x1}. Note that this argument is strictly used for the RAD-REST model and is completely ignored when \code{model = "REST"}.
#' @param station_effort_data A data frame containing information for each camera station. Typically, this is the output of the \code{add_effort} function. Alternatively, a manually prepared data frame may be provided, with the required columns depending on the model:
#'   \itemize{
#'     \item If \code{model = "REST"}, the data frame must contain:
#'       \itemize{
#'         \item \code{Station} (character): Camera station ID.
#'         \item \code{Effort} (numeric): Camera trapping effort (days) for each station.
#'         \item \code{Species} (character): Species name.
#'         \item \code{Y} (numeric): Total number of passes through the focal area for each station-species pair.
#'       }
#'     \item If \code{model = "RAD-REST"}, the data frame must contain:
#'       \itemize{
#'         \item \code{Station} (character): Camera station ID.
#'         \item \code{Effort} (numeric): Camera trapping effort (days) for each station.
#'         \item \code{Species} (character): Species name.
#'         \item \code{N} (numeric): Total number of detected videos for each station-species pair.
#'         \item \code{y_X} columns (\code{y_0}, \code{y_1}, ..., \code{y_max}): Number of videos classified by the number of passes observed.
#'       }
#'   }
#' @param stay_data A data frame returned by the \code{format_stay} function, containing the following columns:
#'   \itemize{
#'     \item \code{Station} (character): Camera station ID.
#'     \item \code{Species} (character): Species name.
#'     \item \code{Stay} (numeric): Staying time (in seconds) within the focal area for each detected pass.
#'     \item \code{Cens} (binary): Indicator for censored staying time (1 = censored, 0 = observed).
#'   }
#' @param random_effect_stay A specification of the random effect structure on staying time (e.g., \code{~ (1 | Station)}). Default is \code{NULL} (no random effects).
#' @param activity_data A data frame containing a \code{time} column, representing detection times transformed into radians. Typically, this is the output of the \code{format_activity} function.
#' @param activity_estimation A character string specifying the method used to estimate activity patterns. Choose \code{"kernel"} for fixed kernel density estimation (Rowcliffe et al. 2014) or \code{"mixture"} for nonparametric Bayesian estimation using von Mises mixture models (Nakashima et al. 2025). Default is \code{"kernel"}.
#' @param C The maximum number of von Mises components used for the mixture model. Required only if \code{activity_estimation = "mixture"}. Default is 10.
#' @param bw_adj A numeric bandwidth adjustment parameter for kernel density estimation. Default is 1.0. See Rowcliffe et al. (2014) for details.
#' @param stay_family A character string specifying the probability distribution used for modeling staying times (e.g., \code{"exponential"}, \code{"gamma"}, \code{"lognormal"}, \code{"weibull"}). This should ideally be selected based on the output of the \code{bayes_stay_selection} function. Default is \code{"lognormal"}.
#' @param focal_area A numeric value representing the area of the camera focal area, in square meters.
#' @param cores An integer specifying the number of CPU cores to use for parallel computation. Default is 3.
#' @param iter An integer specifying the total number of MCMC iterations per chain. Default is 5000.
#' @param warmup An integer specifying the number of warm-up (burn-in) iterations per chain. Default is 1000.
#' @param chains An integer specifying the number of MCMC chains to run. Default is 3.
#' @param thin An integer specifying the thinning interval for MCMC sampling. Default is 2.
#' @param model A character string specifying the model to be used. Choose either \code{"REST"} or \code{"RAD-REST"}. Default is \code{"REST"}.
#' @param all_comb Logical. If \code{TRUE}, all possible combinations of covariates in the density model will be evaluated and compared. If \code{FALSE}, only the specific model provided in \code{formula_density} will be run. Default is \code{FALSE}.
#' @param target_species A character string specifying the target species to be analyzed. Only a single species can be specified at a time.
#' @return A list of class \code{"ResultDensity"}, which includes the following components:
#' \describe{
#'   \item{\code{WAIC}}{An object containing WAIC (Widely Applicable Information Criterion) results for model evaluation.}
#'   \item{\code{summary_result}}{A data frame summarizing the posterior estimates of the parameters, including the mean staying time and density.}
#'   \item{\code{samples}}{A \code{coda::mcmc.list} object containing full MCMC posterior samples for all model parameters.}
#' }
#' The returned object has a custom print method that automatically displays the WAIC and a summary of the estimates.
#' You can access the full MCMC samples via \code{$samples} and visually analyze convergence using the \code{MCMCvis} package.
#'
#' @export
#' @import nimble activity parallel MCMCvis
#' @importFrom dplyr filter pull select mutate rename arrange starts_with all_of
#' @importFrom stats as.formula formula model.frame model.matrix sd var runif median quantile model.response rexp rnorm step dexp pexp dgamma pgamma dlnorm plnorm dweibull pweibull dnbinom na.pass delete.response terms
#' @importFrom extraDistr ddirmnom
#' @importFrom stringr str_remove str_extract str_starts
#' @importFrom tibble rownames_to_column as_tibble
#' @examples
#' \dontrun{
#' station_data_rest <- format_station_data(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_y = "y",
#'   model = "REST"
#' )
#'
#' station_effort_REST <- add_effort(
#'   detection_data = detection_data,
#'   station_data_formatted = station_data_rest,
#'   col_name_station = "Station",
#'   col_name_term = "Term",
#'   col_name_datetime = "DateTime",
#'   plot = TRUE
#' )
#'
#' stay_data <- format_stay(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_stay = "Stay",
#'   col_name_cens = "Cens"
#' )
#'
#' activity_data <- format_activity(
#'   detection_data = detection_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_datetime = "DateTime",
#'   indep_time = 30
#' )
#'
#' rest_model <- bayes_rest(
#'   formula_stay = Stay ~ 1 + x1,
#'   formula_density = ~ 1,
#'   formula_enter = ~ 1,
#'   station_effort_data = station_effort_REST,
#'   stay_data = stay_data,
#'   activity_data = activity_data,
#'   activity_estimation = "kernel",
#'   bw_adj = 1.5,
#'   stay_family = "lognormal",
#'   focal_area = 1.96,
#'   cores = 3,
#'   iter = 5000,
#'   warmup = 1000,
#'   chains = 3,
#'   thin = 2,
#'   model = "REST",
#'   all_comb = FALSE,
#'   target_species = "SP01"
#' )
#' }
bayes_rest_2 <- function(formula_stay,
                       formula_density,
                       formula_enter = ~ 1,   # [FIX 1] デフォルト値を追加（REST時は無視）
                       station_effort_data,
                       stay_data,
                       random_effect_stay = NULL,
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
) {
  ni <- iter
  nt <- thin
  nc <- chains
  nb <- warmup
  nmcmc <- (iter - warmup) * chains / thin

  # ---------------------------------------------------------------------------
  # 入力チェック
  # ---------------------------------------------------------------------------

  if (!stay_family %in% c("lognormal", "gamma", "weibull", "exponential")) {
    stop(paste0("Input stay_family type (", stay_family, ") is incorrect."))
  }

  if (length(target_species) > 1) {
    stop("The 'bayes_rest' function does not support multiple species. Please run the model separately for each species or use 'bayes_rest_multi' for multispecies analysis.")
  }

  if (!activity_estimation %in% c("kernel", "mixture")) {
    stop("`activity_estimation` must be either 'kernel' or 'mixture'.")
  }

  if (!model %in% c("REST", "RAD-REST")) {
    stop("`model` must be either 'REST' or 'RAD-REST'.")
  }

  if (!is.character(target_species) || length(target_species) != 1) {
    stop("`target_species` must be a single character string.")
  }

  if (!is.numeric(focal_area) || length(focal_area) != 1 || focal_area <= 0) {
    stop("`focal_area` must be a positive numeric value representing area in square meters.")
  }

  if (!is.numeric(cores) || cores < 1) {
    stop("`cores` must be a positive integer.")
  }

  if (!is.numeric(iter) || iter <= 0 || length(iter) != 1) {
    stop("`iter` must be a positive integer.")
  }

  if (!is.null(warmup) && (!is.numeric(warmup) || warmup <= 0 || length(warmup) != 1)) {
    stop("`warmup` must be NULL or a positive integer.")
  }

  if (!is.numeric(chains) || chains < 1 || length(chains) != 1) {
    stop("`chains` must be a positive integer.")
  }

  if (!is.numeric(thin) || thin < 1 || length(thin) != 1) {
    stop("`thin` must be a positive integer.")
  }

  if (!is.logical(all_comb) || length(all_comb) != 1) {
    stop("`all_comb` must be either TRUE or FALSE.")
  }

  if (!inherits(formula_stay, "formula")) {
    stop("`formula_stay` must be a valid model formula.")
  }

  if (!inherits(formula_density, "formula")) {
    stop("`formula_density` must be a valid model formula (e.g., `~ 1 + x1`).")
  }

  # [FIX 2] 列名チェックを "Cens" に統一（ドキュメントと実装を一致させる）
  required_cols_stay <- c("Station", "Species", "Stay", "Cens")
  if (!all(required_cols_stay %in% colnames(stay_data))) {
    stop("`stay_data` must contain the following columns: 'Station', 'Species', 'Stay', and 'Cens'.")
  }

  # ---------------------------------------------------------------------------
  # データ準備
  # ---------------------------------------------------------------------------

  act_data <- activity_data %>%
    dplyr::filter(Species == target_species) %>%
    dplyr::pull(time)

  if (activity_estimation == "kernel") {
    model_act <- activity::fitact(act_data,
                                  bw = bw_adj * activity::bwcalc(act_data, K = 3),
                                  reps = 1)
    activity_proportion <- model_act@act
  } else {
    dens.x <- seq(0, 2 * pi, 0.02)
  }

  station_effort_data <- station_effort_data %>%
    dplyr::filter(Species == target_species)

  stay_data_join <- stay_data %>%
    dplyr::filter(Species == target_species) %>%
    dplyr::left_join(station_effort_data,
                     by = intersect(names(stay_data), names(station_effort_data))) %>%
    dplyr::filter(!is.na(Stay))

  N_station  <- nrow(station_effort_data)
  station.id <- unique(station_effort_data$Station)
  stay       <- stay_data_join %>% dplyr::pull(Stay)
  censored   <- stay_data_join %>% dplyr::pull(Cens)   # [FIX 2] "Cens" に統一

  c_time <- stay
  c_time[censored == 0] <- c_time[censored == 0] + 1
  stay[censored == 1] <- NA
  N_stay <- length(stay)

  S        <- focal_area * 10^-6
  N_period <- station_effort_data$Effort[1:N_station] * 60 * 60 * 24

  # stay 共変量
  vars_stay      <- all.vars(formula_stay)
  predictors_stay <- vars_stay[-1]
  model_frame_stay <- stats::model.frame(formula_stay, stay_data_join)
  X_stay          <- stats::model.matrix(stats::as.formula(formula_stay), model_frame_stay)
  nPreds_stay     <- ncol(X_stay)

  # ランダム効果
  if (!is.null(random_effect_stay)) {
    levels_re    <- unique(stay_data_join[[random_effect_stay]])
    nLevels_stay <- length(levels_re)
  } else {
    nLevels_stay <- 0
  }

  # 密度共変量候補
  predictors_density  <- all.vars(formula_density)
  formula_density_all <- list(0)
  if (all_comb) {
    formula_density_all <- full_terms(c(predictors_density))
  } else {
    formula_density_all[[1]] <- formula_density
  }

  # ---------------------------------------------------------------------------
  # 活動割合の推定（mixture モデル）
  # ---------------------------------------------------------------------------

  if (activity_estimation == "mixture") {
    dens.x <- seq(0, 2 * pi, 0.02)
    ndens  <- length(dens.x)
    N      <- length(act_data)

    constants_mix <- list(N = N, C = C, dens.x = dens.x, ndens = ndens)
    data_mix      <- list(act_data = act_data)

    code_mix <- nimble::nimbleCode({
      for (k in 1:(C - 1)) {
        v[k] ~ dbeta(1, alpha)
      }
      alpha ~ dgamma(1, 1)
      w[1:C] <- stick_breaking(v[1:(C - 1)])
      for (k in 1:C) {
        mu_mix[k]    ~ dunif(0, 2 * 3.141592654)
        kappa_mix[k] ~ dgamma(1, 0.01)
      }
      for (n in 1:N) {
        group[n]    ~ dcat(w[1:C])
        act_data[n] ~ dvonMises(mu_mix[group[n]], kappa_mix[group[n]])
        act_data_pred[n] ~ dvonMises(mu_mix[group[n]], kappa_mix[group[n]])
        loglike_obs_act[n] <- dvonMises(act_data[n],      mu_mix[group[n]], kappa_mix[group[n]], log = 1)
        loglike_pred_act[n] <- dvonMises(act_data_pred[n], mu_mix[group[n]], kappa_mix[group[n]], log = 1)
      }
      for (j in 1:ndens) {
        for (i in 1:C) {
          dens.cpt[i, j] <- w[i] * dvonMises(dens.x[j], mu_mix[i], kappa_mix[i], log = 0)
        }
        activity_density[j] <- sum(dens.cpt[1:C, j])
      }
      activity_proportion <- 1.0 / (2 * 3.141592654 * max(activity_density[1:ndens]))
    })

    inits_f_mix <- function() {
      list(
        mu_mix    = stats::runif(C, 0, 2 * pi),
        kappa_mix = stats::rgamma(C, 1, 0.01),
        group     = sample(1:C, size = N, replace = TRUE),
        v         = stats::rbeta(C - 1, 1, 1),
        alpha     = 1
      )
    }

    # von Mises カスタム分布（ワーカー共有用）
    dvonMises_mix <- nimble::nimbleFunction(
      run = function(x = double(0), kappa = double(0), mu = double(0), log = integer(0)) {
        returnType(double(0))
        ccrit <- 1E-6; s <- 1; i <- 1; inc <- 1; satisfied <- FALSE
        while (!satisfied) {
          x_2i <- kappa / (2 * i)
          inc  <- inc * x_2i * x_2i
          s    <- s + inc
          i    <- i + 1
          satisfied <- inc < ccrit
        }
        prob <- exp(kappa * cos(x - mu)) / (2 * pi * s)
        if (log) return(log(prob)) else return(prob)
      }
    )
    rvonMises_mix <- nimble::nimbleFunction(
      run = function(n = integer(0), kappa = double(0), mu = double(0)) {
        returnType(double(0))
        return(0)
      }
    )

    suppressMessages(nimble::registerDistributions(list(
      dvonMises = list(
        BUGSdist = "dvonMises(kappa, mu)",
        types    = c("value = double(0)", "kappa = double(0)", "mu = double(0)"),
        pqAvail  = FALSE
      )
    )))

    run_MCMC_vonMises <- function(info, data, constants, code, params, ni, nt, nb) {
      myModel    <- nimble::nimbleModel(code = code, data = data,
                                        constants = constants, inits = info$inits)
      CmyModel   <- nimble::compileNimble(myModel)
      configModel <- nimble::configureMCMC(myModel, monitors = params)
      myMCMC     <- nimble::buildMCMC(configModel, monitors = params)
      CmyMCMC    <- nimble::compileNimble(myMCMC)
      nimble::runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt,
                      nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
    }

    per_chain_info_mix <- lapply(1:nc, function(i) {
      list(seed = sample(1:9999, 1), inits = inits_f_mix())
    })

    params_mix <- c("activity_density", "activity_proportion",
                    "mu_mix", "kappa_mix", "w", "loglike_obs_act")

    cat("Running MCMC sampling for activity (mixture model). Please wait...\n")

    this_cluster <- parallel::makeCluster(nc)
    on.exit(try(parallel::stopCluster(this_cluster), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(this_cluster, { library(nimble) })
    parallel::clusterExport(this_cluster,
                            c("dvonMises_mix", "rvonMises_mix", "run_MCMC_vonMises"),
                            envir = environment())

    actv_chain_output <- parallel::parLapply(
      cl = this_cluster, X = per_chain_info_mix,
      fun = run_MCMC_vonMises,
      data = data_mix, code = code_mix,
      constants = constants_mix, params = params_mix,
      ni = ni, nt = nt, nb = nb
    )
    parallel::stopCluster(this_cluster)
    on.exit()

    actv_out_trace <- lapply(actv_chain_output, function(ch) {
      ch[, grep("activity_proportion", colnames(ch)), drop = FALSE]
    })

    summary_mix <- MCMCvis::MCMCsummary(actv_chain_output, round = 3)

    activity_density_estimates <- summary_mix %>%
      as.data.frame() %>%
      tibble::rownames_to_column("variable") %>%
      dplyr::filter(stringr::str_starts(variable, "activity_density")) %>%
      dplyr::mutate(x = dens.x)

    loglact <- MCMCvis::MCMCchains(actv_chain_output, params = c("loglike_obs_act"))
  }

  # ===========================================================================
  # REST モデル
  # ===========================================================================

  if (model == "REST") {

    tidy_samples  <- list()
    mcmc_samples  <- list()
    waic          <- numeric(0)

    # von Mises カスタム分布（REST 並列用）
    dvonMises <- nimble::nimbleFunction(
      run = function(x = double(0), kappa = double(0), mu = double(0), log = integer(0)) {
        returnType(double(0))
        ccrit <- 1E-6; s <- 1; i <- 1; inc <- 1; satisfied <- FALSE
        while (!satisfied) {
          x_2i <- kappa / (2 * i)
          inc  <- inc * x_2i * x_2i
          s    <- s + inc
          i    <- i + 1
          satisfied <- inc < ccrit
        }
        prob <- exp(kappa * cos(x - mu)) / (2 * pi * s)
        if (log) return(log(prob)) else return(prob)
      }
    )
    rvonMises <- nimble::nimbleFunction(
      run = function(n = integer(0), kappa = double(0), mu = double(0)) {
        returnType(double(0))
        return(0)
      }
    )

    for (k in seq_along(formula_density_all)) {

      formula_density_k  <- formula_density_all[[k]]
      model_frame_density <- stats::model.frame(formula_density_k, data = station_effort_data)
      X_density           <- stats::model.matrix(stats::as.formula(formula_density_k), model_frame_density)
      nPreds_density      <- ncol(X_density)

      y <- station_effort_data$Y

      data_REST <- list(
        y        = y,
        X_density = X_density,
        stay     = stay,
        censored = censored
      )

      cons_REST <- list(
        N_station          = N_station,
        S                  = S,
        N_period           = N_period,
        nPreds_stay        = nPreds_stay,
        N_stay             = N_stay,
        c_time             = c_time,
        nPreds_density     = nPreds_density,
        stay_family        = stay_family,
        activity_estimation = activity_estimation,
        nLevels_stay       = nLevels_stay
      )

      if (!is.null(random_effect_stay)) {
        cons_REST$group_stay <- as.numeric(factor(stay_data_join[[random_effect_stay]]))
      }
      if (activity_estimation == "kernel") {
        cons_REST$activity_proportion <- activity_proportion
      }
      if (nPreds_stay > 1) {
        data_REST$X_stay <- X_stay
      }

      Model_REST <- nimble::nimbleCode({

        # --- 滞在時間 ---
        for (i in 1:N_stay) {
          censored[i] ~ dinterval(stay[i], c_time[i])

          if (stay_family == "exponential") {
            stay[i]   ~ dexp(rate = 1 / scale[i])
            pred_t[i] ~ dexp(rate = 1 / scale[i])
            loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dexp(stay[i], rate = 1 / scale[i], log = 1) +
              step(censored[i] - 0.5) * log(1 - pexp(c_time[i], rate = 1 / scale[i]))
            loglike_pred_stay[i] <- dexp(pred_t[i], rate = 1 / scale[i], log = 1)
          }
          if (stay_family == "gamma") {
            stay[i]   ~ dgamma(shape = theta_stay, rate = exp(-log(scale[i])))
            pred_t[i] ~ dgamma(shape = theta_stay, rate = exp(-log(scale[i])))
            loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dgamma(stay[i], shape = theta_stay, rate = exp(-log(scale[i])), log = 1) +
              step(censored[i] - 0.5) * log(1 - pgamma(c_time[i], shape = theta_stay, rate = exp(-log(scale[i]))))
            loglike_pred_stay[i] <- dgamma(pred_t[i], shape = theta_stay, rate = exp(-log(scale[i])), log = 1)
          }
          if (stay_family == "lognormal") {
            stay[i]   ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay)
            pred_t[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay)
            loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dlnorm(stay[i], meanlog = log(scale[i]), sdlog = theta_stay, log = 1) +
              step(censored[i] - 0.5) * log(1 - plnorm(c_time[i], meanlog = log(scale[i]), sdlog = theta_stay))
            loglike_pred_stay[i] <- dlnorm(pred_t[i], meanlog = log(scale[i]), sdlog = theta_stay, log = 1)
            meanlog[i]           <- log(scale[i])
          }
          if (stay_family == "weibull") {
            stay[i]   ~ dweibull(shape = theta_stay, scale = scale[i])
            pred_t[i] ~ dweibull(shape = theta_stay, scale = scale[i])
            loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dweibull(stay[i], shape = theta_stay, scale = scale[i], log = 1) +
              step(censored[i] - 0.5) * log(1 - pweibull(c_time[i], shape = theta_stay, scale = scale[i]))
            loglike_pred_stay[i] <- dweibull(pred_t[i], shape = theta_stay, scale = scale[i], log = 1)
          }

          if (nPreds_stay > 1) {
            if (nLevels_stay == 0) {
              log(scale[i]) <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay])
            } else {
              log(scale[i]) <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]]
            }
          } else {
            if (nLevels_stay == 0) {
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

        for (j in 1:nPreds_stay) {
          beta_stay[j] ~ T(dnorm(0, 100), -10, 10)
        }

        if (nLevels_stay > 0) {
          for (lv in 1:nLevels_stay) {
            random_effect_stay[lv] ~ dnorm(0, sd = sigma_stay)
          }
          sigma_stay ~ T(dnorm(0, 100), 0, 10)
        }

        # --- 期待値計算（mean_stay）---
        # [FIX 3] Weibull の平均は scale * Gamma(1 + 1/shape) → 乗算に修正
        if (nPreds_stay == 1) {
          if (stay_family == "exponential") mean_stay <- exp(beta_stay[1])
          if (stay_family == "gamma")       mean_stay <- theta_stay * exp(beta_stay[1])
          if (stay_family == "lognormal")   mean_stay <- exp(beta_stay[1] + theta_stay^2 / 2)
          if (stay_family == "weibull")     mean_stay <- exp(lgamma(1 + 1 / theta_stay) + beta_stay[1])
        }
        if (nPreds_stay > 1) {
          if (nLevels_stay == 0) {
            for (i in 1:N_station) {
              if (stay_family == "exponential") mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
              if (stay_family == "gamma")       mean_stay[i] <- theta_stay * exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
              if (stay_family == "lognormal")   mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + theta_stay^2 / 2)
              if (stay_family == "weibull")     mean_stay[i] <- exp(lgamma(1 + 1 / theta_stay) + inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
            }
          }
          if (nLevels_stay > 0) {
            for (i in 1:N_station) {
              if (stay_family == "exponential") mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]])
              if (stay_family == "gamma")       mean_stay[i] <- theta_stay * exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]])
              if (stay_family == "lognormal")   mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]] + theta_stay^2 / 2)
              if (stay_family == "weibull")     mean_stay[i] <- exp(lgamma(1 + 1 / theta_stay) + inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]])
            }
          }
        }

        # --- 検出回数（y）の尤度 ---
        for (i in 1:N_station) {
          y[i]     ~ dnbinom(size = size, prob = p[i])
          p[i]     <- size / (size + mu[i])
          y_rep[i] ~ dnbinom(size = size, prob = p[i])
          loglike_obs_y[i]  <- dnbinom(y[i],     size, p[i], log = 1)
          loglike_pred_y[i] <- dnbinom(y_rep[i], size, p[i], log = 1)
        }
        size ~ dgamma(1, 1)

        # --- REST 式（mu の計算）---
        # [FIX 5] nPreds_density > 1 のブロック内のスカラー density 定義を削除し、
        #         density[i] をベクトルとして一貫して扱う
        if (nPreds_density == 1) {
          for (i in 1:N_station) {
            if (nPreds_stay == 1) {
              log(mu[i]) <- log(density) + log(S) + log(N_period[i]) - log(mean_stay) + log(activity_proportion)
            } else {
              log(mu[i]) <- log(density) + log(S) + log(N_period[i]) - log(mean_stay[i]) + log(activity_proportion)
            }
          }
          log(density) <- beta_density[1]
          beta_density[1] ~ dnorm(0, sd = 100)
        }
        if (nPreds_density > 1) {
          for (i in 1:N_station) {
            log(density[i]) <- inprod(beta_density[1:nPreds_density], X_density[i, 1:nPreds_density])
            if (nPreds_stay == 1) {
              log(mu[i]) <- log(density[i]) + log(S) + log(N_period[i]) - log(mean_stay) + log(activity_proportion)
            } else {
              log(mu[i]) <- log(density[i]) + log(S) + log(N_period[i]) - log(mean_stay[i]) + log(activity_proportion)
            }
          }
          for (j in 1:nPreds_density) {
            beta_density[j] ~ dnorm(0, sd = 100)
          }
        }
      })

      inits_f <- function() {
        common_inits <- list(
          beta_stay    = stats::runif(nPreds_stay, -1, 1),
          stay         = ifelse(censored == 0, NA, c_time + stats::runif(N_stay, 0.1, 1.0)),
          theta_stay   = stats::runif(1, 0.5, 1.5),
          scale        = stats::runif(N_stay, 0, 2),
          beta_density = stats::rnorm(nPreds_density, 0, 1),
          size         = stats::runif(1, 0.5, 2)
        )
        if (!is.null(random_effect_stay)) {
          common_inits$random_effect_stay <- stats::runif(nLevels_stay, -1, 1)
          common_inits$sigma_stay         <- stats::runif(1, 0.5, 2.5)
        }
        return(common_inits)
      }

      # 並列 MCMC 実行関数
      run_MCMC_REST <- function(info, data, constants, code, params, ni, nt, nb) {
        suppressMessages(nimble::registerDistributions(list(
          dvonMises = list(
            BUGSdist = "dvonMises(kappa, mu)",
            types    = c("value = double(0)", "kappa = double(0)", "mu = double(0)"),
            pqAvail  = FALSE
          )
        )))

        myModel     <- nimble::nimbleModel(code = code, data = data,
                                           constants = constants, inits = info$inits)
        CmyModel    <- nimble::compileNimble(myModel)
        configModel <- nimble::configureMCMC(myModel, monitors = params)

        # [FIX 12] mixture 時のサンプラー登録を control リスト方式に統一
        if (constants$activity_estimation == "mixture") {
          configModel$removeSampler("activity_proportion")
          configModel$addSampler(
            target  = "activity_proportion",
            type    = "prior_samples",
            control = list(samples = info$actv_samples)
          )
        }

        myMCMC  <- nimble::buildMCMC(configModel, monitors = params)
        CmyMCMC <- nimble::compileNimble(myMCMC, project = myModel)

        nimble::runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt,
                        nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
      }

      if (activity_estimation == "mixture") {
        per_chain_info <- lapply(seq_along(actv_out_trace), function(i) {
          list(seed = sample(1:9999, 1), inits = inits_f(),
               actv_samples = as.matrix(actv_out_trace[[i]]))
        })
        nc_run <- length(actv_out_trace)
      } else {
        per_chain_info <- lapply(1:nc, function(i) {
          list(seed = sample(1:9999, 1), inits = inits_f(), actv_samples = NULL)
        })
        nc_run <- nc
      }

      # パラメータリスト
      if (stay_family == "exponential")              prms <- c("scale", "mean_stay")
      if (stay_family %in% c("gamma", "weibull"))    prms <- c("scale", "shape", "mean_stay")
      if (stay_family == "lognormal")                prms <- c("meanlog", "sdlog", "mean_stay")
      prms   <- c(prms, "density", "p", "size")
      params <- c(prms, "loglike_obs_stay", "loglike_obs_y",
                  "loglike_pred_stay", "loglike_pred_y")
      if (activity_estimation == "mixture") params <- c(params, "activity_proportion")

      cat(sprintf("Running MCMC sampling for REST (model %d/%d). Please wait...\n",
                  k, length(formula_density_all)))

      this_cluster <- parallel::makeCluster(nc_run)
      on.exit(try(parallel::stopCluster(this_cluster), silent = TRUE), add = TRUE)

      parallel::clusterEvalQ(this_cluster, { library(nimble) })
      parallel::clusterExport(this_cluster,
                              c("dvonMises", "rvonMises", "run_MCMC_REST"),
                              envir = environment())

      chain_output <- parallel::parLapply(
        cl = this_cluster, X = per_chain_info,
        fun = run_MCMC_REST, data = data_REST, code = Model_REST,
        constants = cons_REST, params = params,
        ni = ni, nt = nt, nb = nb
      )
      parallel::stopCluster(this_cluster)
      on.exit()
      cat("Estimation is finished!\n")

      # WAIC 計算
      loglfstay <- MCMCvis::MCMCchains(chain_output, params = "loglike_obs_stay")
      loglfy    <- MCMCvis::MCMCchains(chain_output, params = "loglike_obs_y")
      loglfall  <- if (activity_estimation == "mixture") cbind(loglfstay, loglfy, loglact) else cbind(loglfstay, loglfy)
      lppd    <- sum(log(colMeans(exp(loglfall))))
      p.waic  <- sum(apply(loglfall, 2, stats::var))
      waic[k] <- (-2) * lppd + 2 * p.waic

      mcmc_samples[[k]] <- MCMCvis::MCMCchains(chain_output, mcmc.list = TRUE, params = prms)
      samples_mat        <- MCMCvis::MCMCchains(chain_output)
      n_iters   <- nrow(samples_mat)
      p_names   <- colnames(samples_mat)
      tidy_samples[[k]] <- data.frame(
        parameter = rep(p_names, each = n_iters),
        value     = as.vector(samples_mat),
        iteration = rep(1:n_iters, times = length(p_names)),
        stringsAsFactors = FALSE
      )
    }
  }

  # ===========================================================================
  # RAD-REST モデル
  # ===========================================================================

  if (model == "RAD-REST") {

    station_data_target <- station_effort_data %>%
      dplyr::filter(Species == target_species)

    model_frame_enter <- stats::model.frame(formula_enter, data = station_data_target)
    X_enter           <- stats::model.matrix(stats::as.formula(formula_enter), model_frame_enter)
    nPreds_enter      <- ncol(X_enter)

    y_enter_cols  <- grep("^y_", names(station_data_target), value = TRUE)
    y_enter       <- station_data_target %>%
      dplyr::select(dplyr::all_of(y_enter_cols)) %>%
      as.matrix()

    N_enter_group <- ncol(y_enter)
    N_det         <- station_data_target %>% dplyr::pull(N)
    N_enter_judge <- apply(y_enter, 1, sum)

    tidy_samples  <- list()
    mcmc_samples  <- list()
    waic          <- numeric(0)

    y          <- station_effort_data %>% dplyr::select(dplyr::starts_with("y_")) %>% as.matrix()
    N_judge    <- apply(y, 1, sum)
    N_group    <- ncol(y)
    N_detection <- station_effort_data %>% dplyr::pull(N)

    if (nPreds_stay > 1) {
      stay_terms           <- stats::delete.response(stats::terms(stats::as.formula(formula_stay)))
      model_frame_stay_st  <- stats::model.frame(stay_terms, data = station_effort_data,
                                                 na.action = stats::na.pass)
      X_stay_station       <- stats::model.matrix(stay_terms, model_frame_stay_st)
    }

    # カスタム分布定義
    ddirchmulti <- nimble::nimbleFunction(
      run = function(x = double(1), alpha = double(1), size = double(0), log = integer(0)) {
        returnType(double(0))
        logProb <- lgamma(size + 1) - sum(lgamma(x + 1)) + lgamma(sum(alpha)) -
          sum(lgamma(alpha)) + sum(lgamma(alpha + x)) -
          lgamma(sum(alpha) + size)
        if (log) return(logProb) else return(exp(logProb))
      }
    )
    rdirchmulti <- nimble::nimbleFunction(
      run = function(n = integer(0), alpha = double(1), size = double(0)) {
        returnType(double(1))
        if (n != 1) print("rdirchmulti only allows n = 1; using n = 1.")
        p <- nimble::rdirch(1, alpha)
        return(nimble::rmulti(1, size = size, prob = p))
      }
    )
    dvonMises <- nimble::nimbleFunction(
      run = function(x = double(0), kappa = double(0), mu = double(0), log = integer(0)) {
        returnType(double(0))
        ccrit <- 1E-6; s <- 1; i <- 1; inc <- 1; satisfied <- FALSE
        while (!satisfied) {
          x_2i <- kappa / (2 * i)
          inc  <- inc * x_2i * x_2i
          s    <- s + inc
          i    <- i + 1
          satisfied <- inc < ccrit
        }
        prob <- exp(kappa * cos(x - mu)) / (2 * pi * s)
        if (log) return(log(prob)) else return(prob)
      }
    )
    rvonMises <- nimble::nimbleFunction(
      run = function(n = integer(0), kappa = double(0), mu = double(0)) {
        returnType(double(0))
        return(0)
      }
    )

    dist_registration_list <- list(
      dvonMises = list(
        BUGSdist = "dvonMises(kappa, mu)",
        types    = c("value = double(0)", "kappa = double(0)", "mu = double(0)"),
        pqAvail  = FALSE
      ),
      ddirchmulti = list(
        BUGSdist = "ddirchmulti(alpha, size)",
        types    = c("value = double(1)", "alpha = double(1)", "size = double(0)"),
        pqAvail  = FALSE
      )
    )
    suppressMessages(nimble::registerDistributions(dist_registration_list))

    for (k in seq_along(formula_density_all)) {

      formula_density_k  <- formula_density_all[[k]]
      model_frame_density <- stats::model.frame(formula_density_k, data = station_effort_data)
      X_density           <- stats::model.matrix(stats::as.formula(formula_density_k), model_frame_density)
      nPreds_density      <- ncol(X_density)

      data_REST <- list(
        y           = y,
        N_detection = N_detection,
        stay        = stay,
        censored    = censored,
        N_judge     = N_judge,
        y_enter     = y_enter
      )

      cons_REST <- list(
        N_station           = N_station,
        S                   = S,
        N_period            = N_period,
        nPreds_stay         = nPreds_stay,
        N_stay              = N_stay,
        c_time              = c_time,
        nPreds_density      = nPreds_density,
        stay_family         = stay_family,
        N_group             = N_group,
        activity_estimation = activity_estimation,
        X_density           = X_density,
        X_enter             = X_enter,
        nPreds_enter        = nPreds_enter,
        N_det               = N_det,
        N_enter_group       = N_enter_group,
        N_enter_judge       = N_enter_judge,
        nLevels_stay        = nLevels_stay
      )

      if (activity_estimation == "kernel") cons_REST$activity_proportion <- activity_proportion
      if (nPreds_stay > 1) {
        data_REST$X_stay          <- X_stay
        cons_REST$X_stay_station  <- X_stay_station
      }

      if (!is.null(random_effect_stay)) {
        re_levels <- levels(factor(station_effort_data[[random_effect_stay]]))
        cons_REST$group_stay         <- as.numeric(factor(stay_data_join[[random_effect_stay]], levels = re_levels))
        cons_REST$group_stay_station <- as.numeric(factor(station_effort_data[[random_effect_stay]], levels = re_levels))
      }

      Model_RAD_REST <- nimble::nimbleCode({

        # --- 滞在時間 ---
        for (i in 1:N_stay) {
          censored[i] ~ dinterval(stay[i], c_time[i])

          if (stay_family == "exponential") {
            stay[i]   ~ dexp(rate = 1 / scale[i])
            pred_t[i] ~ dexp(rate = 1 / scale[i])
            loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dexp(stay[i], rate = 1 / scale[i], log = 1) +
              step(censored[i] - 0.5) * log(1 - pexp(c_time[i], rate = 1 / scale[i]))
            loglike_pred_stay[i] <- dexp(pred_t[i], rate = 1 / scale[i], log = 1)
          }
          if (stay_family == "gamma") {
            stay[i]   ~ dgamma(shape = theta_stay, rate = exp(-log(scale[i])))
            pred_t[i] ~ dgamma(shape = theta_stay, rate = exp(-log(scale[i])))
            loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dgamma(stay[i], shape = theta_stay, rate = exp(-log(scale[i])), log = 1) +
              step(censored[i] - 0.5) * log(1 - pgamma(c_time[i], shape = theta_stay, rate = exp(-log(scale[i]))))
            loglike_pred_stay[i] <- dgamma(pred_t[i], shape = theta_stay, rate = exp(-log(scale[i])), log = 1)
          }
          if (stay_family == "lognormal") {
            stay[i]   ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay)
            pred_t[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay)
            loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dlnorm(stay[i], meanlog = log(scale[i]), sdlog = theta_stay, log = 1) +
              step(censored[i] - 0.5) * log(1 - plnorm(c_time[i], meanlog = log(scale[i]), sdlog = theta_stay))
            loglike_pred_stay[i] <- dlnorm(pred_t[i], meanlog = log(scale[i]), sdlog = theta_stay, log = 1)
            meanlog[i]           <- log(scale[i])
          }
          if (stay_family == "weibull") {
            stay[i]   ~ dweibull(shape = theta_stay, scale = scale[i])
            pred_t[i] ~ dweibull(shape = theta_stay, scale = scale[i])
            loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dweibull(stay[i], shape = theta_stay, scale = scale[i], log = 1) +
              step(censored[i] - 0.5) * log(1 - pweibull(c_time[i], shape = theta_stay, scale = scale[i]))
            loglike_pred_stay[i] <- dweibull(pred_t[i], shape = theta_stay, scale = scale[i], log = 1)
          }

          if (nPreds_stay > 1) {
            if (nLevels_stay == 0) {
              log(scale[i]) <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay])
            } else {
              log(scale[i]) <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]]
            }
          } else {
            if (nLevels_stay == 0) {
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

        for (j in 1:nPreds_stay) {
          beta_stay[j] ~ dnorm(0, sd = 100)
        }

        if (nLevels_stay > 0) {
          for (lv in 1:nLevels_stay) {
            random_effect_stay[lv] ~ dnorm(0, sd = sigma_stay)
          }
          sigma_stay ~ T(dnorm(0, sd = 100), 0, 5)
        }

        # --- 期待値計算（mean_stay）---
        # [FIX 3] Weibull の平均は scale * Gamma(1 + 1/shape) → log スケールで加算
        if (nPreds_stay == 1) {
          if (stay_family == "exponential") mean_stay <- exp(beta_stay[1])
          if (stay_family == "gamma")       mean_stay <- theta_stay * exp(beta_stay[1])
          if (stay_family == "lognormal")   mean_stay <- exp(beta_stay[1] + theta_stay^2 / 2)
          if (stay_family == "weibull")     mean_stay <- exp(lgamma(1 + 1 / theta_stay) + beta_stay[1])
        }
        if (nPreds_stay > 1) {
          if (stay_family == "exponential") {
            for (i in 1:N_station) mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
          }
          if (stay_family == "gamma") {
            for (i in 1:N_station) mean_stay[i] <- theta_stay * exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
          }
          if (stay_family == "lognormal") {
            for (i in 1:N_station) mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + theta_stay^2 / 2)
          }
          if (stay_family == "weibull") {
            for (i in 1:N_station) mean_stay[i] <- exp(lgamma(1 + 1 / theta_stay) + inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]))
          }
        }

        # --- 侵入回数（y / Dirichlet-Multinomial）---
        for (g in 1:N_group) {
          for (jj in 1:nPreds_enter) {
            beta_enter[jj, g] ~ dnorm(0, sd = 100)
          }
        }

        for (i in 1:N_station) {
          for (g in 1:N_group) {
            log(alpha_Dirichlet[i, g]) <- inprod(beta_enter[1:nPreds_enter, g], X_enter[i, 1:nPreds_enter])
          }
          alpha_sum[i] <- sum(alpha_Dirichlet[i, 1:N_group])
          for (g in 1:N_group) {
            p_expected[i, g] <- alpha_Dirichlet[i, g] / alpha_sum[i]
            c_expected[i, g] <- p_expected[i, g] * (g - 1)
          }
          mean_pass[i] <- sum(c_expected[i, 1:N_group])

          y[i, 1:N_group]      ~ ddirchmulti(alpha_Dirichlet[i, 1:N_group], N_judge[i])
          pred_y[i, 1:N_group] ~ ddirchmulti(alpha_Dirichlet[i, 1:N_group], N_judge[i])
          loglike_obs_y[i]     <- ddirchmulti(y[i, 1:N_group],      alpha_Dirichlet[i, 1:N_group], N_judge[i], log = 1)
          loglike_pred_y[i]    <- ddirchmulti(pred_y[i, 1:N_group], alpha_Dirichlet[i, 1:N_group], N_judge[i], log = 1)
        }

        # --- 検出回数（N_detection）---
        for (i in 1:N_station) {
          N_detection[i]     ~ dnbinom(size = size, prob = p[i])
          p[i]               <- size / (size + mu[i])
          N_detection_rep[i] ~ dnbinom(size = size, prob = p[i])
          # [FIX 4] log = 1 を明示して対数尤度を正しく格納
          loglike_obs_detection[i]  <- dnbinom(N_detection[i],     size, p[i], log = 1)
          loglike_pred_detection[i] <- dnbinom(N_detection_rep[i], size, p[i], log = 1)
        }
        size ~ dgamma(1, 1)

        # --- REST 式（mu の計算）---
        # [FIX 5] nPreds_density > 1 時は density[i] をベクトルとして一貫して扱う
        if (nPreds_density == 1) {
          for (i in 1:N_station) {
            if (nPreds_stay == 1) {
              log(mu[i]) <- log(density) + log(S) + log(N_period[i]) - log(mean_stay) + log(activity_proportion) - log(mean_pass[i])
            } else {
              log(mu[i]) <- log(density) + log(S) + log(N_period[i]) - log(mean_stay[i]) + log(activity_proportion) - log(mean_pass[i])
            }
          }
          log(density) <- beta_density[1]
          beta_density[1] ~ dnorm(0, sd = 100)
        }
        if (nPreds_density > 1) {
          for (i in 1:N_station) {
            log(density[i]) <- inprod(beta_density[1:nPreds_density], X_density[i, 1:nPreds_density])
            if (nPreds_stay == 1) {
              log(mu[i]) <- log(density[i]) + log(S) + log(N_period[i]) - log(mean_stay) + log(activity_proportion) - log(mean_pass[i])
            } else {
              log(mu[i]) <- log(density[i]) + log(S) + log(N_period[i]) - log(mean_stay[i]) + log(activity_proportion) - log(mean_pass[i])
            }
          }
          for (j in 1:nPreds_density) {
            beta_density[j] ~ dnorm(0, sd = 100)
          }
        }
      })

      inits_f <- function() {
        # [FIX 8] stay_data_join（フィルタ済み）を参照
        stay_mean_log <- log(mean(stay_data_join$Stay, na.rm = TRUE))
        beta_stay_init    <- stats::rnorm(nPreds_stay, 0, 0.1)
        beta_stay_init[1] <- stay_mean_log

        common_inits <- list(
          beta_stay    = beta_stay_init,
          stay         = ifelse(censored == 0, NA, c_time + stats::runif(N_stay, 0.1, 2.0)),
          theta_stay   = stats::runif(1, 0.8, 1.2),
          beta_density = stats::rnorm(nPreds_density, 0, 0.1),
          beta_enter   = matrix(stats::rnorm(nPreds_enter * N_group, 0, 0.1),
                                nrow = nPreds_enter, ncol = N_group),
          size         = stats::runif(1, 0.8, 1.2),
          # [FIX 7] mean_pass をベクトルで初期化（N_station 次元）
          mean_pass    = stats::runif(N_station, 0.5, 2.0)
        )
        if (!is.null(random_effect_stay)) {
          common_inits$random_effect_stay <- stats::runif(nLevels_stay, -0.1, 0.1)
          common_inits$sigma_stay         <- stats::runif(1, 0.8, 1.5)
        }
        return(common_inits)
      }

      # 並列 MCMC 実行関数
      run_MCMC_RAD <- function(info, data, constants, code, params, ni, nt, nb, is_mixture) {
        myModel     <- nimble::nimbleModel(code = code, data = data,
                                           constants = constants, inits = info$inits)
        CmyModel    <- nimble::compileNimble(myModel)
        configModel <- nimble::configureMCMC(myModel, monitors = params)

        # [FIX 12] control リスト方式で統一
        if (is_mixture) {
          configModel$removeSampler("activity_proportion")
          configModel$addSampler(
            target  = "activity_proportion",
            type    = "prior_samples",
            control = list(samples = info$actv_samples)
          )
        }

        myMCMC  <- nimble::buildMCMC(configModel)
        CmyMCMC <- nimble::compileNimble(myMCMC, project = myModel)

        nimble::runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt,
                        nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
      }

      # パラメータリスト
      if (stay_family == "exponential")           prms <- c("scale", "mean_stay")
      if (stay_family %in% c("gamma", "weibull")) prms <- c("scale", "shape", "mean_stay")
      if (stay_family == "lognormal")             prms <- c("meanlog", "sdlog", "mean_stay")
      prms   <- unique(c(prms, "density", "mean_stay", "mean_pass", "mu",
                         "p", "size", "beta_stay", "beta_density"))
      params <- c(prms, "loglike_obs_stay", "loglike_obs_y",
                  "loglike_pred_stay", "loglike_pred_y")
      if (activity_estimation == "mixture") params <- c(params, "activity_proportion")

      # [FIX 9] nc を上書きせず nc_run を使用
      if (activity_estimation == "mixture") {
        nc_run <- length(actv_out_trace)
        per_chain_info <- lapply(1:nc_run, function(i) {
          list(seed = sample(1:9999, 1), inits = inits_f(),
               actv_samples = as.matrix(actv_out_trace[[i]]))
        })
      } else {
        nc_run <- nc
        per_chain_info <- lapply(1:nc_run, function(i) {
          list(seed = sample(1:9999, 1), inits = inits_f(), actv_samples = NULL)
        })
      }

      cat("Running MCMC sampling for RAD-REST. Please wait...\n")

      this_cluster <- parallel::makeCluster(nc_run)
      # [FIX 10] on.exit でクラスター残留を防止
      on.exit(try(parallel::stopCluster(this_cluster), silent = TRUE), add = TRUE)

      parallel::clusterExport(this_cluster,
                              c("dist_registration_list",
                                "dvonMises", "rvonMises",
                                "ddirchmulti", "rdirchmulti",
                                "run_MCMC_RAD"),
                              envir = environment())
      parallel::clusterEvalQ(this_cluster, {
        library(nimble)
        suppressMessages(nimble::registerDistributions(dist_registration_list))
      })

      chain_output <- parallel::parLapply(
        cl = this_cluster, X = per_chain_info,
        fun = run_MCMC_RAD,
        data = data_REST, code = Model_RAD_REST,
        constants = cons_REST, params = params,
        ni = ni, nt = nt, nb = nb,
        is_mixture = (activity_estimation == "mixture")
      )
      parallel::stopCluster(this_cluster)
      on.exit()
      cat("Estimation is finished!\n")

      # WAIC 計算（Log-Sum-Exp で数値安定化）
      safe_log_mean_exp <- function(x) {
        max_val <- max(x, na.rm = TRUE)
        max_val + log(mean(exp(x - max_val), na.rm = TRUE))
      }

      loglfstay <- MCMCvis::MCMCchains(chain_output, params = "loglike_obs_stay")
      loglfy    <- MCMCvis::MCMCchains(chain_output, params = "loglike_obs_y")
      loglfall  <- if (activity_estimation == "mixture") cbind(loglfstay, loglfy, loglact) else cbind(loglfstay, loglfy)
      lppd    <- sum(apply(loglfall, 2, safe_log_mean_exp))
      p.waic  <- sum(apply(loglfall, 2, stats::var))
      waic[k] <- (-2) * lppd + 2 * p.waic

      mcmc_samples[[k]] <- MCMCvis::MCMCchains(chain_output, mcmc.list = TRUE, params = prms)
      samples_mat        <- MCMCvis::MCMCchains(chain_output)
      n_iters   <- nrow(samples_mat)
      p_names   <- colnames(samples_mat)
      tidy_samples[[k]] <- data.frame(
        parameter = rep(p_names, each = n_iters),
        value     = as.vector(samples_mat),
        iteration = rep(1:n_iters, times = length(p_names)),
        stringsAsFactors = FALSE
      )
    }
  }

  # ===========================================================================
  # 結果集計
  # ===========================================================================

  if (is.null(random_effect_stay)) random_effect_stay <- "NULL"

  check_no_cov <- function(f) {
    if (is.null(f)) return(TRUE)
    f <- stats::as.formula(f)
    if (length(f) == 3) f <- f[-2]
    length(all.vars(f)) == 0
  }

  is_density_global <- check_no_cov(formula_density)
  is_stay_global    <- check_no_cov(formula_stay) && random_effect_stay == "NULL"
  is_enter_global   <- check_no_cov(formula_enter)

  is_pass_global <- if (model == "REST") {
    is_density_global && is_stay_global
  } else {
    is_enter_global
  }

  WAIC <- data.frame(
    Model              = sapply(formula_density_all, function(x) paste(deparse(x), collapse = " ")),
    random_effect_stay = random_effect_stay,
    WAIC               = waic
  ) %>%
    dplyr::arrange(WAIC)

  best.model        <- which.min(waic)
  mcmc_samples_best <- mcmc_samples[[best.model]]
  tidy_samples_best <- tidy_samples[[best.model]]

  if (activity_estimation == "mixture") {
    sample_activity <- MCMCvis::MCMCchains(actv_chain_output, mcmc.list = TRUE,
                                           params = "activity_proportion")
    mcmc_samples_best <- lapply(seq_along(mcmc_samples_best), function(i) {
      coda::as.mcmc(cbind(mcmc_samples_best[[i]], sample_activity[[i]]))
    })
    mcmc_samples_best <- coda::as.mcmc.list(mcmc_samples_best)
  }

  prms <- c("density", "mean_stay")
  if (model == "RAD-REST")              prms <- c(prms, "mean_pass")
  if (activity_estimation == "mixture") prms <- c(prms, "activity_proportion")

  summary_mean_temp <- MCMCvis::MCMCsummary(mcmc_samples_best, params = prms, round = 2) %>%
    tibble::rownames_to_column(var = "Variable") %>%
    tibble::as_tibble() %>%
    dplyr::rename(lower = `2.5%`, median = `50%`, upper = `97.5%`) %>%
    dplyr::mutate(cv = sd / mean)

  summary_mean <- summary_mean_temp %>%
    dplyr::mutate(
      Species = target_species,
      BaseVar = stringr::str_remove(Variable, "\\[.*\\]"),
      idx     = as.integer(stringr::str_extract(Variable, "(?<=\\[)\\d+(?=\\])")),
      Station = ifelse(!is.na(idx), station.id[idx], "All")
    ) %>%
    dplyr::filter(
      !(BaseVar == "density"   & is_density_global & !is.na(idx) & idx != 1),
      !(BaseVar == "mean_stay" & is_stay_global    & !is.na(idx) & idx != 1),
      !(BaseVar == "mean_pass" & is_pass_global    & !is.na(idx) & idx != 1)
    ) %>%
    dplyr::mutate(
      Station = ifelse(
        (BaseVar == "density"   & is_density_global) |
          (BaseVar == "mean_stay" & is_stay_global)    |
          (BaseVar == "mean_pass" & is_pass_global),
        "All", Station
      ),
      Variable = ifelse(
        (BaseVar == "density"   & is_density_global) |
          (BaseVar == "mean_stay" & is_stay_global)    |
          (BaseVar == "mean_pass" & is_pass_global),
        BaseVar, Variable
      )
    ) %>%
    dplyr::select(Species, Station, Variable, mean, sd, lower, median, upper, Rhat, n.eff, cv)

  density_result <- list(
    WAIC           = WAIC,
    summary_result = summary_mean,
    samples        = mcmc_samples_best
  )

  if (activity_estimation == "mixture") {
    density_result$activity_curve <- activity_density_estimates
  }

  class(density_result) <- "ResultDensity"
  density_result
}
