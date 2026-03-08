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
#'     \item \code{Censored} (binary): Indicator for censored staying time (1 = censored, 0 = observed).
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
#' @importFrom stats as.formula formula model.frame model.matrix sd var runif median quantile model.response rexp rnorm step dexp pexp dgamma pgamma dlnorm plnorm dweibull pweibull dnbinom
#' @importFrom dplyr select
#' @importFrom extraDistr ddirmnom
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
#'   formula_enter = ~ 1, # Ignored since model = "REST"
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

bayes_rest_multi <- function(formula_stay,
                             formula_density,
                             formula_enter = NULL, # 今回使用しない場合はNULL許容
                             station_effort_data,
                             stay_data,
                             random_effect_stay = NULL,
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
                             target_species) {

  # ============================================================================
  # 1. エラーチェック（入力検証）
  # ============================================================================
  cat("Running input checks...\n")

  # 1.1 引数の値チェック
  if (!activity_estimation %in% c("kernel", "mixture")) {
    stop("Error: 'activity_estimation' must be either 'kernel' or 'mixture'.")
  }
  if (!stay_family %in% c("exponential", "gamma", "lognormal", "weibull")) {
    stop("Error: 'stay_family' must be 'exponential', 'gamma', 'lognormal', or 'weibull'.")
  }
  if (iter <= warmup) {
    stop("Error: 'iter' must be greater than 'warmup'.")
  }
  if (!is.numeric(focal_area) || focal_area <= 0) {
    stop("Error: 'focal_area' must be a positive number.")
  }

  # 1.2 必須カラムのチェック (想定される標準的なカラム名に基づく)
  required_cols_station <- c("Station", "Y", "N_period", "Species")
  missing_st <- setdiff(required_cols_station, colnames(station_effort_data))
  if (length(missing_st) > 0) stop(paste("Error: station_effort_data is missing columns:", paste(missing_st, collapse=", ")))

  required_cols_stay <- c("Station", "Stay", "Censored", "Species")
  missing_stay <- setdiff(required_cols_stay, colnames(stay_data))
  if (length(missing_stay) > 0) stop(paste("Error: stay_data is missing columns:", paste(missing_stay, collapse=", ")))

  required_cols_act <- c("Time", "Species") # Timeはラジアン(0~2π)を想定
  missing_act <- setdiff(required_cols_act, colnames(activity_data))
  if (length(missing_act) > 0) stop(paste("Error: activity_data is missing columns:", paste(missing_act, collapse=", ")))

  # ============================================================================
  # 2. データの前処理
  # ============================================================================
  cat("Preparing data for target species:", target_species, "...\n")

  # ターゲット種でフィルタリング
  st_data_sp <- station_effort_data[station_effort_data$Species == target_species, ]
  stay_data_sp <- stay_data[stay_data$Species == target_species, ]
  act_data_sp <- activity_data[activity_data$Species == target_species, ]

  if (nrow(st_data_sp) == 0) stop("Error: No data found for target_species in station_effort_data.")
  if (nrow(stay_data_sp) == 0) stop("Error: No data found for target_species in stay_data.")
  if (nrow(act_data_sp) == 0) stop("Error: No data found for target_species in activity_data.")

  # MCMCパラメータの準備
  ni <- iter
  nb <- warmup
  nt <- thin
  nc <- chains

  # モデル用定数・データの切り出し
  S <- focal_area
  N_station <- nrow(st_data_sp)
  N_stay <- nrow(stay_data_sp)
  act_data <- act_data_sp$Time

  # Stayデータの整形 (Censored = 1 なら Stay は NA, c_time を設定)
  censored <- stay_data_sp$Censored
  stay <- stay_data_sp$Stay
  c_time <- stay_data_sp$Stay # 打ち切り時間（観測された最大時間）
  stay[censored == 1] <- NA

  # 共変量行列の作成 (Stay)
  model_frame_stay <- stats::model.frame(formula_stay, data = stay_data_sp)
  X_stay <- stats::model.matrix(formula_stay, model_frame_stay)
  nPreds_stay <- ncol(X_stay)

  nLevels_stay <- 0
  if (!is.null(random_effect_stay)) {
    if (!random_effect_stay %in% colnames(stay_data_sp)) stop("Error: random_effect_stay column not found in stay_data.")
    group_stay <- as.numeric(as.factor(stay_data_sp[[random_effect_stay]]))
    nLevels_stay <- max(group_stay)
  }

  # formula_densityが単一の式かリストか判定し、リストに統一
  if (!is.list(formula_density)) {
    formula_density_all <- list(formula_density)
  } else {
    formula_density_all <- formula_density
  }

  # ============================================================================
  # 3. カスタム分布の定義 (von Mises)
  # ============================================================================
  dvonMises <- nimble::nimbleFunction(
    run = function(x = double(0), kappa = double(0), mu = double(0), log = integer(0)) {
      returnType(double(0))
      ccrit <- 1E-6; s <- 1; i <- 1; inc <- 1; x_2i <- 0; satisfied <- FALSE
      while(!satisfied) {
        x_2i <- kappa / (2 * i)
        inc <- inc * x_2i * x_2i
        s <- s + inc
        i <- i + 1
        satisfied <- inc < ccrit
      }
      prob <- exp(kappa * cos(x - mu)) / (2 * pi * s)
      if (log) return(log(prob)) else return(prob)
    }
  )

  rvonMises <- nimble::nimbleFunction(
    run = function(n = integer(0), kappa = double(0), mu = double(0)) {
      returnType(double(0))
      return(0) # サンプリングは実質的にダミー
    }
  )

  suppressWarnings(suppressMessages(nimble::registerDistributions(list(
    dvonMises = list(
      BUGSdist = "dvonMises(kappa, mu)",
      types = c('value = double(0)', 'kappa = double(0)', 'mu = double(0)'),
      pqAvail = FALSE
    )
  ))))

  # ============================================================================
  # 4. Activity Estimation
  # ============================================================================
  cat("Estimating activity proportion (Method:", activity_estimation, ")...\n")

  activity_proportion <- NA
  actv_out_trace <- NULL
  loglact <- NULL

  if (activity_estimation == "kernel") {
    # 簡易的なカーネル密度推定（activityパッケージ等を使用するのが一般的ですが、ここでは基本実装）
    dens <- stats::density(act_data, bw = "nrd0", adjust = bw_adj, from = 0, to = 2*pi)
    activity_proportion <- 1 / (2 * pi * max(dens$y))
    cat("Kernel activity proportion estimated:", round(activity_proportion, 3), "\n")

  } else if (activity_estimation == "mixture") {
    dens.x <- seq(0, 2 * pi, 0.02)
    ndens <- length(dens.x)
    N_act <- length(act_data)
    constants_act <- list(N = N_act, C = C, dens.x = dens.x, ndens = ndens)
    data_act <- list(act_data = act_data)

    code_act <- nimble::nimbleCode({
      for(k in 1:(C-1)) { v[k] ~ dbeta(1, alpha) }
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
      }
      for (j in 1:ndens) {
        for (i in 1:C) { dens.cpt[i, j] <- w[i] * dvonMises(dens.x[j] , mu_mix[i], kappa_mix[i], log = 0) }
        activity_density[j] <- sum(dens.cpt[1:C, j])
      }
      activity_proportion <- 1.0 / (2 * 3.141592654 * max(activity_density[1:ndens]))
    })

    inits_f_act <- function() {
      list(
        mu_mix = stats::runif(constants_act$C, 0, 2 * pi),
        kappa_mix = stats::rgamma(constants_act$C, 1, 0.01),
        group = sample(1:constants_act$C, size = constants_act$N, replace = TRUE),
        v = stats::rbeta(constants_act$C-1, 1, 1),
        alpha = 1
      )
    }

    run_MCMC_vonMises <- function(info, data, constants, code, params, ni, nt, nb) {
      suppressMessages(nimble::registerDistributions(list(
        dvonMises = list(BUGSdist = "dvonMises(kappa, mu)", types = c('value = double(0)', 'kappa = double(0)', 'mu = double(0)'), pqAvail = FALSE)
      )))
      myModel <- nimble::nimbleModel(code = code, data = data, constants = constants, inits = info$inits)
      CmyModel <- nimble::compileNimble(myModel)
      configModel <- nimble::configureMCMC(myModel, monitors = params)
      myMCMC <- nimble::buildMCMC(configModel, monitors = params)
      CmyMCMC <- nimble::compileNimble(myMCMC)
      nimble::runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
    }

    per_chain_info_act <- lapply(1:nc, function(i) list(seed = sample(1:9999, 1), inits = inits_f_act()))
    params_act <- c("activity_proportion", "loglike_obs_act")

    this_cluster <- parallel::makeCluster(cores)
    parallel::clusterEvalQ(this_cluster, { library(nimble) })
    parallel::clusterExport(this_cluster, c("dvonMises", "rvonMises", "run_MCMC_vonMises"), envir = environment())

    actv_chain_output <- tryCatch({
      parallel::parLapply(cl = this_cluster, X = per_chain_info_act, fun = run_MCMC_vonMises,
                          data = data_act, code = code_act, constants = constants_act, params = params_act, ni = ni, nt = nt, nb = nb)
    }, finally = { parallel::stopCluster(this_cluster) })

    actv_out_trace <- purrr::map(actv_chain_output, ~ .[, grep("activity_proportion", colnames(.))])
    loglact <- MCMCvis::MCMCchains(actv_chain_output, params = c("loglike_obs_act"))
  }

  # ============================================================================
  # 5. REST Model Estimation
  # ============================================================================
  tidy_samples <- mcmc_samples <- list()
  waic <- numeric(length(formula_density_all))

  for (k in 1:length(formula_density_all)) {
    cat(sprintf("\nRunning REST Model MCMC (%d/%d) -> Formula: %s\n", k, length(formula_density_all), deparse(formula_density_all[[k]])))

    model_frame_density <- stats::model.frame(formula_density_all[[k]], data = st_data_sp)
    X_density <- stats::model.matrix(formula_density_all[[k]], model_frame_density)
    nPreds_density <- ncol(X_density)

    data_REST <- list(y = st_data_sp$Y, X_density = X_density, stay = stay, censored = censored)
    cons_REST <- list(
      N_station = N_station, S = S, N_period = st_data_sp$N_period,
      nPreds_stay = nPreds_stay, N_stay = N_stay, c_time = c_time,
      nPreds_density = nPreds_density, stay_family = stay_family,
      activity_estimation = activity_estimation,
      nLevels_stay = nLevels_stay
    )

    if (nLevels_stay > 0) cons_REST$group_stay <- group_stay
    if (activity_estimation == "kernel") cons_REST$activity_proportion <- activity_proportion
    if (nPreds_stay > 1) data_REST$X_stay <- X_stay

    Model_REST <- nimble::nimbleCode({
      if (activity_estimation == "mixture") activity_proportion ~ dunif(0, 1) # Dummy prior

      for(i in 1:N_stay) {
        censored[i] ~ dinterval(stay[i], c_time[i])

        if (stay_family == "lognormal") {
          stay[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay)
          loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dlnorm(stay[i], meanlog = log(scale[i]), sdlog = theta_stay, log = 1) +
            step(censored[i] - 0.5) * log(1 - plnorm(c_time[i], meanlog = log(scale[i]), sdlog = theta_stay))
          meanlog[i] <- log(scale[i])
        } else if (stay_family == "gamma") {
          stay[i] ~ dgamma(shape = theta_stay, rate = exp(-log(scale[i])))
          loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dgamma(stay[i], shape = theta_stay, rate = exp(-log(scale[i])), log = 1) +
            step(censored[i] - 0.5) * log(1 - pgamma(c_time[i], shape = theta_stay, rate = exp(-log(scale[i]))))
        } else if (stay_family == "weibull") {
          stay[i] ~ dweibull(shape = theta_stay, scale = scale[i])
          loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dweibull(stay[i], shape = theta_stay, scale = scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pweibull(c_time[i], shape = theta_stay, scale = scale[i]))
        } else if (stay_family == "exponential") {
          stay[i] ~ dexp(rate = 1/scale[i])
          loglike_obs_stay[i] <- (1 - step(censored[i] - 0.5)) * dexp(stay[i], rate = 1/scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pexp(c_time[i], rate = 1/scale[i]))
        }

        if (nPreds_stay > 1) {
          log(scale[i]) <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + (if(nLevels_stay > 0) random_effect_stay[group_stay[i]] else 0)
        } else {
          log(scale[i]) <- beta_stay[1] + (if(nLevels_stay > 0) random_effect_stay[group_stay[i]] else 0)
        }
      }

      theta_stay ~ dgamma(1, 1)
      if (stay_family == "lognormal") { sdlog <- theta_stay } else { shape <- theta_stay }
      for(j in 1:nPreds_stay) { beta_stay[j] ~ T(dnorm(0, 100), -10, 10) }

      if (nLevels_stay > 0) {
        for(level in 1:nLevels_stay) { random_effect_stay[level] ~ dnorm(0, sd = sigma_stay) }
        sigma_stay ~ T(dnorm(0, 100), 0, 10)
      }

      # mean_stay の計算ブロック
      if (nPreds_stay == 1) {
        if (stay_family == "lognormal")   mean_stay <- exp(beta_stay[1] + theta_stay ^ 2 / 2)
        if (stay_family == "gamma")       mean_stay <- theta_stay * exp(beta_stay[1])
        if (stay_family == "weibull")     mean_stay <- exp(lgamma(1 + 1 / theta_stay)) + exp(beta_stay[1])
        if (stay_family == "exponential") mean_stay <- exp(beta_stay[1])
      } else {
        for(i in 1:N_station){
          val <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + (if(nLevels_stay > 0) random_effect_stay[group_stay[i]] else 0)
          if (stay_family == "lognormal")   mean_stay[i] <- exp(val + theta_stay ^ 2 / 2)
          if (stay_family == "gamma")       mean_stay[i] <- theta_stay * exp(val)
          if (stay_family == "weibull")     mean_stay[i] <- exp(lgamma(1 + 1 / theta_stay)) + exp(val)
          if (stay_family == "exponential") mean_stay[i] <- exp(val)
        }
      }

      for (i in 1:N_station) {
        y[i] ~ dnbinom(size = size, prob = p[i])
        p[i] <- size / (size + mu[i])
        loglike_obs_y[i] <- dnbinom(y[i], size, p[i], log = 1)

        # μの計算
        log_density_val <- if(nPreds_density == 1) beta_density else inprod(beta_density[1:nPreds_density], X_density[i, 1:nPreds_density])
        mean_stay_val <- if(nPreds_stay == 1) mean_stay else mean_stay[i]

        log(mu[i]) <- log_density_val + log(S) + log(N_period[i]) - log(mean_stay_val) + log(activity_proportion)
        if(nPreds_density > 1) density[i] <- exp(log_density_val)
      }

      if(nPreds_density == 1) {
        log(density) <- beta_density
        beta_density ~ dnorm(0, sd = 100)
      } else {
        for(j in 1:nPreds_density) { beta_density[j] ~ dnorm(0, sd = 100) }
      }
      size ~ dgamma(1, 1)
    })

    inits_f_rest <- function() {
      inits <- list(
        beta_stay = stats::runif(nPreds_stay, -1, 1),
        stay = ifelse(censored == 0, NA, c_time + stats::runif(N_stay, 0.1, 1.0)),
        theta_stay = stats::runif(1, 0.5, 1.5),
        scale = stats::runif(N_stay, 0, 2),
        beta_density = stats::rnorm(nPreds_density, 0, 1),
        size = stats::runif(1, 0.5, 2)
      )
      if (nLevels_stay > 0) {
        inits$random_effect_stay <- stats::runif(nLevels_stay, -1, 1)
        inits$sigma_stay <- stats::runif(1, 0.5, 2.5)
      }
      return(inits)
    }

    run_MCMC_REST <- function(info, data, constants, code, params, ni, nt, nb) {
      suppressMessages(nimble::registerDistributions(list(
        dvonMises = list(BUGSdist = "dvonMises(kappa, mu)", types = c('value = double(0)', 'kappa = double(0)', 'mu = double(0)'), pqAvail = FALSE)
      )))
      myModel <- nimble::nimbleModel(code = code, data = data, constants = constants, inits = info$inits)
      CmyModel <- nimble::compileNimble(myModel)
      configModel <- nimble::configureMCMC(myModel, monitors = params)

      if (constants$activity_estimation == "mixture") {
        configModel$removeSampler("activity_proportion")
        configModel$addSampler(target = "activity_proportion", type = 'prior_samples', samples = info$actv_samples)
      }

      myMCMC <- nimble::buildMCMC(configModel, monitors = params)
      CmyMCMC <- nimble::compileNimble(myMCMC)
      nimble::runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
    }

    per_chain_info_rest <- lapply(1:nc, function(i) {
      info <- list(seed = sample(1:9999, 1), inits = inits_f_rest())
      if (activity_estimation == "mixture") info$actv_samples <- as.matrix(actv_out_trace[[i]])
      return(info)
    })

    prms <- c("beta_stay", "beta_density", "density", "p", "size", "scale", "mean_stay")
    if (stay_family %in% c("gamma", "weibull")) prms <- c(prms, "shape")
    if (stay_family == "lognormal") prms <- c(prms, "meanlog", "sdlog")
    params_rest <- c(prms, "loglike_obs_stay", "loglike_obs_y")
    if (activity_estimation == "mixture") params_rest <- c(params_rest, "activity_proportion")

    this_cluster_rest <- parallel::makeCluster(cores)
    parallel::clusterEvalQ(this_cluster_rest, { library(nimble) })
    parallel::clusterExport(this_cluster_rest, c("dvonMises", "rvonMises", "run_MCMC_REST"), envir = environment())

    chain_output <- tryCatch({
      parallel::parLapply(cl = this_cluster_rest, X = per_chain_info_rest, fun = run_MCMC_REST,
                          data = data_REST, code = Model_REST, constants = cons_REST, params = params_rest, ni = ni, nt = nt, nb = nb)
    }, finally = { parallel::stopCluster(this_cluster_rest) })

    # --- WAIC 計算 ---
    loglfstay <- MCMCvis::MCMCchains(chain_output, params = "loglike_obs_stay")
    loglfy <- MCMCvis::MCMCchains(chain_output, params = "loglike_obs_y")
    loglfall <- if (activity_estimation == "mixture") cbind(loglfstay, loglfy, loglact) else cbind(loglfstay, loglfy)

    lppd <- sum(log(colMeans(exp(loglfall))))
    p.waic <- sum(apply(loglfall, 2, stats::var))
    waic[k] <- (-2) * lppd + 2 * p.waic

    # --- サンプル保存 ---
    mcmc_samples[[k]] <- MCMCvis::MCMCchains(chain_output, mcmc.list = TRUE, params = prms)
    samples_mat <- MCMCvis::MCMCchains(chain_output, params = prms)
    tidy_samples[[k]] <- data.frame(
      parameter = rep(colnames(samples_mat), each = nrow(samples_mat)),
      value = as.vector(samples_mat),
      iteration = rep(1:nrow(samples_mat), times = ncol(samples_mat)),
      stringsAsFactors = FALSE
    )
  }

  cat("\nAll models successfully estimated!\n")

  return(list(
    tidy_samples = tidy_samples,
    mcmc_samples = mcmc_samples,
    waic = waic,
    activity_proportion = if(activity_estimation == "kernel") activity_proportion else NULL
  ))
}
