#' Bayesian Parameter Estimation of a Multispecies REST/RAD-REST Model via MCMC Sampling using `nimble`
#'
#' @param formula_stay A model formula for staying times within the focal area (e.g., \code{Stay ~ 1 + x1}). The left-hand side must specify the column name for staying time. Random effects should be specified separately using the \code{random_effect_stay} argument.
#' @param formula_density A model formula for animal density (e.g., \code{~ 1 + x1}). Note that the left-hand side must be omitted as density is a latent parameter.
#' @param formula_enter A model formula for the total number of passes (\code{Y}) entering the focal area (e.g., \code{~ 1 + x1}). Note that the left-hand side must be omitted. This is particularly used for the RAD-REST model or when estimating intrusion rates.
#' @param station_effort_data A data frame containing information for each camera station. Typically, this is the output of the \code{add_effort} function. Alternatively, a manually prepared data frame may be provided with the following required columns:
#'   \itemize{
#'     \item For \code{model = "REST"}:
#'       \itemize{
#'         \item \code{Station} (character): Camera station ID.
#'         \item \code{Effort} (numeric): Camera trapping effort (in days) at each station.
#'         \item \code{Species} (character): Species name.
#'         \item \code{Y} (numeric): Total number of passes through the focal area for each station-species pair.
#'       }
#'     \item For \code{model = "RAD-REST"}:
#'       \itemize{
#'         \item \code{Station} (character): Camera station ID.
#'         \item \code{Effort} (numeric): Camera trapping effort (in days) at each station.
#'         \item \code{Species} (character): Species name.
#'         \item \code{N} (numeric): Total number of detected videos for each station-species pair.
#'         \item \code{y_X} columns (\code{y_0}, \code{y_1}, ..., \code{y_max}): Number of videos categorized by the number of observed passes.
#'       }
#'   }
#' @param stay_data A data frame returned by the \code{format_stay} function, containing the following columns:
#'   \itemize{
#'     \item \code{Station} (character): Camera station ID.
#'     \item \code{Species} (character): Species name.
#'     \item \code{Stay} (numeric): Staying time (in seconds) within the focal area for each detected pass.
#'     \item \code{Censored} (binary): Indicator for censored staying time (1 = censored, 0 = observed).
#'   }
#' @param random_effect_stay A character string specifying a random effect structure on mean staying time (e.g., \code{"~ (1 | Station)"}). Default is \code{NULL}. Note: Species-level random effects are automatically included in this multispecies model and do not need to be manually specified here.
#' @param activity_data A data frame containing a \code{time} column, representing detection times transformed into radians. Typically, this is the output of the \code{format_activity} function.
#' @param activity_estimation A character string specifying the method used to estimate activity patterns. Choose \code{"kernel"} for fixed kernel density estimation (Rowcliffe et al. 2014), or \code{"mixture"} for nonparametric Bayesian estimation using von Mises mixture models (Nakashima et al. 2025). Default is \code{"kernel"}.
#' @param bw_adj A numeric bandwidth adjustment parameter for kernel density estimation. Default is 1.0. See Rowcliffe et al. (2014) for details.
#' @param C An integer specifying the maximum number of von Mises components to use in the mixture model. Required only if \code{activity_estimation = "mixture"}. Default is 10.
#' @param stay_family A character string specifying the probability distribution used to model staying times (e.g., \code{"exponential"}, \code{"gamma"}, \code{"lognormal"}, \code{"weibull"}). This should ideally be selected based on the output of the \code{bayes_stay_selection} function. Default is \code{"lognormal"}.
#' @param focal_area A numeric value representing the size of the focal area, in square meters.
#' @param cores An integer specifying the number of CPU cores to use for parallel computation. Default is 3.
#' @param iter An integer specifying the total number of MCMC iterations per chain. Default is 5000.
#' @param warmup An integer specifying the number of warm-up (burn-in) iterations per chain. Default is 1000.
#' @param chains An integer specifying the number of MCMC chains. Default is 3.
#' @param thin An integer specifying the thinning interval for MCMC sampling. Default is 4 (meaning every 4th sample is kept).
#' @param target_species A character vector specifying the species to be analyzed. Multiple species must be specified for this multispecies model.
#' @return A list of class \code{"ResultDensity"}, which includes the following components:
#' \describe{
#'   \item{\code{WAIC}}{An object containing WAIC (Widely Applicable Information Criterion) results for model evaluation.}
#'   \item{\code{summary_result}}{A data frame summarizing the posterior estimates, including the mean staying time and density across species.}
#'   \item{\code{samples}}{A \code{coda::mcmc.list} object containing full MCMC posterior samples for all model parameters.}
#' }
#' The returned object has a custom print method that displays the WAIC and a summary of the parameter estimates.
#' You can access the full MCMC samples via \code{$samples} and visually analyze convergence using the \code{MCMCvis} package.
#'
#' @export
#' @import nimble activity parallel MCMCvis tibble
#' @importFrom stats as.formula formula model.frame model.matrix sd var runif median quantile model.response rexp rnorm step dexp pexp dgamma pgamma dlnorm plnorm dweibull pweibull dnbinom
#' @importFrom dplyr select
#' @importFrom extraDistr ddirmnom
#' @examples
#' \dontrun{
#' station_data_RAD <- format_station_data(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_y = "y",
#'   model = "RAD-REST"
#' )
#'
#' station_effort_RAD <- add_effort(
#'   detection_data = detection_data,
#'   station_data_formatted = station_data_RAD,
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
#' rest_model <- bayes_rest_multi(
#'   formula_stay = Stay ~ 1,
#'   formula_density = ~ 1,
#'   formula_enter = ~ 1,
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
#'   target_species = c("SP01", "SP02", "SP03", "SP04", "SP05")
#' )
#' }
bayes_rest_multi <- function(formula_stay,
                             formula_density,
                             formula_enter,
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

  # Define functions --------------------------------------------------------
  ni <- iter
  nt <- thin
  nc <- chains
  nb <- warmup
  nmcmc <- (iter - warmup) * chains / thin


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

  # -------------------------------------------------------------------------
  # 共変量を標準化するヘルパー関数
  # 切片列（全要素が1の列）を除く数値列を平均0・標準偏差1に正規化する。
  # 標準偏差が0（定数列）の列はスケーリングをスキップして警告を出す。
  # 戻り値: list(X = 正規化済み行列, center = 各列の平均, scale = 各列のSD)
  # -------------------------------------------------------------------------
  standardize_design_matrix <- function(X) {
    center_vec <- rep(0, ncol(X))
    scale_vec  <- rep(1, ncol(X))

    for (j in seq_len(ncol(X))) {
      col_j <- X[, j]
      # 切片列（全要素が1）はスキップ
      if (all(col_j == 1)) next
      # 数値列のみ正規化
      if (!is.numeric(col_j)) next

      col_mean <- mean(col_j, na.rm = TRUE)
      col_sd   <- sd(col_j,   na.rm = TRUE)

      if (is.na(col_sd) || col_sd == 0) {
        warning(paste0(
          "Column '", colnames(X)[j], "' has zero variance and will not be scaled."
        ))
        next
      }

      X[, j]       <- (col_j - col_mean) / col_sd
      center_vec[j] <- col_mean
      scale_vec[j]  <- col_sd
    }

    list(X = X, center = center_vec, scale = scale_vec)
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

  if (!is.character(random_effect_stay) && !is.null(random_effect_stay)) {
    stop("`random_effect_stay` must be a character vector or NULL.")
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

  if (!(stay_family %in% c("lognormal", "gamma", "weibull", "exponential"))) {
    stop(paste0("Input family type(", stay_family, ") is incorrect."))
  }

  # Extract variable names for stay
  vars_stay <- all.vars(formula_stay)
  response_stay <- vars_stay[1]
  predictors_stay <- vars_stay[-1]

  model_frame_stay <- model.frame(formula_stay, stay_data_join)
  X_stay_raw <- model.matrix(as.formula(formula_stay), model_frame_stay)
  stay <- model.response(model_frame_stay)
  censored <- stay_data_join[["Cens"]]
  is.censored <- stay_data_join[["Cens"]]

  # ----- 滞在時間の共変量を標準化 -----
  scaled_stay   <- standardize_design_matrix(X_stay_raw)
  X_stay        <- scaled_stay$X
  scaling_stay  <- list(center = scaled_stay$center, scale = scaled_stay$scale)
  if (ncol(X_stay) > 1) {
    cat("Covariates in formula_stay have been standardized (mean=0, sd=1).\n")
  }

  # N of random effects
  if (!is.null(random_effect_stay)) {
    levels <- unique(stay_data_join[[random_effect_stay]])
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

  vars_enter <- all.vars(formula_enter)
  predictors_enter <- vars_enter

  # MCMCの計算結果を格納するためのリスト等（ループが無くなるためシンプルな構成でOKです）
  samples <- list()
  density_estimates <- list()

  # === Density data preparation ===
  model_frame_density <- model.frame(formula_density, station_effort_data)
  X_density_raw <- model.matrix(formula_density, model_frame_density)
  nPreds_density_raw <- ncol(X_density_raw)
  X_density_raw <- X_density_raw[1:N_station , , drop = FALSE]

  # ----- 密度の共変量を標準化 -----
  scaled_density   <- standardize_design_matrix(X_density_raw)
  X_density        <- scaled_density$X
  scaling_density  <- list(center = scaled_density$center, scale = scaled_density$scale)
  nPreds_density   <- ncol(X_density)
  if (nPreds_density > 1) {
    cat("Covariates in formula_density have been standardized (mean=0, sd=1).\n")
  }

  # === Enter (Dirichlet Alpha) data preparation ===
  model_frame_enter <- model.frame(formula_enter, station_effort_data)
  X_alpha_raw <- model.matrix(formula_enter, model_frame_enter)
  X_alpha_raw <- X_alpha_raw[1:N_station , , drop = FALSE]

  # ----- 進入率の共変量を標準化 -----
  scaled_alpha   <- standardize_design_matrix(X_alpha_raw)
  X_alpha        <- scaled_alpha$X
  scaling_alpha  <- list(center = scaled_alpha$center, scale = scaled_alpha$scale)
  nPreds_alpha   <- ncol(X_alpha)
  if (nPreds_alpha > 1) {
    cat("Covariates in formula_enter have been standardized (mean=0, sd=1).\n")
  }

  S <- focal_area * 10^-6
  N_period <- station_effort_data$Effort[1:N_station] * 60 * 60 * 24
  N_detection_matrix <- matrix(N_detection, ncol = nSpecies, byrow = FALSE)

  species_id_stay <- stay_data_join %>% pull(Species) %>% factor() %>% as.numeric()
  species_id_ey <- station_effort_data %>% pull(Species) %>% factor() %>% as.numeric()

  # 地点IDの作成 (Station列を使用)
  station_id_ey <- station_effort_data %>% pull(Station) %>% factor() %>% as.numeric()

  data_density <- list(stay = stay, is.censored = is.censored, y = y, N_judge = N_judge, N_detection_matrix = N_detection_matrix)

  if(nPreds_density > 1) data_density$X_density <- X_density
  if(nPreds_alpha > 1) data_density$X_alpha <- X_alpha

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
    station_id_ey = station_id_ey,
    S = S,
    N_period = N_period,
    nPreds_density = nPreds_density,
    nPreds_alpha = nPreds_alpha
  )
  if(activity_estimation == "kernel") {
    cons_density$activity_proportion <- activity_proportion
  }

  if (!is.null(random_effect_stay)) {
    cons_density$group <- as.numeric(factor(stay_data_join[[random_effect_stay]]))
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
      } else {
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

    # Expected value calculation for stay
    if (nPreds_stay == 1) {
      if(stay_family == "exponential") {
        for(m in 1:nSpecies) { mean_stay[m] <- exp(beta_stay[1] + species_effect_stay[m, 1]) }
      }
      if(stay_family == "gamma") {
        for(m in 1:nSpecies) { mean_stay[m] <- theta_stay[m] * exp(beta_stay[1] + species_effect_stay[m, 1]) }
      }
      if(stay_family == "lognormal") {
        for(m in 1:nSpecies) { mean_stay[m] <- exp(beta_stay[1] + species_effect_stay[m, 1] + theta_stay[m] ^ 2 / 2) }
      }
      if(stay_family == "weibull") {
        for(m in 1:nSpecies) { mean_stay[m] <-  lgamma(1 + 1 / theta_stay[m]) + exp(beta_stay[1] + species_effect_stay[m, 1]) }
      }
    }
    if (nPreds_stay > 1) {
      if(stay_family == "exponential") {
        for(m in 1:nSpecies) {
          for(i in 1:N_station) { mean_stay[i, m] <- exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay])) }
        }
      }
      if(stay_family == "gamma") {
        for(m in 1:nSpecies) {
          for(i in 1:N_station) { mean_stay[i, m] <- theta_stay[m] * exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay])) }
        }
      }
      if(stay_family == "lognormal") {
        for(m in 1:nSpecies) {
          for(i in 1:N_station) { mean_stay[i, m] <- exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + theta_stay[m] ^ 2 / 2) }
        }
      }
      if(stay_family == "weibull") {
        for(m in 1:nSpecies) {
          for(i in 1:N_station) { mean_stay[i, m] <-  lgamma(1 + 1 / theta_stay[m]) * exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay[i, 1:nPreds_stay])) }
        }
      }
    }

    # Alpha (Enter) の計算と mean_pass の地点固有化
    if (nPreds_alpha == 1) {
      for (i in 1:N_station) {
        for (m in 1:nSpecies) {
          for (g in 1:N_group) {
            log(alpha_mat[i, m, g]) <- cutpoint[g] + beta_enter[1] + species_effect_alpha[m, 1]
          }
          alpha_sum[i, m] <- sum(alpha_mat[i, m, 1:N_group])
          for (g in 1:N_group) {
            p_expected[i, m, g] <- alpha_mat[i, m, g] / alpha_sum[i, m]
            c_expected[i, m, g] <- p_expected[i, m, g] * (g - 1)
          }
          mean_pass[i, m] <- sum(c_expected[i, m, 1:N_group])
        }
      }
    }
    if (nPreds_alpha > 1) {
      for (i in 1:N_station) {
        for (m in 1:nSpecies) {
          for (g in 1:N_group) {
            log(alpha_mat[i, m, g]) <- cutpoint[g] + inprod(beta_enter[1:nPreds_alpha] + species_effect_alpha[m, 1:nPreds_alpha], X_alpha[i, 1:nPreds_alpha])
          }
          alpha_sum[i, m] <- sum(alpha_mat[i, m, 1:N_group])
          for (g in 1:N_group) {
            p_expected[i, m, g] <- alpha_mat[i, m, g] / alpha_sum[i, m]
            c_expected[i, m, g] <- p_expected[i, m, g] * (g - 1)
          }
          mean_pass[i, m] <- sum(c_expected[i, m, 1:N_group])
        }
      }
    }


    # model for y (Dirichlet-multinomial)
    for (j in 1:N_station_species) {
      y[j, 1:N_group] ~ ddirchmulti(alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j])
      pred_y[j, 1:N_group] ~ ddirchmulti(alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j])
      loglike_obs_y[j] <- ddirchmulti(y[j, 1:N_group], alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j], log = 1)
      loglike_pred_y[j] <- ddirchmulti(pred_y[j, 1:N_group], alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j], log = 1)
    }

    # Priors for Alpha
    for (g in 1:N_group) {
      cutpoint[g] ~ dnorm(0, sd = 100)
    }
    for (k in 1:nPreds_alpha) {
      beta_enter[k] ~ dnorm(0, sd = 100)
      sd_species_alpha[k] ~ T(dnorm(0, sd = 100), 0, 5)
      for (m in 1:nSpecies) {
        species_effect_alpha[m, k] ~ dnorm(0, sd = sd_species_alpha[k])
      }
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
            log(mu[i, m]) <- log(density[m]) + log(S) + log(N_period[i]) - log(mean_stay[m]) + log(activity_proportion[m]) - log(mean_pass[i, m])
          } else {
            log(mu[i, m]) <- log(density[m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m]) - log(mean_pass[i, m])
          }
        }
        log(density[m]) <- beta_density + species_effect_density[m, 1]
      }
      beta_density ~ dnorm(0, sd = 100)
    }

    if(nPreds_density > 1) {
      for(m in 1:nSpecies) {
        for(i in 1:N_station) {
          if(nPreds_stay == 1) {
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[m]) + log(activity_proportion[m]) - log(mean_pass[i, m])
          } else {
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m]) - log(mean_pass[i, m])
          }
          log(density[i, m]) <- inprod(beta_density[1:nPreds_density] + species_effect_density[m, 1:nPreds_density], X_density[i, 1:nPreds_density])
        }
      }
      for(j in 1:nPreds_density) {
        beta_density[j] ~ dnorm(0, sd = 100)
      }
    }

    # Random effects for density
    for(m in 1:nSpecies) {
      for(k in 1:nPreds_density) {
        species_effect_density[m, k] ~ dnorm(0, sd_species_density)
      }
    }
    sd_species_density ~ T(dnorm(0, sd = 100), 0, 5)
  })


  cat("Compiling the model. This may take a moment...\n")
  # 修正後：打ち切りの場合のみ、c_timeより確実に大きい値を初期値として与える
  stay_inits <- rep(NA, N_stay)
  stay_inits[is.censored == 1] <- c_time[is.censored == 1] + 1.0

  inits_f <- function() {
    common_inits <- list(
      beta_stay = runif(nPreds_stay, -1, 1),
      stay = stay_inits,
      theta_stay = runif(nSpecies, 0.5, 4.0),
      species_effect_stay = matrix(runif(nSpecies * nPreds_stay, -1, 1), nrow = nSpecies, ncol = nPreds_stay),

      beta_density = rnorm(nPreds_density, 0, 1),
      species_effect_density = matrix(rnorm(nSpecies * nPreds_density, 0, 0.5), nrow = nSpecies, ncol = nPreds_density),
      sd_density = runif(1, 0.01, 2),
      sd_species_density = runif(1, 0.01, 2),

      beta_enter = rnorm(nPreds_alpha, 0, 0.1),
      cutpoint = rnorm(N_group, 0, 0.1),
      species_effect_alpha = matrix(rnorm(nSpecies * nPreds_alpha, 0, 0.1),
                                    nrow = nSpecies, ncol = nPreds_alpha),
      sd_species_alpha = runif(nPreds_alpha, 0.01, 1)
    )

    if (!is.null(random_effect_stay)) {
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

  # -------------------------------------------------------------------------
  # クラスターの準備と実行
  # -------------------------------------------------------------------------
  this_cluster <- makeCluster(cores)

  clusterEvalQ(this_cluster, {
    library(nimble)
    ddirchmulti <- nimbleFunction(
      run = function(x = double(1), alpha = double(1), size = double(0), log = integer(0)) {
        returnType(double(0))
        logProb <- lgamma(size + 1) - sum(lgamma(x + 1)) + lgamma(sum(alpha)) -
          sum(lgamma(alpha)) + sum(lgamma(alpha + x)) - lgamma(sum(alpha) + size)
        if (log) return(logProb) else return(exp(logProb))
      }
    )
    rdirchmulti <- nimbleFunction(
      run = function(n = integer(0), alpha = double(1), size = double(0)) {
        returnType(double(1))
        if (n != 1) print("rdirchmulti only allows n = 1; using n = 1.")
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

  # --- パラメータ（モニター対象）の設定 ---
  # 既存の基本パラメータ
  if(stay_family == "exponential") prms <- c("scale", "mean_stay")
  if(stay_family == "gamma" | stay_family == "weibull") prms <- c("scale", "shape", "mean_stay")
  if(stay_family == "lognormal") prms <- c("meanlog", "sdlog", "mean_stay")

  prms <- c(prms, "density", "mean_pass")

  # 集約に必要な基本セット
  params <- c(prms, "loglike_obs_stay", "loglike_obs_y", "loglike_obs_detection",
              "loglike_pred_detection", "loglike_pred_stay", "loglike_pred_y")

  # --- ここから重要：すべての回帰係数を追加 ---
  params <- c(params,
              "beta_density", "species_effect_density", # 密度用
              "beta_stay", "species_effect_stay",       # 滞在用
              "beta_enter", "species_effect_alpha")     # 進入用

  # activity_estimation が mixture の場合
  if(activity_estimation == "mixture") {
    params <- c(params, "activity_proportion")
  }

  # ここで run_MCMC_RAD もまとめてクラスターにエクスポート
  clusterExport(this_cluster, c("ddirchmulti", "rdirchmulti", "registerDistributions", "run_MCMC_RAD"), envir = environment())

  cat("Running MCMC sampling. Please wait...\n")
  chain_output <- parLapply(
    cl = this_cluster,
    X = per_chain_info,
    fun = run_MCMC_RAD,
    data = data_density,
    code = code,
    constants = cons_density,
    params = params,
    ni = iter,
    nt = thin,
    nb = warmup
  )

  stopCluster(this_cluster)
  cat("Estimation is finished!\n")
  ## WAIC
  loglfy <- MCMCchains(chain_output, params = c("loglike_obs_y"))
  loglfstay <- MCMCchains(chain_output, params = c("loglike_obs_stay"))
  loglfN <- MCMCchains(chain_output, params = c("loglike_obs_detection"))

  loglfall <- cbind(loglfstay, loglfy, loglfN) # kernel/mixture関係なく結合可能

  lppd <- sum(log(colMeans(exp(loglfall))))
  p.waic <- sum(apply(loglfall, 2, var))
  waic <- (-2) * lppd + 2 * p.waic
  # 結果の集約 -------------------------------------------------------------------

  # --- 共変量やランダム効果がない（全体共通）かの判定 ---
  check_no_cov <- function(f) {
    if (is.null(f)) return(TRUE)
    f <- as.formula(f)
    if (length(f) == 3) f <- f[-2]
    length(all.vars(f)) == 0
  }

  is_density_global <- check_no_cov(formula_density)
  is_stay_global    <- check_no_cov(formula_stay) && (is.null(random_effect_stay) || random_effect_stay == "NULL")
  is_enter_global   <- if (exists("formula_enter")) check_no_cov(formula_enter) else TRUE
  is_pass_global    <- (is_density_global && is_stay_global) || is_enter_global

  # --- mcmc_samples の構築 ---
  mcmc_samples <- coda::as.mcmc.list(chain_output)

  if (activity_estimation == "mixture") {
    sample_activity <- MCMCvis::MCMCchains(actv_chain_output,
                                           mcmc.list = TRUE,
                                           params = "activity_proportion")
    mcmc_samples <- lapply(seq_along(mcmc_samples), function(i) {
      coda::as.mcmc(cbind(mcmc_samples[[i]], sample_activity[[i]]))
    })
    mcmc_samples <- coda::as.mcmc.list(mcmc_samples)
  }

  # 全MCMCサンプルを行列として取得（種ごと係数計算に使用）
  all_samples_mat <- MCMCvis::MCMCchains(mcmc_samples)

  # -------------------------------------------------------------------------
  # ヘルパー1: MCMCsummary を実行して整形するラッパー
  #   exact = TRUE で完全一致検索し、余分なパラメータを拾わないようにする
  # -------------------------------------------------------------------------
  summarize_param <- function(param_name) {
    MCMCvis::MCMCsummary(
      MCMCvis::MCMCchains(mcmc_samples, mcmc.list = TRUE,
                          params = param_name, exact = TRUE),
      round = 4
    ) %>%
      tibble::rownames_to_column(var = "Variable") %>%
      tibble::as_tibble() %>%
      dplyr::rename(lower = `2.5%`, median = `50%`, upper = `97.5%`)
  }

  # -------------------------------------------------------------------------
  # ヘルパー2: beta + species_effect の事後サンプルを合算して種ごとに集約する
  #
  #  【stay / density】
  #    beta_stay[j]  または beta_density[j]  （nPreds==1 なら添字なし）
  #    + species_effect_stay[m, j] または species_effect_density[m, j]
  #    → coef_stay[j] (変数名) / coef_density[j] (変数名)  ×  種
  #
  #  【enter】
  #    beta_enter[j, g] + species_effect_alpha[m, j, g]
  #    → coef_enter[j, g] (変数名, catN)  ×  種
  #    ※ g はディリクレ多項分布のカテゴリ（0-indexed で cat0, cat1, ...）
  #
  #  切片列（j == 1）も出力する：切片の種差は基準密度・滞在時間の種差を意味する
  # -------------------------------------------------------------------------
  make_species_coef_summary <- function(param_type) {

    cfg <- switch(param_type,
                  stay    = list(nPreds = nPreds_stay,    col_names = colnames(X_stay),
                                 beta_pfx = "beta_stay",    eff_pfx = "species_effect_stay"),
                  density = list(nPreds = nPreds_density,  col_names = colnames(X_density),
                                 beta_pfx = "beta_density",  eff_pfx = "species_effect_density"),
                  enter   = list(nPreds = nPreds_alpha,    col_names = colnames(X_alpha),
                                 beta_pfx = "beta_enter",    eff_pfx = "species_effect_alpha")
    )

    rows <- list()

    if (param_type %in% c("stay", "density")) {

      for (m in seq_len(nSpecies)) {
        for (j in seq_len(cfg$nPreds)) {

          # nimble の出力列名：nPreds==1 のとき添字なし、複数のとき [j]
          beta_col <- if (cfg$nPreds == 1) {
            cfg$beta_pfx
          } else {
            paste0(cfg$beta_pfx, "[", j, "]")
          }
          eff_col <- paste0(cfg$eff_pfx, "[", m, ", ", j, "]")

          if (!(beta_col %in% colnames(all_samples_mat))) next
          if (!(eff_col  %in% colnames(all_samples_mat))) next

          samps <- all_samples_mat[, beta_col] + all_samples_mat[, eff_col]

          rows[[length(rows) + 1]] <- tibble::tibble(
            Species  = target_species[m],
            Station  = "All",
            Variable = paste0("coef_", param_type,
                              "[", j, "] (", cfg$col_names[j], ")"),
            mean   = mean(samps),
            sd     = sd(samps),
            lower  = quantile(samps, 0.025),
            median = quantile(samps, 0.500),
            upper  = quantile(samps, 0.975),
            Rhat   = NA_real_,
            n.eff  = NA_real_
          )
        }
      }

    } else if (param_type == "enter") {

      for (m in seq_len(nSpecies)) {
        for (j in seq_len(cfg$nPreds)) {
          for (g in seq_len(N_group)) {

            beta_col <- paste0(cfg$beta_pfx, "[", j, ", ", g, "]")
            eff_col  <- paste0(cfg$eff_pfx,  "[", m, ", ", j, ", ", g, "]")

            if (!(beta_col %in% colnames(all_samples_mat))) next
            if (!(eff_col  %in% colnames(all_samples_mat))) next

            samps <- all_samples_mat[, beta_col] + all_samples_mat[, eff_col]

            rows[[length(rows) + 1]] <- tibble::tibble(
              Species  = target_species[m],
              Station  = "All",
              Variable = paste0("coef_enter[", j, ", ", g, "] (",
                                cfg$col_names[j], ", cat", g - 1, ")"),
              mean   = mean(samps),
              sd     = sd(samps),
              lower  = quantile(samps, 0.025),
              median = quantile(samps, 0.500),
              upper  = quantile(samps, 0.975),
              Rhat   = NA_real_,
              n.eff  = NA_real_
            )
          }
        }
      }
    }

    dplyr::bind_rows(rows)
  }

  # --- 集約の準備 ---
  unique_stations <- station_effort_data$Station[1:N_station]

  # --- density の集約 ---
  raw_density <- summarize_param("density")

  if (nPreds_density == 1) {
    # density[m]：m = 種インデックス
    summary_density <- raw_density %>%
      tidyr::extract(Variable, into = "Species_idx",
                     regex = "\\[(\\d+)\\]", convert = TRUE, remove = FALSE) %>%
      dplyr::mutate(Species = target_species[Species_idx], Station = "All") %>%
      dplyr::select(-Species_idx)
  } else {
    # density[i, m]：i = 地点インデックス, m = 種インデックス
    summary_density <- raw_density %>%
      tidyr::extract(Variable, into = c("Station_idx", "Species_idx"),
                     regex = "\\[(\\d+),\\s*(\\d+)\\]", convert = TRUE, remove = FALSE) %>%
      dplyr::mutate(Species = target_species[Species_idx],
                    Station = unique_stations[Station_idx]) %>%
      dplyr::select(-Station_idx, -Species_idx)
  }

  # --- mean_stay の集約 ---
  raw_stay <- summarize_param("mean_stay")

  if (nPreds_stay == 1) {
    # mean_stay[m]
    summary_stay <- raw_stay %>%
      tidyr::extract(Variable, into = "Species_idx",
                     regex = "\\[(\\d+)\\]", convert = TRUE, remove = FALSE) %>%
      dplyr::mutate(Species = target_species[Species_idx], Station = "All") %>%
      dplyr::select(-Species_idx)
  } else {
    # mean_stay[i, m]
    summary_stay <- raw_stay %>%
      tidyr::extract(Variable, into = c("Station_idx", "Species_idx"),
                     regex = "\\[(\\d+),\\s*(\\d+)\\]", convert = TRUE, remove = FALSE) %>%
      dplyr::mutate(Species = target_species[Species_idx],
                    Station = unique_stations[Station_idx]) %>%
      dplyr::select(-Station_idx, -Species_idx)
  }

  # --- mean_pass の集約（常に [i, m] 形式） ---
  raw_pass <- summarize_param("mean_pass") %>%
    tidyr::extract(Variable, into = c("Station_idx", "Species_idx"),
                   regex = "\\[(\\d+),\\s*(\\d+)\\]", convert = TRUE, remove = FALSE)

  if (is_pass_global) {
    # 地点間で同値のため、地点1のみ採用して Station = "All" とする
    summary_pass <- raw_pass %>%
      dplyr::filter(Station_idx == 1) %>%
      dplyr::mutate(Species  = target_species[Species_idx],
                    Station  = "All",
                    Variable = paste0("mean_pass[", Species_idx, "]")) %>%
      dplyr::select(-Station_idx, -Species_idx)
  } else {
    summary_pass <- raw_pass %>%
      dplyr::mutate(Species = target_species[Species_idx],
                    Station = unique_stations[Station_idx]) %>%
      dplyr::select(-Station_idx, -Species_idx)
  }

  # --- 種ごとの実効係数（beta + species_effect）---
  # nPreds > 1 の場合のみ生成する（切片のみモデルでは不要）
  summary_coef_list <- list(summary_density, summary_stay, summary_pass)

  if (nPreds_stay    > 1) {
    summary_coef_list <- c(summary_coef_list,
                           list(make_species_coef_summary("stay")))
  }
  if (nPreds_density > 1) {
    summary_coef_list <- c(summary_coef_list,
                           list(make_species_coef_summary("density")))
  }
  if (nPreds_alpha   > 1) {
    summary_coef_list <- c(summary_coef_list,
                           list(make_species_coef_summary("enter")))
  }

  # --- 最終結合・CV計算（cv は絶対値） ---
  summary_mean <- dplyr::bind_rows(summary_coef_list) %>%
    dplyr::mutate(cv = abs(sd / mean)) %>%
    dplyr::select(Species, Station, Variable, mean, sd, cv, lower, median, upper, Rhat, n.eff)

  # --- 正規化パラメータの整形 ---
  names(scaling_stay$center)     <- colnames(X_stay)
  names(scaling_stay$scale)      <- colnames(X_stay)
  names(scaling_density$center)  <- colnames(X_density)
  names(scaling_density$scale)   <- colnames(X_density)
  names(scaling_alpha$center)    <- colnames(X_alpha)
  names(scaling_alpha$scale)     <- colnames(X_alpha)

  scaling_params <- list(
    stay    = scaling_stay,
    density = scaling_density,
    enter   = scaling_alpha
  )

  # --- 返り値の構築 ---
  density_result <- list(
    WAIC           = waic,
    summary_result = summary_mean,
    samples        = mcmc_samples,
    scaling_params = scaling_params
  )
  class(density_result) <- "ResultDensity"


  return(density_result)
} # 関数 bayes_rest_multi の終端
time <- Stay <- NULL
