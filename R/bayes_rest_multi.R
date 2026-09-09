#' Bayesian Parameter Estimation of a Multispecies REST/RAD-REST Model via MCMC Sampling using `nimble`
#'
#' @param formula_stay A model formula for staying times within the focal area (e.g., \code{Stay ~ 1 + x1}). The left-hand side must specify the column name for staying time. Random effects should be specified separately using the \code{random_effect_stay} argument.
#' @param formula_density A model formula for animal density (e.g., \code{~ 1 + x1}). Note that the left-hand side must be omitted as density is a latent parameter.
#' @param formula_enter A model formula for the number of passes through the focal area per video (e.g., \code{~ 1 + x1}). Required when \code{model = "RAD-REST"}; ignored when \code{model = "REST"}.
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
#'     \item \code{Cens} (binary): Indicator for censored staying time (1 = censored, 0 = observed).
#'   }
#' @param random_effect_stay A character string specifying the column name in \code{stay_data} to use as a grouping factor for a station-level random intercept on staying time (e.g., \code{"Station"}). Default is \code{NULL}. Note: species-level random effects are automatically included in this multispecies model.
#' @param activity_data A data frame containing a \code{time} column, representing detection times transformed into radians. Typically, this is the output of the \code{format_activity} function.
#' @param activity_estimation A character string specifying the method used to estimate activity patterns. Choose \code{"kernel"} for fixed kernel density estimation (Rowcliffe et al. 2014), or \code{"mixture"} for nonparametric Bayesian estimation using von Mises mixture models (Nakashima et al. 2025). Default is \code{"kernel"}.
#' @param bw_adj A numeric bandwidth adjustment parameter for kernel density estimation. Default is 1.0. See Rowcliffe et al. (2014) for details.
#' @param C An integer specifying the maximum number of von Mises components to use in the mixture model. Required only if \code{activity_estimation = "mixture"}. Default is 10.
#' @param stay_family A character string specifying the probability distribution used to model staying times (e.g., \code{"exponential"}, \code{"gamma"}, \code{"lognormal"}, \code{"weibull"}). Default is \code{"lognormal"}.
#' @param focal_area A numeric value representing the size of the focal area, in square meters.
#' @param cores An integer specifying the number of CPU cores to use for parallel computation. Default is 3.
#' @param iter An integer specifying the total number of MCMC iterations per chain. Default is 5000.
#' @param warmup An integer specifying the number of warm-up (burn-in) iterations per chain. Default is 1000.
#' @param chains An integer specifying the number of MCMC chains. Default is 3.
#' @param thin An integer specifying the thinning interval for MCMC sampling. Default is 4.
#' @param model A character string specifying the model to be used. Choose either \code{"REST"} or \code{"RAD-REST"}. Default is \code{"RAD-REST"}.
#' @param target_species A character vector specifying the species to be analyzed. Multiple species must be specified for this multispecies model.
#' @return A list of class \code{"ResultDensity"} with the following components:
#' \describe{
#'   \item{\code{WAIC}}{A numeric WAIC value for the fitted model.}
#'   \item{\code{summary_result}}{A data frame summarizing posterior estimates (mean, sd, lower, median, upper, Rhat, n.eff, cv) for density, mean_stay, and (for RAD-REST) mean_pass across species and stations.}
#'   \item{\code{samples}}{A \code{coda::mcmc.list} object of full MCMC posterior samples.}
#'   \item{\code{tidy_samples}}{A long-format data frame of all monitored MCMC samples, with columns \code{parameter}, \code{value}, and \code{iteration}.}
#'   \item{\code{scaling_params}}{A list of centering and scaling parameters used to standardize design matrices.}
#' }
#'
#' @export
#' @import nimble activity parallel MCMCvis tibble
#' @importFrom stats as.formula formula model.frame model.matrix sd var runif median quantile model.response rexp rnorm step dexp pexp dgamma pgamma dlnorm plnorm dweibull pweibull dnbinom delete.response terms
#' @importFrom dplyr select filter mutate arrange pull bind_rows rename
#' @importFrom purrr map
#' @examples
#' \dontrun{
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
#'   thin = 4,
#'   model = "RAD-REST",
#'   target_species = c("SP01", "SP02", "SP03")
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
                             model = "RAD-REST",
                             target_species) {

  # Helper functions ----------------------------------------------------------

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

  standardize_design_matrix <- function(X) {
    center_vec <- rep(0, ncol(X))
    scale_vec  <- rep(1, ncol(X))
    for (j in seq_len(ncol(X))) {
      col_j <- X[, j]
      if (all(col_j == 1)) next
      if (!is.numeric(col_j)) next
      col_mean <- mean(col_j, na.rm = TRUE)
      col_sd   <- stats::sd(col_j, na.rm = TRUE)
      if (is.na(col_sd) || col_sd == 0) {
        warning(paste0("Column '", colnames(X)[j], "' has zero variance and will not be scaled."))
        next
      }
      X[, j]        <- (col_j - col_mean) / col_sd
      center_vec[j] <- col_mean
      scale_vec[j]  <- col_sd
    }
    list(X = X, center = center_vec, scale = scale_vec)
  }

  check_no_cov <- function(f) {
    if (is.null(f)) return(TRUE)
    f <- stats::as.formula(f)
    vars <- all.vars(f[[length(f)]])
    return(length(vars) == 0)
  }

  safe_log_mean_exp <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA)
    max_val <- max(x)
    max_val + log(mean(exp(x - max_val)))
  }

  safe_var <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2) return(NA)
    stats::var(x)
  }

  # Input validation ----------------------------------------------------------

  if (!inherits(formula_stay, "formula"))
    stop("`formula_stay` must be a valid formula (e.g., Stay ~ 1 + x1).")
  if (!inherits(formula_density, "formula"))
    stop("`formula_density` must be a valid formula (e.g., ~ 1 + x1).")
  if (!model %in% c("REST", "RAD-REST"))
    stop("`model` must be either 'REST' or 'RAD-REST'.")
  if (model == "RAD-REST" && !inherits(formula_enter, "formula"))
    stop("`formula_enter` must be a valid model formula when model is 'RAD-REST'.")
  if (!is.data.frame(station_effort_data))
    stop("`station_effort_data` must be a data frame.")
  if (!is.data.frame(stay_data))
    stop("`stay_data` must be a data frame.")
  if (!is.character(random_effect_stay) && !is.null(random_effect_stay))
    stop("`random_effect_stay` must be a character vector or NULL.")
  if (!is.data.frame(activity_data))
    stop("`activity_data` must be a data frame.")
  if (!activity_estimation %in% c("kernel", "mixture"))
    stop("`activity_estimation` must be either 'kernel' or 'mixture'.")
  if (!is.numeric(bw_adj) || length(bw_adj) != 1 || bw_adj <= 0)
    stop("`bw_adj` must be a positive number.")
  if (activity_estimation == "mixture" && (missing(C) || !is.numeric(C) || C <= 0))
    stop("`C` must be a positive integer when `activity_estimation = 'mixture'`.")
  if (!stay_family %in% c("lognormal", "gamma", "weibull", "exponential"))
    stop(paste0("Input stay_family type (", stay_family, ") is incorrect."))
  if (!is.numeric(focal_area) || length(focal_area) != 1 || focal_area <= 0)
    stop("`focal_area` must be a positive number.")
  if (!is.numeric(cores)  || cores  < 1) stop("`cores` must be a positive integer.")
  if (!is.numeric(iter)   || iter   <= 0) stop("`iter` must be a positive integer.")
  if (!is.numeric(warmup) || warmup <= 0) stop("`warmup` must be a positive integer.")
  if (!is.numeric(chains) || chains <= 0) stop("`chains` must be a positive integer.")
  if (!is.numeric(thin)   || thin   <= 0) stop("`thin` must be a positive integer.")
  if (!is.character(target_species) || length(target_species) < 2)
    stop("`target_species` must be a character vector with at least two species.")

  # MCMC settings
  ni <- iter
  nt <- thin
  nc <- chains
  nb <- warmup

  # Activity proportion estimation (mixture) ---------------------------------

  actv_out_trace <- list(0)

  if (activity_estimation == "mixture") {
    out_trace_0 <- list()
    for (m in seq_along(target_species)) {
      act_data <- activity_data %>%
        dplyr::filter(Species == target_species[m]) %>%
        dplyr::pull(time)
      dens.x <- seq(0, 2 * pi, 0.02)
      ndens  <- length(dens.x)
      N_act  <- length(act_data)
      constants_act <- list(N = N_act, C = C, dens.x = dens.x, ndens = ndens)
      data_act <- list(act_data = act_data)

      code_act <- nimble::nimbleCode({
        for (k in 1:(C - 1)) { v[k] ~ dbeta(1, alpha) }
        alpha ~ dgamma(1, 1)
        w[1:C] <- stick_breaking(v[1:(C - 1)])
        for (k in 1:C) {
          mu_mix[k]    ~ dunif(0, 2 * 3.141592654)
          kappa_mix[k] ~ dgamma(1, 0.01)
        }
        for (n in 1:N) {
          group[n] ~ dcat(w[1:C])
          act_data[n] ~ dvonMises(mu_mix[group[n]], kappa_mix[group[n]])
          loglike_obs_act[n] <- dvonMises(act_data[n], mu_mix[group[n]], kappa_mix[group[n]], log = 1)
        }
        for (j in 1:ndens) {
          for (i in 1:C) {
            dens.cpt[i, j] <- w[i] * dvonMises(dens.x[j], mu_mix[i], kappa_mix[i], log = 0)
          }
          activity_density[j] <- sum(dens.cpt[1:C, j])
        }
        activity_proportion <- 1.0 / (2 * 3.141592654 * max(activity_density[1:ndens]))
      })

      inits_act <- function() {
        list(mu_mix = stats::runif(C, 0, 2 * pi),
             kappa_mix = stats::rgamma(C, 1, 0.01),
             group = sample(1:C, size = N_act, replace = TRUE),
             v = stats::rbeta(C - 1, 1, 1),
             alpha = 1)
      }

      run_MCMC_vonMises <- function(info, data, constants, code, params, ni, nt, nb) {
        myModel   <- nimble::nimbleModel(code = code, data = data, constants = constants, inits = info$inits)
        CmyModel  <- nimble::compileNimble(myModel)
        configModel <- nimble::configureMCMC(myModel, monitors = params)
        myMCMC    <- nimble::buildMCMC(configModel, monitors = params)
        CmyMCMC   <- nimble::compileNimble(myMCMC)
        nimble::runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt,
                        nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
      }

      per_chain_act <- lapply(seq_len(nc), function(i) list(seed = sample(1:9999, 1), inits = inits_act()))
      params_act <- c("activity_density", "activity_proportion", "mu_mix", "kappa_mix", "w")

      this_cluster_act <- parallel::makeCluster(nc)
      on.exit(try(parallel::stopCluster(this_cluster_act), silent = TRUE), add = TRUE)

      parallel::clusterEvalQ(this_cluster_act, {
        library(nimble)
        dvonMises <- nimble::nimbleFunction(
          run = function(x = double(0), kappa = double(0), mu = double(0), log = integer(0)) {
            returnType(double(0))
            ccrit <- 1E-6; s <- 1; i <- 1; inc <- 1; x_2i <- 0; satisfied <- FALSE
            while (!satisfied) {
              x_2i <- kappa / (2 * i); inc <- inc * x_2i * x_2i; s <- s + inc; i <- i + 1
              satisfied <- inc < ccrit
            }
            prob <- exp(kappa * cos(x - mu)) / (2 * pi * s)
            if (log) return(log(prob)) else return(prob)
          }
        )
        rvonMises <- nimble::nimbleFunction(
          run = function(n = integer(0), kappa = double(0), mu = double(0)) {
            returnType(double(0)); return(0)
          }
        )
        suppressMessages(nimble::registerDistributions(list(
          dvonMises = list(
            BUGSdist = "dvonMises(kappa, mu)",
            types    = c("value = double(0)", "kappa = double(0)", "mu = double(0)"),
            pqAvail  = FALSE
          )
        )))
      })

      parallel::clusterExport(this_cluster_act, c("run_MCMC_vonMises"), envir = environment())
      actv_chain_output <- parallel::parLapply(
        cl = this_cluster_act, X = per_chain_act,
        fun = run_MCMC_vonMises, data = data_act, code = code_act,
        constants = constants_act, params = params_act,
        ni = ni, nt = nt, nb = nb)
      parallel::stopCluster(this_cluster_act)
      on.exit(NULL, add = FALSE)

      out_trace_0[[m]] <- actv_chain_output %>%
        purrr::map(~ .[, grep("activity_proportion", colnames(.), fixed = TRUE), drop = FALSE])
    }

    for (j in 1:nc) {
      actv_list   <- lapply(out_trace_0, function(sp) as.vector(sp[[j]]))
      actv_matrix <- do.call(cbind, actv_list)
      colnames(actv_matrix) <- paste0("activity_proportion[", seq_along(target_species), "]")
      actv_out_trace[[j]] <- coda::mcmc(actv_matrix)
    }
  }

  # Data preparation ----------------------------------------------------------

  target_species <- sort(target_species)

  station_effort_data <- station_effort_data %>%
    dplyr::filter(Species %in% target_species) %>%
    dplyr::arrange(Species, Station)

  stay_data_join <- stay_data %>%
    dplyr::filter(Species %in% target_species) %>%
    dplyr::left_join(station_effort_data,
                     by = intersect(names(stay_data), names(station_effort_data))) %>%
    dplyr::filter(!is.na(Stay)) %>%
    dplyr::arrange(Species, Station)

  N_station_species <- nrow(station_effort_data)
  N_station         <- length(unique(station_effort_data$Station))
  nSpecies          <- length(target_species)
  unique_stations   <- unique(station_effort_data$Station[order(match(station_effort_data$Station,
                                                                       unique(station_effort_data$Station)))])[1:N_station]
  unique_stations   <- station_effort_data$Station[!duplicated(station_effort_data$Station)][1:N_station]

  # Stay data
  model_frame_stay <- stats::model.frame(formula_stay, stay_data_join)
  X_stay_raw       <- stats::model.matrix(stats::as.formula(formula_stay), model_frame_stay)
  stay             <- stats::model.response(model_frame_stay)
  censored         <- stay_data_join[["Cens"]]
  is.censored      <- censored

  scaled_stay  <- standardize_design_matrix(X_stay_raw)
  X_stay       <- scaled_stay$X
  scaling_stay <- list(center = scaled_stay$center, scale = scaled_stay$scale)
  nPreds_stay  <- ncol(X_stay)
  if (nPreds_stay > 1) cat("Covariates in formula_stay have been standardized (mean=0, sd=1).\n")

  names(stay) <- NULL
  c_time <- stay
  c_time[is.censored == 0] <- c_time[is.censored == 0] + 1
  stay[is.censored == 1]   <- NA
  N_stay <- length(stay)

  # X_stay_station: station-level stay covariates (first N_station rows)
  formula_stay_rhs     <- stats::delete.response(stats::terms(stats::as.formula(formula_stay)))
  model_frame_stay_st  <- stats::model.frame(formula_stay_rhs, station_effort_data,
                                              na.action = stats::na.pass)
  X_stay_station_raw   <- stats::model.matrix(formula_stay_rhs, model_frame_stay_st)
  X_stay_station_raw   <- X_stay_station_raw[1:N_station, , drop = FALSE]
  X_stay_station       <- scale(X_stay_station_raw,
                                center = scaling_stay$center,
                                scale  = scaling_stay$scale)
  X_stay_station       <- matrix(as.numeric(X_stay_station),
                                 nrow = N_station, ncol = ncol(X_stay_station_raw))

  # Random effects for stay
  if (!is.null(random_effect_stay)) {
    re_levels_stay  <- unique(stay_data_join[[random_effect_stay]])
    nLevels_stay    <- length(re_levels_stay)
  } else {
    nLevels_stay <- 0
  }

  # Density design matrix
  model_frame_density <- stats::model.frame(formula_density, station_effort_data)
  X_density_raw       <- stats::model.matrix(stats::as.formula(formula_density), model_frame_density)
  X_density_raw       <- X_density_raw[1:N_station, , drop = FALSE]
  scaled_density      <- standardize_design_matrix(X_density_raw)
  X_density           <- scaled_density$X
  scaling_density     <- list(center = scaled_density$center, scale = scaled_density$scale)
  nPreds_density      <- ncol(X_density)
  if (nPreds_density > 1) cat("Covariates in formula_density have been standardized (mean=0, sd=1).\n")

  # Activity proportions (kernel)
  if (activity_estimation == "kernel") {
    activity_proportion <- numeric(nSpecies)
    for (i in seq_len(nSpecies)) {
      act_sp <- activity_data %>% dplyr::filter(Species == target_species[i])
      model_act <- activity::fitact(act_sp %>% dplyr::pull(time),
                                    bw = bw_adj * activity::bwcalc(act_sp %>% dplyr::pull(time), K = 3),
                                    reps = 1)
      activity_proportion[i] <- model_act@act
    }
  }

  S        <- focal_area * 1e-6
  N_period <- station_effort_data$Effort[1:N_station] * 60 * 60 * 24

  species_id_stay    <- as.numeric(factor(stay_data_join$Species, levels = target_species))
  species_id_ey      <- as.numeric(factor(station_effort_data$Species, levels = target_species))
  station_id_ey      <- as.numeric(factor(station_effort_data$Station,
                                          levels = unique(station_effort_data$Station[1:N_station])))

  # Model-specific data
  if (model == "RAD-REST") {
    # Enter design matrix
    model_frame_enter <- stats::model.frame(formula_enter, station_effort_data)
    X_alpha_raw       <- stats::model.matrix(stats::as.formula(formula_enter), model_frame_enter)
    X_alpha_raw       <- X_alpha_raw[1:N_station, , drop = FALSE]
    scaled_alpha      <- standardize_design_matrix(X_alpha_raw)
    X_alpha           <- scaled_alpha$X
    scaling_alpha     <- list(center = scaled_alpha$center, scale = scaled_alpha$scale)
    nPreds_alpha      <- ncol(X_alpha)
    if (nPreds_alpha > 1) cat("Covariates in formula_enter have been standardized (mean=0, sd=1).\n")

    y_mat        <- station_effort_data %>% dplyr::select(dplyr::starts_with("y_")) %>% as.matrix()
    N_judge      <- apply(y_mat, 1, sum)
    N_group      <- ncol(y_mat)
    N_detection  <- station_effort_data %>% dplyr::pull(N)
    N_detection_matrix <- matrix(N_detection, nrow = N_station, ncol = nSpecies, byrow = FALSE)

    data_density <- list(
      stay             = stay,
      is.censored      = is.censored,
      y                = y_mat,
      N_judge          = N_judge,
      N_detection_matrix = N_detection_matrix
    )
    if (nPreds_density > 1) data_density$X_density    <- X_density
    if (nPreds_alpha   > 1) data_density$X_alpha      <- X_alpha
    if (nPreds_stay    > 1) {
      data_density$X_stay         <- X_stay
      data_density$X_stay_station <- X_stay_station
    }

    cons_density <- list(
      N_stay            = N_stay,
      nPreds_stay       = nPreds_stay,
      N_station_species = N_station_species,
      N_station         = N_station,
      c_time            = c_time,
      stay_family       = stay_family,
      N_group           = N_group,
      nSpecies          = nSpecies,
      species_id_stay   = species_id_stay,
      species_id_ey     = species_id_ey,
      station_id_ey     = station_id_ey,
      S                 = S,
      N_period          = N_period,
      nPreds_density    = nPreds_density,
      nPreds_alpha      = nPreds_alpha,
      nLevels_stay      = nLevels_stay
    )
    if (activity_estimation == "kernel") cons_density$activity_proportion <- activity_proportion
    if (!is.null(random_effect_stay)) {
      cons_density$group_stay <- as.numeric(factor(stay_data_join[[random_effect_stay]],
                                                    levels = re_levels_stay))
    }

  } else {
    # REST mode
    scaling_alpha <- list(center = numeric(0), scale = numeric(0))
    nPreds_alpha  <- 0L
    N_group       <- 0L

    Y_vec    <- station_effort_data %>% dplyr::pull(Y)
    Y_matrix <- matrix(Y_vec, nrow = N_station, ncol = nSpecies, byrow = FALSE)

    data_density <- list(
      stay        = stay,
      is.censored = is.censored,
      Y_matrix    = Y_matrix
    )
    if (nPreds_density > 1) data_density$X_density    <- X_density
    if (nPreds_stay    > 1) {
      data_density$X_stay         <- X_stay
      data_density$X_stay_station <- X_stay_station
    }

    cons_density <- list(
      N_stay            = N_stay,
      nPreds_stay       = nPreds_stay,
      N_station_species = N_station_species,
      N_station         = N_station,
      c_time            = c_time,
      stay_family       = stay_family,
      nSpecies          = nSpecies,
      species_id_stay   = species_id_stay,
      S                 = S,
      N_period          = N_period,
      nPreds_density    = nPreds_density,
      nLevels_stay      = nLevels_stay
    )
    if (activity_estimation == "kernel") cons_density$activity_proportion <- activity_proportion
    if (!is.null(random_effect_stay)) {
      cons_density$group_stay <- as.numeric(factor(stay_data_join[[random_effect_stay]],
                                                    levels = re_levels_stay))
    }
  }

  # NIMBLE model code generators ---------------------------------------------

  get_RADREST_code <- function(stay_family) {

    if (stay_family == "exponential") {
      nimble::nimbleCode({
        # [1] Stay model: Exponential
        for (i in 1:N_stay) {
          censored[i] ~ dinterval(stay[i], c_time[i])
          stay[i]   ~ dexp(rate = 1 / scale[i])
          pred_t[i] ~ dexp(rate = 1 / scale[i])
          loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dexp(stay[i], rate = 1 / scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pexp(c_time[i], rate = 1 / scale[i]))
          loglike_pred_stay[i] <- dexp(pred_t[i], rate = 1 / scale[i], log = 1)
          if (nPreds_stay > 1) {
            if (nLevels_stay == 0) { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) }
            else { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]] }
          } else {
            if (nLevels_stay == 0) { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] }
            else { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] + random_effect_stay[group_stay[i]] }
          }
        }
        if (nPreds_stay == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- exp(beta_stay[1] + species_effect_stay[m, 1]) } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay_station[i, 1:nPreds_stay])) } }
        }
        # Stay priors
        for (j in 1:nPreds_stay) {
          beta_stay[j] ~ dnorm(0, sd = 5)
          sigma_species_stay[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_stay[m, j] ~ dnorm(0, sd = sigma_species_stay[j]) }
        }
        if (nLevels_stay > 0) {
          for (k in 1:nLevels_stay) { random_effect_stay[k] ~ dnorm(0, sd = sigma_stay) }
          sigma_stay ~ T(dnorm(0, sd = 2), 0, )
        }
        # [2] Enter model
        for (m in 1:nSpecies) { theta_enter[m] ~ dgamma(2, 2) }
        cutpoint[1] <- 0
        for (g in 2:N_group) { cutpoint[g] ~ dnorm(0, sd = 5) }
        if (nPreds_alpha == 1) {
          for (i in 1:N_station) { for (m in 1:nSpecies) {
            eta[i, m] <- beta_enter[1] + species_effect_alpha[m, 1]
            for (g in 1:N_group) { log_phi[i, m, g] <- cutpoint[g] + (g - 1) * eta[i, m]; phi[i, m, g] <- exp(log_phi[i, m, g]) }
            sum_phi[i, m] <- sum(phi[i, m, 1:N_group])
            for (g in 1:N_group) { p_expected[i, m, g] <- phi[i, m, g] / sum_phi[i, m]; c_expected[i, m, g] <- p_expected[i, m, g] * (g - 1); alpha_mat[i, m, g] <- theta_enter[m] * p_expected[i, m, g] }
            mean_pass[i, m] <- sum(c_expected[i, m, 1:N_group])
          } }
        } else {
          for (i in 1:N_station) { for (m in 1:nSpecies) {
            eta[i, m] <- inprod(beta_enter[1:nPreds_alpha] + species_effect_alpha[m, 1:nPreds_alpha], X_alpha[i, 1:nPreds_alpha])
            for (g in 1:N_group) { log_phi[i, m, g] <- cutpoint[g] + (g - 1) * eta[i, m]; phi[i, m, g] <- exp(log_phi[i, m, g]) }
            sum_phi[i, m] <- sum(phi[i, m, 1:N_group])
            for (g in 1:N_group) { p_expected[i, m, g] <- phi[i, m, g] / sum_phi[i, m]; c_expected[i, m, g] <- p_expected[i, m, g] * (g - 1); alpha_mat[i, m, g] <- theta_enter[m] * p_expected[i, m, g] }
            mean_pass[i, m] <- sum(c_expected[i, m, 1:N_group])
          } }
        }
        for (k in 1:nPreds_alpha) {
          beta_enter[k] ~ dnorm(0, sd = 5)
          sd_species_alpha[k] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_alpha[m, k] ~ dnorm(0, sd = sd_species_alpha[k]) }
        }
        # y model
        for (j in 1:N_station_species) {
          y[j, 1:N_group]      ~ ddirchmulti(alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j])
          pred_y[j, 1:N_group] ~ ddirchmulti(alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j])
          loglike_obs_y[j]  <- ddirchmulti(y[j, 1:N_group], alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j], log = 1)
          loglike_pred_y[j] <- ddirchmulti(pred_y[j, 1:N_group], alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j], log = 1)
        }
        # N_detection model
        for (m in 1:nSpecies) {
          for (i in 1:N_station) {
            N_detection_matrix[i, m] ~ dnbinom(size = size[m], prob = p[i, m])
            p[i, m] <- size[m] / (size[m] + mu[i, m])
            N_detection_rep[i, m]   ~ dnbinom(size = size[m], prob = p[i, m])
            loglike_obs_detection[i, m]  <- dnbinom(N_detection_matrix[i, m], size[m], p[i, m], log = 1)
            loglike_pred_detection[i, m] <- dnbinom(N_detection_rep[i, m],    size[m], p[i, m], log = 1)
          }
          size[m] ~ dgamma(1, 1)
        }
        # Density and REST formula
        if (nPreds_density == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- beta_density[1] + species_effect_density[m, 1]
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m]) - log(mean_pass[i, m])
          } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- inprod(beta_density[1:nPreds_density] + species_effect_density[m, 1:nPreds_density], X_density[i, 1:nPreds_density])
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m]) - log(mean_pass[i, m])
          } }
        }
        for (j in 1:nPreds_density) {
          beta_density[j] ~ dnorm(0, sd = 5)
          sd_species_density[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_density[m, j] ~ dnorm(0, sd = sd_species_density[j]) }
        }
      })

    } else if (stay_family == "gamma") {
      nimble::nimbleCode({
        # [1] Stay model: Gamma
        for (i in 1:N_stay) {
          censored[i] ~ dinterval(stay[i], c_time[i])
          stay[i]   ~ dgamma(shape = theta_stay[species_id_stay[i]], rate = 1 / scale[i])
          pred_t[i] ~ dgamma(shape = theta_stay[species_id_stay[i]], rate = 1 / scale[i])
          loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dgamma(stay[i], shape = theta_stay[species_id_stay[i]], rate = 1 / scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pgamma(c_time[i], shape = theta_stay[species_id_stay[i]], rate = 1 / scale[i]))
          loglike_pred_stay[i] <- dgamma(pred_t[i], shape = theta_stay[species_id_stay[i]], rate = 1 / scale[i], log = 1)
          if (nPreds_stay > 1) {
            if (nLevels_stay == 0) { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) }
            else { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]] }
          } else {
            if (nLevels_stay == 0) { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] }
            else { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] + random_effect_stay[group_stay[i]] }
          }
        }
        for (m in 1:nSpecies) { theta_stay[m] ~ dgamma(shape_stay, rate_stay); shape[m] <- theta_stay[m] }
        shape_stay ~ dgamma(2, 0.5); rate_stay ~ dgamma(2, 0.5)
        if (nPreds_stay == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- theta_stay[m] * exp(beta_stay[1] + species_effect_stay[m, 1]) } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- theta_stay[m] * exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay_station[i, 1:nPreds_stay])) } }
        }
        for (j in 1:nPreds_stay) {
          beta_stay[j] ~ dnorm(0, sd = 5)
          sigma_species_stay[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_stay[m, j] ~ dnorm(0, sd = sigma_species_stay[j]) }
        }
        if (nLevels_stay > 0) {
          for (k in 1:nLevels_stay) { random_effect_stay[k] ~ dnorm(0, sd = sigma_stay) }
          sigma_stay ~ T(dnorm(0, sd = 2), 0, )
        }
        for (m in 1:nSpecies) { theta_enter[m] ~ dgamma(2, 2) }
        cutpoint[1] <- 0
        for (g in 2:N_group) { cutpoint[g] ~ dnorm(0, sd = 5) }
        if (nPreds_alpha == 1) {
          for (i in 1:N_station) { for (m in 1:nSpecies) {
            eta[i, m] <- beta_enter[1] + species_effect_alpha[m, 1]
            for (g in 1:N_group) { log_phi[i, m, g] <- cutpoint[g] + (g - 1) * eta[i, m]; phi[i, m, g] <- exp(log_phi[i, m, g]) }
            sum_phi[i, m] <- sum(phi[i, m, 1:N_group])
            for (g in 1:N_group) { p_expected[i, m, g] <- phi[i, m, g] / sum_phi[i, m]; c_expected[i, m, g] <- p_expected[i, m, g] * (g - 1); alpha_mat[i, m, g] <- theta_enter[m] * p_expected[i, m, g] }
            mean_pass[i, m] <- sum(c_expected[i, m, 1:N_group])
          } }
        } else {
          for (i in 1:N_station) { for (m in 1:nSpecies) {
            eta[i, m] <- inprod(beta_enter[1:nPreds_alpha] + species_effect_alpha[m, 1:nPreds_alpha], X_alpha[i, 1:nPreds_alpha])
            for (g in 1:N_group) { log_phi[i, m, g] <- cutpoint[g] + (g - 1) * eta[i, m]; phi[i, m, g] <- exp(log_phi[i, m, g]) }
            sum_phi[i, m] <- sum(phi[i, m, 1:N_group])
            for (g in 1:N_group) { p_expected[i, m, g] <- phi[i, m, g] / sum_phi[i, m]; c_expected[i, m, g] <- p_expected[i, m, g] * (g - 1); alpha_mat[i, m, g] <- theta_enter[m] * p_expected[i, m, g] }
            mean_pass[i, m] <- sum(c_expected[i, m, 1:N_group])
          } }
        }
        for (k in 1:nPreds_alpha) {
          beta_enter[k] ~ dnorm(0, sd = 5)
          sd_species_alpha[k] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_alpha[m, k] ~ dnorm(0, sd = sd_species_alpha[k]) }
        }
        for (j in 1:N_station_species) {
          y[j, 1:N_group]      ~ ddirchmulti(alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j])
          pred_y[j, 1:N_group] ~ ddirchmulti(alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j])
          loglike_obs_y[j]  <- ddirchmulti(y[j, 1:N_group], alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j], log = 1)
          loglike_pred_y[j] <- ddirchmulti(pred_y[j, 1:N_group], alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j], log = 1)
        }
        for (m in 1:nSpecies) {
          for (i in 1:N_station) {
            N_detection_matrix[i, m] ~ dnbinom(size = size[m], prob = p[i, m])
            p[i, m] <- size[m] / (size[m] + mu[i, m])
            N_detection_rep[i, m]   ~ dnbinom(size = size[m], prob = p[i, m])
            loglike_obs_detection[i, m]  <- dnbinom(N_detection_matrix[i, m], size[m], p[i, m], log = 1)
            loglike_pred_detection[i, m] <- dnbinom(N_detection_rep[i, m],    size[m], p[i, m], log = 1)
          }
          size[m] ~ dgamma(1, 1)
        }
        if (nPreds_density == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- beta_density[1] + species_effect_density[m, 1]
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m]) - log(mean_pass[i, m])
          } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- inprod(beta_density[1:nPreds_density] + species_effect_density[m, 1:nPreds_density], X_density[i, 1:nPreds_density])
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m]) - log(mean_pass[i, m])
          } }
        }
        for (j in 1:nPreds_density) {
          beta_density[j] ~ dnorm(0, sd = 5)
          sd_species_density[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_density[m, j] ~ dnorm(0, sd = sd_species_density[j]) }
        }
      })

    } else if (stay_family == "lognormal") {
      nimble::nimbleCode({
        # [1] Stay model: Lognormal
        for (i in 1:N_stay) {
          censored[i] ~ dinterval(stay[i], c_time[i])
          stay[i]   ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]])
          pred_t[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]])
          loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dlnorm(stay[i], meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]], log = 1) +
            step(censored[i] - 0.5) * log(1 - plnorm(c_time[i], meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]]))
          loglike_pred_stay[i] <- dlnorm(pred_t[i], meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]], log = 1)
          meanlog[i] <- log(scale[i])
          if (nPreds_stay > 1) {
            if (nLevels_stay == 0) { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) }
            else { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]] }
          } else {
            if (nLevels_stay == 0) { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] }
            else { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] + random_effect_stay[group_stay[i]] }
          }
        }
        for (m in 1:nSpecies) { theta_stay[m] ~ dgamma(shape_stay, rate_stay); sdlog[m] <- theta_stay[m] }
        shape_stay ~ dgamma(2, 0.5); rate_stay ~ dgamma(2, 0.5)
        if (nPreds_stay == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- exp(beta_stay[1] + species_effect_stay[m, 1] + theta_stay[m]^2 / 2) } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay_station[i, 1:nPreds_stay]) + theta_stay[m]^2 / 2) } }
        }
        for (j in 1:nPreds_stay) {
          beta_stay[j] ~ dnorm(0, sd = 5)
          sigma_species_stay[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_stay[m, j] ~ dnorm(0, sd = sigma_species_stay[j]) }
        }
        if (nLevels_stay > 0) {
          for (k in 1:nLevels_stay) { random_effect_stay[k] ~ dnorm(0, sd = sigma_stay) }
          sigma_stay ~ T(dnorm(0, sd = 2), 0, )
        }
        # [2] Enter model (cumulative logit)
        for (m in 1:nSpecies) { theta_enter[m] ~ dgamma(2, 2) }
        cutpoint[1] ~ dnorm(0, sd = 3)
        for (g in 2:(N_group - 1)) { delta[g - 1] ~ dgamma(1, 1); cutpoint[g] <- cutpoint[g - 1] + delta[g - 1] }
        for (k in 1:nPreds_alpha) {
          beta_enter[k] ~ dnorm(0, sd = 5)
          sd_species_alpha[k] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_alpha[m, k] ~ dnorm(0, sd = sd_species_alpha[k]) }
        }
        if (nPreds_alpha == 1) {
          for (i in 1:N_station) { for (m in 1:nSpecies) {
            eta[i, m] <- beta_enter[1] + species_effect_alpha[m, 1]
            for (g in 1:(N_group - 1)) { cum_p[i, m, g] <- ilogit(cutpoint[g] - eta[i, m]) }
            p_expected[i, m, 1] <- cum_p[i, m, 1]
            for (g in 2:(N_group - 1)) { p_expected[i, m, g] <- cum_p[i, m, g] - cum_p[i, m, g - 1] }
            p_expected[i, m, N_group] <- 1 - cum_p[i, m, N_group - 1]
            for (g in 1:N_group) { alpha_mat[i, m, g] <- theta_enter[m] * p_expected[i, m, g] }
            mean_pass_raw[i, m] <- sum(p_expected[i, m, 2:N_group] * (1:(N_group - 1)))
            mean_pass[i, m] <- max(mean_pass_raw[i, m], 0.001)
          } }
        } else {
          for (i in 1:N_station) { for (m in 1:nSpecies) {
            eta[i, m] <- inprod(beta_enter[1:nPreds_alpha] + species_effect_alpha[m, 1:nPreds_alpha], X_alpha[i, 1:nPreds_alpha])
            for (g in 1:(N_group - 1)) { cum_p[i, m, g] <- ilogit(cutpoint[g] - eta[i, m]) }
            p_expected[i, m, 1] <- cum_p[i, m, 1]
            for (g in 2:(N_group - 1)) { p_expected[i, m, g] <- cum_p[i, m, g] - cum_p[i, m, g - 1] }
            p_expected[i, m, N_group] <- 1 - cum_p[i, m, N_group - 1]
            for (g in 1:N_group) { alpha_mat[i, m, g] <- theta_enter[m] * p_expected[i, m, g] }
            mean_pass_raw[i, m] <- sum(p_expected[i, m, 2:N_group] * (1:(N_group - 1)))
            mean_pass[i, m] <- max(mean_pass_raw[i, m], 0.001)
          } }
        }
        for (j in 1:N_station_species) {
          y[j, 1:N_group]      ~ ddirchmulti(alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j])
          pred_y[j, 1:N_group] ~ ddirchmulti(alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j])
          loglike_obs_y[j]  <- ddirchmulti(y[j, 1:N_group], alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j], log = 1)
          loglike_pred_y[j] <- ddirchmulti(pred_y[j, 1:N_group], alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j], log = 1)
        }
        for (m in 1:nSpecies) {
          for (i in 1:N_station) {
            N_detection_matrix[i, m] ~ dnbinom(size = size[m], prob = p[i, m])
            p[i, m] <- size[m] / (size[m] + mu[i, m])
            N_detection_rep[i, m]   ~ dnbinom(size = size[m], prob = p[i, m])
            loglike_obs_detection[i, m]  <- dnbinom(N_detection_matrix[i, m], size[m], p[i, m], log = 1)
            loglike_pred_detection[i, m] <- dnbinom(N_detection_rep[i, m],    size[m], p[i, m], log = 1)
          }
          size[m] ~ dgamma(1, 1)
        }
        if (nPreds_density == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- beta_density[1] + species_effect_density[m, 1]
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m]) - log(mean_pass[i, m])
          } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- inprod(beta_density[1:nPreds_density] + species_effect_density[m, 1:nPreds_density], X_density[i, 1:nPreds_density])
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m]) - log(mean_pass[i, m])
          } }
        }
        for (j in 1:nPreds_density) {
          beta_density[j] ~ dnorm(0, sd = 5)
          sd_species_density[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_density[m, j] ~ dnorm(0, sd = sd_species_density[j]) }
        }
      })

    } else { # weibull
      nimble::nimbleCode({
        # [1] Stay model: Weibull
        for (i in 1:N_stay) {
          censored[i] ~ dinterval(stay[i], c_time[i])
          stay[i]   ~ dweibull(shape = theta_stay[species_id_stay[i]], scale = scale[i])
          pred_t[i] ~ dweibull(shape = theta_stay[species_id_stay[i]], scale = scale[i])
          loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dweibull(stay[i], shape = theta_stay[species_id_stay[i]], scale = scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pweibull(c_time[i], shape = theta_stay[species_id_stay[i]], scale = scale[i]))
          loglike_pred_stay[i] <- dweibull(pred_t[i], shape = theta_stay[species_id_stay[i]], scale = scale[i], log = 1)
          if (nPreds_stay > 1) {
            if (nLevels_stay == 0) { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) }
            else { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]] }
          } else {
            if (nLevels_stay == 0) { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] }
            else { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] + random_effect_stay[group_stay[i]] }
          }
        }
        for (m in 1:nSpecies) { theta_stay[m] ~ dgamma(shape_stay, rate_stay); shape[m] <- theta_stay[m] }
        shape_stay ~ dgamma(2, 0.5); rate_stay ~ dgamma(2, 0.5)
        if (nPreds_stay == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- exp(lgamma(1 + 1 / theta_stay[m]) + beta_stay[1] + species_effect_stay[m, 1]) } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- exp(lgamma(1 + 1 / theta_stay[m]) + inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay_station[i, 1:nPreds_stay])) } }
        }
        for (j in 1:nPreds_stay) {
          beta_stay[j] ~ dnorm(0, sd = 5)
          sigma_species_stay[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_stay[m, j] ~ dnorm(0, sd = sigma_species_stay[j]) }
        }
        if (nLevels_stay > 0) {
          for (k in 1:nLevels_stay) { random_effect_stay[k] ~ dnorm(0, sd = sigma_stay) }
          sigma_stay ~ T(dnorm(0, sd = 2), 0, )
        }
        # [2] Enter model (cumulative logit)
        for (m in 1:nSpecies) { theta_enter[m] ~ dgamma(2, 2) }
        cutpoint[1] ~ dnorm(0, sd = 3)
        for (g in 2:(N_group - 1)) { delta[g - 1] ~ dgamma(1, 1); cutpoint[g] <- cutpoint[g - 1] + delta[g - 1] }
        for (k in 1:nPreds_alpha) {
          beta_enter[k] ~ dnorm(0, sd = 5)
          sd_species_alpha[k] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_alpha[m, k] ~ dnorm(0, sd = sd_species_alpha[k]) }
        }
        if (nPreds_alpha == 1) {
          for (i in 1:N_station) { for (m in 1:nSpecies) {
            eta[i, m] <- beta_enter[1] + species_effect_alpha[m, 1]
            for (g in 1:(N_group - 1)) { cum_p[i, m, g] <- ilogit(cutpoint[g] - eta[i, m]) }
            p_expected[i, m, 1] <- cum_p[i, m, 1]
            for (g in 2:(N_group - 1)) { p_expected[i, m, g] <- cum_p[i, m, g] - cum_p[i, m, g - 1] }
            p_expected[i, m, N_group] <- 1 - cum_p[i, m, N_group - 1]
            for (g in 1:N_group) { alpha_mat[i, m, g] <- theta_enter[m] * p_expected[i, m, g] }
            mean_pass_raw[i, m] <- sum(p_expected[i, m, 2:N_group] * (1:(N_group - 1)))
            mean_pass[i, m] <- max(mean_pass_raw[i, m], 0.001)
          } }
        } else {
          for (i in 1:N_station) { for (m in 1:nSpecies) {
            eta[i, m] <- inprod(beta_enter[1:nPreds_alpha] + species_effect_alpha[m, 1:nPreds_alpha], X_alpha[i, 1:nPreds_alpha])
            for (g in 1:(N_group - 1)) { cum_p[i, m, g] <- ilogit(cutpoint[g] - eta[i, m]) }
            p_expected[i, m, 1] <- cum_p[i, m, 1]
            for (g in 2:(N_group - 1)) { p_expected[i, m, g] <- cum_p[i, m, g] - cum_p[i, m, g - 1] }
            p_expected[i, m, N_group] <- 1 - cum_p[i, m, N_group - 1]
            for (g in 1:N_group) { alpha_mat[i, m, g] <- theta_enter[m] * p_expected[i, m, g] }
            mean_pass_raw[i, m] <- sum(p_expected[i, m, 2:N_group] * (1:(N_group - 1)))
            mean_pass[i, m] <- max(mean_pass_raw[i, m], 0.001)
          } }
        }
        for (j in 1:N_station_species) {
          y[j, 1:N_group]      ~ ddirchmulti(alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j])
          pred_y[j, 1:N_group] ~ ddirchmulti(alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j])
          loglike_obs_y[j]  <- ddirchmulti(y[j, 1:N_group], alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j], log = 1)
          loglike_pred_y[j] <- ddirchmulti(pred_y[j, 1:N_group], alpha_mat[station_id_ey[j], species_id_ey[j], 1:N_group], N_judge[j], log = 1)
        }
        for (m in 1:nSpecies) {
          for (i in 1:N_station) {
            N_detection_matrix[i, m] ~ dnbinom(size = size[m], prob = p[i, m])
            p[i, m] <- size[m] / (size[m] + mu[i, m])
            N_detection_rep[i, m]   ~ dnbinom(size = size[m], prob = p[i, m])
            loglike_obs_detection[i, m]  <- dnbinom(N_detection_matrix[i, m], size[m], p[i, m], log = 1)
            loglike_pred_detection[i, m] <- dnbinom(N_detection_rep[i, m],    size[m], p[i, m], log = 1)
          }
          size[m] ~ dgamma(1, 1)
        }
        if (nPreds_density == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- beta_density[1] + species_effect_density[m, 1]
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m]) - log(mean_pass[i, m])
          } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- inprod(beta_density[1:nPreds_density] + species_effect_density[m, 1:nPreds_density], X_density[i, 1:nPreds_density])
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m]) - log(mean_pass[i, m])
          } }
        }
        for (j in 1:nPreds_density) {
          beta_density[j] ~ dnorm(0, sd = 5)
          sd_species_density[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_density[m, j] ~ dnorm(0, sd = sd_species_density[j]) }
        }
      })
    }
  } # end get_RADREST_code

  get_REST_code <- function(stay_family) {

    if (stay_family == "exponential") {
      nimble::nimbleCode({
        # [1] Stay model: Exponential
        for (i in 1:N_stay) {
          censored[i] ~ dinterval(stay[i], c_time[i])
          stay[i]   ~ dexp(rate = 1 / scale[i])
          pred_t[i] ~ dexp(rate = 1 / scale[i])
          loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dexp(stay[i], rate = 1 / scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pexp(c_time[i], rate = 1 / scale[i]))
          loglike_pred_stay[i] <- dexp(pred_t[i], rate = 1 / scale[i], log = 1)
          if (nPreds_stay > 1) {
            if (nLevels_stay == 0) { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) }
            else { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]] }
          } else {
            if (nLevels_stay == 0) { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] }
            else { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] + random_effect_stay[group_stay[i]] }
          }
        }
        if (nPreds_stay == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- exp(beta_stay[1] + species_effect_stay[m, 1]) } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay_station[i, 1:nPreds_stay])) } }
        }
        for (j in 1:nPreds_stay) {
          beta_stay[j] ~ dnorm(0, sd = 5)
          sigma_species_stay[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_stay[m, j] ~ dnorm(0, sd = sigma_species_stay[j]) }
        }
        if (nLevels_stay > 0) {
          for (k in 1:nLevels_stay) { random_effect_stay[k] ~ dnorm(0, sd = sigma_stay) }
          sigma_stay ~ T(dnorm(0, sd = 2), 0, )
        }
        # [2] REST detection model (Y_matrix = total passes)
        for (m in 1:nSpecies) {
          for (i in 1:N_station) {
            Y_matrix[i, m] ~ dnbinom(size = size[m], prob = p[i, m])
            p[i, m] <- size[m] / (size[m] + mu[i, m])
            Y_rep[i, m]    ~ dnbinom(size = size[m], prob = p[i, m])
            loglike_obs_y[i, m]  <- dnbinom(Y_matrix[i, m], size[m], p[i, m], log = 1)
            loglike_pred_y[i, m] <- dnbinom(Y_rep[i, m],    size[m], p[i, m], log = 1)
          }
          size[m] ~ dgamma(1, 1)
        }
        # REST formula (no mean_pass term)
        if (nPreds_density == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- beta_density[1] + species_effect_density[m, 1]
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m])
          } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- inprod(beta_density[1:nPreds_density] + species_effect_density[m, 1:nPreds_density], X_density[i, 1:nPreds_density])
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m])
          } }
        }
        for (j in 1:nPreds_density) {
          beta_density[j] ~ dnorm(0, sd = 5)
          sd_species_density[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_density[m, j] ~ dnorm(0, sd = sd_species_density[j]) }
        }
      })

    } else if (stay_family == "gamma") {
      nimble::nimbleCode({
        for (i in 1:N_stay) {
          censored[i] ~ dinterval(stay[i], c_time[i])
          stay[i]   ~ dgamma(shape = theta_stay[species_id_stay[i]], rate = 1 / scale[i])
          pred_t[i] ~ dgamma(shape = theta_stay[species_id_stay[i]], rate = 1 / scale[i])
          loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dgamma(stay[i], shape = theta_stay[species_id_stay[i]], rate = 1 / scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pgamma(c_time[i], shape = theta_stay[species_id_stay[i]], rate = 1 / scale[i]))
          loglike_pred_stay[i] <- dgamma(pred_t[i], shape = theta_stay[species_id_stay[i]], rate = 1 / scale[i], log = 1)
          if (nPreds_stay > 1) {
            if (nLevels_stay == 0) { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) }
            else { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]] }
          } else {
            if (nLevels_stay == 0) { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] }
            else { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] + random_effect_stay[group_stay[i]] }
          }
        }
        for (m in 1:nSpecies) { theta_stay[m] ~ dgamma(shape_stay, rate_stay); shape[m] <- theta_stay[m] }
        shape_stay ~ dgamma(2, 0.5); rate_stay ~ dgamma(2, 0.5)
        if (nPreds_stay == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- theta_stay[m] * exp(beta_stay[1] + species_effect_stay[m, 1]) } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- theta_stay[m] * exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay_station[i, 1:nPreds_stay])) } }
        }
        for (j in 1:nPreds_stay) {
          beta_stay[j] ~ dnorm(0, sd = 5)
          sigma_species_stay[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_stay[m, j] ~ dnorm(0, sd = sigma_species_stay[j]) }
        }
        if (nLevels_stay > 0) {
          for (k in 1:nLevels_stay) { random_effect_stay[k] ~ dnorm(0, sd = sigma_stay) }
          sigma_stay ~ T(dnorm(0, sd = 2), 0, )
        }
        for (m in 1:nSpecies) {
          for (i in 1:N_station) {
            Y_matrix[i, m] ~ dnbinom(size = size[m], prob = p[i, m])
            p[i, m] <- size[m] / (size[m] + mu[i, m])
            Y_rep[i, m]    ~ dnbinom(size = size[m], prob = p[i, m])
            loglike_obs_y[i, m]  <- dnbinom(Y_matrix[i, m], size[m], p[i, m], log = 1)
            loglike_pred_y[i, m] <- dnbinom(Y_rep[i, m],    size[m], p[i, m], log = 1)
          }
          size[m] ~ dgamma(1, 1)
        }
        if (nPreds_density == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- beta_density[1] + species_effect_density[m, 1]
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m])
          } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- inprod(beta_density[1:nPreds_density] + species_effect_density[m, 1:nPreds_density], X_density[i, 1:nPreds_density])
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m])
          } }
        }
        for (j in 1:nPreds_density) {
          beta_density[j] ~ dnorm(0, sd = 5)
          sd_species_density[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_density[m, j] ~ dnorm(0, sd = sd_species_density[j]) }
        }
      })

    } else if (stay_family == "lognormal") {
      nimble::nimbleCode({
        for (i in 1:N_stay) {
          censored[i] ~ dinterval(stay[i], c_time[i])
          stay[i]   ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]])
          pred_t[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]])
          loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dlnorm(stay[i], meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]], log = 1) +
            step(censored[i] - 0.5) * log(1 - plnorm(c_time[i], meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]]))
          loglike_pred_stay[i] <- dlnorm(pred_t[i], meanlog = log(scale[i]), sdlog = theta_stay[species_id_stay[i]], log = 1)
          meanlog[i] <- log(scale[i])
          if (nPreds_stay > 1) {
            if (nLevels_stay == 0) { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) }
            else { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]] }
          } else {
            if (nLevels_stay == 0) { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] }
            else { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] + random_effect_stay[group_stay[i]] }
          }
        }
        for (m in 1:nSpecies) { theta_stay[m] ~ dgamma(shape_stay, rate_stay); sdlog[m] <- theta_stay[m] }
        shape_stay ~ dgamma(2, 0.5); rate_stay ~ dgamma(2, 0.5)
        if (nPreds_stay == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- exp(beta_stay[1] + species_effect_stay[m, 1] + theta_stay[m]^2 / 2) } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- exp(inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay_station[i, 1:nPreds_stay]) + theta_stay[m]^2 / 2) } }
        }
        for (j in 1:nPreds_stay) {
          beta_stay[j] ~ dnorm(0, sd = 5)
          sigma_species_stay[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_stay[m, j] ~ dnorm(0, sd = sigma_species_stay[j]) }
        }
        if (nLevels_stay > 0) {
          for (k in 1:nLevels_stay) { random_effect_stay[k] ~ dnorm(0, sd = sigma_stay) }
          sigma_stay ~ T(dnorm(0, sd = 2), 0, )
        }
        for (m in 1:nSpecies) {
          for (i in 1:N_station) {
            Y_matrix[i, m] ~ dnbinom(size = size[m], prob = p[i, m])
            p[i, m] <- size[m] / (size[m] + mu[i, m])
            Y_rep[i, m]    ~ dnbinom(size = size[m], prob = p[i, m])
            loglike_obs_y[i, m]  <- dnbinom(Y_matrix[i, m], size[m], p[i, m], log = 1)
            loglike_pred_y[i, m] <- dnbinom(Y_rep[i, m],    size[m], p[i, m], log = 1)
          }
          size[m] ~ dgamma(1, 1)
        }
        if (nPreds_density == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- beta_density[1] + species_effect_density[m, 1]
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m])
          } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- inprod(beta_density[1:nPreds_density] + species_effect_density[m, 1:nPreds_density], X_density[i, 1:nPreds_density])
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m])
          } }
        }
        for (j in 1:nPreds_density) {
          beta_density[j] ~ dnorm(0, sd = 5)
          sd_species_density[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_density[m, j] ~ dnorm(0, sd = sd_species_density[j]) }
        }
      })

    } else { # weibull
      nimble::nimbleCode({
        for (i in 1:N_stay) {
          censored[i] ~ dinterval(stay[i], c_time[i])
          stay[i]   ~ dweibull(shape = theta_stay[species_id_stay[i]], scale = scale[i])
          pred_t[i] ~ dweibull(shape = theta_stay[species_id_stay[i]], scale = scale[i])
          loglike_obs_stay[i]  <- (1 - step(censored[i] - 0.5)) * dweibull(stay[i], shape = theta_stay[species_id_stay[i]], scale = scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pweibull(c_time[i], shape = theta_stay[species_id_stay[i]], scale = scale[i]))
          loglike_pred_stay[i] <- dweibull(pred_t[i], shape = theta_stay[species_id_stay[i]], scale = scale[i], log = 1)
          if (nPreds_stay > 1) {
            if (nLevels_stay == 0) { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) }
            else { log(scale[i]) <- inprod(beta_stay[1:nPreds_stay] + species_effect_stay[species_id_stay[i], 1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]] }
          } else {
            if (nLevels_stay == 0) { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] }
            else { log(scale[i]) <- beta_stay[1] + species_effect_stay[species_id_stay[i], 1] + random_effect_stay[group_stay[i]] }
          }
        }
        for (m in 1:nSpecies) { theta_stay[m] ~ dgamma(shape_stay, rate_stay); shape[m] <- theta_stay[m] }
        shape_stay ~ dgamma(2, 0.5); rate_stay ~ dgamma(2, 0.5)
        if (nPreds_stay == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- exp(lgamma(1 + 1 / theta_stay[m]) + beta_stay[1] + species_effect_stay[m, 1]) } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) { mean_stay[i, m] <- exp(lgamma(1 + 1 / theta_stay[m]) + inprod(beta_stay[1:nPreds_stay] + species_effect_stay[m, 1:nPreds_stay], X_stay_station[i, 1:nPreds_stay])) } }
        }
        for (j in 1:nPreds_stay) {
          beta_stay[j] ~ dnorm(0, sd = 5)
          sigma_species_stay[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_stay[m, j] ~ dnorm(0, sd = sigma_species_stay[j]) }
        }
        if (nLevels_stay > 0) {
          for (k in 1:nLevels_stay) { random_effect_stay[k] ~ dnorm(0, sd = sigma_stay) }
          sigma_stay ~ T(dnorm(0, sd = 2), 0, )
        }
        for (m in 1:nSpecies) {
          for (i in 1:N_station) {
            Y_matrix[i, m] ~ dnbinom(size = size[m], prob = p[i, m])
            p[i, m] <- size[m] / (size[m] + mu[i, m])
            Y_rep[i, m]    ~ dnbinom(size = size[m], prob = p[i, m])
            loglike_obs_y[i, m]  <- dnbinom(Y_matrix[i, m], size[m], p[i, m], log = 1)
            loglike_pred_y[i, m] <- dnbinom(Y_rep[i, m],    size[m], p[i, m], log = 1)
          }
          size[m] ~ dgamma(1, 1)
        }
        if (nPreds_density == 1) {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- beta_density[1] + species_effect_density[m, 1]
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m])
          } }
        } else {
          for (m in 1:nSpecies) { for (i in 1:N_station) {
            log(density[i, m]) <- inprod(beta_density[1:nPreds_density] + species_effect_density[m, 1:nPreds_density], X_density[i, 1:nPreds_density])
            log(mu[i, m]) <- log(density[i, m]) + log(S) + log(N_period[i]) - log(mean_stay[i, m]) + log(activity_proportion[m])
          } }
        }
        for (j in 1:nPreds_density) {
          beta_density[j] ~ dnorm(0, sd = 5)
          sd_species_density[j] ~ T(dnorm(0, sd = 2), 0, )
          for (m in 1:nSpecies) { species_effect_density[m, j] ~ dnorm(0, sd = sd_species_density[j]) }
        }
      })
    }
  } # end get_REST_code

  # Select code
  if (model == "RAD-REST") {
    code <- get_RADREST_code(stay_family)
  } else {
    code <- get_REST_code(stay_family)
  }

  # Initial values ------------------------------------------------------------

  stay_inits <- rep(NA, N_stay)
  stay_inits[is.censored == 1] <- c_time[is.censored == 1] + 1.0

  inits_f <- function() {
    common_inits <- list(
      beta_stay           = stats::runif(nPreds_stay, -1, 1),
      stay                = stay_inits,
      theta_stay          = stats::runif(nSpecies, 0.5, 4.0),
      shape_stay          = stats::runif(1, 0.5, 2.0),
      rate_stay           = stats::runif(1, 0.5, 2.0),
      species_effect_stay = matrix(stats::runif(nSpecies * nPreds_stay, -1, 1),
                                   nrow = nSpecies, ncol = nPreds_stay),
      sigma_species_stay  = stats::runif(nPreds_stay, 0.01, 1),
      beta_density           = stats::rnorm(nPreds_density, 0, 1),
      species_effect_density = matrix(stats::rnorm(nSpecies * nPreds_density, 0, 0.5),
                                      nrow = nSpecies, ncol = nPreds_density),
      sd_species_density     = stats::runif(nPreds_density, 0.01, 2),
      size = stats::rgamma(nSpecies, shape = 1, rate = 1)
    )

    if (model == "RAD-REST") {
      common_inits$theta_enter          <- stats::runif(nSpecies, 1, 5)
      common_inits$cutpoint             <- c(NA, sort(stats::rnorm(N_group - 1, 0, 0.5)))
      common_inits$beta_enter           <- stats::rnorm(nPreds_alpha, 0, 0.1)
      common_inits$species_effect_alpha <- matrix(stats::rnorm(nSpecies * nPreds_alpha, 0, 0.1),
                                                   nrow = nSpecies, ncol = nPreds_alpha)
      common_inits$sd_species_alpha     <- stats::runif(nPreds_alpha, 0.01, 1)
    }

    if (nLevels_stay > 0) {
      common_inits$random_effect_stay <- stats::runif(nLevels_stay, -1, 1)
      common_inits$sigma_stay         <- stats::runif(1, 0.5, 2.5)
    }

    common_inits
  }

  # MCMC parameters -----------------------------------------------------------

  if (stay_family == "exponential") {
    prms <- c("scale", "mean_stay")
  } else if (stay_family %in% c("gamma", "weibull")) {
    prms <- c("scale", "shape", "mean_stay")
  } else {
    prms <- c("meanlog", "sdlog", "mean_stay")
  }

  prms <- c(prms, "density",
            "beta_density", "species_effect_density",
            "beta_stay", "species_effect_stay")

  if (model == "RAD-REST") {
    prms <- c(prms, "mean_pass", "beta_enter", "species_effect_alpha")
    params <- c(prms, "loglike_obs_stay", "loglike_obs_y", "loglike_obs_detection",
                "loglike_pred_stay", "loglike_pred_y", "loglike_pred_detection")
  } else {
    params <- c(prms, "loglike_obs_stay", "loglike_obs_y",
                "loglike_pred_stay", "loglike_pred_y")
  }

  if (activity_estimation == "mixture") params <- c(params, "activity_proportion")

  # Parallel MCMC -------------------------------------------------------------

  if (activity_estimation != "mixture") {
    per_chain_info <- lapply(seq_len(nc), function(i) {
      list(seed = sample(1:9999, 1), inits = inits_f())
    })
  } else {
    per_chain_info <- lapply(seq_along(actv_out_trace), function(i) {
      list(seed = sample(1:9999, 1), inits = inits_f(),
           actv_samples = as.matrix(actv_out_trace[[i]]))
    })
  }

  is_mixture <- activity_estimation == "mixture"

  run_MCMC_multi <- function(info, data, constants, code, params, ni, nt, nb, is_mixture) {
    worker_dir <- file.path(tempdir(), paste0("nimble_worker_", Sys.getpid()))
    dir.create(worker_dir, showWarnings = FALSE)

    myModel  <- nimble::nimbleModel(code = code, data = data, constants = constants, inits = info$inits)
    CmyModel <- nimble::compileNimble(myModel, dirName = worker_dir)
    configModel <- nimble::configureMCMC(myModel, monitors = params)

    if (is_mixture) {
      configModel$removeSampler("activity_proportion")
      configModel$addSampler(target = "activity_proportion", type = 'prior_samples',
                             control = list(samples = info$actv_samples))
    }

    myMCMC  <- nimble::buildMCMC(configModel)
    CmyMCMC <- nimble::compileNimble(myMCMC, project = myModel, dirName = worker_dir)
    nimble::runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, nchains = 1,
                    setSeed = info$seed, samplesAsCodaMCMC = TRUE)
  }

  cat("Compiling the model. This may take a moment...\n")

  this_cluster <- parallel::makeCluster(nc)
  on.exit(try(parallel::stopCluster(this_cluster), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(this_cluster, {
    library(nimble)
    ddirchmulti <- nimble::nimbleFunction(
      run = function(x = double(1), alpha = double(1), size = double(0), log = integer(0)) {
        returnType(double(0))
        logProb <- lgamma(size + 1) - sum(lgamma(x + 1)) + lgamma(sum(alpha)) -
          sum(lgamma(alpha)) + sum(lgamma(alpha + x)) - lgamma(sum(alpha) + size)
        if (log) return(logProb) else return(exp(logProb))
      }
    )
    rdirchmulti <- nimble::nimbleFunction(
      run = function(n = integer(0), alpha = double(1), size = double(0)) {
        returnType(double(1))
        if (n != 1) print("rdirchmulti only allows n = 1; using n = 1.")
        p <- rdirch(1, alpha)
        return(rmulti(1, size = size, prob = p))
      }
    )
    suppressMessages(nimble::registerDistributions(list(
      ddirchmulti = list(
        BUGSdist = "ddirchmulti(alpha, size)",
        types = c('value = double(1)', 'alpha = double(1)', 'size = double(0)'),
        pqAvail = FALSE
      )
    )))
  })
  parallel::clusterExport(this_cluster, c("run_MCMC_multi"), envir = environment())

  cat("Running MCMC sampling. Please wait...\n")
  chain_output <- parallel::parLapply(
    cl = this_cluster,
    X  = per_chain_info,
    fun = run_MCMC_multi,
    data      = data_density,
    code      = code,
    constants = cons_density,
    params    = params,
    ni         = ni,
    nt         = nt,
    nb         = nb,
    is_mixture = is_mixture
  )

  parallel::stopCluster(this_cluster)
  cat("Estimation is finished!\n")

  # WAIC computation ----------------------------------------------------------

  loglfstay <- MCMCvis::MCMCchains(chain_output, params = "loglike_obs_stay")
  loglfy    <- MCMCvis::MCMCchains(chain_output, params = "loglike_obs_y")

  if (model == "RAD-REST") {
    loglfN   <- MCMCvis::MCMCchains(chain_output, params = "loglike_obs_detection")
    loglfall <- cbind(loglfstay, loglfy, loglfN)
  } else {
    loglfall <- cbind(loglfstay, loglfy)
  }

  lppd   <- sum(apply(loglfall, 2, safe_log_mean_exp), na.rm = TRUE)
  p.waic <- sum(apply(loglfall, 2, safe_var),          na.rm = TRUE)
  waic   <- (-2) * lppd + 2 * p.waic

  # tidy_samples --------------------------------------------------------------

  samples_mat <- MCMCvis::MCMCchains(chain_output)
  n_iters     <- nrow(samples_mat)
  p_names     <- colnames(samples_mat)
  tidy_samples <- data.frame(
    parameter = rep(p_names, each = n_iters),
    value     = as.vector(samples_mat),
    iteration = rep(seq_len(n_iters), times = length(p_names)),
    stringsAsFactors = FALSE
  )

  # mcmc.list
  mcmc_samples <- coda::as.mcmc.list(chain_output)

  # Summary -------------------------------------------------------------------

  check_no_cov <- function(f) {
    if (is.null(f)) return(TRUE)
    f <- stats::as.formula(f)
    vars <- all.vars(f[[length(f)]])
    length(vars) == 0
  }

  is_density_global <- check_no_cov(formula_density)
  is_stay_global    <- check_no_cov(formula_stay) && is.null(random_effect_stay)
  is_enter_global   <- if (model == "RAD-REST") check_no_cov(formula_enter) else TRUE
  is_pass_global    <- is_enter_global

  all_samples_mat <- MCMCvis::MCMCchains(mcmc_samples)

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
    for (m in seq_len(nSpecies)) {
      for (j in seq_len(cfg$nPreds)) {
        beta_col <- if (cfg$nPreds == 1) cfg$beta_pfx else paste0(cfg$beta_pfx, "[", j, "]")
        eff_col  <- paste0(cfg$eff_pfx, "[", m, ", ", j, "]")
        if (!(beta_col %in% colnames(all_samples_mat))) next
        if (!(eff_col  %in% colnames(all_samples_mat))) next
        samps <- all_samples_mat[, beta_col] + all_samples_mat[, eff_col]
        rows[[length(rows) + 1]] <- tibble::tibble(
          Species  = target_species[m],
          Station  = "All",
          Variable = paste0("coef_", param_type, "[", j, "] (", cfg$col_names[j], ")"),
          mean   = mean(samps),
          sd     = stats::sd(samps),
          lower  = stats::quantile(samps, 0.025),
          median = stats::quantile(samps, 0.500),
          upper  = stats::quantile(samps, 0.975),
          Rhat   = NA_real_,
          n.eff  = NA_real_
        )
      }
    }
    dplyr::bind_rows(rows)
  }

  # density summary
  raw_density <- summarize_param("density")
  if (nPreds_density == 1) {
    summary_density <- raw_density %>%
      tidyr::extract(Variable, into = "Species_idx",
                     regex = "\\[(\\d+)\\]", convert = TRUE, remove = FALSE) %>%
      dplyr::mutate(Species = target_species[Species_idx], Station = "All") %>%
      dplyr::select(-Species_idx)
  } else {
    summary_density <- raw_density %>%
      tidyr::extract(Variable, into = c("Station_idx", "Species_idx"),
                     regex = "\\[(\\d+),\\s*(\\d+)\\]", convert = TRUE, remove = FALSE) %>%
      dplyr::mutate(Species = target_species[Species_idx],
                    Station = unique_stations[Station_idx]) %>%
      dplyr::select(-Station_idx, -Species_idx)
  }

  # mean_stay summary
  raw_stay <- summarize_param("mean_stay") %>%
    tidyr::extract(Variable, into = c("Station_idx", "Species_idx"),
                   regex = "\\[(\\d+),\\s*(\\d+)\\]", convert = TRUE, remove = FALSE)
  if (is_stay_global) {
    summary_stay <- raw_stay %>%
      dplyr::filter(Station_idx == 1) %>%
      dplyr::mutate(Species  = target_species[Species_idx],
                    Station  = "All",
                    Variable = paste0("mean_stay[", Species_idx, "]")) %>%
      dplyr::select(-Station_idx, -Species_idx)
  } else {
    summary_stay <- raw_stay %>%
      dplyr::mutate(Species = target_species[Species_idx],
                    Station = unique_stations[Station_idx]) %>%
      dplyr::select(-Station_idx, -Species_idx)
  }

  summary_coef_list <- list(summary_density, summary_stay)

  # mean_pass summary (RAD-REST only)
  if (model == "RAD-REST") {
    raw_pass <- summarize_param("mean_pass") %>%
      tidyr::extract(Variable, into = c("Station_idx", "Species_idx"),
                     regex = "\\[(\\d+),\\s*(\\d+)\\]", convert = TRUE, remove = FALSE)
    if (is_pass_global) {
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
    summary_coef_list <- c(summary_coef_list, list(summary_pass))
  }

  if (nPreds_stay    > 1) summary_coef_list <- c(summary_coef_list, list(make_species_coef_summary("stay")))
  if (nPreds_density > 1) summary_coef_list <- c(summary_coef_list, list(make_species_coef_summary("density")))
  if (model == "RAD-REST" && nPreds_alpha > 1)
    summary_coef_list <- c(summary_coef_list, list(make_species_coef_summary("enter")))

  summary_mean <- dplyr::bind_rows(summary_coef_list) %>%
    dplyr::mutate(cv = abs(sd / mean)) %>%
    dplyr::select(Species, Station, Variable, mean, sd, lower, median, upper, Rhat, n.eff, cv)

  # Scaling params
  names(scaling_stay$center)    <- colnames(X_stay)
  names(scaling_stay$scale)     <- colnames(X_stay)
  names(scaling_density$center) <- colnames(X_density)
  names(scaling_density$scale)  <- colnames(X_density)

  scaling_params <- list(stay = scaling_stay, density = scaling_density)

  if (model == "RAD-REST") {
    names(scaling_alpha$center) <- colnames(X_alpha)
    names(scaling_alpha$scale)  <- colnames(X_alpha)
    scaling_params$enter <- scaling_alpha
  }

  # Return --------------------------------------------------------------------

  density_result <- list(
    WAIC           = waic,
    summary_result = summary_mean,
    samples        = mcmc_samples,
    tidy_samples   = tidy_samples,
    scaling_params = scaling_params
  )
  class(density_result) <- "ResultDensity"

  return(density_result)
}
