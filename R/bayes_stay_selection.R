#' Model selection for staying time analysis using WAIC for REST/RAD-REST models based on MCMC samplings with nimble
#'
#' @param formula_stay A model formula for staying time within a focal area. For example, `Stay ~ 1 + x1`.
#' @param random_effect_stay A character string specifying the column name in `stay_data` to use as a grouping factor for a station-level random intercept on staying time (e.g., `"Station"`). Default is `NULL` (no random effects).
#' @param stay_data A data frame returned by the `format_stay` function, containing the processed staying time data.
#' @param col_name_cens A string specifying the column name in `stay_data` that indicates whether the observation is censored (1) or not (0).
#' @param stay_family A character string specifying the probability distribution of staying time. Choose from `"exponential"`, `"gamma"`, `"weibull"`, or `"lognormal"`.
#' @param cores The number of CPU cores to use for parallel computation. Default is 3.
#' @param iter The total number of MCMC iterations per chain. Default is 5000.
#' @param warmup The number of warm-up (burn-in) iterations per chain. Default is 1000.
#' @param thin The thinning interval for MCMC sampling. Default is 4.
#' @param chains The number of MCMC chains. Default is 3.
#' @param all_comb A logical value indicating whether to compare models with all possible combinations of covariates. If `FALSE`, only the specified model in `formula_stay` is evaluated. Default is FALSE.
#' @param target_species A character string specifying the species of interest. Only a single species can be specified.
#' @return A list of class \code{"ResultStay"}, which includes the following components:
#' \describe{
#'   \item{\code{WAIC}}{A data frame of WAIC values for each candidate stay model, sorted from best to worst.}
#'   \item{\code{Bayesian_p_value}}{Bayesian p-value for the best model, used to assess model fit.}
#'   \item{\code{summary_result}}{A data frame summarizing posterior estimates of the mean staying time.}
#'   \item{\code{samples}}{A \code{coda::mcmc.list} object containing MCMC samples for all parameters.}
#'   \item{\code{tidy_samples}}{A long-format data frame of all monitored MCMC samples from the best model, with columns \code{parameter}, \code{value}, and \code{iteration}.}
#' }
#' @export
#' @import dplyr nimble parallel MCMCvis
#' @importFrom tidyr unite extract
#' @importFrom purrr map
#' @importFrom stringr str_extract str_detect
#' @importFrom stats as.formula formula model.frame model.matrix sd var runif median quantile model.response rexp rnorm step dexp pexp dgamma pgamma dlnorm plnorm dweibull pweibull dnbinom
#' @examples
#' stay_data <- format_stay(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_stay = "Stay",
#'   col_name_cens = "Cens"
#' )
#' bayes_stay_selection(
#'   formula_stay   = Stay ~ 1 + x1,
#'   stay_data      = stay_data,
#'   col_name_cens  = "Cens",
#'   stay_family    = "lognormal",
#'   cores          = 2,
#'   iter           = 5000,
#'   warmup         = 1000,
#'   chains         = 2,
#'   thin           = 4,
#'   all_comb       = FALSE,
#'   target_species = "SP01"
#' )

bayes_stay_selection <- function(
    formula_stay      = Stay ~ 1,
    random_effect_stay = NULL,
    stay_data         = stay_data,
    col_name_cens     = "Cens",
    stay_family       = "lognormal",
    cores             = 3,
    iter              = 5000,
    warmup            = 1000,
    chains            = 3,
    thin              = 4,
    all_comb          = FALSE,
    target_species    = NULL
) {

  # Input checks --------------------------------------------------------------

  if (!stay_family %in% c("lognormal", "gamma", "weibull", "exponential")) {
    stop(paste0("Input stay_family (", stay_family, ") is incorrect."))
  }
  if (!inherits(formula_stay, "formula")) {
    stop("formula_stay must be a formula.")
  }
  if (!is.null(random_effect_stay) && !(random_effect_stay %in% colnames(stay_data))) {
    stop(paste0("random_effect_stay column '", random_effect_stay, "' not found in stay_data."))
  }
  if (!is.data.frame(stay_data)) {
    stop("stay_data must be a data frame.")
  }
  if (!(col_name_cens %in% colnames(stay_data))) {
    stop(paste0("Column '", col_name_cens, "' not found in stay_data."))
  }
  if (!is.numeric(cores)   || cores   < 1 || (cores   %% 1 != 0)) stop("cores must be a positive integer.")
  if (!is.numeric(iter)    || iter    < 1 || (iter    %% 1 != 0)) stop("iter must be a positive integer.")
  if (!is.null(warmup) && (!is.numeric(warmup) || warmup < 1 || (warmup %% 1 != 0))) stop("warmup must be a positive integer or NULL.")
  if (!is.numeric(chains)  || chains  < 1 || (chains  %% 1 != 0)) stop("chains must be a positive integer.")
  if (!is.numeric(thin)    || thin    < 1 || (thin    %% 1 != 0)) stop("thin must be a positive integer.")
  if (!is.null(target_species) && !all(target_species %in% stay_data$Species)) {
    stop("Some values in target_species are not found in stay_data$Species.")
  }

  # Helper: all covariate combinations ----------------------------------------

  bit.test <- function(number, n) {
    (number %/% (2^n)) %% 2
  }

  full_terms <- function(x) {
    lapply(seq_len(2^length(x)), function(i) {
      r <- ""
      for (j in seq_len(length(x))) {
        r <- paste(r, ifelse(bit.test((i - 1), (j - 1)), paste0(x[j], " + "), ""), sep = "")
      }
      r <- paste0("Stay ~ 1 + ", r, "1")
      strsplit(r, " \\+ 1$")[[1]][1]
    })
  }

  # Prepare data ---------------------------------------------------------------

  stay_data <- stay_data %>%
    dplyr::filter(Species == target_species) %>%
    dplyr::arrange(Station)

  station.id <- stay_data %>% dplyr::pull(Station) %>% unique()

  vars_stay       <- all.vars(formula_stay)
  predictors_stay <- vars_stay[-1]

  formula_stay_all <- list()
  if (all_comb) {
    formula_stay_all <- full_terms(c(predictors_stay))
  } else {
    formula_stay_all[[1]] <- formula_stay
  }

  mcmc_samples <- tidy_samples <- list()
  waic <- p_value <- numeric(0)

  # MCMC loop over candidate models -------------------------------------------

  for (k in seq_along(formula_stay_all)) {

    formula_stay_k <- formula_stay_all[[k]]

    model_frame_stay <- model.frame(formula_stay_k, stay_data)
    X_stay           <- model.matrix(stats::as.formula(formula_stay_k), model_frame_stay)
    stay             <- model.response(model_frame_stay)
    censored         <- stay_data[[col_name_cens]]

    nPreds_stay <- ncol(X_stay)
    N_station   <- length(unique(stay_data$Station))

    names(stay) <- NULL
    c_time <- stay
    c_time[censored == 0] <- c_time[censored == 0] + 0.1
    stay[censored == 1]   <- NA
    N_stay <- length(stay)

    if (!is.null(random_effect_stay)) {
      levels_stay  <- unique(stay_data[[random_effect_stay]])
      nLevels_stay <- length(levels_stay)
    } else {
      nLevels_stay <- 0
    }

    data_stay <- list(stay = stay, censored = censored)
    cons_stay <- list(
      N_stay      = N_stay,
      nPreds_stay = nPreds_stay,
      N_station   = N_station,
      c_time      = c_time
    )

    if (nPreds_stay > 1) {
      data_stay$X_stay <- X_stay

      # Station-level design matrix for mean_stay (one row per station)
      stay_data_station <- stay_data %>%
        dplyr::arrange(Station) %>%
        dplyr::distinct(Station, .keep_all = TRUE)
      model_frame_station <- model.frame(formula_stay_k, stay_data_station)
      X_stay_station      <- model.matrix(stats::as.formula(formula_stay_k), model_frame_station)
      cons_stay$X_stay_station <- X_stay_station
    }

    if (!is.null(random_effect_stay)) {
      cons_stay$group_stay  <- as.numeric(factor(stay_data[[random_effect_stay]]))
      cons_stay$nLevels_stay <- nLevels_stay
    } else {
      cons_stay$nLevels_stay <- 0
    }

    # NIMBLE model code --------------------------------------------------------

    code <- nimbleCode({
      for (i in 1:N_stay) {
        censored[i] ~ dinterval(stay[i], c_time[i])

        if (stay_family == "exponential") {
          stay[i]   ~ dexp(rate = 1/scale[i])
          pred_t[i] ~ dexp(rate = 1/scale[i])
          loglike_obs[i] <- (1 - step(censored[i] - 0.5)) * dexp(stay[i], rate = 1/scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pexp(c_time[i], rate = 1/scale[i]))
          loglike_pred[i] <- (1 - step(censored[i] - 0.5)) *
            dexp(pred_t[i], rate = 1/scale[i], log = 1) +
            step(censored[i] - 0.5) *
            log(1 - pexp(c_time[i], rate = 1/scale[i]))
        }

        if (stay_family == "gamma") {
          stay[i]   ~ dgamma(shape = theta_stay, rate = exp(-log(scale[i])))
          pred_t[i] ~ dgamma(shape = theta_stay, rate = exp(-log(scale[i])))
          loglike_obs[i] <- (1 - step(censored[i] - 0.5)) * dgamma(stay[i], shape = theta_stay, rate = exp(-log(scale[i])), log = 1) +
            step(censored[i] - 0.5) * log(1 - pgamma(c_time[i], shape = theta_stay, rate = exp(-log(scale[i]))))
          loglike_pred[i] <- (1 - step(censored[i] - 0.5)) *
            dgamma(pred_t[i], shape = theta_stay, rate = exp(-log(scale[i])), log = 1) +
            step(censored[i] - 0.5) *
            log(1 - pgamma(c_time[i], shape = theta_stay, rate = exp(-log(scale[i]))))
        }

        if (stay_family == "lognormal") {
          stay[i]   ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay)
          pred_t[i] ~ dlnorm(meanlog = log(scale[i]), sdlog = theta_stay)
          loglike_obs[i] <- (1 - step(censored[i] - 0.5)) * dlnorm(stay[i], meanlog = log(scale[i]), sdlog = theta_stay, log = 1) +
            step(censored[i] - 0.5) * log(1 - plnorm(c_time[i], meanlog = log(scale[i]), sdlog = theta_stay))
          loglike_pred[i] <- (1 - step(censored[i] - 0.5)) *
            dlnorm(pred_t[i], meanlog = log(scale[i]), sdlog = theta_stay, log = 1) +
            step(censored[i] - 0.5) *
            log(1 - plnorm(c_time[i], meanlog = log(scale[i]), sdlog = theta_stay))
          meanlog[i] <- log(scale[i])
        }

        if (stay_family == "weibull") {
          stay[i]   ~ dweibull(shape = theta_stay, scale = scale[i])
          pred_t[i] ~ dweibull(shape = theta_stay, scale = scale[i])
          loglike_obs[i] <- (1 - step(censored[i] - 0.5)) * dweibull(stay[i], shape = theta_stay, scale = scale[i], log = 1) +
            step(censored[i] - 0.5) * log(1 - pweibull(c_time[i], shape = theta_stay, scale = scale[i]))
          loglike_pred[i] <- (1 - step(censored[i] - 0.5)) *
            dweibull(pred_t[i], shape = theta_stay, scale = scale[i], log = 1) +
            step(censored[i] - 0.5) *
            log(1 - pweibull(c_time[i], shape = theta_stay, scale = scale[i]))
        }

        # Linear predictor for log(scale[i])
        if (nPreds_stay == 1) {
          if (nLevels_stay > 0) {
            log(scale[i]) <- beta_stay[1] + random_effect_stay[group_stay[i]]
          } else {
            log(scale[i]) <- beta_stay[1]
          }
        } else {
          if (nLevels_stay > 0) {
            log(scale[i]) <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay]) + random_effect_stay[group_stay[i]]
          } else {
            log(scale[i]) <- inprod(beta_stay[1:nPreds_stay], X_stay[i, 1:nPreds_stay])
          }
        }
      }

      if (stay_family == "lognormal") {
        sdlog <- theta_stay
      } else {
        shape <- theta_stay
      }

      # Priors
      for (j in 1:nPreds_stay) {
        beta_stay[j] ~ dnorm(0, sd = 5)
      }
      theta_stay ~ dgamma(1, 1)

      if (nLevels_stay > 0) {
        for (k in 1:nLevels_stay) {
          random_effect_stay[k] ~ dnorm(0, sd = sigma_stay)
        }
        sigma_stay ~ T(dnorm(0, sd = 5), 0, 10)
      }

      # Expected staying time
      if (nPreds_stay == 1) {
        if (stay_family == "exponential") mean_stay <- exp(beta_stay[1])
        if (stay_family == "gamma")       mean_stay <- theta_stay * exp(beta_stay[1])
        if (stay_family == "lognormal")   mean_stay <- exp(beta_stay[1] + theta_stay ^ 2 / 2)
        if (stay_family == "weibull")     mean_stay <- exp(beta_stay[1] + lgamma(1 + 1 / theta_stay))
      }
      if (nPreds_stay > 1) {
        for (i in 1:N_station) {
          if (stay_family == "exponential") mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay_station[i, 1:nPreds_stay]))
          if (stay_family == "gamma")       mean_stay[i] <- theta_stay * exp(inprod(beta_stay[1:nPreds_stay], X_stay_station[i, 1:nPreds_stay]))
          if (stay_family == "lognormal")   mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay_station[i, 1:nPreds_stay]) + theta_stay ^ 2 / 2)
          if (stay_family == "weibull")     mean_stay[i] <- exp(inprod(beta_stay[1:nPreds_stay], X_stay_station[i, 1:nPreds_stay]) + lgamma(1 + 1 / theta_stay))
        }
      }

      sum_loglike_obs  <- sum(loglike_obs[1:N_stay])
      sum_loglike_pred <- sum(loglike_pred[1:N_stay])
      deviance_obs     <- -2 * sum_loglike_obs
      deviance_pred    <- -2 * sum_loglike_pred
    })

    cat("Compiling the model. This may take a moment...\n")

    this_cluster <- parallel::makeCluster(chains)
    on.exit(try(parallel::stopCluster(this_cluster), silent = TRUE), add = TRUE)

    run_MCMC_stay <- function(info, data, constants, code, params, iter, thin, warmup) {
      myModel  <- nimble::nimbleModel(code = code, data = data, constants = constants, inits = info$inits)
      CmyModel <- nimble::compileNimble(myModel)
      myMCMC   <- nimble::buildMCMC(CmyModel, monitors = params)
      CmyMCMC  <- nimble::compileNimble(myMCMC)
      nimble::runMCMC(
        CmyMCMC,
        niter              = iter,
        nburnin            = warmup,
        thin               = thin,
        nchains            = 1,
        setSeed            = info$seed,
        samplesAsCodaMCMC  = TRUE
      )
    }

    inits_f <- function() {
      inits <- list(
        beta_stay  = stats::runif(nPreds_stay, -1, 1),
        stay       = ifelse(censored == 0, NA, c_time + stats::runif(N_stay, 0.1, 1.0)),
        theta_stay = stats::runif(1, 0.5, 1.5),
        scale      = stats::runif(N_stay, 0, 2)
      )
      if (!is.null(random_effect_stay)) {
        inits$random_effect_stay <- stats::runif(nLevels_stay, -1, 1)
        inits$sigma_stay         <- stats::runif(1, 0.5, 2.5)
      }
      inits
    }

    per_chain_info <- lapply(seq_len(chains), function(i) {
      list(seed = sample(1:9999, 1), inits = inits_f())
    })

    if (stay_family == "exponential")              prms <- c("scale", "mean_stay")
    if (stay_family %in% c("gamma", "weibull"))    prms <- c("scale", "shape", "mean_stay")
    if (stay_family == "lognormal")                prms <- c("meanlog", "sdlog", "mean_stay")

    params <- c(prms, "beta_stay", "sum_loglike_obs", "sum_loglike_pred",
                "deviance_obs", "deviance_pred", "loglike_obs")

    cat("Running MCMC sampling. Please wait...\n")

    parallel::clusterEvalQ(this_cluster, { library(nimble) })
    parallel::clusterExport(this_cluster, varlist = c("run_MCMC_stay"), envir = environment())

    chain_output <- parallel::parLapply(
      cl        = this_cluster,
      X         = per_chain_info,
      fun       = run_MCMC_stay,
      data      = data_stay,
      code      = code,
      constants = cons_stay,
      params    = params,
      iter      = iter,
      thin      = thin,
      warmup    = warmup
    )
    parallel::stopCluster(this_cluster)
    cat("Estimation is finished!\n")

    # Summarize results --------------------------------------------------------

    mcmc_samples[[k]] <- MCMCvis::MCMCchains(chain_output, mcmc.list = TRUE, params = prms)

    samples_mat <- MCMCvis::MCMCchains(chain_output)
    n_iters     <- nrow(samples_mat)
    p_names     <- colnames(samples_mat)
    tidy_samples[[k]] <- data.frame(
      parameter = rep(p_names, each = n_iters),
      value     = as.vector(samples_mat),
      iteration = rep(seq_len(n_iters), times = length(p_names)),
      stringsAsFactors = FALSE
    )

    # Bayesian p-value
    deviance_obs_samples  <- samples_mat[, "deviance_obs"]
    deviance_pred_samples <- samples_mat[, "deviance_pred"]
    p_value[k] <- mean(deviance_pred_samples > deviance_obs_samples)

    # WAIC (log-sum-exp for numerical stability)
    loglf <- MCMCvis::MCMCchains(chain_output, params = "loglike_obs")
    safe_log_mean_exp <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) return(NA)
      m <- max(x)
      m + log(mean(exp(x - m)))
    }
    safe_var <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 2) return(NA)
      stats::var(x)
    }
    lppd   <- sum(apply(loglf, 2, safe_log_mean_exp), na.rm = TRUE)
    p.waic <- sum(apply(loglf, 2, safe_var), na.rm = TRUE)
    waic[k] <- (-2) * lppd + 2 * p.waic
  }

  # WAIC comparison ------------------------------------------------------------

  re_label <- if (is.null(random_effect_stay)) "NULL" else random_effect_stay
  WAIC <- data.frame(
    Model         = as.character(unlist(formula_stay_all)),
    Family        = rep(stay_family, length(waic)),
    Random_effect = re_label,
    WAIC          = waic,
    stringsAsFactors = FALSE
  ) %>% dplyr::arrange(WAIC)

  best.model <- which(WAIC[1, 1] == unlist(formula_stay_all))

  mcmc_samples_best <- mcmc_samples[[best.model]]
  tidy_samples_best <- tidy_samples[[best.model]]
  p_value_best      <- p_value[best.model]

  # Summary of mean_stay for best model ----------------------------------------

  # Recompute nPreds for best model
  formula_stay_best   <- formula_stay_all[[best.model]]
  model_frame_best    <- model.frame(formula_stay_best, stay_data)
  X_stay_best         <- model.matrix(stats::as.formula(formula_stay_best), model_frame_best)
  nPreds_stay_best    <- ncol(X_stay_best)

  summary_mean_temp <- MCMCvis::MCMCsummary(mcmc_samples_best, params = "mean_stay", round = 2) %>%
    tibble::rownames_to_column(var = "Variable") %>%
    tibble::as_tibble() %>%
    dplyr::rename(lower = `2.5%`, median = `50%`, upper = `97.5%`) %>%
    dplyr::mutate(cv = sd / mean)

  if (nPreds_stay_best > 1) {
    summary_mean <- dplyr::bind_cols(
      tibble::tibble(Station = station.id),
      summary_mean_temp
    )
  } else {
    summary_mean <- summary_mean_temp %>%
      dplyr::mutate(Station = "All") %>%
      dplyr::select(Station, dplyr::everything())
  }

  stay_result <- list(
    WAIC             = WAIC,
    Bayesian_p_value = p_value_best,
    summary_result   = summary_mean,
    samples          = mcmc_samples_best,
    tidy_samples     = tidy_samples_best,
    target_species   = target_species,
    stay_family      = stay_family
  )
  class(stay_result) <- "ResultStay"

  stay_result
}
