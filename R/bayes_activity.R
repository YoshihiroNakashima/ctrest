#' Activity level estimation using a von Mises mixture distribution with a stick-breaking prior
#'
#' Estimates the proportion of time an animal is active via Bayesian MCMC using a
#' nonparametric von Mises mixture model (Nakashima et al. 2025). This is a standalone
#' version of the activity estimation that is also embedded within \code{bayes_rest}
#' when \code{activity_estimation = "mixture"}.
#'
#' @param activity_data A data frame containing a \code{time} column (detection times
#'   in radians, 0–2π) and a \code{Species} column. Typically the output of
#'   \code{format_activity}.
#' @param C The maximum number of von Mises components. Default is 10. Must be >= 2.
#' @param cores An integer specifying the number of CPU cores (currently unused; the
#'   number of parallel workers equals \code{chains}). Default is 3.
#' @param iter An integer specifying the total number of MCMC iterations per chain.
#'   Default is 5000.
#' @param warmup An integer specifying the number of warm-up (burn-in) iterations per
#'   chain. Default is 1000.
#' @param chains An integer specifying the number of MCMC chains. Default is 3.
#' @param thin An integer specifying the thinning interval. Default is 2.
#' @param target_species A single character string specifying the species to analyse.
#'
#' @return A list of class \code{"ResultActivity"} with the following components:
#' \describe{
#'   \item{\code{summary_result}}{A data frame with the posterior summary of
#'     \code{activity_proportion} (mean, sd, lower, median, upper, Rhat, n.eff, cv).}
#'   \item{\code{activity_curve}}{A data frame with the posterior summary of the
#'     activity density function at each time point (columns: variable, mean, sd,
#'     lower, median, upper, Rhat, n.eff, x), aligned with the output of
#'     \code{bayes_rest} when \code{activity_estimation = "mixture"}.}
#'   \item{\code{WAIC}}{WAIC for the fitted mixture model.}
#'   \item{\code{Bayesian_p_value}}{Bayesian p-value (proportion of posterior draws
#'     where replicated deviance exceeds observed deviance).}
#'   \item{\code{samples}}{A \code{coda::mcmc.list} of the main monitored parameters
#'     (\code{activity_proportion}, \code{activity_density}).}
#'   \item{\code{tidy_samples}}{A long-format data frame of all monitored MCMC
#'     samples, with columns \code{parameter}, \code{value}, and \code{iteration}.}
#'   \item{\code{act_data}}{Numeric vector of raw detection times (radians) for the
#'     target species. Used by \code{printResultActivity} to generate plots.}
#'   \item{\code{target_species}}{Character string of the analysed species name.}
#' }
#' Use \code{printResultActivity} to display a summary and activity curve plot.
#'
#' @export
#' @import dplyr nimble parallel MCMCvis
#' @importFrom stats rbeta rgamma runif median quantile var
#' @examples
#' \dontrun{
#' activity_data <- format_activity(
#'   detection_data    = detection_data,
#'   col_name_station  = "Station",
#'   col_name_species  = "Species",
#'   col_name_datetime = "DateTime",
#'   indep_time        = 30
#' )
#' result <- bayes_activity(
#'   activity_data  = activity_data,
#'   C              = 10,
#'   iter           = 5000,
#'   warmup         = 1000,
#'   chains         = 3,
#'   thin           = 2,
#'   target_species = "SP01"
#' )
#' printResultActivity(result)
#' }
bayes_activity <- function(
    activity_data,
    C              = 10,
    cores          = 3,
    iter           = 5000,
    warmup         = 1000,
    chains         = 3,
    thin           = 2,
    target_species
) {

  # --- Input validation -------------------------------------------------------

  if (!is.data.frame(activity_data))
    stop("'activity_data' must be a data frame. Did you run format_activity() first?", call. = FALSE)

  req_cols  <- c("time", "Species")
  miss_cols <- setdiff(req_cols, colnames(activity_data))
  if (length(miss_cols) > 0)
    stop(
      sprintf(
        "Column(s) not found in 'activity_data': %s\n  Did you run format_activity() first?\n  Available columns: %s",
        paste(miss_cols, collapse = ", "),
        paste(colnames(activity_data), collapse = ", ")
      ),
      call. = FALSE
    )

  if (missing(target_species) || is.null(target_species))
    stop("'target_species' must be specified.", call. = FALSE)
  if (!is.character(target_species) || length(target_species) != 1)
    stop("'target_species' must be a single character string.", call. = FALSE)
  if (!target_species %in% activity_data$Species)
    stop(
      sprintf(
        "Species '%s' not found in 'activity_data$Species'.\n  Available species: %s",
        target_species,
        paste(sort(unique(activity_data$Species)), collapse = ", ")
      ),
      call. = FALSE
    )

  if (!is.numeric(C) || length(C) != 1 || C < 2 || C != floor(C))
    stop("'C' must be an integer >= 2.", call. = FALSE)
  if (!is.numeric(iter)   || length(iter)   != 1 || iter   <= 0 || iter   != floor(iter))
    stop("'iter' must be a positive integer.", call. = FALSE)
  if (!is.numeric(warmup) || length(warmup) != 1 || warmup <= 0 || warmup != floor(warmup))
    stop("'warmup' must be a positive integer.", call. = FALSE)
  if (warmup >= iter)
    stop(sprintf("'warmup' (%d) must be less than 'iter' (%d).", warmup, iter), call. = FALSE)
  if (!is.numeric(chains) || length(chains) != 1 || chains <= 0 || chains != floor(chains))
    stop("'chains' must be a positive integer.", call. = FALSE)
  if (!is.numeric(thin)   || length(thin)   != 1 || thin   <= 0 || thin   != floor(thin))
    stop("'thin' must be a positive integer.", call. = FALSE)

  # --- Prepare data -----------------------------------------------------------

  act_data_vec <- activity_data %>%
    dplyr::filter(Species == target_species) %>%
    dplyr::pull(time)

  N <- length(act_data_vec)
  if (N == 0)
    stop(
      sprintf("No detection records found for species '%s' in 'activity_data'.", target_species),
      call. = FALSE
    )
  if (N < C)
    warning(
      sprintf(
        "Number of detections (N = %d) is less than C = %d. The model is overparameterised; consider reducing C.",
        N, C
      ),
      call. = FALSE
    )

  dens.x <- seq(0, 2 * pi, 0.02)
  ndens  <- length(dens.x)

  constants <- list(N = N, C = C, dens.x = dens.x, ndens = ndens)
  data      <- list(act_data = act_data_vec)

  # --- NIMBLE model -----------------------------------------------------------
  # Node names aligned with bayes_rest (activity_estimation = "mixture"):
  #   act_data, activity_density, activity_proportion, loglike_obs_act

  code <- nimbleCode({
    for (k in 1:(C-1)) {
      v[k] ~ dbeta(1, alpha)
    }
    alpha ~ dgamma(1, 1)
    w[1:C] <- stick_breaking(v[1:(C-1)])
    for (k in 1:C) {
      mu_mix[k]    ~ dunif(0, 2 * 3.141592654)
      kappa_mix[k] ~ dgamma(1, 0.01)
    }
    for (n in 1:N) {
      group[n]          ~ dcat(w[1:C])
      act_data[n]       ~ dvonMises(mu_mix[group[n]], kappa_mix[group[n]])
      act_data_pred[n]  ~ dvonMises(mu_mix[group[n]], kappa_mix[group[n]])
      loglike_obs_act[n]  <- dvonMises(act_data[n],      mu_mix[group[n]], kappa_mix[group[n]], log = 1)
      loglike_pred_act[n] <- dvonMises(act_data_pred[n], mu_mix[group[n]], kappa_mix[group[n]], log = 1)
    }
    for (j in 1:ndens) {
      for (i in 1:C) {
        dens.cpt[i, j] <- w[i] * dvonMises(dens.x[j], mu_mix[i], kappa_mix[i], log = 0)
      }
      activity_density[j] <- sum(dens.cpt[1:C, j])
    }
    activity_proportion <- 1.0 / (2 * 3.141592654 * max(activity_density[1:ndens]))
    sum_loglike_obs  <- sum(loglike_obs_act[1:N])
    sum_loglike_pred <- sum(loglike_pred_act[1:N])
    deviance_obs     <- -2 * sum_loglike_obs
    deviance_pred    <- -2 * sum_loglike_pred
  })

  # --- Inits ------------------------------------------------------------------

  inits_f <- function() {
    list(
      mu_mix    = stats::runif(C, 0, 2 * pi),
      kappa_mix = stats::rgamma(C, 1, 0.01),
      group     = sample(seq_len(C), size = N, replace = TRUE),
      v         = stats::rbeta(C - 1, 1, 1),
      alpha     = 1
    )
  }

  per_chain_info <- lapply(seq_len(chains), function(i) {
    list(seed = sample(1:9999, 1), inits = inits_f())
  })

  params <- c(
    "activity_density", "activity_proportion",
    "mu_mix", "kappa_mix", "w",
    "loglike_obs_act", "loglike_pred_act",
    "deviance_obs", "deviance_pred"
  )

  # --- Parallel MCMC ----------------------------------------------------------

  cat("Compiling the model. This may take a moment...\n")

  this_cluster <- parallel::makeCluster(chains)
  on.exit(try(parallel::stopCluster(this_cluster), silent = TRUE), add = TRUE)

  # Define dvonMises/rvonMises in each worker's global environment
  parallel::clusterEvalQ(this_cluster, {
    library(nimble)
    dvonMises <- nimble::nimbleFunction(
      run = function(x = double(0), kappa = double(0), mu = double(0), log = integer(0)) {
        returnType(double(0))
        ccrit <- 1E-6; s <- 1; i <- 1; inc <- 1; x_2i <- 0; satisfied <- FALSE
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
    suppressMessages(nimble::registerDistributions(list(
      dvonMises = list(
        BUGSdist = "dvonMises(kappa, mu)",
        types    = c("value = double(0)", "kappa = double(0)", "mu = double(0)"),
        pqAvail  = FALSE
      )
    )))
  })

  run_MCMC_vonMises <- function(info, data, constants, code, params, iter, thin, warmup) {
    myModel     <- nimble::nimbleModel(code = code, data = data, constants = constants, inits = info$inits)
    CmyModel    <- nimble::compileNimble(myModel)
    configModel <- nimble::configureMCMC(myModel, monitors = params)
    myMCMC      <- nimble::buildMCMC(configModel)
    CmyMCMC     <- nimble::compileNimble(myMCMC)
    nimble::runMCMC(
      CmyMCMC,
      niter             = iter,
      nburnin           = warmup,
      thin              = thin,
      nchains           = 1,
      setSeed           = info$seed,
      samplesAsCodaMCMC = TRUE
    )
  }

  parallel::clusterExport(this_cluster, varlist = c("run_MCMC_vonMises"), envir = environment())

  cat("Running MCMC sampling. Please wait...\n")

  actv_chain_output <- parallel::parLapply(
    cl        = this_cluster,
    X         = per_chain_info,
    fun       = run_MCMC_vonMises,
    data      = data,
    code      = code,
    constants = constants,
    params    = params,
    iter      = iter,
    thin      = thin,
    warmup    = warmup
  )
  parallel::stopCluster(this_cluster)
  cat("Estimation is finished!\n")

  # --- Summarize results ------------------------------------------------------

  # activity_proportion posterior summary
  summary_result <- MCMCvis::MCMCsummary(
    MCMCvis::MCMCchains(actv_chain_output, mcmc.list = TRUE, params = "activity_proportion"),
    round = 4
  ) %>%
    tibble::rownames_to_column(var = "Variable") %>%
    tibble::as_tibble() %>%
    dplyr::rename(lower = `2.5%`, median = `50%`, upper = `97.5%`) %>%
    dplyr::mutate(cv = sd / mean)

  # Activity density curve (aligned with bayes_rest$activity_curve)
  activity_curve <- MCMCvis::MCMCsummary(
    MCMCvis::MCMCchains(actv_chain_output, mcmc.list = TRUE, params = "activity_density"),
    round = 5
  ) %>%
    tibble::rownames_to_column(var = "variable") %>%
    tibble::as_tibble() %>%
    dplyr::rename(lower = `2.5%`, median = `50%`, upper = `97.5%`) %>%
    dplyr::mutate(x = dens.x)

  # coda mcmc.list for key parameters
  samples <- MCMCvis::MCMCchains(
    actv_chain_output, mcmc.list = TRUE,
    params = c("activity_proportion", "activity_density")
  )

  # Long-format tidy samples (all monitored parameters)
  samples_mat <- MCMCvis::MCMCchains(actv_chain_output)
  n_iters     <- nrow(samples_mat)
  p_names     <- colnames(samples_mat)
  tidy_samples <- data.frame(
    parameter = rep(p_names, each = n_iters),
    value     = as.vector(samples_mat),
    iteration = rep(seq_len(n_iters), times = length(p_names)),
    stringsAsFactors = FALSE
  )

  # Bayesian p-value
  deviance_obs_samples  <- samples_mat[, "deviance_obs"]
  deviance_pred_samples <- samples_mat[, "deviance_pred"]
  bayesian_p_value      <- mean(deviance_pred_samples > deviance_obs_samples)

  # WAIC (log-sum-exp for numerical stability)
  loglf <- MCMCvis::MCMCchains(actv_chain_output, params = "loglike_obs_act")
  safe_log_mean_exp <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    m <- max(x)
    m + log(mean(exp(x - m)))
  }
  safe_var <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2) return(NA_real_)
    stats::var(x)
  }
  lppd   <- sum(apply(loglf, 2, safe_log_mean_exp), na.rm = TRUE)
  p.waic <- sum(apply(loglf, 2, safe_var),          na.rm = TRUE)
  waic   <- (-2) * lppd + 2 * p.waic

  # --- Return -----------------------------------------------------------------

  result <- list(
    summary_result   = summary_result,
    activity_curve   = activity_curve,
    WAIC             = waic,
    Bayesian_p_value = bayesian_p_value,
    samples          = samples,
    tidy_samples     = tidy_samples,
    act_data         = act_data_vec,
    target_species   = target_species
  )
  class(result) <- "ResultActivity"
  result
}

utils::globalVariables(c("v", "mu_mix", "kappa_mix", "activity_density", "x", "variable", "sd", "mean"))
