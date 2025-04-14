#' Activity level estimation by von Mises mixture distribution using stick-breaking prior
#'
#' @param activity_data A data frame containing a "time" column, representing detection times transformed into radians. Typically, this is the output of the format_activity function.
#' @param C The maximum number of components in the von Mises mixture distribution. Defaults to 10.
#' @param cores The number of CPU cores to use for parallel computation. Default is 3.
#' @param iter The total number of MCMC iterations per chain. Default is 5000
#' @param warmup The number of warm-up (burn-in) iterations per chain. Default is 1000.
#' @param chains The number of MCMC chains. Default is 2.
#' @param target_species The species name of interest.  Only one species can be specified.
#' @return the estimated proportion of activity time and A figure showing the temporal variation in activity level.
#' @import dplyr ggplot2 nimble
#' @importFrom tidyr unite extract
#' @importFrom purrr map
#' @importFrom stringr str_extract
#' @importFrom stats rbeta rgamma
#' @importFrom stats median model.response quantile
#' @export
#' @examples
#' activity_data <- format_activity(
#'   detection_data = detection_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_datetime = "DateTime",
#'   target_species = "SP01",
#'   indep_time = 30
#' )
#' bayes_activity(
#'   activity_data = activity_data,
#'   C = 10,
#'   cores = 3,
#'   iter = 2000,
#'   warmup = 1000,
#'   chains = 3,
#'   thin = 2,
#'   target_species = "SP01")

bayes_activity <- function(
    activity_data,
    C = 10,
    cores = 3,
    iter = 5000,
    warmup = 1000,
    chains = 3,
    thin = 2,
    target_species = NULL
) {

  dens.x <- seq(0, 2 * pi, 0.02)
  ndens <- length(dens.x)

  if (is.null(target_species)) {
    stop("target_species must be specified.")
  }

  activity_data <- activity_data %>%
    filter(Species == target_species)

  time <- activity_data$time
  N <- length(time)
  constants <- list(N = N, C = C, dens.x = dens.x, ndens = ndens)
  data <- list(time = time)

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
      time[n] ~ dvonMises(mu_mix[group[n]], kappa_mix[group[n]])
    }
    for (j in 1:ndens) {
      for (i in 1:C) {
        dens.cpt[i, j] <- w[i] * dvonMises(dens.x[j] , mu_mix[i], kappa_mix[i], log = 0)
      }
      activity_dens[j] <- sum(dens.cpt[1:C, j])
    }
    activity_proportion <- 1.0 / (2 * 3.141592654 * max(activity_dens[1:ndens]));
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

  run_MCMC_vonMises <- function(info, data, constants, code, params, iter, thin, warmup) {
    myModel <- nimbleModel(code = code,
                           data = data,
                           constants = constants,
                           inits = info$inits)

    CmyModel <- compileNimble(myModel)

    configModel <- configureMCMC(myModel, monitors = params)
    myMCMC <- buildMCMC(configModel, monitors = params)
    CmyMCMC <- compileNimble(myMCMC)

    results <- runMCMC(CmyMCMC, niter = iter, nburnin = warmup, thin = thin, nchains = 1, setSeed = info$seed, samplesAsCodaMCMC = TRUE)
  }

  per_chain_info <- lapply(1:chains, function(i) {
    list(
      seed = sample(1:9999, 1),
      inits = inits_f()
    )
  })

  params <- c("activity_dens", "activity_proportion", "mu_mix", "kappa_mix", "w")
  cat("Running MCMC sampling. Please wait...\n")

  this_cluster <- makeCluster(chains)
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
                            iter = iter, thin = thin, warmup = warmup)
  # on.exit(stopCluster(this_cluster), add = TRUE)
  stopCluster(this_cluster)
  cat("Estimation is finished!\n")

  # Summarize results -------------------------------------------------
  summary_activity <- MCMCvis::MCMCchains(actv_chain_output,
                                          mcmc.list = TRUE,
                                          params = "activity_proportion") %>%
    MCMCsummary(., round = 2) %>%
    tibble::rownames_to_column(., var = "Variable") %>%
    tibble(.) %>%
    rename(lower = `2.5%`, median = `50%`, upper = `97.5%`)

  mcmc_samples <- MCMCvis::MCMCchains(actv_chain_output,
                                      mcmc.list = TRUE,
                                      params = "activity_dens")

  summary_dens <- MCMCsummary(mcmc_samples, round = 5) %>%
    tibble::rownames_to_column(., var = "Variable") %>%
    tibble(.) %>%
    rename(lower = `2.5%`, median = `50%`, upper = `97.5%`) %>%
    mutate(x = dens.x)

  model_act <- fitact(time, bw = 1.5 * bwcalc(time, K = 3), reps = 1)
  dens_est_rw <- data.frame(model_act@pdf)

  g <- ggplot() +
  geom_histogram(data = data.frame(time = time),
                 mapping = aes(x = time, y = after_stat(density)),
                 bins = 30, fill = "gray") +
  geom_line(data = summary_dens,
            mapping = aes(x = x, y = mean, colour = "vonMises mixture"),
            linewidth = 1) +
  geom_line(data = dens_est_rw,
            mapping = aes(x = x, y = y, colour = "Kernel density"),
            linewidth = 1) +
  geom_ribbon(data = summary_dens,
              mapping = aes(x = x, ymin = lower, ymax = upper),
              fill = "red", alpha = 0.3) +
  scale_color_manual(name = "Density Estimations",
                     values = c("vonMises mixture" = "red",
                                "Kernel density" = "blue")) +
  labs(x = "Radian time",
       y = "Probability density")

  return(
    list(summary_result = summary_activity,
         plot = g
    )
  )
}
v <- mu_mix <- kappa_mix <- variable <- cv <- mad <- q5 <- q95 <- rhat <- ess_bulk <- ess_tail <- density <- x <- NULL
