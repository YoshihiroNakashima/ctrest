#' Activity level estimation by von Misex mixture distribution using stick-breaking prior
#'
#' @param activity_data A vector of detection times transformed into radians
#' @param C ""
#' @param plot If TRUE, plots the expected values of residence times.
#' @param iter The number of iterations. The default value is 2000.
#' @param cores The number of cores used for calculations.
#' @param warmup The number of warmup iterations. The default value is NULL (half of the iterations will be used for warmup).
#' @param chains The number of chains. The default value is 2.
#' @param thin The thinning interval. The default value is 1.
#' @param all_comb If TRUE, models with all combinations of covariates are compared. If FALSE, only the designated model in model_formula is run.
#' @param target_species Species name of interest.
#' @return Activity level
#' @import dplyr ggplot2 nimble
#' @importFrom tidyr unite extract
#' @importFrom purrr map
#' @importFrom stringr str_extract
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
#' library(nimble)
#' bayes_activity(
#'   activity_data = activity_data,
#'   C = 10,
#'   plot = TRUE,
#'   iter = 2000,
#'   warmup = 1000,
#'   chains = 2,
#'   thin = 1,
#'   all_comb = FALSE,
#'   target_species = "SP01")

bayes_activity <- function(
    activity_data = activity_data,
    C = 10,
    plot = TRUE,
    cores = 3,
    iter = 3000,
    warmup = 1000,
    chains = 3,
    thin = 5,
    target_species = NULL
) {
  ni <- iter
  nt <- thin
  nc <- chains
  nb <- warmup

  dens.x <- seq(0, 2 * pi, 0.02)
  ndens <- length(dens.x)

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
      mu_mix[k] ~ dunif(0, 2 * 3.141592654)  # von Misesの平均
      kappa_mix[k] ~ dgamma(1, 0.01)  # von Misesの集中度パラメータ
    }
    for(n in 1:N) {
      group[n] ~ dcat(w[1:C])
      time[n] ~ dvonMises(mu_mix[group[n]], kappa_mix[group[n]])
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

  # 結果の処理
  out_vonMises <- coda::mcmc.list(actv_chain_output) %>% posterior::as_draws()

  model_summary <- summary(out_vonMises)　%>%
    mutate(cv = sd / mean) %>%
    select(variable, mean, median, sd, cv, mad, q5, q95, rhat, ess_bulk, ess_tail) %>%
    mutate(Species = target_species)

  activity_density_estimates <- model_summary %>%
    dplyr::filter(str_starts(variable, "activity_dens")) %>%
    dplyr::mutate(x = dens.x)

  model_act <- fitact(time, bw = 1.5 * bwcalc(time, K = 3), reps = 1)
  dens_est_rw <- data.frame(model_act@pdf)

  tidy_samples <- MCMCvis::MCMCchains(actv_chain_output) %>%
    as_tibble() %>%
    pivot_longer(cols = everything(),
                 names_to = "parameter",
                 values_to = "value") %>%
    group_by(parameter) %>%
    mutate(iteration = row_number()) %>%
    ungroup()

  if (plot == TRUE) {
    g <- ggplot() +
      geom_histogram(data = data.frame(time = time),
                     mapping = aes(x = time, y = after_stat(density)),
                     bins = 30, fill = "gray") +
      geom_line(data = activity_density_estimates,
                mapping = aes(x = x, y = mean, colour = "vonMises mixture"),
                linewidth = 1) +
      geom_line(data = dens_est_rw,
                mapping = aes(x = x, y = y, colour = "Kernel density"),
                linewidth = 1) +
      geom_ribbon(data = activity_density_estimates,
                  mapping = aes(x = x, ymin = q5, ymax = q95),
                  fill = "red", alpha = 0.3) +
      scale_color_manual(name = "Density Estimations",
                         values = c("vonMises mixture" = "red",
                                    "Kernel density" = "blue")) +
      labs(x = "Radian time",
           y = "Probability density")

    print(g)
  }
  # actv_out_trace <- actv_chain_output %>%
  #   map(~ .[, grep(paste("actv", collapse = "|"), colnames(.))])

  return(list(model_summary = model_summary,
              plot = g,
              tidy_samples = tidy_samples,
              # actv_out_trace = actv_out_trace,
              actv_chain_output = actv_chain_output)
  )
}

