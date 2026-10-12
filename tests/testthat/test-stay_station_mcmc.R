# Slow tests: the mean staying time of each station used in the density formula.
# It must be calculated from the covariates of the station and the random effect of the group of
# the station (not from the covariates of the staying time records).
# Run with: Sys.setenv(CTREST_SLOW_TESTS = "true")

skip_on_cran()
skip_if_not(nzchar(Sys.getenv("CTREST_SLOW_TESTS")), "set CTREST_SLOW_TESTS=true to run MCMC tests")

make_sim_stay <- function() {
  set.seed(202)
  species <- c("SP1", "SP2", "SP3")
  cam <- expand.grid(Transect = paste0("T", 1:4), cam = 1:5, stringsAsFactors = FALSE)
  cam$Station <- sprintf("%s_%02d", cam$Transect, cam$cam)
  cam$x <- round(rnorm(nrow(cam)), 3)
  cam$Effort <- 60
  cam <- cam[order(cam$Station, method = "radix"), ]
  u <- c(T1 = 0.6, T2 = -0.4, T3 = 0, T4 = 0.3)               # transect effect on log staying time
  S <- 2.67e-6
  st <- list(); sy <- list()
  for (m in seq_along(species)) {
    ml <- log(6) + 0.8 * cam$x + u[cam$Transect]             # log scale of the staying time of each camera
    mean_stay <- exp(ml + 0.5^2 / 2)
    N <- rpois(nrow(cam), 10 * S * cam$Effort * 86400 / mean_stay)
    y <- t(vapply(N, function(n) as.vector(rmultinom(1, n, c(0.1, 0.8, 0.1))), numeric(3)))
    colnames(y) <- paste0("y_", 0:2)
    st[[m]] <- data.frame(Station = cam$Station, Species = species[m], N = N, y,
                          Effort = cam$Effort, x = cam$x, Transect = cam$Transect)
    sy[[m]] <- data.frame(Station = rep(cam$Station, N), Species = species[m],
                          Stay = rlnorm(sum(N), rep(ml, N), 0.5), Cens = 0L,
                          Transect = rep(cam$Transect, N))
  }
  list(station = do.call(rbind, st), stay = do.call(rbind, sy),
       act = data.frame(Species = rep(species, each = 100), Station = "x", time = runif(300, 0, 2 * pi)),
       species = species, cam = cam)
}

sim <- make_sim_stay()
mcmc <- list(iter = 300, warmup = 150, chains = 2, thin = 1)
draws_of <- function(fit) do.call(rbind, lapply(fit$samples, as.matrix))
rest_data <- sim$station[, c("Station", "Species", "Effort", "x", "Transect")]
rest_data$Y <- sim$station$y_1 + 2 * sim$station$y_2          # total number of passes (REST)

fit_single <- function(model, formula_stay, random_effect_stay = NULL) {
  do.call(bayes_rest, c(list(
    formula_stay = formula_stay, formula_density = ~ 1, formula_enter = ~ 1,
    station_effort_data = if (model == "REST") rest_data else sim$station,
    stay_data = sim$stay, activity_data = sim$act, random_effect_stay = random_effect_stay,
    activity_estimation = "kernel", stay_family = "lognormal", focal_area = 2.67,
    model = model, target_species = "SP1"), mcmc))
}

# expected mean staying time of each camera (draws x cameras), lognormal family
expected_mean_stay <- function(S, x = NULL, tr_idx = NULL) {
  eta <- S[, "beta_stay[1]"]
  eta <- matrix(eta, nrow(S), nrow(sim$cam))
  if (!is.null(x)) eta <- eta + outer(S[, "beta_stay[2]"], x)
  if (!is.null(tr_idx)) eta <- eta + S[, sprintf("random_effect_stay[%d]", 1:4), drop = FALSE][, tr_idx]
  exp(eta + S[, "sdlog"]^2 / 2)
}

tr_idx <- as.integer(factor(sim$cam$Transect))     # stations are in the order of sim$cam

for (model in c("REST", "RAD-REST")) {
  test_that(paste0("bayes_rest (", model, "): mean staying time uses the covariates of each station"), {
    fit <- fit_single(model, Stay ~ 1 + x)
    S <- draws_of(fit)
    ms <- S[, sprintf("mean_stay[%d]", seq_len(nrow(sim$cam)))]
    expect_equal(ms, expected_mean_stay(S, x = sim$cam$x), tolerance = 1e-10, ignore_attr = TRUE)
    # the effect of the covariate is recovered (true: 0.8)
    expect_gt(mean(S[, "beta_stay[2]"]), 0.6)
    expect_lt(mean(S[, "beta_stay[2]"]), 1.0)
  })

  test_that(paste0("bayes_rest (", model, "): mean staying time includes the random effect of the group of each station"), {
    fit <- fit_single(model, Stay ~ 1 + x, random_effect_stay = "Transect")
    S <- draws_of(fit)
    ms <- S[, sprintf("mean_stay[%d]", seq_len(nrow(sim$cam)))]
    expect_equal(ms, expected_mean_stay(S, x = sim$cam$x, tr_idx = tr_idx), tolerance = 1e-10, ignore_attr = TRUE)
    # random effect only (no covariate): still station-specific
    fit0 <- fit_single(model, Stay ~ 1, random_effect_stay = "Transect")
    S0 <- draws_of(fit0)
    ms0 <- S0[, sprintf("mean_stay[%d]", seq_len(nrow(sim$cam)))]
    expect_equal(ms0, expected_mean_stay(S0, tr_idx = tr_idx), tolerance = 1e-10, ignore_attr = TRUE)
    expect_gt(max(apply(ms0, 2, mean)) / min(apply(ms0, 2, mean)), 1.1)    # differs among the groups
  })
}

test_that("without covariates and random effects the mean staying time is a single value", {
  fit <- fit_single("REST", Stay ~ 1)
  expect_true("mean_stay" %in% colnames(draws_of(fit)))
  expect_false("mean_stay[1]" %in% colnames(draws_of(fit)))
})

test_that("random_effect_stay must be a column of the station data", {
  stay_obs <- sim$stay; stay_obs$Observer <- "a"
  expect_error(
    bayes_rest(formula_stay = Stay ~ 1, formula_density = ~ 1, formula_enter = ~ 1,
               station_effort_data = sim$station, stay_data = stay_obs, activity_data = sim$act,
               random_effect_stay = "Observer", activity_estimation = "kernel", focal_area = 2.67,
               model = "RAD-REST", target_species = "SP1"),
    "must be a column of both")
})

test_that("bayes_rest_multi: mean staying time uses the covariates and the random effect of each station", {
  fit <- do.call(bayes_rest_multi, c(list(
    formula_stay = Stay ~ 1 + x, formula_density = ~ 1, formula_enter = ~ 1,
    station_effort_data = sim$station, stay_data = sim$stay, activity_data = sim$act,
    random_effect_stay = "Transect", activity_estimation = "kernel", stay_family = "lognormal",
    focal_area = 2.67, model = "RAD-REST", target_species = sim$species), mcmc))
  S <- draws_of(fit)
  cs <- fit$scaling_params$stay
  z <- (sim$cam$x - cs$center[["x"]]) / cs$scale[["x"]]
  re <- S[, sprintf("random_effect_stay[%d]", 1:4)][, tr_idx]
  for (m in seq_along(sim$species)) {
    ms <- S[, sprintf("mean_stay[%d, %d]", seq_len(nrow(sim$cam)), m)]
    b0 <- S[, "beta_stay[1]"] + S[, sprintf("species_effect_stay[%d, 1]", m)]
    b1 <- S[, "beta_stay[2]"] + S[, sprintf("species_effect_stay[%d, 2]", m)]
    expected <- exp(b0 + outer(b1, z) + re + S[, sprintf("sdlog[%d]", m)]^2 / 2)
    expect_equal(ms, expected, tolerance = 1e-10, ignore_attr = TRUE)
  }
})
