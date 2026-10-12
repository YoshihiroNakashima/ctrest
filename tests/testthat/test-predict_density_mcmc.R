# Slow tests: bayes_rest() and bayes_rest_multi() are run on a small simulated data set
# (short MCMC; compiles nimble models) and predict() is compared with the density nodes
# of the fitted model itself. Run with: Sys.setenv(CTREST_SLOW_TESTS = "true")

skip_on_cran()
skip_if_not(nzchar(Sys.getenv("CTREST_SLOW_TESTS")), "set CTREST_SLOW_TESTS=true to run MCMC tests")

make_sim <- function() {
  set.seed(101)
  species <- c("SP1", "SP2", "SP3")
  cam <- expand.grid(Transect = c("A", "B"), cam = 1:4, stringsAsFactors = FALSE)
  cam$x <- round(rnorm(nrow(cam), 500, 100))          # covariate of a camera (same in both years)
  st <- merge(cam, data.frame(Year = c(0, 1)))
  st$Station <- sprintf("%s%02d_%d", st$Transect, st$cam, st$Year)   # one unique ID per camera x year
  st$Effort <- round(runif(nrow(st), 40, 90))
  st <- st[order(st$Station, method = "radix"), ]
  rows <- lapply(species, function(sp) {
    mu <- 8 * exp(0.5 * (sp == "SP1") * st$Year)
    N <- rpois(nrow(st), mu)
    Nj <- rbinom(nrow(st), N, 0.8)
    y <- t(vapply(Nj, function(n) as.vector(rmultinom(1, n, c(0.2, 0.6, 0.2))), numeric(3)))
    colnames(y) <- paste0("y_", 0:2)
    data.frame(Station = st$Station, Species = sp, N = N, y, Effort = st$Effort, x = st$x,
               Year = st$Year, Transect = st$Transect, Y = rpois(nrow(st), mu))
  })
  list(
    station = do.call(rbind, rows),
    stay = do.call(rbind, lapply(species, function(sp) {
      data.frame(Station = sample(st$Station, 40, TRUE), Species = sp,
                 Stay = rlnorm(40, log(8), 0.5), Cens = 0L)
    })),
    act = do.call(rbind, lapply(species, function(sp) {
      data.frame(Species = sp, Station = sample(st$Station, 120, TRUE), time = rnorm(120, pi, 0.8) %% (2 * pi))
    })),
    species = species
  )
}

draws_of <- function(fit) do.call(rbind, lapply(fit$samples, as.matrix))

sim <- make_sim()
mcmc <- list(iter = 300, warmup = 150, chains = 2, thin = 1)

test_that("bayes_rest returns density predicted from covariates", {
  fit <- do.call(bayes_rest, c(list(
    formula_stay = Stay ~ 1, formula_density = ~ 1 + Year, formula_enter = ~ 1,
    station_effort_data = sim$station, stay_data = sim$stay, activity_data = sim$act,
    random_effect_density = "Transect", activity_estimation = "kernel",
    stay_family = "lognormal", focal_area = 2.67, model = "RAD-REST",
    target_species = "SP1"), mcmc))

  expect_s3_class(fit, "ResultDensity")
  expect_s3_class(fit$summary_prediction, "tbl_df")
  expect_equal(nrow(fit$summary_prediction), 4)          # Year = 0, 1; types "typical" and "average"
  expect_setequal(fit$summary_prediction$Type, c("typical", "average"))
  expect_true(all(fit$summary_prediction$mean > 0))

  # predict() reproduces the density[i] nodes of the model at the station covariates
  dat <- sim$station[sim$station$Species == "SP1", ]
  S <- draws_of(fit)
  d <- S[, paste0("density[", seq_len(nrow(dat)), "]")]
  expect_equal(predict(fit, dat["Year"], type = "typical")$mean, unname(colMeans(d)), tolerance = 1e-8)

  # transect-level prediction equals density[i] * exp(transect_effect[t])
  te <- S[, paste0("transect_effect_density[", 1:2, "]")]
  tr <- as.numeric(factor(dat$Transect))
  expect_equal(predict(fit, dat[c("Year", "Transect")], type = "typical")$mean,
               unname(colMeans(d * exp(te[, tr]))), tolerance = 1e-8)

  # average density = typical * exp(sigma^2 / 2) (draw by draw)
  sc <- S[, "sigma_camera_density"]; st <- S[, "sigma_transect_density"]
  expect_equal(predict(fit, dat["Year"], type = "average")$mean,
               unname(colMeans(d * exp(sc^2 / 2 + st^2 / 2))), tolerance = 1e-8)
  expect_equal(predict(fit, dat[c("Year", "Transect")], type = "average")$mean,
               unname(colMeans(d * exp(te[, tr]) * exp(sc^2 / 2))), tolerance = 1e-8)

  # mean over the cameras of the estimation (default newdata) and over a group
  av <- predict(fit, average = TRUE, type = "average")
  expect_equal(av$n, nrow(dat))
  expect_equal(av$mean, mean(colMeans(d * exp(te[, tr]) * exp(sc^2 / 2))), tolerance = 1e-8)
  av_year <- predict(fit, average = TRUE, by = "Year", type = "typical")
  expect_equal(av_year$Year, c(0, 1))
  expect_equal(av_year$mean[1], mean(colMeans(d[, dat$Year == 0] * exp(te[, tr[dat$Year == 0]]))),
               tolerance = 1e-8)

  # the global model keeps summary_transect and returns no prediction table
  fit0 <- do.call(bayes_rest, c(list(
    formula_stay = Stay ~ 1, formula_density = ~ 1, formula_enter = ~ 1,
    station_effort_data = sim$station, stay_data = sim$stay, activity_data = sim$act,
    random_effect_density = "Transect", activity_estimation = "kernel",
    stay_family = "lognormal", focal_area = 2.67, model = "RAD-REST",
    target_species = "SP1"), mcmc))
  expect_null(fit0$summary_prediction)
  expect_false(is.null(fit0$summary_transect))
  expect_equal(predict(fit0, data.frame(Transect = c("A", "B")), type = "typical")$mean,
               fit0$summary_transect$mean, tolerance = 1e-8)
  # both types are available for the global model as well
  expect_equal(nrow(predict(fit0)), 2)
})

test_that("bayes_rest checks newdata_density before the MCMC run", {
  expect_error(
    bayes_rest(formula_stay = Stay ~ 1, formula_density = ~ 1 + Year, formula_enter = ~ 1,
               station_effort_data = sim$station, stay_data = sim$stay, activity_data = sim$act,
               focal_area = 2.67, model = "REST", target_species = "SP1",
               newdata_density = data.frame(Yr = 1)),
    "`newdata_density` must contain the covariate")
})

test_that("bayes_rest_multi returns species-specific density predicted from covariates", {
  fit <- do.call(bayes_rest_multi, c(list(
    formula_stay = Stay ~ 1, formula_density = ~ 1 + Year + x, formula_enter = ~ 1,
    station_effort_data = sim$station, stay_data = sim$stay, activity_data = sim$act,
    random_effect_density = "Transect", activity_estimation = "kernel",
    stay_family = "lognormal", focal_area = 2.67, model = "RAD-REST",
    target_species = sim$species), mcmc))

  expect_s3_class(fit$summary_prediction, "tbl_df")
  expect_setequal(unique(fit$summary_prediction$Species), sim$species)
  expect_true(all(c("Profile", "Year", "x") %in% names(fit$summary_prediction)))

  # predict() reproduces density[i, m] of the model at the station covariates
  dat <- sim$station[sim$station$Species == "SP1", ]
  dat <- dat[order(dat$Station, method = "radix"), ]
  S <- draws_of(fit)
  sc <- S[, "sigma_camera_density"]; st <- S[, "sigma_transect_density"]

  # the log-likelihoods are used for WAIC only and are not returned
  expect_false(any(grepl("^loglike_", colnames(S))))
  expect_true(is.finite(fit$WAIC))

  # posterior draws: one column per row of newdata, all draws, in the order of the draws
  d1 <- S[, sprintf("density[%d, %d]", seq_len(nrow(dat)), 1)]
  dd <- predict(fit, dat[c("Year", "x")], species = sim$species[1], type = "typical", draws = TRUE)
  expect_equal(nrow(dd), nrow(S) * nrow(dat))
  expect_equal(matrix(dd$density, nrow = nrow(S)), unname(d1), tolerance = 1e-10)
  expect_equal(sort(unique(dd$draw)), seq_len(nrow(S)))

  for (m in seq_along(sim$species)) {
    d <- S[, sprintf("density[%d, %d]", seq_len(nrow(dat)), m)]
    expect_equal(predict(fit, dat[c("Year", "x")], species = sim$species[m], type = "typical")$mean,
                 unname(colMeans(d)), tolerance = 1e-8)
    te <- S[, sprintf("transect_effect_density[%d, %d]", 1:2, m)]
    tr <- as.numeric(factor(dat$Transect))
    expect_equal(predict(fit, dat[c("Year", "x", "Transect")], species = sim$species[m], type = "typical")$mean,
                 unname(colMeans(d * exp(te[, tr]))), tolerance = 1e-8)
    # average density over cameras, and its mean over the cameras of the estimation
    expect_equal(predict(fit, dat[c("Year", "x")], species = sim$species[m], type = "average")$mean,
                 unname(colMeans(d * exp(sc^2 / 2 + st^2 / 2))), tolerance = 1e-8)
    av <- predict(fit, average = TRUE, type = "average", species = sim$species[m])
    expect_equal(av$n, nrow(dat))
    expect_equal(av$mean, mean(colMeans(d * exp(te[, tr]) * exp(sc^2 / 2))), tolerance = 1e-8)
  }
})
