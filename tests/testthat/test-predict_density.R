# Tests of predict.ResultDensity() with hand-made posterior draws (no MCMC is run).
# Expected values are computed by explicit loops, independently of the package code.

make_mcmc_list <- function(draws, n_chains = 2) {
  # draws: matrix (rows = draws) split evenly into chains
  idx <- split(seq_len(nrow(draws)), rep(seq_len(n_chains), length.out = nrow(draws)))
  coda::mcmc.list(lapply(idx, function(i) coda::mcmc(draws[i, , drop = FALSE])))
}

# ---- single-species model: continuous + factor covariates, transect effects ----

set.seed(11)
n_draw <- 60
dat1 <- data.frame(
  x  = seq(0, 100, length.out = 30),
  g  = rep(c("a", "b", "c"), 10),
  Tr = rep(c("t1", "t2", "t3"), each = 10),
  stringsAsFactors = FALSE
)
dm1 <- .make_density_model(~ 1 + x + g, dat1, species = "SP", multi = FALSE,
                           random_effect_density = "Tr",
                           transect_levels = levels(factor(dat1$Tr)),
                           camera_data = dat1)
draws1 <- cbind(
  matrix(rnorm(n_draw * 4, 0, 0.3), n_draw, 4,
         dimnames = list(NULL, paste0("beta_density[", 1:4, "]"))),
  matrix(rnorm(n_draw * 3, 0, 0.2), n_draw, 3,
         dimnames = list(NULL, paste0("transect_effect_density[", 1:3, "]"))),
  sigma_camera_density   = runif(n_draw, 0.4, 1.0),
  sigma_transect_density = runif(n_draw, 0.1, 0.5),
  junk = rnorm(n_draw)
)
obj1 <- structure(list(samples = make_mcmc_list(draws1), density_model = dm1,
                       target_species = "SP"), class = "ResultDensity")

# draws x rows matrix of the density. type "typical": exp(eta [+ transect effect]);
# type "average": additionally exp(sigma_camera^2 / 2) (and exp(sigma_transect^2 / 2) if the
# transect is not given)
expected_single <- function(draws, x, g, tr = NULL, type = "typical") {
  # model.matrix columns: (Intercept), x, gb, gc
  out <- matrix(NA_real_, nrow(draws), length(x))
  for (r in seq_len(nrow(draws))) {
    for (k in seq_along(x)) {
      eta <- draws[r, "beta_density[1]"] + draws[r, "beta_density[2]"] * x[k] +
        draws[r, "beta_density[3]"] * (g[k] == "b") + draws[r, "beta_density[4]"] * (g[k] == "c")
      if (!is.null(tr)) eta <- eta + draws[r, paste0("transect_effect_density[", tr[k], "]")]
      if (type == "average") {
        eta <- eta + draws[r, "sigma_camera_density"]^2 / 2
        if (is.null(tr)) eta <- eta + draws[r, "sigma_transect_density"]^2 / 2
      }
      out[r, k] <- exp(eta)
    }
  }
  out
}
# the chains are stacked in order, so reorder the hand-made draws the same way
chain_order <- unlist(split(seq_len(n_draw), rep(1:2, length.out = n_draw)), use.names = FALSE)
draws1_ord <- draws1[chain_order, ]

summ_of <- function(e) {
  data.frame(mean = colMeans(e), sd = apply(e, 2, sd),
             lower = apply(e, 2, quantile, 0.025, names = FALSE),
             median = apply(e, 2, median),
             upper = apply(e, 2, quantile, 0.975, names = FALSE))
}

test_that("typical density equals exp(linear predictor), original covariate scale", {
  nd <- data.frame(x = c(10, 55, 90), g = c("a", "b", "c"))
  res <- predict(obj1, newdata = nd, type = "typical")
  e <- expected_single(draws1_ord, nd$x, nd$g)
  expect_equal(res$mean,   colMeans(e), tolerance = 1e-10)
  expect_equal(res$median, apply(e, 2, median), tolerance = 1e-10)
  expect_equal(res$sd,     apply(e, 2, sd), tolerance = 1e-10)
  expect_equal(res$lower,  apply(e, 2, quantile, 0.025, names = FALSE), tolerance = 1e-10)
  expect_equal(res$upper,  apply(e, 2, quantile, 0.975, names = FALSE), tolerance = 1e-10)
  expect_equal(res$Species, rep("SP", 3))
  expect_equal(res$Type, rep("typical", 3))
  expect_equal(res$x, nd$x)
})

test_that("average density = typical * exp(sigma^2 / 2); both types are returned next to each other", {
  nd <- data.frame(x = c(10, 55, 90), g = c("a", "b", "c"))
  res <- predict(obj1, newdata = nd)                      # type = "both" (default)
  expect_equal(nrow(res), 6)
  expect_equal(res$Type, rep(c("typical", "average"), 3))
  expect_equal(res$x, rep(nd$x, each = 2))
  for (ty in c("typical", "average")) {
    e <- expected_single(draws1_ord, nd$x, nd$g, type = ty)    # no transect given
    r <- res[res$Type == ty, ]
    expect_equal(r$mean, colMeans(e), tolerance = 1e-10)
    expect_equal(r$upper, apply(e, 2, quantile, 0.975, names = FALSE), tolerance = 1e-10)
  }
  expect_true(all(res$mean[res$Type == "average"] > res$mean[res$Type == "typical"]))
})

test_that("transect effects are added only when the transect column is supplied", {
  nd <- data.frame(x = c(10, 55, 90), g = c("a", "b", "c"), Tr = c("t3", "t1", "t2"))
  # typical: transect effect added
  res <- predict(obj1, newdata = nd, type = "typical")
  e <- expected_single(draws1_ord, nd$x, nd$g, tr = c(3, 1, 2))
  expect_equal(res$mean, colMeans(e), tolerance = 1e-10)
  res0 <- predict(obj1, newdata = nd[, c("x", "g")], type = "typical")
  expect_false(isTRUE(all.equal(res$mean, res0$mean)))
  # average with a given transect: only the camera-level variance is added
  ra <- predict(obj1, newdata = nd, type = "average")
  ea <- expected_single(draws1_ord, nd$x, nd$g, tr = c(3, 1, 2), type = "average")
  expect_equal(ra$mean, colMeans(ea), tolerance = 1e-10)
})

test_that("prob changes the credible interval", {
  nd <- data.frame(x = 10, g = "a")
  r90 <- predict(obj1, nd, prob = 0.90, type = "typical")
  e <- expected_single(draws1_ord, 10, "a")
  expect_equal(r90$lower, unname(quantile(e, 0.05)), tolerance = 1e-10)
  expect_equal(r90$upper, unname(quantile(e, 0.95)), tolerance = 1e-10)
})

test_that("average = TRUE averages the density over rows for each posterior draw", {
  nd <- data.frame(x = c(10, 55, 90), g = c("a", "b", "c"), Tr = c("t3", "t1", "t2"))
  for (ty in c("typical", "average")) {
    res <- predict(obj1, nd, type = ty, average = TRUE)
    e <- expected_single(draws1_ord, nd$x, nd$g, tr = c(3, 1, 2), type = ty)
    s <- summ_of(matrix(rowMeans(e), ncol = 1))
    expect_equal(nrow(res), 1)
    expect_equal(res$n, 3L)
    expect_equal(res$Type, ty)
    expect_equal(res$mean, s$mean, tolerance = 1e-10)
    expect_equal(res$sd, s$sd, tolerance = 1e-10)
    expect_equal(res$lower, s$lower, tolerance = 1e-10)
    expect_equal(res$upper, s$upper, tolerance = 1e-10)
  }
  # the mean of the densities is not the density at the mean covariate (nonlinear)
  mean_cov <- predict(obj1, data.frame(x = mean(nd$x), g = "b", Tr = "t1"), type = "typical")
  expect_false(isTRUE(all.equal(predict(obj1, nd, type = "typical", average = TRUE)$mean, mean_cov$mean)))
})

test_that("by averages within groups", {
  nd <- data.frame(x = c(10, 55, 90, 20), g = c("a", "b", "c", "b"), Tr = c("t3", "t1", "t2", "t1"),
                   grp = c("u", "v", "u", "v"))
  res <- predict(obj1, nd, type = "typical", average = TRUE, by = "grp")
  expect_equal(res$grp, c("u", "v"))
  expect_equal(res$n, c(2L, 2L))
  e <- expected_single(draws1_ord, nd$x, nd$g, tr = c(3, 1, 2, 1))
  expect_equal(res$mean[1], mean(rowMeans(e[, c(1, 3)])), tolerance = 1e-10)
  expect_equal(res$mean[2], mean(rowMeans(e[, c(2, 4)])), tolerance = 1e-10)
  both <- predict(obj1, nd, average = TRUE, by = "grp")
  expect_equal(both$grp, c("u", "u", "v", "v"))
  expect_equal(both$Type, rep(c("typical", "average"), 2))
})

test_that("average = TRUE without newdata uses the cameras of the estimation", {
  res <- predict(obj1, average = TRUE)
  expect_equal(res$n, c(30L, 30L))
  for (ty in c("typical", "average")) {
    e <- expected_single(draws1_ord, dat1$x, dat1$g, tr = as.integer(factor(dat1$Tr)), type = ty)
    expect_equal(res$mean[res$Type == ty], mean(colMeans(e)), tolerance = 1e-10)
  }
  no_cam <- obj1; no_cam$density_model$camera_data <- NULL
  expect_error(predict(no_cam, average = TRUE), "re-estimate")
})

test_that("invalid arguments give informative errors", {
  expect_error(predict(obj1, data.frame(x = 1)), "must contain the covariate")
  expect_error(predict(obj1, data.frame(x = 1, g = "z")), "new level")
  expect_error(predict(obj1, data.frame(x = NA_real_, g = "a")), "missing values")
  expect_error(predict(obj1, data.frame(x = 1, g = "a", Tr = "t9")), "not found in the fitted model")
  expect_error(predict(obj1, data.frame(x = 1, g = "a", Species = "SP")), "Species")
  expect_error(predict(obj1, data.frame(x = 1, g = "a", Type = "x")), "Type")
  expect_error(predict(obj1, "a"), "data frame")
  expect_error(predict(obj1, data.frame(x = 1, g = "a"), prob = 1.2), "prob")
  expect_error(predict(obj1, data.frame(x = 1, g = "a"), species = "XX"))
  expect_error(predict(obj1, data.frame(x = 1, g = "a"), type = "mean"))
  expect_error(predict(obj1, data.frame(x = 1, g = "a"), by = "x"), "average = TRUE")
  expect_error(predict(obj1, data.frame(x = 1, g = "a"), average = TRUE, by = "zzz"), "not found")
  expect_error(predict(obj1, data.frame(x = 1, g = "a"), average = NA), "average")
})

test_that("results without the variances still give the typical density", {
  no_sigma <- draws1[, !grepl("^sigma_", colnames(draws1))]
  obj <- obj1; obj$samples <- make_mcmc_list(no_sigma)
  nd <- data.frame(x = 10, g = "a")
  expect_equal(nrow(predict(obj, nd, type = "typical")), 1)
  expect_error(predict(obj, nd, type = "average"), "sigma_camera_density")
  expect_error(predict(obj, nd), "re-estimate")
})

test_that("old results without density_model give a clear error", {
  old <- structure(list(samples = obj1$samples), class = "ResultDensity")
  expect_error(predict(old), "re-estimate")
})

test_that("default newdata gives one profile per covariate", {
  res <- predict(obj1, type = "typical")
  # x is continuous (25-point grid), g has 3 levels
  expect_equal(nrow(res), 25 + 3)
  expect_true(all(c("Profile", "x", "g") %in% names(res)))
  px <- res[res$Profile == "x", ]
  expect_equal(range(px$x), range(dat1$x))
  expect_equal(length(unique(px$g)), 1)         # the other covariate is held at one value
  pg <- res[res$Profile == "g", ]
  expect_setequal(as.character(pg$g), c("a", "b", "c"))
  expect_equal(length(unique(pg$x)), 1)
  expect_equal(unique(pg$x), mean(dat1$x))
  expect_equal(nrow(predict(obj1)), 2 * (25 + 3))      # both types
})

test_that("the printed table explains the column Type", {
  out <- capture.output(print(predict(obj1, data.frame(x = 10, g = "a"))))
  expect_true(any(grepl("Type: typical", out)))
  expect_true(any(grepl("average = mean density", out)))
  expect_s3_class(predict(obj1, data.frame(x = 10, g = "a")), "tbl_df")
})

test_that("draws = TRUE returns the posterior draws behind the summaries", {
  nd <- data.frame(x = c(10, 55), g = c("a", "b"), Tr = c("t3", "t1"))
  d <- predict(obj1, nd, type = "average", draws = TRUE)
  expect_equal(nrow(d), 2 * n_draw)
  expect_true(all(c("Species", "x", "g", "Tr", "Type", "draw", "density") %in% names(d)))
  expect_false(any(c("mean", "sd", "lower", "median", "upper") %in% names(d)))
  e <- expected_single(draws1_ord, nd$x, nd$g, tr = c(3, 1), type = "average")
  d1 <- d[d$x == 10, ]
  expect_equal(d1$density[order(d1$draw)], e[, 1], tolerance = 1e-10)       # draw = index of the posterior draw
  expect_equal(sort(unique(d$draw)), seq_len(n_draw))
  # the summaries are the summaries of the draws
  s <- predict(obj1, nd, type = "average")
  expect_equal(s$mean[1], mean(d1$density), tolerance = 1e-10)
  expect_equal(s$upper[1], unname(quantile(d1$density, 0.975)), tolerance = 1e-10)
  expect_false(inherits(d, "ctrest_prediction"))
})

test_that("draws with average = TRUE averages draw by draw within groups", {
  nd <- data.frame(x = c(10, 55, 90, 20), g = c("a", "b", "c", "b"), Tr = c("t3", "t1", "t2", "t1"),
                   grp = c("u", "v", "u", "v"))
  d <- predict(obj1, nd, type = "typical", average = TRUE, by = "grp", draws = TRUE, n_draws = 10)
  expect_equal(nrow(d), 2 * 10)
  expect_equal(d$n[d$grp == "u"][1], 2L)
  e <- expected_single(draws1_ord, nd$x, nd$g, tr = c(3, 1, 2, 1))
  du <- d[d$grp == "u", ]
  expect_equal(du$density, rowMeans(e[du$draw, c(1, 3)]), tolerance = 1e-10)
})

test_that("invalid draws arguments give errors", {
  nd <- data.frame(x = 1, g = "a")
  expect_error(predict(obj1, nd, n_draws = 5), "draws = TRUE")
  expect_error(predict(obj1, nd, draws = NA), "draws")
  expect_error(predict(obj1, nd, draws = TRUE, n_draws = 0), "positive")
})

# ---- multispecies model: standardized covariate, species effects, transect effects ----

set.seed(12)
dat2 <- data.frame(Year = rep(c(0, 1), each = 12), Tr = rep(c("A", "B", "C", "D"), 6),
                   stringsAsFactors = FALSE)
cen <- c("(Intercept)" = 0, Year = mean(dat2$Year))
sca <- c("(Intercept)" = 1, Year = sd(dat2$Year))
dm2 <- .make_density_model(~ 1 + Year, dat2, center = unname(cen), scale = unname(sca),
                           species = c("S1", "S2", "S3"), multi = TRUE,
                           random_effect_density = "Tr", transect_levels = c("A", "B", "C", "D"),
                           camera_data = dat2)
nm <- c(paste0("beta_density[", 1:2, "]"),
        paste0("species_effect_density[", rep(1:3, each = 2), ", ", rep(1:2, 3), "]"),
        paste0("transect_effect_density[", rep(1:4, 3), ", ", rep(1:3, each = 4), "]"))
draws2 <- cbind(
  matrix(rnorm(n_draw * length(nm), 0, 0.3), n_draw, length(nm), dimnames = list(NULL, nm)),
  sigma_camera_density = runif(n_draw, 0.4, 1.0),
  sigma_transect_density = runif(n_draw, 0.1, 0.5)
)
obj2 <- structure(list(samples = make_mcmc_list(draws2), density_model = dm2,
                       target_species = c("S1", "S2", "S3")), class = "ResultDensity")
draws2_ord <- draws2[chain_order, ]

expected_multi <- function(draws, m, year, tr = NULL, type = "typical") {
  z <- (year - cen[["Year"]]) / sca[["Year"]]
  vapply(seq_len(nrow(draws)), function(r) {
    b0 <- draws[r, "beta_density[1]"] + draws[r, sprintf("species_effect_density[%d, 1]", m)]
    b1 <- draws[r, "beta_density[2]"] + draws[r, sprintf("species_effect_density[%d, 2]", m)]
    eta <- b0 + b1 * z
    if (!is.null(tr)) eta <- eta + draws[r, sprintf("transect_effect_density[%d, %d]", tr, m)]
    if (type == "average") {
      eta <- eta + draws[r, "sigma_camera_density"]^2 / 2
      if (is.null(tr)) eta <- eta + draws[r, "sigma_transect_density"]^2 / 2
    }
    exp(eta)
  }, numeric(1))
}

test_that("multispecies prediction uses species-specific coefficients and the stored scaling", {
  res <- predict(obj2, data.frame(Year = c(0, 1)))
  expect_equal(nrow(res), 12)
  expect_equal(res$Species, rep(c("S1", "S2", "S3"), each = 4))
  for (m in 1:3) for (yr in 0:1) for (ty in c("typical", "average")) {
    row <- res[res$Species == paste0("S", m) & res$Year == yr & res$Type == ty, ]
    expect_equal(row$mean, mean(expected_multi(draws2_ord, m, yr, type = ty)), tolerance = 1e-10)
  }
})

test_that("multispecies transect-level prediction picks the species x transect effect", {
  nd <- expand.grid(Tr = c("D", "B"), Year = c(0, 1), stringsAsFactors = FALSE)
  res <- predict(obj2, nd, species = c("S3", "S1"), type = "typical")
  expect_equal(unique(res$Species), c("S3", "S1"))
  tr_index <- c(A = 1, B = 2, C = 3, D = 4)
  for (k in seq_len(nrow(res))) {
    m <- as.integer(substr(res$Species[k], 2, 2))
    expect_equal(res$mean[k],
                 mean(expected_multi(draws2_ord, m, res$Year[k], tr = tr_index[[res$Tr[k]]])),
                 tolerance = 1e-10)
  }
})

test_that("multispecies averages are species-specific", {
  nd <- data.frame(Tr = c("A", "C"), Year = c(0, 1))
  res <- predict(obj2, nd, type = "average", average = TRUE)
  expect_equal(res$Species, c("S1", "S2", "S3"))
  tr_index <- c(A = 1, C = 3)
  for (m in 1:3) {
    e <- cbind(expected_multi(draws2_ord, m, 0, tr = 1, type = "average"),
               expected_multi(draws2_ord, m, 1, tr = 3, type = "average"))
    expect_equal(res$mean[m], mean(rowMeans(e)), tolerance = 1e-10)
  }
  by_year <- predict(obj2, average = TRUE, by = "Year", type = "typical")   # the cameras of the estimation
  expect_equal(nrow(by_year), 6)
  expect_equal(by_year$n, rep(12L, 6))
})

test_that("n_draws returns the same random draws for every row, species and type", {
  nd <- data.frame(Year = c(0, 1))
  set.seed(5); d <- predict(obj2, nd, draws = TRUE, n_draws = 20)
  expect_equal(length(unique(d$draw)), 20)
  expect_equal(nrow(d), 3 * 2 * 2 * 20)                    # species x rows x types x draws
  sets <- tapply(d$draw, list(d$Species, d$Year, d$Type), function(z) paste(sort(z), collapse = ","))
  expect_equal(length(unique(as.vector(sets))), 1)           # joint draws
  set.seed(5); expect_equal(d, predict(obj2, nd, draws = TRUE, n_draws = 20))   # reproducible
  # the density in a draw is the density of that posterior draw
  full <- predict(obj2, nd, draws = TRUE)
  key <- function(z) paste(z$Species, z$Year, z$Type, z$draw)
  expect_equal(d$density, full$density[match(key(d), key(full))], tolerance = 1e-12)
  # more draws than available returns all of them
  expect_equal(nrow(predict(obj2, nd, draws = TRUE, n_draws = 1000, type = "typical")), 3 * 2 * n_draw)
})

test_that("default newdata for a discrete covariate returns each value", {
  res <- predict(obj2, type = "typical")
  expect_equal(nrow(res), 3 * 2)
  expect_false("Profile" %in% names(res))      # a single covariate needs no Profile column
  expect_setequal(res$Year, c(0, 1))
})

# ---- model without covariates ----

dm0 <- .make_density_model(~ 1, dat2, species = c("S1", "S2", "S3"), multi = TRUE,
                           random_effect_density = "Tr", transect_levels = c("A", "B", "C", "D"),
                           camera_data = dat2)
obj0 <- structure(list(samples = make_mcmc_list(draws2), density_model = dm0,
                       target_species = c("S1", "S2", "S3")), class = "ResultDensity")

test_that("a model without covariates predicts one density per species", {
  res <- predict(obj0)
  expect_equal(nrow(res), 6)
  for (m in 1:3) {
    eta <- draws2_ord[, "beta_density[1]"] + draws2_ord[, sprintf("species_effect_density[%d, 1]", m)]
    expect_equal(res$mean[res$Species == paste0("S", m) & res$Type == "typical"], mean(exp(eta)), tolerance = 1e-10)
    adj <- draws2_ord[, "sigma_camera_density"]^2 / 2 + draws2_ord[, "sigma_transect_density"]^2 / 2
    expect_equal(res$mean[res$Species == paste0("S", m) & res$Type == "average"], mean(exp(eta + adj)), tolerance = 1e-10)
  }
  tr <- predict(obj0, data.frame(Tr = c("A", "C")), species = "S2", type = "typical")
  expect_equal(nrow(tr), 2)
  e <- exp(draws2_ord[, "beta_density[1]"] + draws2_ord[, "species_effect_density[2, 1]"] +
             draws2_ord[, "transect_effect_density[3, 2]"])
  expect_equal(tr$mean[2], mean(e), tolerance = 1e-10)
  # mean over the cameras of the estimation
  av <- predict(obj0, average = TRUE, type = "typical", species = "S2")
  expect_equal(av$n, nrow(dat2))
})

# ---- helpers used by bayes_rest() / bayes_rest_multi() ----

test_that(".summary_prediction returns NULL without covariates and warns (not fails) on errors", {
  expect_null(.summary_prediction(obj0))
  expect_s3_class(.summary_prediction(obj2), "tbl_df")
  expect_true(all(c("typical", "average") %in% .summary_prediction(obj2)$Type))
  expect_warning(res <- .summary_prediction(obj2, data.frame(Z = 1)), "could not be computed")
  expect_null(res)
})

test_that(".check_newdata_density validates before estimation", {
  dat <- data.frame(Year = c(0, 1, 0, 1), Transect = c("A", "A", "B", "B"))
  expect_silent(.check_newdata_density(NULL, ~ 1 + Year, dat, "Transect"))
  expect_silent(.check_newdata_density(data.frame(Year = 1, Transect = "B"), ~ 1 + Year, dat, "Transect"))
  expect_error(.check_newdata_density(data.frame(Yr = 1), ~ 1 + Year, dat, "Transect"),
               "`newdata_density` must contain the covariate")
  expect_error(.check_newdata_density(data.frame(Year = 1, Transect = "Z"), ~ 1 + Year, dat, "Transect"),
               "not found")
  expect_error(.check_newdata_density(data.frame(Year = 1, Type = "a"), ~ 1 + Year, dat), "Type")
  expect_error(.check_newdata_density(NULL, ~ 1 + Year, dat, "Nope"), "not a column")
  expect_error(.check_newdata_density(1:3, ~ 1 + Year, dat), "data frame")
})

test_that("a failure to store prediction information warns instead of discarding the fit", {
  expect_warning(
    res <- .make_density_model_safe(~ 1 + Year, dat2, center = c(0, 0, 0), scale = c(1, 1, 1),
                                    species = "S", multi = TRUE),
    "could not be stored")
  expect_null(res)
  obj_null <- structure(list(samples = obj2$samples, density_model = NULL), class = "ResultDensity")
  expect_null(.summary_prediction(obj_null))
  expect_warning(.summary_prediction(obj_null, data.frame(Year = 1)), "could not be computed")
})

test_that("the stored formula is a string, so the result does not carry an environment", {
  expect_type(dm1$formula, "character")
  make_dm <- function() {
    big <- numeric(1e6)
    f <- stats::as.formula("~ 1 + x")
    .make_density_model(f, dat1, species = "SP", multi = FALSE, camera_data = dat1)
  }
  expect_lt(length(serialize(make_dm(), NULL)), 1e5)
})
