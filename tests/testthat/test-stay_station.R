# Tests of the station-level design matrix of the staying time model (no MCMC is run).
# The mean staying time of each station used in the density formula must be calculated from the
# covariates of the station, not from those of the staying time records.

station <- data.frame(
  Station = sprintf("S%02d", 1:6),
  x  = c(-1, 0, 1, 2, -2, 0.5),
  g  = c("a", "b", "c", "a", "b", "c"),
  stringsAsFactors = FALSE
)
# staying time records: several per station, in an order that differs from the station order
stay <- data.frame(
  Station = rep(c("S03", "S01", "S06", "S02", "S05", "S04"), times = c(2, 3, 1, 4, 2, 2)),
  Stay    = c(1:14)
)
stay <- merge(stay, station, by = "Station", sort = FALSE)

make_obs <- function(formula, data) {
  mf <- stats::model.frame(formula, data)
  list(mf = mf, X = stats::model.matrix(stats::as.formula(formula), mf))
}

test_that("the matrix has one row per station with the covariates of that station", {
  o <- make_obs(Stay ~ 1 + x, stay)
  X <- .station_design_matrix(Stay ~ 1 + x, o$mf, o$X, station)
  expect_equal(dim(X), c(6, 2))
  expect_equal(colnames(X), colnames(o$X))
  expect_equal(X[, "x"], station$x)
  expect_equal(X[, "(Intercept)"], rep(1, 6))
  # the i-th staying time record is NOT the i-th station
  expect_false(isTRUE(all.equal(o$X[seq_len(6), "x"], station$x)))
})

test_that("factor covariates give the same columns, levels and contrasts as the records", {
  o <- make_obs(Stay ~ 1 + x + g, stay)
  X <- .station_design_matrix(Stay ~ 1 + x + g, o$mf, o$X, station)
  expect_equal(colnames(X), colnames(o$X))                    # (Intercept), x, gb, gc
  expect_equal(X[, "gb"], as.numeric(station$g == "b"))
  expect_equal(X[, "gc"], as.numeric(station$g == "c"))
  # a level that has no staying time record cannot be used
  stay2 <- stay[stay$g != "c", ]
  o2 <- make_obs(Stay ~ 1 + g, stay2)
  expect_error(.station_design_matrix(Stay ~ 1 + g, o2$mf, o2$X, station), "do not match")
})

test_that("the intercept-only model gives one row per station", {
  o <- make_obs(Stay ~ 1, stay)
  X <- .station_design_matrix(Stay ~ 1, o$mf, o$X, station)
  expect_equal(dim(X), c(6, 1))
})

test_that("missing covariates give informative errors", {
  o <- make_obs(Stay ~ 1 + x, stay)
  expect_error(.station_design_matrix(Stay ~ 1 + x, o$mf, o$X, station[, c("Station", "g")]),
               "must also be columns of `station_effort_data`")
  station_na <- station; station_na$x[2] <- NA
  expect_error(.station_design_matrix(Stay ~ 1 + x, o$mf, o$X, station_na), "missing values")
})
