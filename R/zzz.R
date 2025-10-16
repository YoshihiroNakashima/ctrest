.onLoad <- function(libname, pkgname) {
  utils::globalVariables(c(
    "alpha_dir",
    "r",
    "eps",
    "Y",
    ".",
    "WAIC",
    "group",
    "index",
    "lower",
    "parameter",
    "theta",
    "u",
    "upper",
    "value",
    "returnType",
    "2.5%",
    "50%",
    "97.5%",
    "Indep",
    "N_detection_rep",
    "Variable",
    "act_data_pred",
    "alpha_Dirichlet",
    "beta_dens",
    "beta_density",
    "beta_stay",
    "censored",
    "group_stay",
    "loglact",
    "n.eff",
    "pred_t",
    "pred_y",
    "random_effect_stay",
    "size",
    "species_effect_density",
    "species_effect_stay",
    "theta_stay",
    "y_rep"
  ))
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("log<-"))
}

ddirchmulti <- nimbleFunction(
  run = function(x = double(1),
                 alpha = double(1),
                 size = double(0),
                 log = integer(0)) {
    returnType(double(0))
    logProb <- lgamma(size + 1) - sum(lgamma(x + 1)) + lgamma(sum(alpha)) -
      sum(lgamma(alpha)) + sum(lgamma(alpha + x)) -
      lgamma(sum(alpha) + size)
    if (log)
      return(logProb)
    else
      return(exp(logProb))
  }
)

rdirchmulti <- nimbleFunction(
  run = function(n = integer(0),
                 alpha = double(1),
                 size = double(0)) {
    returnType(double(1))
    if (n != 1)
      print("rdirchmulti only allows n = 1; using n = 1.")
    p <- rdirch(1, alpha)
    return(rmulti(1, size = size, prob = p))
  }
)
