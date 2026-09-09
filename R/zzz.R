# Package-level global variable declarations (suppresses R CMD CHECK NOTEs)
utils::globalVariables(c(
  # NIMBLE node names used inside nimbleCode({}) blocks
  "alpha_dir", "r", "eps", "Y", ".", "WAIC",
  "group", "index", "lower", "parameter", "theta", "u", "upper",
  "returnType", "2.5%", "50%", "97.5%",
  "Indep", "N_detection_rep", "Variable",
  "act_data_pred", "alpha_Dirichlet",
  "beta_dens", "beta_density", "beta_stay",
  "censored", "group_stay", "loglact",
  "n.eff", "pred_t", "pred_y", "random_effect_stay",
  "size", "species_effect_density", "species_effect_stay",
  "theta_stay", "y_rep",
  # dplyr / rlang pronouns used inside mutate/filter
  ".parsed_dt", ".diff_min", ".indep", ".pick_order"
))

# The log<- replacement function is used by NIMBLE's log(x) <- ... syntax.
if (getRversion() >= "2.15.1") {
  utils::globalVariables("log<-")
}

# Custom distribution: Dirichlet-Multinomial
# Used by bayes_rest_multi for the RAD-REST pass-count model.
# Workers redefine these via clusterEvalQ; this definition makes them available
# in the main process and for roxygen/devtools::check() to find.
ddirchmulti <- nimble::nimbleFunction(
  run = function(x     = double(1),
                 alpha = double(1),
                 size  = double(0),
                 log   = integer(0)) {
    returnType(double(0))
    logProb <- lgamma(size + 1) - sum(lgamma(x + 1)) +
               lgamma(sum(alpha)) - sum(lgamma(alpha)) +
               sum(lgamma(alpha + x)) - lgamma(sum(alpha) + size)
    if (log) return(logProb) else return(exp(logProb))
  }
)

rdirchmulti <- nimble::nimbleFunction(
  run = function(n     = integer(0),
                 alpha = double(1),
                 size  = double(0)) {
    returnType(double(1))
    if (n != 1) print("rdirchmulti only allows n = 1; using n = 1.")
    p <- rdirch(1, alpha)
    return(rmulti(1, size = size, prob = p))
  }
)
