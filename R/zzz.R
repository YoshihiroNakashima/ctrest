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
    "returnType"
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

