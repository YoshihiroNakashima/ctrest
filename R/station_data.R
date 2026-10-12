#' Virtual dataset of camera station data
#'
#' An example of dataset containing info on camera stations. The density of the
#' species in \code{\link{detection_data}} was simulated to depend on \code{x1},
#' \code{x2} and \code{Transect} (see \code{\link{detection_data}}), so the data
#' can be used to try covariates and transect random effects in the density
#' model (\code{formula_density}, \code{random_effect_density}).
#'
#' @format ## `station_data`
#' A data frame with 100 rows and 4 columns:
#' \describe{
#'   \item{Station}{station id detecting an animal pass}
#'   \item{x1}{a continious covairate}
#'   \item{x2}{a categorical covairate ("A", "B" or "C")}
#'   \item{Transect}{transect id ("TR01" to "TR10"; 10 cameras per transect)}
#' }
#' @source {simulated data}
"station_data"
