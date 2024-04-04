#' Virtual dataset of camera trapping detection data
#'
#' An example of dataset used for estimating animal density using REST model
#'
#' @format ## `detection_data`
#' A data frame with 5 columns:
#' \describe{
#'   \item{station}{station id detecting an animal pass}
#'   \item{datetime}{datatime infomation of detections}
#'   \item{species}{species name detected}
#'   \item{pass_number}{the number of times an animal passed through the focal area in a single video}
#'   \item{stay}{staying time of animals within a focal area}
#' }
#' @source {simulated data}
"detection_data"
