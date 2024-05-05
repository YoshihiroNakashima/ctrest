#' Virtual dataset of camera trapping detection data
#'
#' An example of dataset used for estimating animal density using REST/REST-RAD model
#'
#' @format ## `detection_data`
#' A data frame with 5 columns:
#' \describe{
#'   \item{Station}{station id detecting an animal pass}
#'   \item{DateTime}{datatime infomation of detections}
#'   \item{Term}{survey terms}
#'   \item{Species}{species name detected}
#'   \item{Nfocal}{the number of times an animal passed through the focal area in a single video}
#'   \item{Stay}{staying time of animals within a focal area}
#'   \item{Cens}{censored timing for staying time}
#' }
#' @source {simulated data}
"detection_data"
