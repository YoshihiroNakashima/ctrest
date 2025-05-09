% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_station_data.R
\name{format_station_data}
\alias{format_station_data}
\title{Aggregate the number of animal passes within focal areas for each station from detection data.}
\usage{
format_station_data(
  detection_data,
  station_data,
  col_name_station,
  col_name_species,
  col_name_y,
  model
)
}
\arguments{
\item{detection_data}{A data frame containing information for individual videos, with six columns:
\itemize{
\item Camera station ID (character)
\item Capture datetime (character, not used in this function)
\item Species name (character)
\item Number of passes through the focal area (numeric, not used in this function)
\item Staying time (seconds) within the focal area for each pass in a video (numeric)
\item Whether the observation of the pass was censored (1 = censored, 0 = observed)
}}

\item{station_data}{A data frame containing information for each camera station, with one row per station.
It must include at least the camera station IDs (character) and may contain additional columns such as
location coordinates, covariates, and other relevant information. Unlike \code{detection_data}, this data frame
lists all the camera stations that have been set up.}

\item{col_name_station}{A string specifying the column name containing station ID information.}

\item{col_name_species}{A string specifying the column name containing detected species names.}

\item{col_name_y}{A string specifying the column name containing the number of animal passes within a focal area.
\itemize{
\item If \code{model = "REST"}, this column must have values for all detections. Any \code{NA} values will be replaced with \code{1}.
\item If \code{model = "RAD-REST"}, only measured detections have count values, and unmeasured ones should be \code{NA}.
}}

\item{model}{A string specifying the model type, either \code{"REST"} or \code{"RAD-REST"}.}
}
\value{
A data frame with aggregated detection counts. The structure varies depending on the model:
\itemize{
\item If \code{model = "REST"}, the data frame contains:
\itemize{
\item \code{Station} (character): Camera station ID.
\item \code{Species} (character): Species name.
\item \code{Y} (numeric): The total number of passes through the focal area for each station-species pair.
\item All columns from \code{station_data} merged by \code{Station}.
}
\item If \code{model = "RAD-REST"}, the data frame contains:
\itemize{
\item \code{Station} (character): Camera station ID.
\item \code{Species} (character): Species name.
\item \code{N} (numeric): The total number of detected videos for each station-species pair.
\item \code{y_X} columns (\code{y_0}, \code{y_1}, ..., \code{y_max}): Number of detections classified by pass count.
\item All columns from \code{station_data} merged by \code{Station}.
If no matching \code{Station} is found in \code{station_data}, the corresponding values will be \code{NA}.
}
}
}
\description{
Aggregate the number of animal passes within focal areas for each station from detection data.
}
\examples{
format_station_data(
  detection_data = detection_data,
  station_data = station_data,
  col_name_station = "Station",
  col_name_species = "Species",
  col_name_y = "y",
  model = "REST"
)
}
