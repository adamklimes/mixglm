#' Tree cover data for South America
#'
#' Tree cover, precipitation, temperature, and predicted precipitation data for
#' a random set of 5000 locations in the South America.
#'
#' @format ## `treeCover`
#' A data frame with 5,000 rows and 6 columns:
#' \describe{
#'   \item{x}{Longitude in degrees}
#'   \item{y}{Latitude in degrees}
#'   \item{treeCover}{Percentage of tree cover in a 250m large pixel}
#'   \item{precip}{Mean yearly precipitation computed from 1981-2010 in millimeters}
#'   \item{temp}{Mean temperature computed from 1981-2010 in °C}
#'   \item{precip2100}{Predicted yearly precipitation for 2081-2100 based on RCP4.5 scenatio in millimeters}
#' }
#' @source
#' \describe{
#'   \item{treeCover}{<https://doi.org/10.5067/MODIS/MOD44B.061/>}
#'   \item{precip, temp}{<https://chelsa-climate.org/exchelsa-extended-bioclim/>}
#'   \item{precip2100}{Gutiérrez et al., 2021; Iturbide et al., 2021}
#' }
"treeCover"
