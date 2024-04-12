#' Violent crime data in United States
#'
#' A subset of data of violent crime per 100,000 people in western states from 1960 to 2019. 
#'
#' @format ## `Violent`
#' A data frame with 60 rows and 11 columns:
#' \describe{
#'   \item{AZ}{Arizona}
#'   \item{CA}{California}
#'   \item{CO}{Colorado}
#'   \item{ID}{Idaho}
#'   \item{MT}{Montana}
#'   \item{NV}{Nevada}
#'   \item{NM}{New Mexico}
#'   \item{OR}{Oregon}
#'   \item{UT}{Utah}
#'   \item{WA}{Washington}
#'   \item{WY}{Wyoming}
#'   ...
#' }
#' @source <https://www.disastercenter.com/crime/>
"Violent"


#' Property crime in United States
#'
#' A subset of data of property crime per 100,000 people in western states from 1960 to 2019. 
#'
#' @format ## `Property`
#' A data frame with 60 rows and 11 columns:
#' \describe{
#'   \item{AZ}{Arizona}
#'   \item{CA}{California}
#'   \item{CO}{Colorado}
#'   \item{ID}{Idaho}
#'   \item{MT}{Montana}
#'   \item{NV}{Nevada}
#'   \item{NM}{New Mexico}
#'   \item{OR}{Oregon}
#'   \item{UT}{Utah}
#'   \item{WA}{Washington}
#'   \item{WY}{Wyoming}
#'   ...
#' }
#' @source <https://www.disastercenter.com/crime/>
"Property"


#' Westen states in United States
#'
#' A sp map data of the western states in United States
#'
#' @format ## `WestStates`
#' A SpatialPolygonsDataFrame data of the western states in United States
#' \describe{
#'   \item{FID}{The number ID of the western states}
#'   \item{State_Code}{Abbreviations of the state names}
#'   \item{State_Name}{Names of the states}
#'   A SpatialPolygonsDataFrame data of the western states in United States
#' }
#' @source <https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html>
"WestStates"

