#' Example fossil dataset
#'
#' A \emph{PackedSpatVector} object with locations and ages of fossils
#' for a virtual species. Use \code{\link[terra]{wrap}} to convert back to a
#' \emph{SpatVector}.
#'
#' @format A \code{\link[terra:vect]{SpatVect}} with columns:
#' \describe{
#'   \item{date_estimate}{ Age estimate of fossil record}
#'   \item{ci_lower, ci_upper}{ Lower and upper CI of fossil age}
#'   \item{occ}{ Binary identifier of whether an occurence orf background point}
#' }
#' @source data-raw/ex_foss.R
#' @usage data(ex_foss)
#' @format A \code{\link[terra:vect]{PackedSpatVector}}
#' @examples
#' \dontrun{
#' library(terra)
#' ex_foss
#' ex_foss <- terra::unwrap(ex_foss)
#' ex_foss
#' }
"ex_foss"

#' Example habitat suitability
#'
#' A \emph{PackedSpatRaster} object with "true" contemporary habitat suitability
#' for a virtual species. Use \code{\link[terra]{wrap}} to convert back to a
#' \code{\link[terra:rast]{SpatRast}}.
#'
#' @source data-raw/ex_foss.R
#' @usage data(true_suit)
#' @format A \code{\link[terra:rast]{PackedSpatRaster}}
#' @examples
#' \dontrun{
#' library(terra)
#' true_suit
#' true_suit <- terra::unwrap(true_suit)
#' true_suit
#' }
"true_suit"

#' Timesteps for TraCE-Sahul
#'
#' A data.table object containing the timesteps for each layer in the
#' TraCE-Sahul downscaled climate data.
#'
#' @source data-raw/ex_foss.R
#' @usage data(TraCESahul_timesteps)
#' @format A \code{\link[data.table:data.table]{DataTable}} of 31,740 rows and
#' 11 columns
#' @examples
#' \dontrun{
#' library(data.table)
#' TraCESahul_timesteps
#' }
"TraCESahul_timesteps"
