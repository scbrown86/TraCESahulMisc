#' Set time index for TraCE-Sahul SpatRaster
#'
#' @description
#' This function sets the \code{\link[terra:time]{time index}} of a SpatRaster
#' representing the full (or a subset) from the \emph{TraCE-Sahul} decadal
#' dataset. Internally the function reads a CSV of timesteps from
#' \code{inst/extdata/TraCE-Sahul_timesteps.csv}.
#'
#' @details
#' The function will stop if \code{x} is not a SpatRaster.
#'
#' The function \strong{should} work with the subsets of the decadal data
#' e.g \emph{TraCE_22ka_downscaled_pr_decadal_21k_1500CE_biascorr_000001.nc}
#' provided that they haven't been renamed!
#'
#' @note
#' There are no checks done to ensure that the dataset passed is decadal data
#' from \emph{TraCE-Sahul}.
#'
#' @param x A \code{\link[terra:rast]{SpatRast}} object from the \emph{TraCE-Sahul} dataset
#' @return A \emph{SpatRaster} with updated time index
#'
#' @examples
#' \dontrun{
#' library(terra)
#' climData <- rast("TraCE_22ka_downscaled_pr_decadal_21k_1500CE_biascorr_000001.nc")
#' climData <- setTime(climData)
#' climData
#' }
#' @export
setTime <- function(x) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for setTime()")
  }
  stopifnot("`x` must be a `terra::SpatRaster`" = inherits(x, "SpatRaster"))
  # get the number of layers
  nt <- terra::nlyr(x)
  # read in the data
  path <- system.file("extdata", "TraCE-Sahul_timesteps.csv", package = "TraCESahulMisc")
  data <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (nt == 25860) { # if full file, then replace with full time index
    time_index <- data$DecYear
    stopifnot(length(unique(time_index)) ==  length(time_index))
    stopifnot(nt == length(time_index))
    # Overwrite existing time variable
    terra::time(x, "yearmonths") <- time_index
    return(x)
  } else { # else use the dates for the step only
    # grab the "step" from the filename
    step <- as.integer(gsub(".nc", "", sapply(strsplit(basename(x@pntr$filenames()), "_"), utils::tail, 1)))
    time_index <- data[data$file_step == step, "DecYear"]
    terra::time(x, "yearmonths") <- time_index
    return(x)
  }
}
