#' Set time index for TraCE-Sahul SpatRaster
#'
#' @description
#' This function sets the \code{\link[terra:time]{time index}} of a SpatRaster
#' representing the full (or a subset) from the \emph{TraCE-Sahul} decadal
#' dataset. Internally the function reads a CSV of timesteps from
#' \code{inst/extdata/TraCE-Sahul_timesteps.csv}.
#'
#' @details
#' The function will stop if \code{x} is not a file.path.
#'
#' The function \strong{should} will with the subsets of the decadal data
#' e.g \emph{TraCE_22ka_downscaled_pr_decadal_21k_1500CE_biascorr_000001.nc}
#' provided that they haven't been renamed!
#'
#' @import terra
#' @import ncdfCF
#'
#' @note
#' There are no checks done to ensure that the dataset passed is decadal data
#' from \emph{TraCE-Sahul}.
#'
#' @param x A filepath pointing to the \emph{TraCE-Sahul} dataset(s)
#' @return A \emph{SpatRaster} with updated time index showing the relevant year
#'
#' A \emph{data.table} in the \code{.GlobalEnv} named \code{TraCESahul_time_steps} with the relevant details for all timesteps for \code{x}
#'
#' @examples
#' \dontrun{
#' library(terra)
#' fn <- "TraCE_22ka_downscaled_pr_decadal_21k_1500CE_biascorr_000001.nc"
#' climData <- setTime(fn)
#' climData
#' }
#'
setTime <- function(x) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for setTime()")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for setTime()")
  }
  if (!requireNamespace("ncdfCF", quietly = TRUE)) {
    stop("Package 'ncdfDF' is required for setTime()")
  }
  stopifnot("`x` must be a valid `file.path` to the TraCE-Sahul data" = file.exists(x))
  # get the number of layers
  nt <- terra::nlyr(x)
  # read in the data
  path <- system.file("extdata", "TraCE-Sahul_timesteps.csv", package = "TraCESahulMisc")
  data <- data.table::fread(path, stringsAsFactors = FALSE)
  # grab time index based on x
  if (nt == 25860) {
    time_index <- data[["DecYear"]]
    stopifnot(length(unique(time_index)) == length(time_index))
    stopifnot(nt == length(time_index))
    assign("TraCESahul_time_steps", data, envir = .GlobalEnv)
  } else {
    step <- as.integer(gsub(".nc", "", sapply(strsplit(basename(x@pntr$filenames()), "_"), utils::tail, 1)))
    time_index <- data[file_step == step, ][["DecYear"]]
    assign("TraCESahul_time_steps", data[file_step == step, ], envir = .GlobalEnv)
  }
  # assign the relevant time steps to global environment
  #assign("TraCESahul_time_index", time_index, envir = .GlobalEnv)
  # setTime on raster
  terra::time(x, tstep = "") <- time_index
  return(x)
}
