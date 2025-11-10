#' Classify the input files based on whats provided
#'
#' @description
#' This function classifies the input files into a number of different categories
#' depending on whats passed which is then used to switch the function used for
#' importing the data. This function is not exported, and should not ever need to
#' be called directly by the user.
#'
classify_TraCESahul <- function(files) {
  stopifnot(length(files) > 0)
  bfiles <- basename(files)
  # Extract variable name
  vars <- sub("^TraCE_22ka_downscaled_([^_]+)_.*$", "\\1", bfiles)
  stopifnot("Only a single variable can be processed at a time. Check your input files." = length(unique(vars)) == 1)
  is_1500_1990 <- endsWith(bfiles, "1500_1990_biascorr.nc")
  is_decadal_single <- grepl("21k_1500CE_biascorr\\.nc$|_biascorr_[0-9]{2}\\.nc$", bfiles)
  # Chunk pattern (01–06)
  chunk_pattern <- "_biascorr_([0-9]{2})\\.nc$"
  chunk_extract <- function(x) {
    m <- regexec(chunk_pattern, x)
    mm <- regmatches(x, m)
    if (length(mm[[1]]) == 2) as.integer(mm[[1]][2]) else NA_integer_
  }
  chunk_nums <- vapply(bfiles, chunk_extract, integer(1))
  is_chunk <- !is.na(chunk_nums)
  n <- length(files)
  # Monthly data from 1500–1990
  if (n == 1 && is_1500_1990) {
    return(list(case = "1500_1990_single", var = unique(vars), n = 1))
  }
  # Single decadal chunk
  if (n == 1 && is_decadal_single) {
    return(list(case = "decadal_single", var = unique(vars), n = 1, chunk_num = chunk_nums))
  }
  # Six chunked files in order
  if (n == 6 &&
      all(is_chunk) &&
      all.equal(unname(chunk_nums), 1:6)) {
    return(list(case = "decadal_chunks_6", var = unique(vars), n = 6))
  }
  # Seven-file case: six decadal chunks + monhthly data
  if (n == 7) {
    # first six must be chunks 1..6 *in order*
    first6 <- bfiles[1:6]
    last1  <- bfiles[7]
    chunk_nums_first6 <- vapply(first6, chunk_extract, integer(1))
    if (all.equal(unname(chunk_nums_first6), 1:6) &&
        ends_with(last1, "1500_1990_biascorr.nc")) {
      return(list(case = "chunks_plus_1500_1990", var = unique(vars), n = 7))
    }
  }
  # Two-files: full decadal file then monthly data
  if (n == 2 &&
      is_decadal_single[1] &&
      is_1500_1990[2]) {
    return(list(case = "decadal_single_plus_1500_1990", var = unique(vars), n = 2))
  }
  stop("Filename set does not match any expected ordered pattern. Have the filenames been altered?")
}

#' Imports the TraCE-Sahul dataset
#'
#' @description
#' This function returns a \code{\link[terra:rast]{SpatRast}} representing the
#' full (or a subset) \emph{TraCE-Sahul} dataset. Internally the
#' function reads a CSV of timesteps from \code{inst/extdata/TraCE-Sahul_timesteps.csv},
#' and exports it to the \code{.GlobalEnv}. This data.table is \emph{only} relevant
#' to the decadal data.
#'
#' If the monthly data from 1500-1989 is also provided, it will be appended.
#'
#' @details
#' The function will stop if \code{x} is not a file.path or series of file.paths.
#'
#' The function will work with the subsets of the decadal data
#' e.g \emph{TraCE_22ka_downscaled_pr_decadal_21k_1500CE_biascorr_000001.nc}
#' provided that they haven't been renamed!
#'
#' @note
#' The function will fail if the \emph{TraCE-Sahul} datasets have been renamed.
#'
#' filepaths must be given in order (e.g. 01.nc, 02.nc, 03.nc, 04.nc, 05.nc, 06.nc, 1500_1990_biascorr.nc)
#'
#' @import ncdfCF
#' @import terra
#' @importFrom pbapply pblapply
#'
#' @param files A character vector of filepaths pointing to the \emph{TraCE-Sahul} dataset(s).
#' @return A \emph{SpatRaster} with updated \code{\link[terra]{time}} index showing the relevant year, and \code{\link[terra]{depth}} representing the month
#'
#' A \emph{data.table} in the \code{.GlobalEnv} named \code{TraCESahul_time_steps} with the relevant details for all timesteps for \code{files}
#'
#' @examples
#' \dontrun{
#' fn <- "TraCE_22ka_downscaled_pr_decadal_21k_1500CE_biascorr_000001.nc"
#' pr_chunk01 <- import_TraCESahul(fn)
#' pr_chunk01
#' }
#' @export
#'
import_TraCESahul <- function(files) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for import_TraCESahul()")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for import_TraCESahul()")
  }
  # check the files
  out <- classify_TraCESahul(files)

  handle_decadal_chunks <- function(files) {
    path <- system.file("extdata", "TraCE-Sahul_timesteps.csv",
                        package = "TraCESahulMisc")
    data <- data.table::fread(path, stringsAsFactors = FALSE)
    ds <- terra::rast(files)
    terra::depth(ds) <- rep(1:12, times = terra::nlyr(ds)/12)
    terra::depthName(ds) <- "Month"
    terra::depthUnit(ds) <- "calendar month"
    terra::time(ds) <- data[["Year"]]
    assign("TraCESahul_time_steps", data, envir = .GlobalEnv)
    return(ds)
  }
  handle_two_file_combo <- function(files) {
    path <- system.file("extdata", "TraCE-Sahul_timesteps.csv",
                        package = "TraCESahulMisc")
    data <- data.table::fread(path, stringsAsFactors = FALSE)
    ds <- terra::rast(files[1])
    terra::depth(ds) <- rep(1:12, times = terra::nlyr(ds)/12)
    terra::depthName(ds) <- "Month"
    terra::depthUnit(ds) <- "calendar month"
    terra::time(ds) <- data[["Year"]]
    ds_m <- terra::rast(files[2])
    terra::depth(ds_m) <- rep(1:12, times = terra::nlyr(ds_m)/12)
    terra::depthName(ds_m) <- "Month"
    terra::depthUnit(ds_m) <- "calendar month"
    terra::time(ds_m) <- rep(1500:1989, each = 12)
    ds_c <- c(ds, ds_m)
    assign("TraCESahul_time_steps", data, envir = .GlobalEnv)
    return(ds_c)
  }
  handle_1500_1990 <- function(files) {
    ds <- terra::rast(files)
    terra::depth(ds) <- rep(1:12, times = terra::nlyr(ds)/12)
    terra::depthName(ds) <- "Month"
    terra::depthUnit(ds) <- "calendar month"
    terra::time(ds) <- rep(1500:1989, each = 12)
    return(ds)
  }
  handle_single_decadal <- function(files,...) {
    path <- system.file("extdata", "TraCE-Sahul_timesteps.csv",
                        package = "TraCESahulMisc")
    data <- data.table::fread(path, stringsAsFactors = FALSE)
    chunk <- unname(out$chunk_num)
    data <- data.table::copy(data)[file_step == chunk, ]
    ds <- terra::rast(files)
    terra::depth(ds) <- rep(1:12, times = terra::nlyr(ds)/12)
    terra::depthName(ds) <- "Month"
    terra::depthUnit(ds) <- "calendar month"
    terra::time(ds) <- data[["Year"]]
    assign("TraCESahul_time_steps", data, envir = .GlobalEnv)
    return(ds)
  }
  handle_chunks_plus <- function(files) {
    path <- system.file("extdata", "TraCE-Sahul_timesteps.csv",
                        package = "TraCESahulMisc")
    data <- data.table::fread(path, stringsAsFactors = FALSE)
    ds <- terra::rast(files[-length(files)])
    terra::depth(ds) <- rep(1:12, times = terra::nlyr(ds)/12)
    terra::depthName(ds) <- "Month"
    terra::depthUnit(ds) <- "calendar month"
    terra::time(ds) <- data[["Year"]]
    ds_m <- terra::rast(files[length(files)])
    terra::depth(ds_m) <- rep(1:12, times = terra::nlyr(ds_m)/12)
    terra::depthName(ds_m) <- "Month"
    terra::depthUnit(ds_m) <- "calendar month"
    terra::time(ds_m) <- rep(1500:1989, each = 12)
    ds_c <- c(ds, ds_m)
    assign("TraCESahul_time_steps", data, envir = .GlobalEnv)
    return(ds_c)
  }
  # processing dispatch
  switch(
    out$case,
    "decadal_chunks_6" = handle_decadal_chunks(files),
    "decadal_single_plus_1500_1990" = handle_two_file_combo(files),
    "1500_1990_single" = handle_1500_1990(files),
    "decadal_single" = handle_single_decadal(files),
    "chunks_plus_1500_1990" = handle_chunks_plus(files),
    stop("Unhandled case")
  )
}

