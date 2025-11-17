#' Classify the input files based on whats provided
#'
#' @description
#' This function classifies the input files into a number of different categories
#' depending on whats passed which is then used to switch the function used for
#' importing the data. This function is not exported, and should not ever need to
#' be called directly by the user.
#'
#' @noRd
classify_TraCESahul <- function(files) {
  stopifnot(length(files) > 0)
  bfiles <- basename(files)
  vars <- sub("^TraCE_22ka_downscaled_([^_]+)_.*$", "\\1", bfiles)
  stopifnot("Only a single variable can be processed at a time. Check your input files." = length(unique(vars)) == 1)
  is_1500_1990 <- endsWith(bfiles, "1500_1990_biascorr.nc")
  chunk_pattern <- "_biascorr_([0-9]{2})\\.nc$"
  chunk_extract <- function(x) {
    m <- regexec(chunk_pattern, x)
    mm <- regmatches(x, m)
    if (length(mm[[1]]) == 2) as.integer(mm[[1]][2]) else NA_integer_
  }
  chunk_nums <- vapply(bfiles, chunk_extract, integer(1))
  is_chunk <- !is.na(chunk_nums)
  n <- length(files)
  if (n == 1 && is_1500_1990) {
    return(list(case = "1500_1990_single", var = unique(vars), n = 1))
  }
  if (n == 1 && (endsWith(bfiles, "21k_1500CE_biascorr.nc") || is_chunk)) {
    return(list(case = "decadal_single", var = unique(vars), n = 1, chunk_num = chunk_nums))
  }
  if (n == 6 && all(is_chunk) && identical(unname(chunk_nums), 1:6)) {
    return(list(case = "decadal_chunks_6", var = unique(vars), n = 6))
  }
  if (n == 7) {
    first6 <- bfiles[1:6]
    last1  <- bfiles[7]
    nums6 <- vapply(first6, chunk_extract, integer(1))
    if (identical(unname(nums6), 1:6) &&
        endsWith(last1, "1500_1990_biascorr.nc"))
      return(list(case = "chunks_plus_1500_1990", var = unique(vars), n = 7))
  }
  if (n == 2 &&
      (endsWith(bfiles[1], "21k_1500CE_biascorr.nc") || is_chunk[1]) &&
      is_1500_1990[2]) {
    return(list(case = "decadal_single_plus_1500_1990", var = unique(vars), n = 2, chunk_num = chunk_nums[1]))
  }
  stop("Filename set does not match any expected ordered pattern. Have you changed the filenames?")
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
#' e.g \emph{TraCE_22ka_downscaled_pr_decadal_21k_1500CE_biascorr_01.nc}
#' provided that they haven't been renamed!
#'
#' @note
#' The function will fail if the \emph{TraCE-Sahul} datasets have been renamed.
#'
#' filepaths must be given in order (e.g. 01.nc, 02.nc, 03.nc, 04.nc, 05.nc, 06.nc, 1500_1990_biascorr.nc)
#'
#' @import ncdfCF
#' @import data.table
#' @importFrom terra rast time depth depthName depthUnit
#' @importFrom pbapply pblapply
#'
#' @param files A character vector of filepaths pointing to the \emph{TraCE-Sahul} dataset(s).
#' @return A \emph{SpatRaster} with updated \code{\link[terra]{time}} index showing the relevant year, and \code{\link[terra]{depth}} representing the month
#'
#' A \emph{data.table} in the \code{.GlobalEnv} named \code{TraCESahul_time_steps} with the relevant details for all timesteps for \code{files}
#'
#' @examples
#' \dontrun{
#' fn <- c("TraCE_22ka_downscaled_pr_decadal_21k_1500CE_biascorr_06.nc",
#'         "TraCE_22ka_downscaled_pr_1500_1990_biascorr.nc")
#'         sahul_pr <- import_TraCESahul(fn)
#' sahul_pr
#' fn <- c("TraCE_22ka_downscaled_tasmax_decadal_21k_1500CE_biascorr_06.nc",
#'         "TraCE_22ka_downscaled_tasmax_1500_1990_biascorr.nc")
#' sahul_tasmax <- import_TraCESahul(fn)
#' sahul_tasmax
#' fn <- c("TraCE_22ka_downscaled_tasmin_decadal_21k_1500CE_biascorr_06.nc",
#'         "TraCE_22ka_downscaled_tasmin_1500_1990_biascorr.nc")
#' sahul_tasmin <- import_TraCESahul(fn)
#' sahul_tasmin
#' }
#' @export
#'
import_TraCESahul <- function(files) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is required.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
  out <- classify_TraCESahul(files)
  set_monthly_metadata <- function(ds, years) {
    terra::depth(ds) <- rep(1:12, times = terra::nlyr(ds) / 12)
    terra::depthName(ds) <- "Month"
    terra::depthUnit(ds) <- "calendar month"
    terra::time(ds) <- years
    attr(ds, "TraCESahul") <- TRUE
    ds
  }
  load_timestep_table <- function() {
    path <- system.file("extdata", "TraCE-Sahul_timesteps.csv", package = "TraCESahulMisc")
    data.table::fread(path)
  }
  get_monthly_years <- function() rep(1500:1989, each = 12)
  handle_decadal_core <- function(f) {
    data <- load_timestep_table()
    ds   <- terra::rast(f)
    ds   <- set_monthly_metadata(ds, data[["Year"]])
    assign("TraCESahul_time_steps", data, envir = .GlobalEnv)
    ds
  }
  handle_decadal_single <- function(f,...) {
    data <- load_timestep_table()
    chunk <- unname(out$chunk_num)
    if (!is.na(chunk)) {
      data <- data[data$file_step == chunk, ]
    }
    ds <- terra::rast(f)
    ds <- set_monthly_metadata(ds, data$Year)
    assign("TraCESahul_time_steps", data, envir = .GlobalEnv)
    ds
  }
  handle_monthly_only <- function(f,...) {
    ds <- terra::rast(f)
    set_monthly_metadata(ds, get_monthly_years())
  }
  handle_two_sets <- function(f_decadal, f_monthly,...) {
    data <- load_timestep_table()
    chunk <- unname(out$chunk_num)
    data <- data[data$file_step == chunk, ]
    ds1  <- set_monthly_metadata(terra::rast(f_decadal), data$Year)
    ds2  <- set_monthly_metadata(terra::rast(f_monthly), get_monthly_years())
    assign("TraCESahul_time_steps", data, envir = .GlobalEnv)
    c(ds1, ds2)
  }
  switch(out$case,
         "decadal_chunks_6" = handle_decadal_core(files),
         "decadal_single_plus_1500_1990" = handle_two_sets(files[1], files[2]),
         "1500_1990_single" = handle_monthly_only(files),
         "decadal_single" = handle_decadal_single(files),
         "chunks_plus_1500_1990" = handle_two_sets(files[-length(files)], files[length(files)]),
    stop("Unhandled case.")
  )
}


