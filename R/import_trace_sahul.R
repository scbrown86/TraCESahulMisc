#' Classify the input files based on whats provided
#'
#' @description
#' This function classifies the input files into a number of different categories
#' depending on whats passed which is then used to switch the function used for
#' importing the data. This function is not exported, and should not ever need to
#' be called directly by the user.
#'
#' @keywords internal
#'
#' @noRd
classify_TraCESahul <- function(files) {
  stopifnot(length(files) > 0)
  bfiles <- basename(files)
  pat_chunk   <- "^TraCE-Sahul_decadal_22k_1500CE_([a-z]+)_([0-9]{2})\\.nc$"
  pat_single  <- "^TraCE-Sahul_decadal_22k_1500CE_([a-z]+)\\.nc$"
  pat_annual  <- "^TraCE-Sahul_annual_1500_1990_([a-z]+)\\.nc$"
  m_chunk  <- regexec(pat_chunk, bfiles)
  r_chunk  <- regmatches(bfiles, m_chunk)
  m_single <- regexec(pat_single, bfiles)
  r_single <- regmatches(bfiles, m_single)
  m_annual <- regexec(pat_annual, bfiles)
  r_annual <- regmatches(bfiles, m_annual)
  is_chunk  <- lengths(r_chunk) == 3
  is_single <- lengths(r_single) == 2
  is_annual <- lengths(r_annual) == 2
  if (!all(is_chunk | is_single | is_annual)) {
    stop("Filenames do not match expected TraCE-Sahul naming convention.")
  }
  vars <- character(length(bfiles))
  vars[is_chunk]  <- vapply(r_chunk[is_chunk],  `[`, character(1), 2)
  vars[is_single] <- vapply(r_single[is_single], `[`, character(1), 2)
  vars[is_annual] <- vapply(r_annual[is_annual], `[`, character(1), 2)
  if (length(unique(vars)) != 1) {
    stop("All files must use the same variable.")
  }
  chunk_nums <- rep(NA_integer_, length(bfiles))
  chunk_nums[is_chunk] <- as.integer(vapply(r_chunk[is_chunk],
                                            `[`, character(1), 3))
  n <- length(files)
  if (n == 1 && is_annual) {
    return(list(case = "1500_1990_single", var = vars[1], n = 1))
  }
  if (n == 1 && (is_single || is_chunk)) {
    return(list(case = "decadal_single", var = vars[1], n = 1, chunk_num = chunk_nums))
  }
  if (n == 6 && all(is_chunk) && identical(chunk_nums, 1:6)) {
    return(list(case = "decadal_chunks_6", var = vars[1], n = 6))
  }
  if (n == 7 && all(is_chunk[1:6]) && is_annual[7] && identical(chunk_nums[1:6], 1:6)) {
    return(list(case = "chunks_plus_1500_1990", var = vars[1], n = 7))
  }
  if (n == 2 && (is_single[1] || is_chunk[1]) && is_annual[2]) {
    return(list(case = "decadal_single_plus_1500_1990", var = vars[1],
                n = 2, chunk_num = chunk_nums[1]))
  }
  stop("Filename set is valid but not a recognised combination or order.")
}

#' Imports the TraCE-Sahul dataset
#'
#' @description
#' This function returns a \code{\link[terra:rast]{SpatRast}} representing the
#' full (or a subset) \emph{TraCE-Sahul} dataset. Internally the
#' function reads a CSV of timesteps from \code{inst/extdata/TraCE-Sahul_timesteps.csv},
#' and exports it to the \code{.GlobalEnv}. This data.table contains all the timesteps
#' for every layer in the TraCE-Sahul dataset.
#'
#' If the monthly data from 1500-1989 is also provided, it will be appended.
#'
#' @details
#' The function will stop if \code{x} is not a file.path or series of file.paths.
#'
#' The function will work with the subsets of the decadal data
#' e.g \emph{TraCE-Sahul_decadal_22k_1500CE_pr_01.nc} provided that they haven't
#' been renamed!
#'
#' @note
#' The function will fail if the \emph{TraCE-Sahul} datasets have been renamed.
#'
#' filepaths must be given in order (e.g. 01.nc, 02.nc, 03.nc, 04.nc, 05.nc, 06.nc),
#' with the monthly data from 1500 to 1989 provided last
#'
#' The file can be read in directly with \code{terra::rast}, but note that due to the way
#' \pkg{terra} handles the time index the dates will appear incorrect.
#' See Example 5 using the test dataset below.
#'
#' @import data.table
#' @importFrom terra rast time depth depthName depthUnit
#' @importFrom pbapply pblapply
#'
#' @param files A character vector of filepaths pointing to the
#' \emph{TraCE-Sahul} dataset(s).
#' @param aoi an \code{\link[terra]{ext}} object or four parameter vector
#' (\code{c(xmin,xmax,ymin,ymax)}) passed to \code{\link[terra]{rast}} to
#' spatially subset the \emph{TraCE-Sahul} data. Default = \code{NULL} which will
#' return the entire spatial domain.
#'
#' @return A \emph{SpatRaster} with updated \code{\link[terra]{time}} index
#' showing the relevant year, and \code{\link[terra]{depth}} representing the month
#'
#' A \emph{data.table} in the \code{.GlobalEnv} named \code{TraCESahul_time_steps} with the relevant details for all timesteps for \code{files}
#'
#' @examples
#' \dontrun{
#' base_dir <- "~/Documents/TraCE-Sahul"
#' vars <- c("pr", "tasmax", "tasmin")
#'
#' # Example 1: full decadal data (no chunk) + annual, all variables
#' load_full <- function(base, var) {
#'   f1 <- file.path(base, var,
#'                   sprintf("TraCE-Sahul_decadal_22k_1500CE_%s.nc", var))
#'   f2 <- file.path(base, var,
#'                   sprintf("TraCE-Sahul_annual_1500_1990_%s.nc", var))
#'   import_TraCESahul(files = c(f1, f2))
#' }
#' sahul_full <- lapply(vars, function(v) load_full(base_dir, v))
#' names(sahul_full) <- vars
#' sahul_full
#'
#' # Example 2: chunk 6 only + annual, all variables
#' load_chunk6 <- function(base, var) {
#'   f1 <- file.path(base, var,
#'                   sprintf("TraCE-Sahul_decadal_22k_1500CE_%s_06.nc", var))
#'   f2 <- file.path(base, var,
#'                   sprintf("TraCE-Sahul_annual_1500_1990_%s.nc", var))
#'   import_TraCESahul(files = c(f1, f2))
#' }
#' sahul_chunk6 <- lapply(vars, function(v) load_chunk6(base_dir, v))
#' names(sahul_chunk6) <- vars
#' sahul_chunk6
#'
#' # Example 3: chunk 3 only, all variables
#' load_chunk3 <- function(base, var) {
#'   f <- file.path(base, var,
#'                  sprintf("TraCE-Sahul_decadal_22k_1500CE_%s_03.nc", var))
#'   import_TraCESahul(files = f)
#' }
#' sahul_chunk3 <- lapply(vars, function(v) load_chunk3(base_dir, v))
#' names(sahul_chunk3) <- vars
#' sahul_chunk3
#'
#' # Example 4: monthly only, all variables
#' load_monthly <- function(base, var) {
#'   f <- file.path(base, var,
#'                  sprintf("TraCE-Sahul_annual_1500_1990_%s.nc", var))
#'  import_TraCESahul(files = f)
#'  }
#' sahul_monthly <- lapply(vars, function(v) load_monthly(base_dir, v))
#' names(sahul_monthly) <- vars
#' sahul_monthly
#'
#' # Example 5: Read in directly with terra
#' pr_sahul <- terra::rast("~/Documents/TraCE-Sahul/pr/TraCE-Sahul_annual_1500_1990_pr.nc")
#' pr_sahul
#' # class       : SpatRaster
#' # size        : 100, 200, 5880  (nrow, ncol, nlyr)
#' # resolution  : 0.05, 0.05  (x, y)
#' # extent      : 125, 135, -13, -8  (xmin, xmax, ymin, ymax)
#' # coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84)
#' # source      : TraCE-Sahul_annual_1500_1990_pr.nc
#' # varname     : pr (Precipitation)
#' # names       :     pr_1,     pr_2,     pr_3,     pr_4,     pr_5,     pr_6,      ...
#' # unit        : mm/month
#' # time (ymnts): 1970-Jan to 2459-Dec (5880 steps)
#' }
#' @export
#'
import_TraCESahul <- function(files, aoi = NULL) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is required.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
  out <- classify_TraCESahul(files)
  set_monthly_metadata <- function(ds, years) {
    terra::depth(ds) <- rep(1:12, times = terra::nlyr(ds) / 12)
    terra::depthName(ds) <- "Month"
    terra::depthUnit(ds) <- "calendar month"
    terra::time(ds) <- years
    terra::crs(ds) <- "EPSG:4326"
    attr(ds, "TraCESahul") <- TRUE
    ds
  }
  load_timestep_table <- function() {
    path <- system.file("extdata", "TraCE-Sahul_timesteps.csv", package = "TraCESahulMisc")
    data.table::fread(path)
  }
  get_monthly_years <- function() rep(1500:1989, each = 12)
  handle_decadal_core <- function(f, aoi) {
    data <- load_timestep_table()
    ds   <- terra::rast(f, win = aoi)
    terra::crs(ds) <- "EPSG:4326"
    ds   <- set_monthly_metadata(ds, data[["YearsCE"]])
    assign("TraCESahul_time_steps", data, envir = .GlobalEnv)
    ds
  }
  handle_decadal_single <- function(f, aoi) {
    data <- load_timestep_table()
    chunk <- unname(out$chunk_num)
    if (!is.na(chunk)) {
      data <- data.table::copy(data)[data$file_step == chunk, ]
    }
    ds <- terra::rast(f, win = aoi)
    ds <- set_monthly_metadata(ds, data[["YearsCE"]])
    assign("TraCESahul_time_steps", data, envir = .GlobalEnv)
    ds
  }
  handle_monthly_only <- function(f, aoi) {
    ds <- terra::rast(f, win = aoi)
    set_monthly_metadata(ds, get_monthly_years())
  }
  handle_two_sets <- function(f_decadal, f_monthly,aoi) {
    data <- load_timestep_table()
    chunk <- unname(out$chunk_num)
    if (is.null(chunk) || is.na(chunk)) {
      data <- data.table::copy(data)[!is.na(data$file_step), ]
    } else if (!is.na(chunk)) {
      data <- data.table::copy(data)[data$file_step == chunk, ]
    }
    ds1  <- set_monthly_metadata(terra::rast(f_decadal, win = aoi), data$YearsCE)
    ds2  <- set_monthly_metadata(terra::rast(f_monthly, win = aoi), get_monthly_years())
    assign("TraCESahul_time_steps", data, envir = .GlobalEnv)
    c(ds1, ds2)
  }
  switch(out$case,
         "decadal_chunks_6" = handle_decadal_core(files, aoi),
         "decadal_single_plus_1500_1990" = handle_two_sets(files[1], files[2], aoi),
         "1500_1990_single" = handle_monthly_only(files, aoi),
         "decadal_single" = handle_decadal_single(files, aoi),
         "chunks_plus_1500_1990" = handle_two_sets(files[-length(files)],
                                                   files[length(files)],
                                                   aoi),
    stop("Unhandled case.")
  )
}


