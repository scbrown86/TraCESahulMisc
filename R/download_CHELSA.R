#' Source daily CHELSA climate rasters
#'
#' @description
#' Downloads a daily CHELSA v2.1 GeoTIFF for a given date and variable via GDAL's
#' virtual file system (\href{https://gdal.org/en/stable/user/virtual_file_systems.html#network-based-file-systems}{\emph{vsicurl}}),
#' reads only a buffered window around \code{template}, resamples to the template grid, optionally converts units, and
#' writes the result to disk.
#'
#' @details
#'
#' The crop window is hard coded to the extent of \code{template}, buffered by 50 km.
#' Output files are written to \code{file.path(dir, var)} and named using the source
#' layer name (for example, \code{"CHELSA_pr_01_01_1980_V.2.1.tif"}).
#'
#' Note that the final raster is resampled to match the exact extent of the template raster,
#' but data is read/downloaded for the extent of the template + 50km.
#'
#' If \code{convert = TRUE}:
#'
#' \code{tasmin} and \code{tasmax} are converted from \emph{Kelvin} to \emph{Celsius}, \code{pr} is
#' converted from \emph{mm/day} to \emph{kg m-2 s-1}.
#'
#' @importFrom terra rast ext buffer vect crs resample setValues values units varnames writeRaster
#'
#' @param x Date. Must be a \code{\link[base]{Date}} class object/vector
#' @param var Character. CHELSA variable name. Must be one of \code{pr}, \code{tasmax}, or \code{tasmin}.
#' @param dir Character. Base output directory. A subdirectory named \code{var} is
#'   created if needed.
#' @param template \code{\link[terra:rast]{SpatRast}} grid for resampling and cropping. Should be a \emph{TraCE-Sahul} grid.
#' @param algo Character. Resampling method passed to \code{\link[terra]{resample}}.
#' @param mask Logical. If \code{TRUE} mask to the non-NA values of the template raster
#' @param convert Logical. If \code{TRUE}, apply variable-specific unit conversion.
#' @param ... Additional arguments passed to \code{\link[terra]{writeRaster}}.
#'
#' @return \code{file.path} to the cropped and resampled raster
#'
#' @examples
#' \dontrun{
#' library(terra)
#'
#' template <- rast(
#'   ext(105, 161.25, -50, 10),
#'   resolution = 0.05,
#'   crs = "EPSG:4326"
#' )
#'
#' r <- download_CHELSA(
#'   x = as.Date("1980-01-01"),
#'   var = "pr",
#'   dir = tempdir(),
#'   template = template,
#'   algo = "cubicspline",
#'   convert = TRUE,
#'   overwrite = TRUE
#' )
#'
#' r
#' }
#' @export
download_CHELSA <- function(x, var, dir, template,
                            algo = "cubicspline", mask = TRUE,
                            convert = TRUE, ...) {
  if (!inherits(x, "Date")) {
    stop("'x' must be a Date")
  }
  match.arg(var, c("pr", "tasmax", "tasmin"), several.ok = FALSE)
  if (!inherits(template, "SpatRaster")) {
    stop("'template' must be a SpatRaster")
  }
  out_dir <- file.path(dir, var)
  if (!dir.exists(out_dir)) {
    message(sprintf("Creating download directory: %s", out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)
  }
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  base_url <- "/vsicurl/https://os.unil.cloud.switch.ch/chelsa02/"
  model <- "chelsa"
  extent <- "global"
  url <- paste0(base_url, model, "/", extent,
                "/daily/", var, "/", format(x, "%Y"), "/")
  file <- sprintf("CHELSA_%s_%s_%s_%s_V.2.1.tif",
                  var, format(x, "%d"), format(x, "%m"), format(x, "%Y"))
  src <- paste0(url, file)
  dl_ext <- terra::ext(
    terra::buffer(
      terra::vect(terra::ext(template), crs = terra::crs(template)),
      width = 50 * 1000))
  r <- terra::rast(src, win = dl_ext)
  r_res <- terra::resample(r, template, method = algo)
  if (mask) {
    r_res <- terra::mask(r_res, template)
  }
  if (convert) {
    if (var %in% c("tasmax", "tasmin")) {
      r_res <- terra::setValues(x = r_res,
                                values = round(terra::values(r_res) - 273.15, 2))
      terra::units(r_res) <- "deg_C"
      terra::varnames(r_res) <- var
    } else {
      r_res <- terra::setValues(x = r_res,values = terra::values(r_res) / 86400)
      terra::units(r_res) <- "kg m-2 s-1"
      terra::varnames(r_res) <- var
    }
  }
  out_file <- file.path(out_dir, paste0(terra::names(r), ".tif"))
  dots <- list(...)
  if (!is.null(dots$gdal)) {
    gdal_opts <- dots$gdal
    dots$gdal <- NULL
  } else {
    gdal_opts <- c("COMPRESS=LZW", "TFW=NO", "PREDICTOR=3")
  }
  do.call(terra::writeRaster,
          c(list(x = r_res, filename = out_file, gdal = gdal_opts), dots))
  return(out_file)
}
