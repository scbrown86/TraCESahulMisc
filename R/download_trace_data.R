#' Download and extract the example TraCE-Sahul dataset
#'
#' @description
#' Downloads and extracts a file from box (~1 GB) which contains a spatial
#' subset of the \emph{TraCE-Sahul} dataset for the Arafura–Timor region
#' (125–135°E, 8–13°S) from 21k BP to 1990 CE. Will extract to a folder called
#' \emph{Trace-Sahul} in the specified directory.
#'
#' @param destdir Character. Path to the directory where the data should be downloaded and extracted.
#' @param ... Optional parameters to \code{\link[utils]{download.file}}. \code{mode = "wb"} is fixed.
#' @return Invisibly returns the path to the extracted directory.
#' @examples
#' \dontrun{
#' download_trace_data("~/Downloads")
#' }
#' @note The timeout for the download is set at 30 minutes.
#' @export
download_trace_data <- function(destdir,...) {
  url <- "https://adelaideuniversity.box.com/shared/static/nv7x31cvzjax1q9f9qc1rd8oxjff38e7.gz"
  stopifnot(is.character(destdir), is.character(url))
  destdir <- path.expand(destdir)
  if (!dir.exists(destdir)) {
    dir.create(destdir, recursive = TRUE)
  }
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(1800, old_timeout)) #0.5 hour timeout
  destfile <- file.path(destdir, "TraCE-Sahul.tar.gz")
  message("Downloading: ", url)
  utils::download.file(url, destfile, mode = "wb", quiet = FALSE)
  message("Extracting: ", destfile)
  utils::untar(destfile, exdir = destdir)
  invisible(destdir)
}


#' Download and extract the example TraCE-Sahul CMIP6 dataset
#'
#' @description
#' Downloads and extracts a file from box (~150 MB) which contains a spatial
#' subset of the \emph{TraCE-Sahul} CMIP6 dataset for the Arafura–Timor region
#' (125–135°E, 8–13°S) from 1900 to 2100 for 4 SSP scenarios. Will extract to a
#' folder called \emph{CMIP6} in the specified directory.
#'
#' @param destdir Character. Path to the directory where the data should be downloaded and extracted.
#' @param ... Optional parameters to \code{\link[utils]{download.file}}. \code{mode = "wb"} is fixed.
#' @return Invisibly returns the path to the extracted directory.
#' @examples
#' \dontrun{
#' download_trace_data("~/Downloads")
#' }
#' @note The timeout for the download is set at 30 minutes.
#' @export
download_trace_cmip6_data <- function(destdir,...) {
  url <- "https://adelaideuniversity.box.com/shared/static/utid2eudvfivqyznq2poqql845isirkq.gz"
  stopifnot(is.character(destdir), is.character(url))
  destdir <- path.expand(destdir)
  if (!dir.exists(destdir)) {
    dir.create(destdir, recursive = TRUE)
  }
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(1800, old_timeout)) #0.5 hour timeout
  destfile <- file.path(destdir, "TraCE-Sahul-CMIP6.tar.gz")
  message("Downloading: ", url)
  utils::download.file(url, destfile, mode = "wb", quiet = FALSE)
  message("Extracting: ", destfile)
  utils::untar(destfile, exdir = destdir)
  invisible(destdir)
}
