#' @title summarise_TraCESahul
#' @description Summarise TraCESahul climate data into annual, monthly, or seasonal periods.
#'
#' This function aggregates \emph{TraCESahul} SpatRaster data by year, month, or
#' austral season. Monthly and seasonal summaries can be calculated across the
#' entire record or within user defined windows based on the time dimension.
#' Seasonal grouping follows standard austral seasons, with DJF treated
#' correctly when time intervals are uneven. The summarising function may be
#' any supported \code{\link[terra:tapp]{terra::tapp}} method or a custom
#' function.
#'
#' Data before and after 1500 CE are kept separate for all summary types and
#' are \emph{never} combined in the same aggregation.
#'
#' @param x a \code{\link[terra:rast]{SpatRast}}. \strong{Must} be created using \code{\link{import_TraCESahul}}.
#' @param type which timesteps to summarise over. Can be seasonal or annual. Default: "annual".
#' @param sumfun summary function to be applied. See \code{\link[terra]{tapp}} for details. Default: "mean".
#' @param window the size of the window to aggregate across. Default = NULL.
#' \emph{e.g.} setting a value of 30 means that data will be averaged over 30 years (\emph{i.e.} 3 decadal timesteps, or 30 years for post-1500 data)
#' @param ... additional argument to \code{\link[terra]{tapp}}
#' @return a \code{\link[terra:rast]{SpatRast}} summarised over the chosen value
#' @examples
#' \dontrun{
#' pr_summary <- summarise_TraCEShaul(
#' x = pr_chunks[[1:120]], # first 120 layers only
#' type = "seasonal",
#' sumfun = "mean")
#' }
#' @export
#' @importFrom terra rast tapp time depth names nlyr depthName depthUnit
#'
summarise_TraCESahul <- function(x, type = "annual", sumfun = "mean", window = NULL, ...) {
  stopifnot("x must be output from `import_TraCESahul`" = isTRUE(attr(x, "TraCESahul")))
  stopifnot(
    "window must be NULL or numeric > 10" =
      (is.null(window)) || (is.numeric(window) && window > 10)
  )
  type <- match.arg(type, choices = c("annual", "monthly", "seasonal"),
                    several.ok = FALSE)

  if (is.character(sumfun) && !sumfun %in%
      c("sum", "mean", "median", "modal", "which", "which.min", "which.max",
        "min", "max", "prod", "any", "all", "sd", "std", "first")) {
    sumfun <- match.fun(sumfun)
  }
  years <- terra::time(x)
  months <- terra::depth(x)
  n <- terra::nlyr(x)
  # annual
  if (type == "annual") {
    out <- terra::tapp(x, index = years, fun = sumfun, ...)
    terra::time(out) <- unique(years)
    names(out) <- paste0("Year_", terra::time(out))
    return(out)
  }
  # compute pre-post 1500
  zone <- ifelse(years < 1500, "pre", "post")
  # monthly no window
  if (type == "monthly" && is.null(window)) {
    grp <- paste0(zone, "_", months) # ensure pre-post 1500 split
    out <- terra::tapp(x, index = grp, fun = sumfun, ...)
    rep_time <- tapply(years, grp, median)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    names(out) <- names(out)
    # use terra::depth for months
    n_out <- terra::nlyr(out)
    terra::depth(out) <- rep(1:12, length.out = n_out)
    terra::depthName(out) <- "Month"
    terra::depthUnit(out) <- "month"
    return(out)
  }
  # monthly windowed
  if (type == "monthly" && !is.null(window)) {
    grp <- character(n)
    for (z in c("pre", "post")) {
      idx <- which(zone == z)
      if (length(idx) == 0) next
      nz <- length(idx)
      if (nz %% window != 0) {
        warning(
          sprintf("Last window incomplete in %s1500 zone: %s layers.",
                  ifelse(z == "pre", "<", ">="),
                  nz - floor(nz / window) * window),
          call. = TRUE, immediate. = TRUE)
      }
      window_id <- rep(seq_len(ceiling(nz / window)), each = window)[1:nz]
      grp[idx] <- paste0("M", months[idx], "_", z, "_", window_id)
    }
    out <- terra::tapp(x, index = grp, fun = sumfun, ...)
    rep_time <- tapply(years, grp, max)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    names(out) <- names(out)
    terra::depth(out) <- rep(1:12, length.out = terra::nlyr(out))
    terra::depthName(out) <- "Month"
    terra::depthUnit(out) <- "month"
    return(out)
  }
  # seasonal assignment
  season_name <- character(n)
  season_name[months %in% c(12, 1, 2)] <- "DJF"
  season_name[months %in% c(3, 4, 5)] <- "MAM"
  season_name[months %in% c(6, 7, 8)] <- "JJA"
  season_name[months %in% c(9, 10, 11)] <- "SON"
  season_year_index <- integer(n)
  si <- 1L
  season_year_index[1] <- si
  for (i in 2:n) {
    if (months[i] == 12) si <- si + 1L
    season_year_index[i] <- si
  }
  rep_year <- years
  djf_dec <- which(season_name == "DJF" & months == 12)
  rep_year[djf_dec] <- ifelse(djf_dec < n, years[djf_dec + 1], years[djf_dec])
  zone <- ifelse(rep_year < 1500, "pre", "post") # pre-post 1500 split
  # seasonal no window
  if (type == "seasonal" && is.null(window)) {
    grp <- paste0(zone, "_", season_name) # ensure pre-post 1500 split
    out <- terra::tapp(x, grp, fun = sumfun, ...)
    rep_time <- tapply(rep_year, grp, median)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    names(out) <- names(out)
    terra::depthName(out) <- "Season"
    terra::depthUnit(out) <- ""
    return(out)
  }
  # seasonal with window
  if (type == "seasonal" && !is.null(window)) {
    grp <- character(n)
    for (z in c("pre", "post")) {
      idx <- which(zone == z)
      if (length(idx) == 0) next
      yr_min_z <- min(rep_year[idx])
      yr_max_z <- max(rep_year[idx])
      breaks_z <- seq(yr_min_z, yr_max_z + window, by = window)
      window_id_z <- cut(rep_year[idx], breaks = breaks_z,
                         include.lowest = TRUE, labels = FALSE)
      grp[idx] <- paste0(season_name[idx], "_", z, "_", window_id_z)
    }
    out <- terra::tapp(x, grp, fun = sumfun, ...)
    yr_lookup <- tapply(rep_year, grp, max)
    terra::time(out) <- as.numeric(yr_lookup[names(out)])
    names(out) <- names(out)
    return(out)
  }
}
