#' @title Summarise \emph{TraCESahul} climate data
#'
#' @description Aggregate an imported \emph{TraCESahul} SpatRaster to annual, monthly, or seasonal summaries.
#' Monthly and seasonal summaries can either be taken across the full record or grouped into user-defined
#' windows of years. Windows are right-aligned, and any incomplete trailing windows are excluded.
#' Seasonal summaries use austral seasons, with DJF December handled as part of the following year.
#' Data from before and after 1500 CE are always processed separately and never mixed.
#'
#' @param x A SpatRaster imported using \code{\link{import_TraCESahul}} with a \code{TraCESahul} attribute.
#' @param type Character; one of \code{"annual"}, \code{"monthly"}, or \code{"seasonal"} specifying the summary type.
#' @param sumfun A function or character string indicating the summarising function to use. Can be any
#' function supported by \code{\link[terra:tapp]{terra::tapp}} or a custom function.
#' @param window Numeric; optional number of years to group together for monthly or seasonal summaries.
#' Must be greater than or equal to 10. Windows are right-aligned and incomplete windows are excluded.
#' @param ... Additional arguments passed to \code{\link[terra:tapp]{terra::tapp}}.
#'
#' @return A SpatRaster summarised according to the specified type and window. Monthly outputs will have
#' 12 layers per timezone (pre- or post-1500), and seasonal outputs will have four layers per timezone when
#' \code{window = NULL}.
#'
#' @details Right-aligned windows ensure that the last window ends on the maximum year in the data,
#' without including incomplete periods. For DJF December is considered part of the following year,
#' ensuring correct austral season alignment.

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
      (is.null(window)) || (is.numeric(window) && window >= 10)
  )
  type <- match.arg(type, choices = c("annual", "monthly", "seasonal"), several.ok = FALSE)
  # make sure sumfun is suitable
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
  # define timezones for pre/post 1500
  zone_monthly <- ifelse(years < 1500, "pre", "post")
  # monthly, separate pre/post, 12 or 24 layers
  if (type == "monthly" && is.null(window)) {
    grp <- paste0(zone_monthly, "_", months)
    out <- terra::tapp(x, index = grp, fun = sumfun, ...)
    rep_time <- tapply(years, grp, median)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    names(out) <- names(out)
    terra::depth(out) <- rep(1:12, length.out = terra::nlyr(out))
    terra::depthName(out) <- "Month"
    terra::depthUnit(out) <- "month"
    return(out)
  }
  # monthly window; bin by right-aligned years within each zone
  if (type == "monthly" && !is.null(window)) {
    grp <- character(n)
    for (z in c("pre", "post")) {
      idx <- which(zone_monthly == z)
      if (length(idx) == 0) next
      yrs_z <- years[idx]
      ymin <- min(yrs_z)
      ymax <- max(yrs_z)
      nwin <- floor((ymax - ymin + 1)/window)
      if (nwin == 0) next
      win_ends <- ymax - ((nwin:1 - 1) * window)
      win_starts <- win_ends - window + 1
      win_map <- mapply(function(s, e) s:e, win_starts, win_ends, SIMPLIFY = FALSE)
      win_id <- integer(length(yrs_z))
      for (i in seq_along(win_map)) {
        win_id[yrs_z %in% win_map[[i]]] <- i
      }
      keep <- win_id > 0
      grp[idx[keep]] <- paste0(month.abb[months[idx[keep]]], "_", z, "_", win_id[keep])
    }
    # only keep layers in a group
    keep_layers <- grp != ""
    out <- terra::tapp(x[[keep_layers]], index = grp[keep_layers], fun = sumfun, ...)
    rep_time <- tapply(years[keep_layers], grp[keep_layers], max)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    terra::depth(out) <- rep(1:12, length.out = terra::nlyr(out))
    terra::depthName(out) <- "Month"
    terra::depthUnit(out) <- "month"
    names(out) <- names(out)
    return(out)
  }
  # seasonal assignment
  season_name <- character(n)
  season_name[months %in% c(12, 1, 2)] <- "DJF"
  season_name[months %in% c(3, 4, 5)] <- "MAM"
  season_name[months %in% c(6, 7, 8)] <- "JJA"
  season_name[months %in% c(9, 10, 11)] <- "SON"
  # representative year for seasonal (DJF Dec takes following year if present)
  rep_year <- years
  djf_dec <- which(season_name == "DJF" & months == 12)
  rep_year[djf_dec] <- ifelse(djf_dec < n, years[djf_dec + 1], years[djf_dec])
  zone_seasonal <- ifelse(rep_year < 1500, "pre", "post")
  # seasonal, separate pre/post, up to 8 layers
  if (type == "seasonal" && is.null(window)) {
    grp <- paste0(zone_seasonal, "_", season_name)
    out <- terra::tapp(x, grp, fun = sumfun, ...)
    rep_time <- tapply(rep_year, grp, median)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    names(out) <- names(out)
    terra::depthName(out) <- "Season"
    terra::depthUnit(out) <- ""
    return(out)
  }
  # seasonal window; bin by right-aligned years within each zone
  if (type == "seasonal" && !is.null(window)) {
    grp <- character(n)
    for (z in c("pre", "post")) {
      idx <- which(zone_seasonal == z)
      if (length(idx) == 0) next
      yrs_z <- rep_year[idx]
      ymin <- min(yrs_z)
      ymax <- max(yrs_z)
      nwin <- floor((ymax - ymin + 1)/window)
      if (nwin == 0) next
      win_ends <- ymax - ((nwin:1 - 1) * window)
      win_starts <- win_ends - window + 1
      win_map <- mapply(function(s, e) s:e, win_starts, win_ends, SIMPLIFY = FALSE)
      win_id <- integer(length(yrs_z))
      for (i in seq_along(win_map)) {
        win_id[yrs_z %in% win_map[[i]]] <- i
      }
      keep <- win_id > 0
      grp[idx[keep]] <- paste0(season_name[idx[keep]], "_", z, "_", win_id[keep])
    }
    keep_layers <- grp != ""
    out <- terra::tapp(x[[keep_layers]], index = grp[keep_layers], fun = sumfun, ...)
    rep_time <- tapply(rep_year[keep_layers], grp[keep_layers], max)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    names(out) <- names(out)
    return(out)
  }
}
