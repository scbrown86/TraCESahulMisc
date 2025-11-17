#' @title Summarise \emph{TraCESahul} climate data
#'
#' @description
#' Summarise an imported \emph{TraCESahul} SpatRaster to annual, monthly, or
#' seasonal climatologies. Monthly and seasonal summaries can be taken across the
#' full record or grouped into right-aligned windows of years. Windows are
#' guaranteed not to mix data from before and after 1500 CE.
#'
#' Seasonal summaries follow austral seasons. December is treated as belonging to
#' the following year/timestep to ensure correct DJF alignment, including when
#' time steps are spaced irregularly or when windowed summaries are requested.
#'
#' @param x A SpatRaster created by \code{\link{import_TraCESahul}} containing a
#'   \code{TraCESahul} attribute.
#' @param type Character; one of \code{"annual"}, \code{"monthly"}, or
#'   \code{"seasonal"}, specifying the temporal summary to compute.
#' @param sumfun A character string or function giving the aggregation method.
#'   May be any function supported by \code{\link[terra:tapp]{terra::tapp}} or a
#'   user-defined function.
#' @param window Numeric; optional number of years to group together when
#'   producing monthly or seasonal summaries. Must be greater than or equal to
#'   10. Windows are right-aligned and incomplete trailing windows are excluded.
#'   When \code{window = 10}, pre-1500 data are treated as already representing a
#'   decadal climatology and handled accordingly.
#' @param ... Additional arguments passed to \code{\link[terra:tapp]{terra::tapp}}.
#'
#' @return A SpatRaster summarised according to the chosen \code{type} and
#'   \code{window}. Monthly outputs contain one layer per month for each
#'   timezone (pre- and post-1500). Seasonal outputs contain one layer per
#'   austral season for each timezone when \code{window = NULL}. Windowed
#'   outputs contain one layer per window–season or window–month combination.
#'
#' @details
#' Windows are constructed using right-aligned year ranges so the final window
#' ends on the latest year in the data. Pre- and post-1500 periods are always
#' windowed separately. For seasonal summaries, December is shifted into the
#' following year before window assignment to preserve correct DJF grouping.
#' When \code{window = 10}, pre-1500 data are already decadal climatologies and
#' are treated accordingly when deriving annual, monthly, or seasonal outputs.
#'
#' @examples
#' \dontrun{
#' # summarise the data across 30 year intervals
#' pr_monthly_win <- summarise_TraCESahul(x = sahul_pr, type = "monthly",
#'                                        sumfun = "mean", window = 30)
#' tasmax_monthly_win <- summarise_TraCESahul(x = sahul_tasmax, type = "monthly",
#'                                            sumfun = "mean", window = 30)
#' tasmin_monthly_win <- summarise_TraCESahul(sahul_tasmin, type = "monthly",
#'                                            sumfun = "mean", window = 30)
#' }
#'
#' @export
#' @importFrom terra rast tapp time depth names nlyr depthName depthUnit
#' @importFrom stats median

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
  u <- terra::units(x)[1]
  v <- terra::varnames(x)[1]
  ln <- terra::longnames(x)[1]
  # annual
  if (type == "annual" && is.null(window)) {
    out <- terra::tapp(x, index = years, fun = sumfun, ...)
    terra::time(out) <- unique(years)
    names(out) <- paste0("Year_", terra::time(out))
    terra::depthName(out) <- "Annual"
    terra::depthUnit(out) <- "year"
    terra::units(out) <- u
    terra::varnames(out) <- v
    terra::longnames(out) <- ln
    attr(out, "TraCESahul") <- TRUE
    return(out)
  }
  # annual window; pre-1500: treat window==10 as already 10-yr climatology, otherwise apply windows
  if (type == "annual" && !is.null(window)) {
    zone_annual <- ifelse(years < 1500, "pre", "post")
    grp <- character(length(years))
    # pre-1500 handling
    pre_idx <- which(zone_annual == "pre")
    if (length(pre_idx) > 0) {
      if (window == 10) {
        # window == 10, collapse each 10-yr climatology to one annual layer
        grp[pre_idx] <- paste0("pre_", years[pre_idx])
      } else {
        # window > 10, apply right-aligned windows within pre-1500
        yrs_z <- years[pre_idx]
        ymin <- min(yrs_z)
        ymax <- max(yrs_z)
        nwin <- floor((ymax - ymin + 1L) / window)
        if (nwin > 0) {
          win_ends <- ymax - ((nwin:1 - 1) * window)
          win_starts <- win_ends - window + 1L
          win_map <- mapply(function(s, e) s:e, win_starts, win_ends, SIMPLIFY = FALSE)
          win_id <- integer(length(yrs_z))
          for (i in seq_along(win_map)) win_id[yrs_z %in% win_map[[i]]] <- i
          keep <- win_id > 0
          grp[pre_idx[keep]] <- paste0("pre_", win_id[keep])
        }
      }
    }
    # post-1500 handling (always windowed)
    post_idx <- which(zone_annual == "post")
    if (length(post_idx) > 0) {
      yrs_z <- years[post_idx]
      ymin <- min(yrs_z)
      ymax <- max(yrs_z)
      nwin <- floor((ymax - ymin + 1L) / window)
      if (nwin > 0) {
        win_ends <- ymax - ((nwin:1 - 1) * window)
        win_starts <- win_ends - window + 1L
        win_map <- mapply(function(s, e) s:e, win_starts, win_ends, SIMPLIFY = FALSE)
        win_id <- integer(length(yrs_z))
        for (i in seq_along(win_map)) win_id[yrs_z %in% win_map[[i]]] <- i
        keep <- win_id > 0
        grp[post_idx[keep]] <- paste0("post_", win_id[keep])
      }
    }
    keep_layers <- grp != ""
    out <- terra::tapp(x[[keep_layers]], index = grp[keep_layers], fun = sumfun, ...)
    rep_time <- tapply(years[keep_layers], grp[keep_layers], max)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    terra::depthName(out) <- "Annual"
    terra::depthUnit(out) <- "year"
    names(out) <- names(out)
    terra::units(out) <- u
    terra::varnames(out) <- v
    terra::longnames(out) <- ln
    attr(out, "TraCESahul") <- TRUE
    return(out)
  }
  # define timezones for pre/post 1500
  zone_monthly <- ifelse(years < 1500, "pre", "post")
  # monthly, separate pre/post, 12 or 24 layers
  if (type == "monthly" && is.null(window)) {
    grp <- paste0(zone_monthly, "_", months)
    out <- terra::tapp(x, index = grp, fun = sumfun, ...)
    rep_time <- tapply(years, grp, stats::median)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    terra::depth(out) <- rep(1:12, length.out = terra::nlyr(out))
    terra::depthName(out) <- "Month"
    terra::depthUnit(out) <- "month"
    names(out) <- paste0(month.abb, "_", terra::time(out))
    terra::units(out) <- u
    terra::varnames(out) <- v
    terra::longnames(out) <- ln
    attr(out, "TraCESahul") <- TRUE
    return(out)
  }
  # monthly window; bin by right-aligned years within each zone
  if (type == "monthly" && !is.null(window)) {
    grp <- character(n)
    # pre-1500 handling
    pre_idx <- which(zone_monthly == "pre")
    if (length(pre_idx) > 0) {
      if (window == 10) {
        # pre-1500 are already 10-yr monthly climatologies so keep as-is
        grp[pre_idx] <- paste0("pre_", month.abb[months[pre_idx]], "_", years[pre_idx])
      } else {
        # window > 10, must apply right-aligned windows to pre-1500 too
        yrs_z <- years[pre_idx]
        ymin <- min(yrs_z)
        ymax <- max(yrs_z)
        nwin <- floor((ymax - ymin + 1) / window)
        if (nwin > 0) {
          win_ends <- ymax - ((nwin:1 - 1) * window)
          win_starts <- win_ends - window + 1
          win_map <- mapply(function(s, e) s:e, win_starts, win_ends, SIMPLIFY = FALSE)
          win_id <- integer(length(yrs_z))
          for (i in seq_along(win_map)) {
            win_id[yrs_z %in% win_map[[i]]] <- i
          }
          keep <- win_id > 0
          grp[pre_idx[keep]] <- paste0(month.abb[months[pre_idx[keep]]], "_pre_", win_id[keep])
        }
      }
    }
    # post-1500 handling
    post_idx <- which(zone_monthly == "post")
    if (length(post_idx) > 0) {
      yrs_z <- years[post_idx]
      ymin <- min(yrs_z)
      ymax <- max(yrs_z)
      nwin <- floor((ymax - ymin + 1) / window)
      if (nwin > 0) {
        win_ends <- ymax - ((nwin:1 - 1) * window)
        win_starts <- win_ends - window + 1
        win_map <- mapply(function(s, e) s:e, win_starts, win_ends, SIMPLIFY = FALSE)
        win_id <- integer(length(yrs_z))
        for (i in seq_along(win_map)) {
          win_id[yrs_z %in% win_map[[i]]] <- i
        }
        keep <- win_id > 0
        grp[post_idx[keep]] <- paste0(month.abb[months[post_idx[keep]]], "_post_", win_id[keep])
      }
    }
    # only keep valid layers
    keep_layers <- grp != ""
    out <- terra::tapp(x[[keep_layers]], index = grp[keep_layers], fun = sumfun, ...)
    rep_time <- tapply(years[keep_layers], grp[keep_layers], max)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    terra::depth(out) <- rep(1:12, length.out = terra::nlyr(out))
    terra::depthName(out) <- "Month"
    terra::depthUnit(out) <- "month"
    names(out) <- paste0(names(out))
    terra::units(out) <- u
    terra::varnames(out) <- v
    terra::longnames(out) <- ln
    attr(out, "TraCESahul") <- TRUE
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
    rep_time <- tapply(rep_year, grp, stats::median)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    names(out) <- names(out)
    terra::depthName(out) <- "Season"
    terra::depthUnit(out) <- ""
    terra::units(out) <- u
    terra::varnames(out) <- v
    terra::longnames(out) <- ln
    attr(out, "TraCESahul") <- TRUE
    return(out)
  }
  # seasonal window; handle pre-1500 and post-1500 separately and preserve DJF logic
  if (type == "seasonal" && !is.null(window)) {
    grp <- character(n)
    # pre-1500 handling
    pre_idx <- which(zone_seasonal == "pre")
    if (length(pre_idx) > 0) {
      if (window == 10) {
        # window == 10, recompute seasonal summaries from existing decadal climatology
        # DJF requires: Dec(previous), Jan/Feb(current). Shift December forward
        yrs_z <- rep_year[pre_idx]
        months_z <- months[pre_idx]
        # assign base window ids (each decadal climatology)
        base_id <- match(yrs_z, sort(unique(yrs_z)))
        # shift December forward so DJF aligns correctly
        shifted_id <- base_id
        dec_pos <- which(months_z == 12)
        if (length(dec_pos) > 0) shifted_id[dec_pos] <- base_id[dec_pos] + 1L
        # drop first shifted window (no December available before it)
        keep <- shifted_id > 1
        if (any(keep)) {
          grp[pre_idx[keep]] <- paste0(season_name[pre_idx[keep]], "_pre_", shifted_id[keep])
        }
      } else {
        # window > 10, apply right-aligned windows to pre-1500
        yrs_z <- rep_year[pre_idx]
        ymin <- min(yrs_z)
        ymax <- max(yrs_z)
        nwin <- floor((ymax - ymin + 1)/window)
        if (nwin > 0) {
          win_ends <- ymax - ((nwin:1 - 1) * window)
          win_starts <- win_ends - window + 1
          win_map <- mapply(function(s, e) s:e, win_starts, win_ends, SIMPLIFY = FALSE)
          win_id <- integer(length(yrs_z))
          for (i in seq_along(win_map)) win_id[yrs_z %in% win_map[[i]]] <- i
          # shift December into following window so DJF is correct
          dec_pos <- which(months[pre_idx] == 12)
          if (length(dec_pos) > 0) win_id[dec_pos] <- win_id[dec_pos] + 1L
          keep <- win_id > 0 & win_id <= nwin
          grp[pre_idx[keep]] <- paste0(season_name[pre_idx[keep]], "_pre_", win_id[keep])
        }
      }
    }
    # post-1500 handling (always windowed)
    post_idx <- which(zone_seasonal == "post")
    if (length(post_idx) > 0) {
      yrs_z <- rep_year[post_idx]
      ymin <- min(yrs_z)
      ymax <- max(yrs_z)
      nwin <- floor((ymax - ymin + 1)/window)
      if (nwin > 0) {
        win_ends <- ymax - ((nwin:1 - 1) * window)
        win_starts <- win_ends - window + 1
        win_map <- mapply(function(s, e) s:e, win_starts, win_ends, SIMPLIFY = FALSE)
        win_id <- integer(length(yrs_z))
        for (i in seq_along(win_map)) win_id[yrs_z %in% win_map[[i]]] <- i
        # shift December into following window for DJF
        dec_pos <- which(months[post_idx] == 12)
        if (length(dec_pos) > 0) win_id[dec_pos] <- win_id[dec_pos] + 1L
        keep <- win_id > 0 & win_id <= nwin
        grp[post_idx[keep]] <- paste0(season_name[post_idx[keep]], "_post_", win_id[keep])
      }
    }
    keep_layers <- grp != ""
    out <- terra::tapp(x[[keep_layers]], index = grp[keep_layers], fun = sumfun, ...)
    rep_time <- tapply(rep_year[keep_layers], grp[keep_layers], max)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    names(out) <- names(out)
    terra::depthName(out) <- "Season"
    terra::depthUnit(out) <- ""
    terra::units(out) <- u
    terra::varnames(out) <- v
    terra::longnames(out) <- ln
    attr(out, "TraCESahul") <- TRUE
    return(out)
  }
}
