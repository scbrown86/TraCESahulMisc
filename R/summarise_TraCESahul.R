#' @title Summarise \emph{TraCESahul} climate data
#'
#' @description
#' Summarise an imported \emph{TraCESahul} SpatRaster to annual, monthly, or
#' seasonal climatologies. Monthly and seasonal summaries can be taken across the
#' full record or grouped into right-aligned windows of years.
#'
#' Windows are always right-aligned to the most recent data, with older,
#' incomplete leading years dropped as necessary, regardless of whether the
#' underlying data are annual (post-1500) or decadal (pre-1500).
#'
#' Windows are guaranteed not to mix data from before and after 1500 CE.
#'
#' Seasonal summaries follow austral seasons. December is treated as belonging to
#' the following year/timestep to ensure correct DJF alignment, including when
#' time steps are spaced irregularly or when windowed summaries are requested.
#' A DJF window is only returned if it contains a complete December + January +
#' February triple. Any DJF window that would otherwise need December data from
#' before the record began, or from across the 1500 CE boundary, is dropped in
#' its entirety. This necessarily always removes the earliest DJF window from
#' each time period (pre-, post-1500).
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
#' \code{window}. Monthly outputs contain one layer per month for each
#' period (pre- and post-1500). Seasonal outputs contain one layer per austral
#' season for each period when \code{window = NULL}. Windowed outputs contain
#' one layer per window-season, window-month, or window combination.
#'
#' @details
#' Windows are constructed separately for the pre- and post-1500 periods, from
#' whatever unique years are actually present, so they work correctly whether
#' the underlying data are annual (post-1500) or decadal (pre-1500). Each zone
#' keeps as many complete, right-aligned windows as its data allows. Remaining
#' years at the older end that don't fill a complete window are dropped.
#'
#' For seasonal summaries, December is shifted into the following year before
#' window assignment to preserve correct DJF grouping. December is only used
#' if the following timestep sits in the same period; if not (or if there
#' is no following timestep at all), that December cannot be paired with a
#' Jan/Feb without crossing the decadal/annual boundary, so it's dropped rather
#' than reassigned to the other period. After windows are assigned, any DJF window
#' that doesn't contain equal counts of December, January, and February layers
#' is dropped entirely, since a partial DJF window can only occur at the very
#' start of a periods record.
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
#' @importFrom utils tail

summarise_TraCESahul <- function(x, type = "annual", sumfun = "mean", window = NULL, ...) {
  stopifnot("x must be output from `import_TraCESahul`" = isTRUE(attr(x, "TraCESahul")))
  stopifnot(
    "window must be NULL or numeric >= 10" =
      is.null(window) || (is.numeric(window) && window >= 10),
    "window must be NULL or divisible by 10" =
      is.null(window) || (window %% 10 == 0)
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
  zone <- ifelse(years < 1500, "pre", "post")
  # right-aligned window lookup for a set of "years", dropping whatever incomplete
  # leading years don't fill a complete window. Works for any step size (e.g.
  # annual post-1500 data, decadal pre-1500 data) since it counts actual unique
  # years present rather than assuming a fixed year-to-year step
  make_window_lookup <- function(yrs, window) {
    uy <- sort(unique(yrs))
    step <- if (length(uy) > 1) min(diff(uy)) else window
    steps_per_window <- window / step
    nwin <- length(uy) %/% steps_per_window
    if (nwin == 0) return(stats::setNames(integer(0), character(0)))
    keep_uy <- utils::tail(uy, nwin * steps_per_window)
    return(stats::setNames(rep(seq_len(nwin), each = steps_per_window), keep_uy))
  }
  apply_window_lookup <- function(lookup, yrs) {
    w <- unname(lookup[as.character(yrs)])
    w[is.na(w)] <- 0L
    return(w)
  }
  apply_common_attrs <- function(out) {
    terra::units(out) <- u
    terra::varnames(out) <- v
    terra::longnames(out) <- ln
    attr(out, "TraCESahul") <- TRUE
    return(out)
  }
  # annual
  if (type == "annual") {
    if (is.null(window)) {
      grp <- as.character(years)
    } else {
      grp <- character(n)
      for (z in c("pre", "post")) {
        idx <- which(zone == z)
        if (length(idx) == 0) next
        lookup <- make_window_lookup(years[idx], window)
        w <- apply_window_lookup(lookup, years[idx])
        keep <- w > 0
        grp[idx[keep]] <- paste0(z, "_", w[keep])
      }
    }
    keep_layers <- grp != ""
    out <- terra::tapp(x[[keep_layers]], index = grp[keep_layers], fun = sumfun, ...)
    rep_time <- tapply(years[keep_layers], grp[keep_layers], max)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    names(out) <- paste0("Year_", terra::time(out))
    terra::depthName(out) <- "Annual"
    terra::depthUnit(out) <- "year"
    out <- apply_common_attrs(out)
    return(out)
  }
  # monthly, separate pre/post, one layer per month per zone (or per window)
  if (type == "monthly") {
    if (is.null(window)) {
      grp <- paste0(zone, "_", months)
      time_fun <- function(y) round(stats::median(y))
    } else {
      grp <- character(n)
      time_fun <- max
      for (z in c("pre", "post")) {
        idx <- which(zone == z)
        if (length(idx) == 0) next
        lookup <- make_window_lookup(years[idx], window)
        w <- apply_window_lookup(lookup, years[idx])
        keep <- w > 0
        grp[idx[keep]] <- paste0(month.abb[months[idx[keep]]], "_", z, "_", w[keep])
      }
    }
    keep_layers <- grp != ""
    out <- terra::tapp(x[[keep_layers]], index = grp[keep_layers], fun = sumfun, ...)
    rep_time <- tapply(years[keep_layers], grp[keep_layers], time_fun)
    rep_month <- tapply(months[keep_layers], grp[keep_layers], unique)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    terra::depth(out) <- as.integer(rep_month[names(out)])
    terra::depthName(out) <- "Month"
    terra::depthUnit(out) <- "month"
    names(out) <- paste0(month.abb[terra::depth(out)], "_", terra::time(out))
    out <- apply_common_attrs(out)
    return(out)
  }
  # seasonal assignment (austral seasons)
  season_name <- character(n)
  season_name[months %in% c(12, 1, 2)] <- "DJF"
  season_name[months %in% c(3, 4, 5)] <- "MAM"
  season_name[months %in% c(6, 7, 8)] <- "JJA"
  season_name[months %in% c(9, 10, 11)] <- "SON"
  # December belongs to the following timestep's year for correct DJF alignment.
  # a December is only usable if the following timestep sits in the same period,
  # otherwise there's no valid Jan/Feb to pair it with without crossing the
  # boundary, so it's dropped rather than reassigned to the other period.
  rep_year <- years
  dec_idx <- which(months == 12)
  nxt <- pmin(dec_idx + 1, n)
  rep_year[dec_idx] <- years[nxt]
  dec_usable <- rep(TRUE, n)
  dec_usable[dec_idx] <- (dec_idx < n) & (zone[dec_idx] == zone[nxt])
  effective_year <- ifelse(months == 12, rep_year, years)
  # seasonal, no window: full-record climatology, separate pre/post
  if (type == "seasonal" && is.null(window)) {
    grp <- character(n)
    usable <- months != 12 | dec_usable
    grp[usable] <- paste0(ifelse(effective_year[usable] < 1500, "pre", "post"),
                          "_", season_name[usable])
    keep_layers <- grp != ""
    out <- terra::tapp(x[[keep_layers]], index = grp[keep_layers], fun = sumfun, ...)
    rep_time <- tapply(effective_year[keep_layers], grp[keep_layers], stats::median)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    names(out) <- names(out)
    terra::depthName(out) <- "Season"
    terra::depthUnit(out) <- ""
    out <- apply_common_attrs(out)
    return(out)
  }
  # seasonal window; right-aligned per zone, December shifted into the following
  # window, any DJF window without a complete Dec + Jan + Feb triple dropped
  if (type == "seasonal" && !is.null(window)) {
    grp <- character(n)
    for (z in c("pre", "post")) {
      idx_all <- which(zone == z)
      if (length(idx_all) == 0) {
        next
      }
      lookup <- make_window_lookup(years[idx_all], window)
      idx <- idx_all[months[idx_all] != 12 | dec_usable[idx_all]]
      w <- apply_window_lookup(lookup, effective_year[idx])
      keep <- w > 0
      grp[idx[keep]] <- paste0(season_name[idx[keep]], "_", z, "_", w[keep])
    }
    djf_groups <- unique(grp[grp != "" & grepl("^DJF_", grp)])
    for (g in djf_groups) {
      in_g <- grp == g
      m <- months[in_g]
      if (sum(m == 12) != sum(m == 1) || sum(m == 1) != sum(m == 2)) {
        grp[in_g] <- ""
      }
    }
    keep_layers <- grp != ""
    out <- terra::tapp(x[[keep_layers]], index = grp[keep_layers], fun = sumfun, ...)
    rep_time <- tapply(effective_year[keep_layers], grp[keep_layers], max)
    terra::time(out) <- as.numeric(rep_time[names(out)])
    names(out) <- names(out)
    terra::depthName(out) <- "Season"
    terra::depthUnit(out) <- ""
    out <- apply_common_attrs(out)
    return(out)
  }
}
