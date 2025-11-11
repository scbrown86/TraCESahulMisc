#' @title summarise_TraCESahul
#' @description FUNCTION_DESCRIPTION
#' @param x a \code{\link[terra:rast]{SpatRast}}. \strong{Must} be created using \code{\link{import_TraCESahul}}.
#' @param type which timesteps to summarise over. Can be seasonal or annual. Default: "annual".
#' @param sumfun summary function to be applied. See \code{\link[terra]{tapp}} for details. Default: "mean".
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
summarise_TraCESahul <- function(x, type = "annual", sumfun = "mean",
                                 window = NULL,
                                 ...) {
  stopifnot("x must be output from `import_TraCESahul`" = isTRUE(attr(x, "TraCESahul")))
  stopifnot(
    "window must be NULL or numeric > 10" =
      (is.null(window)) ||
      (is.numeric(window) && window > 10)
  )
  type <- match.arg(type, choices = c("annual", "monthly", "seasonal"),
                    several.ok = FALSE)
  # Check sumfun against terra::tapp fun methods
  if (is.character(sumfun) && !sumfun %in%
      c("sum", "mean", "median", "modal", "which", "which.min", "which.max",
        "min", "max", "prod", "any", "all", "sd", "std", "first")) {
    sumfun <- match.fun(sumfun)
  }
  years  <- terra::time(x)
  months <- terra::depth(x)
  # annual
  if (type == "annual") {
    ann <- terra::tapp(x, index = yrs, fun = sumfun, ...)
    terra::time(ann) <- unique(years)
    return(ann)
  }
  # monthly
  if (type == "monthly" && is.null(window)) {
    mon <- terra::tapp(x, index = months, fun = sumfun, ...)
    terra::time(mon) <- rep(median(unique(years)), 12)
    names(mon) <- month.abb
    terra::depth(mon) <- 1:12
    terra::depthName(mon) <- "Month"
    terra::depthUnit(mon) <- "month"
    return(mon)
  } else if (type == "monthly" && is.numeric(window)) {
    n <- terra::nlyr(x)
    breaks <- rep(seq(1, n + window, by = window), each = window)
    breaks <- breaks[1:n]
    if (n %% window != 0) {
      mx <- min(table(breaks))
      warning(sprintf("Final window not divisible by number of layers. Final window calculated on %s layers.", mx),
              call. = TRUE, immediate. = TRUE)
    }
    grp <- paste0(months, "_", breaks)
    stopifnot("Too many groups?" = terra::nlyr(x) == length(grp))
    mon <- terra::tapp(x, index = grp, fun = sumfun, ...)
    terra::time(mon) <- rep(tapply(years, breaks, function(z) max(z)), each = 12)
    names(mon) <- paste0(month.abb,"_", terra::time(mon))
    terra::depth(mon) <- rep(1:12, each = terra::nlyr(mon)/12)
    terra::depthName(mon) <- "Month"
    terra::depthUnit(mon) <- "month"
    return(mon)
  }
  # seasonal
  n <- length(months)
  season_name <- character(n)
  season_name[months %in% c(12, 1, 2)]  <- "DJF"
  season_name[months %in% c(3, 4, 5)]   <- "MAM"
  season_name[months %in% c(6, 7, 8)]   <- "JJA"
  season_name[months %in% c(9, 10, 11)] <- "SON"
  season_year_index <- integer(n)
  si <- 1L
  season_year_index[1] <- si
  for (i in 2:n) {
    if (months[i] == 12) si <- si + 1L
    season_year_index[i] <- si
  }
  # for DJF: Dec uses following year if it exists
  rep_year <- numeric(n)
  for (i in seq_len(n)) {
    if (season_name[i] == "DJF" && months[i] == 12) {
      rep_year[i] <- if (i < n) years[i + 1] else years[i]
    } else {
      rep_year[i] <- years[i]
    }
  }
  if (is.null(window)) {
    # summarise across all seasons (so only 4 returned values)
    grp <- season_name
    out <- terra::tapp(x, grp, fun = sumfun, ...)
    terra::time(out) <- rep(median(years), 4)
    terra::depthName(out) <- "austral season"
    terra::depthUnit(out) <- ""
    return(out)
  }
  # Compute number of seasonal windows
  yr_min <- min(rep_year)
  yr_max <- max(rep_year)
  breaks <- seq(yr_min, yr_max + window, by = window)
  window_id <- cut(rep_year, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  # Final group = season + window
  grp <- paste0(season_name, "_", window_id)
  out <- terra::tapp(x, index = grp, fun = sumfun, ...)
  # Representative time for each group = max rep_year in that group
  out_names <- unique(grp)
  yr_lookup <- tapply(rep_year, grp, function(z) max(z))
  terra::time(out) <- as.numeric(yr_lookup[out_names])
  names(out) <- paste0(c("DJF", "MAM", "JJA", "SON"), "_", abs(terra::time(out)))
  out
}
