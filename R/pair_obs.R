#' Pair fossil observations with environmental raster data
#'
#' @description
#' Extract environmental values for fossil or archaeological point records by
#' matching each record to one or more raster layers based on its calibrated
#' age range. Temporal matching uses either exact time slices or right-aligned
#' averaging windows (when \code{window} is supplied). Spatial extraction can
#' use either a fixed-radius buffer or a nearest-neighbour search. The function
#' supports parallel processing via \pkg{\link{future}} and returns a long-format
#' \code{\link[data.table:data.table]{data.table}} containing extracted values
#' for all environmental rasters and mask layers.
#'
#' @details
#' The input \code{data} must contain the columns \code{ID}, \code{Lat},
#' \code{Lon}, \code{Age}, \code{AgeMin}, and \code{AgeMax}. \code{AgeMin}
#' represents the youngest bound of the calibrated confidence interval, and
#' \code{AgeMax} represents the oldest bound. These are matched to raster
#' times \code{ras_time}, either as exact years (\code{window = NULL}) or as
#' right-aligned averaging windows of length \code{window} years
#' (\code{window >= 10}).
#'
#' Spatial extraction is performed using either:
#' \itemize{
#'   \item{\strong{Buffer-based extraction}}{ if \code{buff_width} is supplied;
#'   extraction is performed within a circular buffer of \code{buff_width}
#'   kilometres around the point.}
#'   \item{\strong{Neighbour extraction}}{ if \code{neigh} is supplied; the
#'   \code{neigh} nearest raster cells are identified (4 or 8 neighbours).}
#' }
#'
#' If \code{dist_cut} is supplied, points may be snapped to the nearest valid
#' raster cell within \code{dist_cut} kilometres. Records with no temporal
#' match or exceeding the snapping threshold are removed (and reported) from the
#' output.
#'
#' All rasters in \code{ras_list} must have identical geometry and number of
#' layers. \code{mask_layer} must have the same geometry and temporal depth.
#' The function automatically assigns \code{ras_time} to all rasters, and checks
#' for spatial consistency.
#'
#' When running in parallel, rasters are wrapped using
#' \code{\link[terra:wrap]{terra::wrap}} to safely pass them to parallel workers.
#'
#' @param data A \code{data.frame} or \code{data.table} containing fossil or
#' observational records. Must include columns \code{ID}, \code{Lat}, \code{Lon},
#' \code{Age}, \code{AgeMin}, and \code{AgeMax}. Column order does not matter.
#' @param ras_list A named list of \code{\link[terra:rast]{SpatRaster}} objects
#' containing environmental variables. All rasters must be aligned, have identical
#' geometry, and share the same temporal dimension.
#' @param mask_layer A \code{SpatRaster} containing categorical mask layers
#' (for example land–sea masks). Must have the same number of layers and geometry
#' as each element of \code{ras_list}.
#' @param ras_time A numeric vector giving the times or time-step labels for each
#' raster layer. Length must equal the number of layers in each raster.
#' @param buff_width Optional numeric. Buffer width in kilometres for
#' buffer-based extraction. Cannot be supplied together with \code{neigh}.
#' @param neigh Optional integer. Number of nearest neighbours to use for
#' neighbour-based extraction. Must be either \code{4} (rook) or \code{8}
#' (queen). Cannot be supplied together with \code{buff_width}.
#' @param window Optional numeric. Length of right-aligned temporal averaging
#' window in years. Must be \code{>= 10}. If \code{NULL}, raster times are treated
#' as instantaneous values. Use the same value from
#' \code{\link{summarise_TraCESahul}} if used.
#' @param prec Integer. Number of decimal places to round extracted variables.
#' @param summ_stat A summary function (for example \code{"mean"} or
#' \code{"median"}) applied when aggregating across neighbours or buffer areas.
#' @param dist_cut Optional numeric. Maximum snapping distance in kilometres.
#' Points further than \code{dist_cut} km from the nearest valid raster cell
#' are removed from the output.
#' @param cores Integer. Number of CPU cores to use for parallel extraction.
#' If \code{cores = 1L}, processing is sequential.
#' @param ... Additional arguments passed to internal methods.
#'
#' @return
#' A \code{data.table} where each row corresponds to one ID–Year combination.
#' Columns include:
#' \itemize{
#'   \item{\code{ID}}{ Record identifier.}
#'   \item{\code{Lon}, \code{Lat}}{ Final coordinates used after optional snapping.}
#'   \item{\code{Year}}{ Temporal index matched to the fossil’s age range.}
#'   \item{Environmental variables}{ One column per raster in \code{ras_list}.}
#'   \item{\code{LandSea}}{ Modal value from \code{mask_layer}.}
#' }
#'
#' Returned values are rounded to \code{prec} decimal places and duplicates are
#' removed.
#'
#' @examples
#' \dontrun{
#' # input rasters
#' ## from bioclim_TraCESahul
#' bio_sel <- names(bioclims) %in% c("bio01", "bio05", "bio06", "bio07",
#'                                   "bio12", "bio13", "bio14", "bio17",
#'                                   "bio18")
#' bioclim_subset <- bioclims[bio_sel]
#' mask_ras <- bioclim_subset[[1]]
#' mask_ras <- ifel(is.na(mask_ras), NA, 1)
#' mask_ras
#' # match the fossils to the data
#' data(ex_foss)
#' ex_foss <- terra::unwrap(ex_foss)
#' ex_foss
#' ex_foss <- data.frame(ID = 1:nrow(ex_foss),
#'                       Lat = terra::crds(ex_foss)[, 2],
#'                       Lon = terra::crds(ex_foss)[, 1],
#'                       Age = ex_foss$date_estimate,
#'                       AgeMin = ex_foss$ci_upper,
#'                       AgeMax = ex_foss$ci_lower)
#' foss_matched <- pair_obs(data = ex_foss,
#'                          ras_list = bioclim_subset,
#'                          mask_layer = mask_ras,
#'                          ras_time = terra::time(bioclim_subset$bio01),
#'                          window = 30, neigh = 8,
#'                          prec = 3, summ_stat = "mean",
#'                          dist_cut = 50, cores = 4L)
#' }
#'
#' @seealso
#' \code{\link{terra}}, \code{\link{future}}, \code{\link{future.apply}}
#'
#' @export
#'
pair_obs <- function(data, ras_list, mask_layer, ras_time, buff_width = NULL,
                     window = NULL, neigh = 8, prec = 2, summ_stat = "mean",
                     dist_cut = NULL, cores = 4L, ...) {
  # Required argument checks
  if (missing(data) || is.null(data)) {
    stop("'data' must be provided and must not be NULL.")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.")
  }
  req_cols <- c("ID", "Lat", "Lon", "Age", "AgeMin", "AgeMax")
  if (!all(req_cols %in% names(data))) {
    stop("Input 'data' must contain columns: ",
         paste(req_cols, collapse = ", "))
  }
  if (missing(ras_list) || is.null(ras_list)) {
    stop("'ras_list' must be provided and must not be NULL.")
  }
  if (!is.list(ras_list) || length(ras_list) == 0L) {
    stop("'ras_list' must be a non-empty list of SpatRaster objects.")
  }
  if (!all(vapply(ras_list, inherits, logical(1), "SpatRaster"))) {
    stop("'ras_list' must contain only terra::SpatRaster objects.")
  }
  if (missing(mask_layer) || is.null(mask_layer)) {
    stop("'mask_layer' must be provided and must not be NULL.")
  }
  if (!inherits(mask_layer, "SpatRaster")) {
    stop("'mask_layer' must be a terra::SpatRaster.")
  }
  if (missing(ras_time) || is.null(ras_time)) {
    stop("'ras_time' must be provided and must not be NULL.")
  }
  if (!is.numeric(ras_time)) {
    stop("'ras_time' must be numeric.")
  }
  # Raster / time alignment
  nl <- terra::nlyr(ras_list[[1]])
  if (length(ras_time) != nl) {
    stop("Length of 'ras_time' (", length(ras_time),
         ") must equal the number of layers in ras_list (", nl, ").")
  }
  if (terra::nlyr(mask_layer) != nl) {
    stop("'mask_layer' must have the same number of layers as ras_list.")
  }
  # Exactly one of neigh or buff_width must be supplied
  if (is.null(neigh) && is.null(buff_width)) {
    stop("Either 'neigh' or 'buff_width' must be supplied.")
  }
  if (!is.null(neigh) && !is.null(buff_width)) {
    stop("Supply only one of 'neigh' or 'buff_width'.")
  }
  # neigh validation
  if (!is.null(neigh)) {
    if (!is.numeric(neigh) || length(neigh) != 1L) {
      stop("'neigh' must be a single integer (4 or 8).")
    }
    if (!neigh %in% c(4, 8)) {
      stop("'neigh' must be 4 or 8.")
    }
  }
  # buff_width validation
  if (!is.null(buff_width)) {
    if (!is.numeric(buff_width) || length(buff_width) != 1L || buff_width <= 0) {
      stop("'buff_width' must be a single positive numeric value.")
    }
  }
  # prec validation
  if (!is.numeric(prec) || length(prec) != 1L || prec %% 1 != 0) {
    stop("'prec' must be a single integer.")
  }
  # optional warning about snapping
  if (!is.null(dist_cut)) {
    warning(
      sprintf("Lon/Lat points will be snapped to nearest points on rasters. Cutoff distance = %d km.", dist_cut),
      call. = FALSE, immediate. = TRUE
    )
  }
  # check raster geometry consistency
  stopifnot(check_geom(ras_list))
  stopifnot(terra::compareGeom(ras_list[[1]], mask_layer))
  # Apply the time to all the rasters
  ras_list <- lapply(ras_list, function(x) {
    terra::time(x) <- ras_time
    x
  })
  # make sure the mask has time info
  terra::time(mask_layer) <- ras_time
  # extract
  if (cores > 1L) {
    message(sprintf("\nExtracting data from rasters in parallel using %s cores and future.apply", cores))
    warning("all rasters must be wrapped before parallel extraction. Wrapping now.",
            call. = FALSE, immediate. = TRUE)
    ras_names <- names(ras_list)
    ras_list <- pbapply::pblapply(seq_along(ras_list), function(x) {
      terra::wrap(ras_list[[x]])
    })
    names(ras_list) <- ras_names
    mask_layer <- terra::wrap(mask_layer)
  } else {
    message("\nExtracting data from rasters sequentially")
  }
  wkt_proj <- 'PROJCS["Sahul_Lambert_Azimuthal",
  GEOGCS["GCS_WGS_1984",
  DATUM["D_WGS_1984",
  SPHEROID["WGS_1984",6378137.0,298.257223563]],
  PRIMEM["Greenwich",0.0],
  UNIT["Degree",0.0174532925199433]],
  PROJECTION["Lambert_Azimuthal_Equal_Area"],
  PARAMETER["False_Easting",0.0],
  PARAMETER["False_Northing",0.0],
  PARAMETER["Central_Meridian",135],
  PARAMETER["Latitude_Of_Origin",-20],
  UNIT["Meter",1.0]]'
  env_pairing <- parallel_env_match(data, window = window, n_cores = cores,
                                    ...)
  # Send warning if any records were removed.
  if (length(unique(env_pairing[["ID"]])) != length(unique(data[["ID"]]))) {
    warning(sprintf("\nThere were %s samples removed. They were > cutoff distance.",
                    length(unique(data[["ID"]])) - length(unique(env_pairing[["ID"]]))),
            immediate. = TRUE)
    warning("\nSamples removed: ", paste(
      unique(data[["ID"]])[!unique(data[["ID"]]) %in%
                             unique(env_pairing[["ID"]])], collapse = ", "),
      immediate. = TRUE)
  }
  # Round the extracted env variables
  cols <- names(env_pairing)[-c(1:4)] # exc ID, Lon, Lat, Year
  env_pairing[,(cols) := round(.SD, prec), .SDcols = cols]
  env_pairing <- unique(env_pairing, by = c("ID", cols))
  message("\nDone!")
  return(env_pairing)
}


#' check_geom
#'
#' @param x list of rasters
#'
#' @returns BOOL
#' @noRd
#'
check_geom <- function(x) {
  ok <- vapply(
    x[-1],
    function(r) terra::compareGeom(x[[1]], r, stopOnError = FALSE),
    logical(1)
  )
  if (!all(ok)) {
    stop("Rasters do not share identical geometry.")
  }
  invisible(TRUE)
}

#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom terra unwrap rast vect crs crds as.points project nearest distance buffer cellFromXY xyFromCell extract
#' @importFrom data.table copy data.table merge.data.table setcolorder rbindlist
#' @importFrom sf st_as_sf st_as_sfc st_geometry
#' @importFrom FNN get.knnx
#'
#' @noRd
#'
#' @keywords internal
parallel_env_match <- function(data, window, n_cores = 1L,...) {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = n_cores)
  idx <- seq_len(nrow(data))
  res_list <- future.apply::future_lapply(
    idx,
    FUN = function(i,...) {
      tryCatch(
        {
          sub <- data.table::copy(data)[i, ]

          # Extract the ras_time which fit in the CI of the fossil age
          ## Always rounds 'down' to oldest interval
          ras_sub <- get_time_indices(ras_time, sub$AgeMin, sub$AgeMax, win = window)
          if (is.null(ras_sub)) {
            return(NULL) # if no matching time index, return NULL
          }
          ## remove duplicates if narrow CI results in same layers being returned
          ras_sub <- ras_sub[!duplicated(ras_sub)]
          # If length of ras_sub is 1, will return index for age
          # if (length(ras_sub) ==  1) {
          #   ras_sub <- which.min(abs(ras_time - sub$Age))
          # }
          ## template raster for point alignment
          if (inherits(ras_list[[1]], "PackedSpatRaster")) {
            template_rast <- terra::unwrap(ras_list[[1]])[[ras_sub]][[1]]
            mask_layer <- terra::unwrap(mask_layer) # unwrap mask
          } else {
            template_rast <- ras_list[[1]][[ras_sub]][[1]]
          }
          ## convert sub to vect
          coords <- terra::vect(sub, geom = c("Lon", "Lat"), crs = "EPSG:4326")
          ## Snap to nearest point if snap is requested
          if (!is.null(dist_cut)) {
            ## project to requested projection
            coords <- terra::project(x = coords, y = wkt_proj)
            ## convert first layer of temporal subset raster to points
            ras_points <- terra::as.points(template_rast, values = FALSE, na.rm = TRUE)
            ras_points <- terra::project(ras_points, y = wkt_proj)
            snap_idx <- terra::nearest(coords, ras_points)[1, ] # only ever the first point
            # distance in km
            snap_dist <- as.vector(terra::distance(coords, snap_idx, unit = "km"))
            ## if snapping, replace the geometry of the coords
            if (snap_dist <= dist_cut) {
              ## cant directly edit a terra::vect, need to convert to sf first?
              coords <- sf::st_as_sf(coords)
              sf::st_geometry(coords) <- sf::st_as_sfc(sf::st_as_sf(ras_points[snap_idx$to_id, ]))
              coords <- terra::vect(coords)
              coords <- terra::project(coords, y = "EPSG:4326")
            } else {
              return(NULL)
            }
          }
          ## Extract either the neighbours or calculate a buffer
          if (!is.null(buff_width)) {
            coords <- terra::project(coords, y = wkt_proj)
            ## buffer in km
            nei <- terra::buffer(coords, width = buff_width*1000)
            ## back to WGS84
            nei <- terra::project(nei, "EPSG:4326")
            nr <- FALSE
            coords <- terra::project(coords, "EPSG:4326")
          } else if (!is.null(neigh)) {
            ## find the n nearest neighbours
            ## use this instead of terra::adjacent as raster method includes NA cells.
            ras_points <- terra::as.points(template_rast,
                                           values = FALSE, na.rm = TRUE)
            ## do the matching in projected coordinates
            nei <- c(FNN::get.knnx(
              data = terra::crds(terra::project(ras_points, y = wkt_proj)),
              query = terra::crds(terra::project(coords, y = wkt_proj)),
              k = neigh)$nn.index)
            nei <- terra::cellFromXY(template_rast,
                                     xy = terra::crds(ras_points[nei, ]))
            nr <- TRUE
          }
          # Extract the env variables and summarise according to summstat
          ext <- lapply(seq_along(ras_list), function(x, ...) {
            if (inherits(ras_list[[x]], "PackedSpatRaster")) {
              r <- terra::unwrap(ras_list[[x]])
            } else {
              r <- ras_list[[x]]
            }
            ex <- tryCatch(
              if (!nr) {
                # extract the mean value in the buffer
                t(terra::extract(r[[ras_sub]], y = nei, fun = "mean",
                               weights = TRUE, na.rm = TRUE, raw = FALSE,
                               ID = FALSE))
              } else if (nr) {
                # if only 1 layer, extract cells directly
                if (length(ras_sub) == 1L) {
                  apply(as.matrix(r[[ras_sub]][nei]), 2, summ_stat, na.rm = TRUE)
                } else {
                  # extract values from all neighbours, then calc mean at each time point
                  apply(terra::extract(r[[ras_sub]], nei), 2, summ_stat, na.rm = TRUE)
                }
              } else {
                stop()
              },
              error = function(e) e)
            if (inherits(ex, "error")) {
              ex <- data.table::data.table("ID" = sub$ID,
                               "Year" = as.vector(ras_time[ras_sub]),
                               "DT2" = rep(NA, length(ras_sub))
              )
              colnames(ex)[3] <- names(ras_list)[x]
              return(ex)
            } else {
              ex <- data.table::data.table("ID" = sub$ID,
                                           "Year" = as.vector(ras_time[ras_sub]),
                                           "DT2" = as.vector(ex))
              colnames(ex)[3] <- names(ras_list)[x]
              return(ex)
            }
          })
          # Extract modal value from mask
          m_ext <- tryCatch(
            if (!nr) {
              # extract the modal value in the buffer
              t(terra::extract(mask_layer[[ras_sub]], y = nei, fun = terra::modal,
                               weights = FALSE, na.rm = TRUE, raw = FALSE,
                               ID = FALSE))
            } else if (nr) {
              # extract values from all neighbours, then calc modal value at each time point
              apply(terra::extract(mask_layer[[ras_sub]], nei), 2, terra::modal, na.rm = TRUE)
            } else {
              stop()
            },
            error = function(e) e)
          if (inherits(m_ext, "error")) {
            m_ext <- data.table::data.table("ID" = sub$ID,
                                "Year" = as.vector(ras_time[ras_sub]),
                                "DT2" = rep(NA, length(ras_sub))
            )
            colnames(m_ext)[3] <- "LandSea"
          } else {
            m_ext <- data.table("ID" = sub$ID,
                                "Year" = as.vector(ras_time[ras_sub]),
                                "DT2" = as.vector(m_ext))
            colnames(m_ext)[3] <- "LandSea"
          }
          ext[[length(ext) + 1]] <- m_ext
          mergedDT <- Reduce(
            function(x, y)
              data.table::merge.data.table(x = x, y = y, by = c("ID", "Year"),
                    all.x = TRUE, all.y = TRUE),
            ext)
          ## add in the snapped coordinates
          mergedDT[, Lon := terra::crds(coords)[, 1]][, Lat := terra::crds(coords)[, 2]]
          ## re-order columns
          neworder <- c("ID", "Lon", "Lat", "Year", names(ras_list), "LandSea")
          data.table::setcolorder(mergedDT, neworder)
          return(mergedDT)
        },
        error = function(e) NULL
      )
    },
    future.seed = TRUE,
    future.globals = TRUE,
    future.packages = c("terra", "data.table", "sf", "FNN", "exactextractr")
  )
  res_list <- Filter(Negate(is.null), res_list)
  if (length(res_list) == 0L) {
    return(data.table::data.table())
  }
  out <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
  return(out)
}

#' @param ras_time ras_time
#' @param AgeMin minAge for record
#' @param AgeMax maxAge for record
#' @param win the window used for the averaging
#'
#' @returns vector
#' @noRd
#'
get_time_indices <- function(ras_time, AgeMin, AgeMax, win = NULL) {
  # no temporal smoothing, treat as exact years
  if (is.null(win)) {
    idx <- which(ras_time >= AgeMax & ras_time <= AgeMin)
    if (length(idx) == 0L) return(NULL)
    return(idx)  }
  # right-aligned window averaging
  if (win < 10) {
    stop("Window length must be >= 10 when smoothing is used.")
  }
  # window is right-aligned: [ras_time - win + 1, ras_time]
  layer_start <- ras_time - win + 1
  layer_end   <- ras_time
  # standard interval-overlap test
  idx <- which(pmax(layer_start, AgeMax) <= pmin(layer_end, AgeMin))
  if (length(idx) == 0L) {
    return(NULL)
  } else {
    return(idx)
  }
}
