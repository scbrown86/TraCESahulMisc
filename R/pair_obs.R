pair_obs <- function(data, ras_list, mask_layer, ras_time, buff_width = NULL,
                     neigh = 8, prec = 2, summ_stat = "mean", dist_cut = NULL,
                     cores = 4L, ...) {
  #### Inputs:
  ####   * data = data.frame/data.table with columns: ID, Lat, Lon, Age, AgeMin, AgeMax.
  ####   Columns can be in any order.
  ####   * ras_list = a list of the environmental rasters
  ####   * mask_layer = raster stack of masks. Needs to have same number of layers as the rasters in ras_list
  ####   * ras_time = seq of years/steps for the environmental rasters
  ####   * buff_width = buffer width in km (e.g. 75 = 75,000m)
  ####   * neigh = Number of neighbours to extract values from. 8 = queens criterion, 4 = rook.
  ####   * cores = number of cores to run in parallel
  ####   * prec = precision of rounded outputs
  ####   * summstat = summary stat for outputs
  ####   * dist_cut = cutoff distance for snapping tolerance in kilometres
  ####   * EPSG = EPSG code/CRS/WKT used for projecting locations before snapping points. Default = Plate Carree
  require(raster)
  require(FNN)
  require(exactextractr)
  require(sf)
  require(lwgeom)
  require(data.table)
  require(doSNOW)
  if(c(missing(neigh) & missing(buff_width)) | c(is.null(neigh) & is.null(buff_width))) {
    stop("neigh or buff_width are not supplied!")
  }
  if(!is.null(neigh) & !is.null(buff_width)) {
    stop("supply only one of neigh or buff_width")
  }
  if(!is.null(neigh) & !neigh %in% c(4,8)) {
    stop("neigh must be 4 or 8 for the meantime...")
  }
  message("Preparing data...")
  # Check that data has columns ID, Lat, Long, Age, AgeMin, AgeMax
  if(
    sum(sapply(colnames(data), function(x) x %in% c("ID", "Lat", "Lon", "Age", "AgeMin", "AgeMax")))
    != 6) {
    stop("Data has the wrong column names")
  }
  # make sure precision is int/numeric
  if (!prec %% 1 == 0) {
    stop("'Prec' must be supplied as an integer")
  }
  # warn that points will be snapped if passed
  if (!is.null(dist_cut)) {
    warning(
      sprintf("Lon/Lat points will be snapped to nearest points on rasters. Cutoff dist = %d km. Points further than the cutoff will be removed from output", dist_cut),
      call. = FALSE, immediate. = TRUE)
  }
  # Check the rasters in the list are identical
  stopifnot(check_geom(ras_list))
  # Check the mask is identical to the rasters
  stopifnot(terra::compareGeom(ras_list[[1]], mask_layer))
  stopifnot(terra::nlyr(mask_layer) == terra::nlyr(ras_list[[1]]))
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
    ras_list <- pbapply::pblapply(seq_along(ras_list), function(x) {
      terra::wrap(ras_list[[x]])
    })
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
  env_pairing <- parallel_env_match()


  env_pairing <- foreach(i = seq_len(nrow(data)), .options.snow = opts,
                         .packages = c("raster", "data.table", "sp", "sf", "FNN", "exactextractr"),
                         .combine = function(...) rbindlist(list(...)),
                         #.verbose = TRUE,
                         .multicombine = TRUE,
                         .inorder = TRUE,
                         .errorhandling = "remove") %dopar% {






                           # Extract modal value from mask
                           m_ext <- tryCatch(
                             if (!nr) {
                               # extract the modal value in the buffer
                               t(exact_extract(mask_layer[[ras_sub]], y = nei,
                                               fun = "mode"))
                             } else if (nr) {
                               # extract values from all neighbours, then calc modal value at each time point
                               apply(t(extract(mask_layer[[ras_sub]], y = nei, df = TRUE))[-1,],
                                     1, raster::modal, na.rm = TRUE)
                             } else {
                               stop()
                             },
                             error = function(e) e)
                           if(inherits(m_ext, "error")) {
                             m_ext <- data.table("ID" = sub$ID,
                                                 "Year" = as.vector(ras_time[ras_sub]),
                                                 "DT2" = rep(NA, length(ras_sub))
                             )
                             colnames(m_ext)[3] <- "SeaLandIce"
                           } else {
                             m_ext <- data.table("ID" = sub$ID,
                                                 "Year" = as.vector(ras_time[ras_sub]),
                                                 "DT2" = as.vector(m_ext))
                             colnames(m_ext)[3] <- "SeaLandIce"
                           }
                           ext[[length(ext) + 1]] <- m_ext
                           mergedDT <- Reduce(
                             function(x, y)
                               merge(x = x, y = y, by = c("ID", "Year"),
                                     all.x = TRUE, all.y = TRUE),
                             ext)
                           ## add in the snapped coordinates
                           mergedDT[, Lon := st_coordinates(coords)[1]][, Lat := st_coordinates(coords)[2]]
                           ## re-order columns
                           neworder <- c("ID", "Lon", "Lat", "Year", names(ras_list), "SeaLandIce")
                           setcolorder(mergedDT, neworder)
                           return(mergedDT)
                         }
  registerDoSEQ(); stopCluster(cls); invisible(gc())
  # Send warning if any records were removed.
  if (length(unique(env_pairing$ID)) != length(unique(data$ID))) {
    warning(sprintf("\nThere were %s samples removed. They were > cutoff distance.", length(unique(data$ID)) - length(unique(env_pairing$ID))),
            immediate. = TRUE)
    warning("\nSamples removed: ", paste(
      unique(data$ID)[!unique(data$ID) %in% unique(env_pairing$ID)], collapse = ", "),
      immediate. = TRUE)
  }
  # Round the extracted env variables
  cols <- names(env_pairing)[-c(1:2)]
  env_pairing[,(cols) := round(.SD, prec), .SDcols = cols]
  env_pairing <- unique(env_pairing, by = c("ID", cols))
  return(env_pairing)
  message("\nDone!")
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
#' @import data.table
#' @keywords internal
parallel_env_match <- function(data, n_cores = cores) {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future::multisession, workers = n_cores)

  idx <- seq_len(nrow(data))

  res_list <- future.apply::future_lapply(
    idx,
    FUN = function(i) {
      tryCatch(
        {
          sub <- data.table::copy(data)[i, ]

          # Extract the ras_time which fit in the CI of the fossil age
          ## Always rounds 'down' to oldest interval
          ras_sub <- seq(which.min(abs(ras_time - sub$AgeMax)),
                         which.min(abs(ras_time - sub$AgeMin)), 1)
          ## remove duplicates if narrow CI results in same layers being returned
          ras_sub <- ras_sub[!duplicated(ras_sub)]
          # If length of ras_sub is 1, will return index for age
          if (length(ras_sub) ==  1) {
            ras_sub <- which.min(abs(ras_time - sub$Age))
          }
          ## template raster for point alignment
          if (inherits(ras_list[[1]], "PackedSpatRaster")) {
            template_rast <- terra::unwrap(ras_list[[1]])[[ras_sub]][[1]]
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

          #### FROM HERE ####

          # Extract the env variables and summarise according to summstat
          ext <- lapply(seq_along(ras_list), function(x, ...) {
            r <- ras_list[[x]]
            ex <- tryCatch(
              if (!nr) {
                # extract the mean value in the buffer
                t(exact_extract(r[[ras_sub]], y = nei,
                                weights = area(r[[ras_sub]][[1]]),
                                fun = "weighted_mean"))
              } else if (nr) {
                # if only 1 layer, extract cells directly
                if (length(ras_sub == 1)) {
                  apply(as.matrix(r[[ras_sub]][nei]), 2, summ_stat, na.rm = TRUE)
                } else {
                  # extract values from all neighbours, then calc mean at each time point
                  apply(t(raster::extract(r[[ras_sub]], y = nei, df = TRUE))[-1, ],
                        1, summ_stat, na.rm = TRUE)
                }
              } else {
                stop()
              },
              error = function(e) e)
            if(inherits(ex, "error")) {
              ex <- data.table("ID" = sub$ID,
                               "Year" = as.vector(ras_time[ras_sub]),
                               "DT2" = rep(NA, length(ras_sub))
              )
              colnames(ex)[3] <- names(ras_list)[x]
              return(ex)
            } else {
              ex <- data.table("ID" = sub$ID,
                               "Year" = as.vector(ras_time[ras_sub]),
                               "DT2" = as.vector(ex))
              colnames(ex)[3] <- names(ras_list)[x]
              return(ex)
            }
          })
        },
        error = function(e) NULL
      )
    },
    future.packages = c("terra", "data.table", "sf", "FNN", "exactextractr")
  )

  res_list <- Filter(Negate(is.null), res_list)

  if (length(res_list) == 0L) {
    return(data.table::data.table())
  }

  out <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
  return(out)
}

