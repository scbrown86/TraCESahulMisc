## code to prepare `ex_foss` dataset goes here

# this code assumes you have downloaded the example TraCE-Sahul data
library(terra)
library(virtualspecies)

dmon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

ex_data <- lapply(c("pr", "tasmax", "tasmin"), function(i) {
  r <- rast(list.files(path = "D:/TraCE-Sahul/",
                       recursive = TRUE,
                       pattern = sprintf("TraCE_22ka_downscaled_%s_1500_1990_biascorr.nc", i),
                       full.names = TRUE),
            lyrs = 5761:5880)*1
  time(r) <- seq(as.Date("1980-01-16"), by ="month", l = nlyr(r))
  r <- tapp(r, "month", mean, na.rm = TRUE)
  names(r) <- paste0(month.abb,"_", i)
  r})
ex_data[[1]] <- ex_data[[1]]*(86400) #kg m-2 s-1 --> mm/day
ex_data
# check of rainfall
# panel(ex_data[[1]], range = c(0, 150),
#       col = hcl.colors(100, "Roma"),
#       fill_range = TRUE)

ex_data <- rast(ex_data)

# Generate species
{set.seed(9621);
  realistic.sp <- generateRandomSp(ex_data,
                                   rescale = TRUE,
                                   plot = TRUE,
                                   realistic.sp = TRUE,
                                   convert.to.PA = TRUE)}
# realistic.sp

# Sampling of 'presence only' occurrences


{set.seed(8945); sp.locs <- terra::spatSample(x = realistic.sp$suitab.raster,
                                              size = 75, method = "weights",
                                              replace = TRUE, na.rm = TRUE,
                                              as.raster = FALSE, as.df = FALSE,
                                              as.points = TRUE,
                                              values = TRUE, xy = TRUE,
                                              cells = FALSE)
  }
sp.locs
# p.dens <- density(sp.locs$lyr.1, bw = "SJ",
#                   kernel = "biweight",
#                   from = 0, to = 1)
# plot(p.dens, xlim = c(0, 1))
# rug(sp.locs$lyr.1)

# plot(realistic.sp$suitab.raster, range = c(0, 1),
#      col = hcl.colors(100, "Batlow"),
#      fill_range = TRUE,
#      fun = function() points(sp.locs, pch = 19))

generate_random_dates <- function(n) {
  # n = number of samples
  # Generate smooth "time trajectory" using cumulative random walk
  base <- cumsum(rnorm(n, mean = 0, sd = 500))  # smooth variation
  base <- (base - min(base)) / (max(base) - min(base))  # rescale 0–1
  dates <- -20900 + base * (1850 + 20900)  # map to range -21000 to 1989

  # Generate small normally distributed confidence intervals (CI)
  ci_width <- abs(rnorm(n, mean = 100, sd = 50))  # average CI width ±100 yrs
  lower <- dates - ci_width / 2
  upper <- dates + ci_width / 2

  # Assemble output
  data.frame(
    sample_id = seq_len(n),
    date_estimate = round(dates, 0),
    ci_lower = round(lower, 0),
    ci_upper = round(upper, 0)
  )
}

{set.seed(123)
  sp.sample.years <- generate_random_dates(nrow(sp.locs))
  }
sp.sample.years <- sp.sample.years[order(sp.sample.years$date_estimate), ]
ex_foss <- cbind(sp.locs, sp.sample.years[, -1])
ex_foss$lyr.1 <- NULL
ex_foss
ex_foss <- terra::wrap(ex_foss)

# ex_foss
usethis::use_data(ex_foss, version = 3, overwrite = TRUE)
true_suit <- terra::wrap(realistic.sp$suitab.raster)
usethis::use_data(true_suit, version = 3, overwrite = TRUE)
