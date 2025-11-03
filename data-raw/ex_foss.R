## code to prepare `ex_foss` dataset goes here
library(virtualspecies)
library(geodata)

dmon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

worldclim <- lapply(c("prec", "tmax", "tmin"), function(i) {
  r <- rast(list.files(path = "C:/Users/Stu/Downloads/wc2.1_country/",
                       pattern = sprintf("%s.*\\.tif$", i), full.names = TRUE),
            win = ext(112.5, 154, -45, -10)
  )*1
  names(r) <- paste0(month.abb,"_", i)
  r})
worldclim[[1]] <- worldclim[[1]]/dmon #mm/day
# worldclim

# check of rainfall
# panel(worldclim[[1]], range = c(0, 10),
#       col = hcl.colors(100, "Roma"),
#       fill_range = TRUE)

aus_stack <- rast(worldclim)

# Generate species
{set.seed(9621);
  realistic.sp <- generateRandomSp(aus_stack,
                                   rescale = TRUE,
                                   plot = TRUE,
                                   realistic.sp = TRUE,
                                   convert.to.PA = FALSE)}
# realistic.sp

{set.seed(8945); sp.locs <- terra::spatSample(x = realistic.sp$suitab.raster,
                                              size = 150, method = "weights",
                                              replace = TRUE, na.rm = TRUE,
                                              as.raster = FALSE, as.df = FALSE,
                                              as.points = TRUE,
                                              values = TRUE, xy = TRUE,
                                              cells = FALSE)}
# sp.locs
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
  sp.sample.years <- generate_random_dates(nrow(sp.locs))}
sp.sample.years <- sp.sample.years[order(sp.sample.years$date_estimate), ]
# sp.sample.years

ex_foss <- cbind(sp.locs, sp.sample.years[, -1])
ex_foss$lyr.1 <- NULL
ex_foss <- terra::wrap(ex_foss)
# ex_foss
usethis::use_data(ex_foss, version = 3, overwrite = TRUE)
