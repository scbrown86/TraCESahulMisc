## code to prepare `ex_foss` dataset goes here

# this code assumes you have downloaded the example TraCE-Sahul data
library(terra)
library(virtualspecies)
library(fastbioclim)

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
ex_data <- fastbioclim::derive_bioclim(
  bios = 1:19,
  tmin = ex_data[[3]],
  tmax = ex_data[[2]],
  prcp = ex_data[[1]],
  overwrite = TRUE)
plot(ex_data)

# check of rainfall
# panel(ex_data[[1]], range = c(0, 150),
#       col = hcl.colors(100, "Roma"),
#       fill_range = TRUE)

# ex_data <- rast(ex_data)*1

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

# Timesteps
library(data.table)

# Create the paleo timesteps
n <- 25860
n_decades <- n / 12  # 2155
steps <- seq(1, by = 120, length.out = n_decades)
ends <- steps + 119
mids_months <- (steps + ends) / 2
mids_years <- mids_months / 12
time_bp <- 22000 - mids_years
time_steps <- data.table(
  layerID = 1:n,
  dec = rep(1:n_decades, each = 12),
  Month = rep(1:12, times = n_decades),
  YearsBP = rep(ceiling(time_bp), each = 12),
  file_step = c(rep(1:5, each = 4812), rep(6, times = 1800)))
# time_steps
time_steps[, YearsCE := 1950 - YearsBP]
time_steps[, dec_year := YearsCE + (Month - 1 + 15.5 / 30.4375) / 12]
time_steps[, `:=`(StartYearCE = YearsCE - 5,
                  EndYearCE = YearsCE + 4,
                  StartYearBP = YearsBP + 5,
                  EndYearBP = YearsBP - 4)]
# Post 1500 timesteps
monthly_years <- rep(1500:1989, each = 12)
monthly_months <- rep(1:12, times = 490)
monthly_rows <- data.table(
  layerID = (n + 1):(n + 490 * 12),
  dec = NA_integer_,
  Month = monthly_months,
  YearsBP = 1950L - monthly_years,
  YearsCE = monthly_years,
  file_step = NA_integer_)
monthly_rows[, dec_year := YearsCE + (Month - 1 + 15.5 / 30.4375) / 12]
monthly_rows[, `:=`(StartYearCE = YearsCE,
                    EndYearCE = YearsCE,
                    StartYearBP = YearsBP,
                    EndYearBP = YearsBP)]

# combine
TraCESahul_timesteps <- rbind(time_steps, monthly_rows)
usethis::use_data(TraCESahul_timesteps, version = 3, overwrite = TRUE)
