## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----loadLibs-----------------------------------------------------------------
library(TraCESahulMisc)
library(terra)
library(randomForest)

# only for plotting
library(rnaturalearth)
conts <- vect(rnaturalearth::ne_coastline(scale = 50))
conts <- crop(conts, ext(125, 135, -13, -8))

## ----downloadData-------------------------------------------------------------
base_dir <- "D:/TraCE-Sahul"
# download_trace_data(base_dir) # Run once only

## ----importData---------------------------------------------------------------
vars <- c("pr", "tasmax", "tasmin")

load_tracesahul <- function(base, var) {
  f1 <- file.path(base, var,
                  sprintf("TraCE_22ka_downscaled_%s_decadal_21k_1500CE_biascorr.nc", var))
  f2 <- file.path(base, var,
                  sprintf("TraCE_22ka_downscaled_%s_1500_1990_biascorr.nc", var))
  import_TraCESahul(c(f1, f2))
}
sahul <- lapply(vars, function(v) load_tracesahul(base_dir, v))
names(sahul) <- vars
sahul

## ----timeAgg------------------------------------------------------------------
pr_monthly_win <- summarise_TraCESahul(sahul$pr, type = "monthly",
                                       sumfun = "mean", window = 30)
tasmax_monthly_win <- summarise_TraCESahul(sahul$tasmax, type = "monthly",
                                           sumfun = "mean", window = 30)
tasmin_monthly_win <- summarise_TraCESahul(sahul$tasmin, type = "monthly",
                                           sumfun = "mean", window = 30)

lyr_names <- rep(month.abb, nlyr(pr_monthly_win)/12)
names(pr_monthly_win) <- paste0(lyr_names, "_pr")
names(tasmax_monthly_win) <- paste0(lyr_names, "_tasmax")
names(tasmin_monthly_win) <- paste0(lyr_names, "_tasmin")

## ----bioclim, eval = FALSE, echo=TRUE-----------------------------------------
# bioclims <- bioclim_TraCESahul(
#   tasmax = tasmax_monthly_win,
#   tasmin = tasmin_monthly_win,
#   pr = pr_monthly_win,
#   # store in a folder in basedir. Not recommended
#   outdir = paste0(base_dir,"/bioclims"),
#   bioclims = c(1, 4, 5, 6, 7, 12, 13, 14, 17, 18), # default is to use 1:19
#   collate = TRUE)
# bioclims

## ----bioClimLoad, echo = FALSE, eval = TRUE-----------------------------------
# This section will run the bioclim extraction only if the output sds doesnt already exist
if (!file.exists(file.path(paste0(base_dir,"/bioclims"), "bioclims_TraCESahul.RDS"))) {
  bioclims <- bioclim_TraCESahul(
  tasmax = tasmax_monthly_win,
  tasmin = tasmin_monthly_win,
  pr = pr_monthly_win,
  # store in a folder in basedir. Not recommended
  outdir = paste0(base_dir,"/bioclims"),
  bioclims = c(1, 4, 5, 6, 7, 12, 13, 14, 17, 18), # default is to use 1:19
  collate = TRUE)
} else {
  bioclims <- terra::unwrap(readRDS(file.path(paste0(base_dir,"/bioclims"), "bioclims_TraCESahul.RDS")))
}

## ----plotBio, echo=FALSE,fig.align='center',fig.dpi=320, fig.height=8, fig.width=9, out.width = "100%"----
n <- nlyr(bioclims$bio01)
panel(bioclims$bio01[[round(quantile(1:n, probs = seq(0, 1, length.out = 6)))]],
      fill_range = TRUE, range = c(10, 30),
      col = map.pal("bgyr", n = 100),
      fun = function() lines(conts, lwd = 1, col  ="#000000"),
      box = TRUE,
      background = "grey90",
      nc = 2,
      plg = list(
        legend = "Suitability",
        title = "BIO01 (°C)",
        title.x = 136.75,
        title.y = -10.7,
        title.srt = 90,
        title.cex = 1.1,
        digits = 0,
        bty = "n",
        size = c(1, 1),
        at = seq(10, 30, l = 5),
        tick = "through"
      ),
      pax = list(
        retro = TRUE
      ),
      loc.main = "bottomleft",
      grid = TRUE)

## -----------------------------------------------------------------------------
# make a mask using bio01
mask_ras <- bioclims$bio01
mask_ras <- ifel(is.na(mask_ras), NA, 1)

## ----exFoss-------------------------------------------------------------------
data(ex_foss)
ex_foss <- terra::unwrap(ex_foss)

ex_foss <- data.table::data.table(
  ID = 1:nrow(ex_foss),
  Lat = crds(ex_foss)[,2],
  Lon = crds(ex_foss)[,1],
  Age = ex_foss$date_estimate,
  AgeMin = ex_foss$ci_upper,
  AgeMax = ex_foss$ci_lower
)


## ----showexFoss, echo=FALSE---------------------------------------------------
head(ex_foss)

## ----fossMatch----------------------------------------------------------------
foss_matched_bioclim <- pair_obs(
  data = ex_foss,
  ras_list = bioclims,
  mask_layer = mask_ras,
  ras_time = terra::time(bioclims$bio01),
  window = 30, # 30-year window
  neigh = 8, # 8-neighbours (rook)
  prec = 3, # 3 decimal places
  summ_stat = "mean",
  dist_cut = 50, # km
  # parallel is disabled and will be set to 1 if any other value is given
  cores = 1L, 
  buff_width = NULL
)

## ----extractBG----------------------------------------------------------------
# which timesteps from the TraCE-Sahul data are included in our fossil matched data
tidx <- which(terra::time(bioclims$bio01) %in% foss_matched_bioclim$Year)

# extract bg points.
## seed is set for reproducability
{set.seed(84564)
bg_points <- lapply(tidx, function(t) {
  nsamps <- 50
  preds <- rast(lapply(seq_along(bioclims), function(i,...) bioclims[[i]][[t]]))
  names(preds) <- sapply(strsplit(names(preds), "_"), "[", 1)
  data.table::setDT(terra::spatSample(preds, nsamps, na.rm = TRUE))
})}

# collate the sample and round
bg_points <- data.table::rbindlist(bg_points)
cols <- names(bg_points)
bg_points[,(cols) := round(.SD, 3), .SDcols = cols]

## ----quickBG_look, echo = FALSE-----------------------------------------------
head(bg_points[,1:5])

## ----buildRF------------------------------------------------------------------
training <- data.table::rbindlist(list(
  data.table::copy(foss_matched_bioclim)[, -c("ID", "Lon", "Lat", "Year", "LandSea")][, Occ := "1"],
  data.table::copy(bg_points)[, Occ := "0"]
))

training[, Occ := factor(Occ, levels = c("0","1"))]

prNum <- as.numeric(table(training$Occ)["1"])
spsize <- c("0" = prNum, "1" = prNum)

rf_dws <- randomForest(
  Occ ~ .,
  mtry = 5,
  data = training,
  ntree = 1000,
  importance = TRUE,
  sampsize = spsize,
  replace = TRUE)

## ----rfSummary, echo=FALSE----------------------------------------------------
rf_dws

## ----makePreds----------------------------------------------------------------
# time steps to extract at
ts <- c(1, 65, 422, 507, 730)

# HS predictions
preds_stack <- rast(lapply(ts, function(t) {
  r <- rast(lapply(seq_along(bioclims), function(i) {
    bioclims[[i]][[t]]
  }))
  # Give the correct names to the layers in each stack
  names(r) <- sapply(strsplit(names(r), "_"), "[", 1)
  # predict prob of occurence
  pr <- predict(r, rf_dws, type = "prob", index = 2)
  names(pr) <- terra::time(bioclims$bio01)[t]
  return(pr)
}))

preds_stack

## ----getTrueSuit--------------------------------------------------------------
# get the true_suit raster
data(true_suit)
true_suit <- unwrap(true_suit)
names(true_suit) <- "true suitability"

# add the true_suit to the preds_stack
preds_stack <- c(preds_stack, true_suit)
preds_stack

## ----plotPreds, echo=FALSE, fig.dpi=320, fig.align='center',fig.height=8, fig.width=9,out.width = "100%"----
panel(preds_stack, 
      range = c(0, 1), fill_range = TRUE,
      fun = function() lines(conts, lwd = 1, col  ="#000000"),
      box = TRUE,
      background = "grey90",
      nc = 2,
      plg = list(
        legend = "Suitability",
        title = "Suitability",
        title.x = 136.75,
        title.y = -10.7,
        title.srt = 90,
        title.cex = 1.1,
        digits = 2,
        bty = "n",
        size = c(1, 1),
        at = seq(0, 1, l = 6),
        tick = "through"
      ),
      pax = list(
        retro = TRUE
      ),
      loc.main = "bottomleft",
      grid = TRUE)

