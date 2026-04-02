# TraCESahulMisc

**TraCESahulMisc** provides helper functions, workflows, and datasets for downloading, importing, processing, and analysing downscaled TraCE-Sahul palaeoclimate data in R. It supports a complete pipeline: acquiring (an example of) raw TraCE-Sahul files, importing them as `terra::SpatRaster` objects with correct metadata, generating monthly, seasonal, or annual summaries, deriving BIOCLIM variables, and pairing environmental rasters with fossil or observational point data.

![downscaling comparison](overview_image.png) 
_A comparison between the downscaled TraCE-Sahul and raw TraCE-21 data. The top row shows the effect of the downscaling on mean annual temperature, while the bottom row shows the effect on mean annual precipitation._

The package was designed specifically for researchers working with the TraCE-Sahul climate reconstructions, with the aim of automating the aggregation of palaeo-environmental data, and time-series environmental analyses across the Sahul region to be used in species distribution modelling.

## Installation

The easiest way to install the package is to use `remotes` as below.

``` r
remotes::install_github("scbrown86/TraCESahulMisc", build_vignettes = FALSE)
```

The package contains a pre-built vignette that *should* not be re-built on install so please make sure you set `build_vignettes = FALSE` when installing.

## Vignette

The package constains a small vignette showing some of the functionality. 

It can be found oneline [here](https://scbrown86.github.io/TraCESahulMisc/TraCESahulMisc_workflow.html), or can be viewed in RStudio as below

```r
vignette("TraCESahulMisc_workflow")
```

## Key Features

The package helps to automate some basic tasks that are common with palaeo climate reconstructions. It has been built specifically to work the TraCE-Sahul dataset, but *may* work with other datasets provided they have a time attribute (see the [Terra package](https://rspatial.github.io/terra/reference/time.html) for details)

-   Download an example TraCE-Sahul climate dataset.
-   Import multi-layer NetCDF files as annotated `SpatRaster` objects.
-   Summarise TraCE datasets to monthly, seasonal, or annual climatologies.
-   Compute BIOCLIM variables.
-   Pair climate rasters with fossil or observational point data.
-   Includes two example datasets.

See the package vignette for a full worked example on using the dataset

## License

MIT License.
