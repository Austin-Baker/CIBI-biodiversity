library(raster)
library(dplyr)
library(stringr)
library(units)
library(terra)
library(parallel)

vegtype_path <- "vegtype_ecoregion_maps/"
veg_files <- list.files(vegtype_path)

map_aggregator <- function(f){
test <- raster(paste0(vegtype_path, f))
test_2 <- aggregate(test, fact=33)

writeRaster(test_2, paste0("veg_maps_agg/", f), overwrite = TRUE)
}

completed_jobs <- list.files("veg_maps_agg")
bins_to_do <- veg_files[!veg_files %in% completed_jobs]

mclapply(bins_to_do, FUN = map_aggregator, mc.cores = 10)