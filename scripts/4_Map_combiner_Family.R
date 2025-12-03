library(sf)
library(raster)
library(dplyr)
library(stringr)
library(units)
library(terra)
library(tigris)
library(data.table)

## Load BOLD data with ecoregions and landfire added
BOLD_data_lf <- readRDS("data/BOLD_data_lf.rds")

# Load projected California plot
CA_lf_p <- rast("data/CA_lf_p.tif")

CA_lf_p_2 <- aggregate(CA_lf_p, fact=33)

## list BINs to group together
BIN_Order <- BOLD_data_lf %>%
  select(BIN, Order, Family) %>%
  filter(!is.na(BIN)) %>%
  unique() %>%
  arrange(BIN) %>%
  group_by(BIN) %>%
  summarize_all(~paste(unique(na.omit(.)), collapse = ',')) %>%
  mutate(Order_Family = paste0(Order, "_", Family)) %>%
  filter(!str_detect(Order_Family, ","))

# fix file names
vegtype_path <- "veg_maps_agg/"
veg_files <- list.files(vegtype_path, full.names=TRUE)
for (file in veg_files) {
  file_name <- basename(file)
  new_name <- gsub(":", "_", file_name, fixed = TRUE)
  new_file_path <- file.path(dirname(file), new_name)
  file.rename(file, new_file_path)
}
# loop for every family to combine BINs
Family_loop <- unique(BIN_Order$Order_Family)

for(f in Family_loop){
  BIN_Family <- BIN_Order %>%
    filter(Order_Family==f) %>%
    mutate(BIN = str_replace(BIN, ":", "\\_"))
  
  f_paths <- paste0("veg_maps_agg/", BIN_Family$BIN, ".tif")
  
  # Load, crop, and extend rasters to match CA_lf_p's extent
  rasters_fixed <- lapply(f_paths, function(f) {
    r <- rast(f)
    #r <- crop(r, CA_lf_p)    # Crop to common extent
    #r <- extend(r, CA_lf_p)  # Extend to ensure exact same extent
    r <- terra::project(r, crs(CA_lf_p_2))  # Ensure same CRS
    r <- terra::resample(r, CA_lf_p_2, method = "near")
    return(r)
  })
  
  # Stack the fixed rasters
  raster_stack <- rast(rasters_fixed)
  
  # Sum the rasters, treating NAs as 0
  combined_raster <- app(raster_stack, sum, na.rm = TRUE)
  
  # Save the final raster
  writeRaster(combined_raster, paste0("Family_maps/", f, ".tif"), overwrite = TRUE)
  
}
  