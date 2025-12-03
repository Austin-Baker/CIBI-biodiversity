library(sf)
library(dplyr)
library(stringr)
library(units)
library(terra)
library(tigris)
library(data.table)
library(parallel)

#### Load Data ####

# Load EPA ecoregions
EPA4 <- st_read("data/ca_eco_l4")

EPA4_tf <- EPA4 %>% 
  dplyr::select(US_L4NAME) %>% 
  st_transform(3310) %>%
  rename(shape_geometry = geometry)

# Load projected California plot
CA_lf_p <- rast("data/CA_lf_p.tif")

CA_levels <- readRDS("data/CA_levels.rds")

## Load BOLD data with ecoregions and landfire added
BOLD_data_lf <- readRDS("data/BOLD_data_lf.rds")


### filter data for unique species+lat+lon combinations, filter out missing data ###
BOLD_data_unique <- unique(BOLD_data_lf[,.(BIN,Lat,Lon,US_L3NAME,US_L4NAME,EVT_NAME)])
species_list <- BOLD_data_unique[,.N, by=.(BIN)]
species_list_clean <- filter(species_list, BIN!="NA")
species_list_clean_big <- species_list_clean %>%
  filter(N>29)

#### Loop for species distributions ####
# Beginning of loop
map_creator <- function(species_name){
  
  print(species_name)
  # filter for a species name
  sp_data <- BOLD_data_unique[BIN == species_name, .(Lat, Lon,US_L3NAME,US_L4NAME,EVT_NAME)]
  # create range based on L4 ecoregion
  L4_range_sp <- EPA4_tf %>% 
    filter(US_L4NAME %in% sp_data$US_L4NAME) %>%
    st_union()
  
  saveRDS(L4_range_sp, paste0("./ecoregion_maps/", species_name, ".rds"))
  
  # crop landfire dataset by polygon
  L4_polygon_sp <- EPA4_tf %>% 
    filter(US_L4NAME %in% sp_data$US_L4NAME) %>%
    rename(geometry = shape_geometry) %>%
    dplyr::select(-US_L4NAME) %>%
    st_transform(crs(CA_lf_p))
  
  sp_landfire_raster <- terra::crop(CA_lf_p, L4_polygon_sp)
  
  sp_landfire_raster2 <- terra::mask(sp_landfire_raster, L4_polygon_sp)
  
  lf_names_sp <- CA_levels[[1]]
  
  lf_not_sp <- lf_names_sp$Value[!lf_names_sp$EVT_NAME %in% sp_data$EVT_NAME]
  
  sp_landfire_raster3 <- subst(sp_landfire_raster2, lf_not_sp, NaN) 
  
  sp_raster_values <- unique(sp_landfire_raster3)
  
  reclassify_matrix <- sp_raster_values %>%
    as.data.frame() %>%
    rename(Value = EVT_NAME) %>%
    left_join(lf_names_sp) %>%
    mutate(new_val = ifelse(!is.na(Value), 1, NA)) %>%
    dplyr::select(Value, new_val) %>%
    bind_rows(data.frame(Value = 9995, new_val = NA)) %>%
    as.matrix()
  
  masked_polygon_sp <- terra::classify(sp_landfire_raster3, reclassify_matrix)
  
  writeRaster(masked_polygon_sp, file=paste0("vegtype_ecoregion_maps/", species_name,".tif"))
  
}
completed_jobs <- list.files("vegtype_ecoregion_maps")
completed_jobs_names <- str_replace(completed_jobs, ".tif", "")
bins_to_do <- species_list_clean_big$BIN[!species_list_clean_big$BIN %in% completed_jobs_names]

mclapply(bins_to_do, FUN = map_creator, mc.cores = 1)

# end of loop

# fix file names
#ecogregion_path <- "ecoregion_maps/"
#eco_files <- list.files(ecogregion_path, full.names=TRUE)
#for (file in eco_files) {
#  file_name <- basename(file)
#  new_name <- gsub(":", "_", file_name, fixed = TRUE)
#  new_file_path <- file.path(dirname(file), new_name)
#  file.rename(file, new_file_path)
#}
#vegtype_path <- "vegtype_ecoregion_maps/"
#veg_files <- list.files(vegtype_path, full.names=TRUE)
#for (file in veg_files) {
#  file_name <- basename(file)
#  new_name <- gsub(":", "_", file_name, fixed = TRUE)
#  new_file_path <- file.path(dirname(file), new_name)
#  file.rename(file, new_file_path)
#}


### Plot map by stacking rasters together

# List all raster files in the folder
#veg_files2 <- list.files(vegtype_path, pattern = "\\.tif$", full.names = TRUE)

# Select only the first 3 rasters for testing
#test_files <- veg_files2[1:10]

# Load, crop, and extend rasters to match CA_lf_p's extent
#rasters_fixed <- lapply(veg_files2, function(f) {
#  r <- rast(f)
#  r <- crop(r, CA_lf_p)    # Crop to common extent
#  r <- extend(r, CA_lf_p)  # Extend to ensure exact same extent
#  return(r)
#})

# Stack the fixed rasters
#raster_stack <- rast(rasters_fixed)

# Sum the rasters, treating NAs as 0
#combined_raster <- app(raster_stack, sum, na.rm = TRUE)

# Save the final raster
#writeRaster(combined_raster, "combined_vegtype_ecoregion_raster.tif", overwrite = TRUE)

# Plot the result
#plot(combined_raster)

