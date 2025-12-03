library(sf)
library(raster)
library(dplyr)
library(stringr)
library(units)
library(terra)
library(tigris)
library(data.table)


# Load projected California plot
CA_lf_p <- rast("data/CA_lf_p.tif")

CA_lf_p_2 <- aggregate(CA_lf_p, fact=33)

  f_paths <- list.files("Order_maps/", full.names = TRUE)
  
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
  writeRaster(combined_raster, "Richness_map.tif", overwrite = TRUE)
  
  # Load and plot map
  richness_map <- rast("Richness_map.tif")
  
  total_BIN_richness <- ggplot() +
    geom_spatraster(data=richness_map, aes(fill=sum))+
    theme_cowplot()+
    scale_fill_viridis_c(name="BIN \nRichness", na.value = NA)+
    theme(axis.line= element_blank(), 
          axis.ticks = element_blank(),
          axis.text = element_blank())+
    annotation_scale()+
    annotation_north_arrow(location= "tl", height=unit(1, "cm"), width=unit(1, "cm"))
  
  
  
  