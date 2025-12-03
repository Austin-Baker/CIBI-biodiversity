### Data Preparation ###

library(sf)
library(raster)
library(dplyr)
library(stringr)
library(ggplot2)
library(units)
library(cowplot)
library(terra)
library(tigris)

setwd("~/Desktop/NHMLAC/Analyses/Insects_California")
#### Load Data ####

## load landfire raster
lf_raster <- rast("data/raw_data/Land_fire/LF2023_EVT_240_CONUS_2/Tif/LC23_EVT_240.tif") 
plot(lf_raster) # this appears to be a better format than raster::raster

## California borders

#california <- read_sf("data/raw_data/ca-state-boundary/") %>%
#  st_transform(crs(lf))
california <- tigris::states() %>%
  filter(STUSPS == "CA")

CA_lf <- crop(lf_raster, extent(st_transform(california, crs(lf_raster))))
plot(CA_lf)
writeRaster(CA_lf, "data/CA_landfire_raster.tif") # open from here after running first time

## land fire categories
#main_land_type_df <- readRDS("data/land_fire_1km.rds")

## load California grid
#california_grid_df <- readRDS("data/clean_data/california_grid_1km.rds") %>%
  #st_transform(3310)

# load insect data
BOLD_data <- fread("data/250104_all_insects.csv")

# add ecoregions to BOLD datasheet
EPA3 <- st_read("lands/EPA_Ecoregions_California/ca_eco_l3")

EPA3_tf <- EPA3 %>% 
  dplyr::select(US_L3NAME) %>% 
  st_transform(3310) %>%
  rename(shape_geometry = geometry)

EPA4 <- st_read("lands/EPA_Ecoregions_California/ca_eco_l4")

EPA4_tf <- EPA4 %>% 
  dplyr::select(US_L4NAME) %>% 
  st_transform(3310) %>%
  rename(shape_geometry = geometry)

CA_lf_p <- terra::project(CA_lf,crs(EPA4_tf)) 

# Save raster that will be used in analysis
writeRaster(CA_lf_p, "data/CA_lf_p.tif")


unique_field_data <- BOLD_data %>%
  distinct(`Field ID`, Lat, Lon) %>%
  filter(!is.na(Lat))

lat_long <- unique_field_data %>% 
  st_as_sf(
    coords = c("Lon", "Lat"),
    agr = "constant",
    crs = 4326,  ##WGS84   
    stringsAsFactors = FALSE,
    remove = FALSE) %>%
  st_transform(3310)

landfire_points <- terra::extract(CA_lf, lat_long)

obs_L3 <- st_join(lat_long, EPA3_tf, join = st_within) %>%
  rename(L3_point_geometry = geometry)

obs_L4 <- st_join(lat_long, EPA4_tf, join = st_within) %>%
  rename(L4_point_geometry = geometry)

last_col_L4 <- obs_L4[, ncol(obs_L4)]

obs_L3_df <- obs_L3 %>% st_drop_geometry()
obs_L4_df <- obs_L4 %>% st_drop_geometry()

obs_L3_L4 <- obs_L3_df %>%
  bind_cols(obs_L4_df %>% select(US_L4NAME)) %>%
  cbind(landfire_points)

# manual fix
obs_L3_L4$US_L3NAME[obs_L3_L4$`Field ID` == "CIBI01950"] <- "Coast Range"
obs_L3_L4$US_L4NAME[obs_L3_L4$`Field ID` == "CIBI01950"] <- "Northern Franciscan Redwood Forest"
obs_L3_L4$US_L3NAME[obs_L3_L4$`Field ID` == "CIBI01826"] <- "Coast Range"
obs_L3_L4$US_L4NAME[obs_L3_L4$`Field ID` == "CIBI01826"] <- "Fort Bragg/Fort Ross Terraces"

BOLD_data_lf <- BOLD_data %>%
  left_join(obs_L3_L4 %>% select(`Field ID`, US_L3NAME, US_L4NAME, EVT_NAME), by = "Field ID") %>%
  filter(!is.na(US_L3NAME) & US_L3NAME != "")

# Save insect data that will be used in analysis
saveRDS(BOLD_data_lf, "data/BOLD_data_lf.rds")
