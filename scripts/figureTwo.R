library(tidyverse)
library(scales)
library(lubridate)
library(janitor)
library(rgbif)
library(sf)
library(spData)
library(terra)
library(tidycensus)
library(patchwork)
library(USAboundaries)
library(gggibbous)
library(vegan)
library(brms)
library(parameters)
library(tidybayes) 
library(iNEXT)
library(ComplexUpset)

#### Read in and process data: ####
dataJan04 <- readxl::read_xlsx("./250626_BOLD_data_clean.xlsx")

dataJan04 <- dataJan04 %>%
  filter(BIN != "NA")

# Fix some messed up names in the observatory sites:
dataJan04 <- dataJan04 %>% 
  mutate(correctedSiteName = case_when(`Exact Site` == "site 5" &
                                         Lat == "41.095" &
                                         Lon == "-122.115" ~ "Kerry Landreth McCloud River Preserve",
                                       `Exact Site` == "site 6" &
                                         Lat == "40.101" &
                                         Lon == "-122.054" ~ "Dye Creek Preserve",
                                       `Exact Site` == "site 7" &
                                         Lat == "39.625" &
                                         Lon == "-120.574" ~ "San Francisco State University Sierra Nevada Field Campus",
                                       `Exact Site` == "site 4" &
                                         Lat == "38.864" &
                                         Lon == "-122.402" ~ "Donald and Sylvia McLaughlin Natural Reserve",
                                       `Exact Site` == "NA" &
                                         Lat == "38.552307" &
                                         Lon == "-121.729164" ~ "Davis",
                                       `Exact Site` == "Point Reyes" ~ "Point Reyes Field Station",
                                       `Exact Site` == "NA" &
                                         Lat == "33.323868" &
                                         Lon == "-118.438835" ~ "Catalina Island, Wild Boar Gully",
                                       `Exact Site` == "Pinon Flats" &
                                         Lat == "33.6057" &
                                         Lon == "-116.4506" ~ "Boyd Deep Canyon Desert Research Center (Pinon Flats)",
                                       `Exact Site` == "Motte Rimrick Reserve" ~ "Motte Rimrock Reserve",
                                       `Exact Site` == "Oasis de los Osos" &
                                         Lat == "33.8894" &
                                         Lon == "-116.6864" ~ "James Reserve Low Elevation (Oasis de los Osos)",
                                       `Exact Site` == "site 1" &
                                         Lat == "34.149" &
                                         Lon == "-116.456" ~ "Burns Pinon Ridge Reserve",
                                       `Exact Site` == "site 3" &
                                         Lat == "34.526" &
                                         Lon == "-120.409" ~ "Jack and Laura Dangermond Preserve",
                                       `Exact Site` == "Granite Mtns Res. MT1" &
                                         Lat == "34.78278" &
                                         Lon == "-115.65083" ~ "Granite Mountains Reserve",
                                       `Exact Site` == "site 2" &
                                         Lat == "35.203" &
                                         Lon == "-118.622" ~ "Frank and Joan Randall Preserve",
                                       `Exact Site` == "Oak woodlands at Reinhardt Redwood Regional Park" ~ "Reinhardt Redwood Regional Park (Oak Woodland)",
                                       `Exact Site` == "Reinhart Redwood Regional Park (Oak Woodland)" ~ "Reinhardt Redwood Regional Park (Oak Woodland)",
                                       `Exact Site` == "Redwoods at Reinhardt Redwood Regional Park" ~ "Reinhardt Redwood Regional Park (Redwood Grove)",
                                       `Exact Site` == "Reinhart Redwood Regional Park (Redwood Grove)" ~ "Reinhardt Redwood Regional Park (Redwood Grove)",
                                       `Exact Site` == "Trap B" ~ "Santa Cruz Island, Trap B",
                                       `Exact Site` == "Trap A" ~ "Santa Cruz Island, Trap A",
                                       TRUE ~ `Exact Site`)) 

# Info about observatory/snapshot is in dataJan04$`Project Code`; if it starts with "CIO" or CIBIP, and in Field ID, if it starts with GMP:
dataJan04 <- dataJan04 %>%
  mutate(projectType = case_when(grepl(pattern = "CIO", x = dataJan04$`Project Code`) ~ "Observatory",
                                 grepl(pattern = "CIBIP", x = dataJan04$`Project Code`) ~ "Observatory",
                                 grepl(pattern = "GMP", x = dataJan04$`Field ID`) ~ "Observatory",
                                 TRUE ~ "Snapshot"))

# Group by correctedSiteName so we can fix the lat/long of varying precisions:
dataJan04 <- dataJan04 %>%
  group_by(correctedSiteName) %>%
  mutate(preciseLat = max(as.character(Lat)),
         preciseLon = max(as.character(Lon))) %>%
  ungroup()

#### Figure 2: Geographic results ####
# Panel A, map of all sampling sites####
# Read in the ecoregion data:
ecoregions <- st_read("./objects/us_eco_l3_state_boundaries/us_eco_l3_state_boundaries.shp") 
ecoregions <- st_transform(ecoregions, 
                           src = st_crs(ecoregions),
                           crs = 3857)
# Get a shapefile for California:
california <- tigris::states() %>%
  filter(NAME == "California")
california <- st_transform(california, 
                           src = st_crs(california),
                           crs = 3857)

# Get just CA ecoregions
caEcoregions <- sf::st_intersection(ecoregions,
                                    california)

# Generate an sf object for the sampling locations:
locations <- dataJan04 %>%
  mutate(samplingType = case_when(`Sampling Protocol` == "Malaise Trap" ~ "Malaise Trap",
                                  TRUE ~ "Other")) %>%
  select(c("Region",
           "Sector",
           "Exact Site",
           "preciseLat",
           "preciseLon",
           "Elev",
           "Collection Date Accuracy",
           "Habitat",
           "samplingType",
           "projectType")) %>%
  distinct()

samplingLocations <- st_as_sf(locations, 
                              coords = c("preciseLon", "preciseLat"))
samplingLocations <- st_set_crs(samplingLocations, 
                                "+proj=longlat +datum=WGS84") 
samplingLocations <- st_transform(samplingLocations, 
                                  src = st_crs(samplingLocations),
                                  crs = 3857)

# Color palette for the ecoregions:
ecoregionColors <- c("#f57d03",
                     "#ffeb05",
                     "#415dac",
                     "#6e96f0",
                     "#c5a4cf",
                     "#afb32c",
                     "#794f47",
                     "#c3195d",
                     "#f4b400",
                     "#0a743d",
                     "#db4336",
                     "#7c3692",
                     "#69dff0")
names(ecoregionColors) <- c("Coast Range",
                            "Central Basin and Range",
                            "Mojave Basin and Range",
                            "Cascades",
                            "Sierra Nevada",
                            "Central California Foothills and Coastal Mountains",
                            "Central California Valley",
                            "Klamath Mountains/California High North Coast Range",
                            "Southern California Mountains",
                            "Northern Basin and Range",
                            "Sonoran Basin and Range",
                            "Southern California/Northern Baja Coast",
                            "Eastern Cascades Slopes and Foothills")
allSitesMap <- ggplot() + 
  geom_sf(data = caEcoregions, 
          mapping = aes(fill = US_L3NAME),
          alpha = 0.4,
          linewidth = 0) + 
  scale_fill_manual(values = ecoregionColors) +
  geom_sf(data = filter(samplingLocations,
                        projectType == "Snapshot"),
          color = "black",
          fill = "#de3a3a",
          shape = 21) +
  geom_sf(data = filter(samplingLocations,
                        projectType == "Observatory"),
          color = "black",
          fill = "#0dd6a0",
          shape = 24) +
  theme_void() +
  theme(legend.position = "none",
        legend.direction = "vertical") +
  ggspatial::annotation_scale() +
  ggspatial::annotation_north_arrow(location = "tr",
                                    height = unit(1, "cm"),
                                    width = unit(1, "cm"))

allSitesMap

# Panel B, Number of BINs per ecoregion ####
# Generate an sf object for the samples:
sampleLocations <- dataJan04 %>%
  select(c("BIN",
           "preciseLat",
           "preciseLon",
           "Exact Site",
           "Collection Date.y",
           "projectType")) 

samplingLocations <- st_as_sf(sampleLocations, 
                              coords = c("preciseLon", "preciseLat"))
samplingLocations <- st_set_crs(samplingLocations, 
                                "+proj=longlat +datum=WGS84") 
samplingLocations <- st_transform(samplingLocations, 
                                  src = st_crs(samplingLocations),
                                  crs = 3857)

# Associate each sample with an ecoregion:
matches <- st_intersects(samplingLocations, 
                         caEcoregions)
matches <- matches %>% 
  lapply(function(x) ifelse(is_empty(x), NA, x)) %>% 
  unlist() 

matchingNames <- caEcoregions$US_L3NAME[matches]

# Add those ecoregion names back in to the sample dataframe:
samplingLocations$ecoregion <- matchingNames

# Now get distinct BINs for each region:
sampleEcoregions <- samplingLocations %>%
  as.data.frame() %>%
  dplyr::select("BIN",
                "ecoregion") %>%
  distinct()

# Add a column with ecoregion abbreviations for plotting:
sampleEcoregions$abbreviation <- sapply(stringr::str_extract_all(sampleEcoregions$ecoregion,
                                                                 '[A-Z]'),
                                        paste0, 
                                        collapse = '')

# Create a light color palette for the ecoregions, to match the alpha = 0.4 map:
ecoregionColorsLight <- c("#fbcdad",
                          "#fff7ae",
                          "#bbc3df",
                          "#c9d7f9",
                          "#e8dcec",
                          "#e0e1b6",
                          "#cdbfbd",
                          "#e8acc4",
                          "#fbe2ac",
                          "#afcbba",
                          "#f1bbb8",
                          "#ceb7d6",
                          "#c8f2f9")
names(ecoregionColorsLight) <- c("Coast Range",
                                 "Central Basin and Range",
                                 "Mojave Basin and Range",
                                 "Cascades",
                                 "Sierra Nevada",
                                 "Central California Foothills and Coastal Mountains",
                                 "Central California Valley",
                                 "Klamath Mountains/California High North Coast Range",
                                 "Southern California Mountains",
                                 "Northern Basin and Range",
                                 "Sonoran Basin and Range",
                                 "Southern California/Northern Baja Coast",
                                 "Eastern Cascades Slopes and Foothills")

# Plot:
samplesPerEcoregion <- ggplot(data = filter(sampleEcoregions,
                                            !is.na(ecoregion)), 
                              mapping = aes(x = fct_infreq(abbreviation),
                                            fill = ecoregion,
                                            color = ecoregion)) +
  geom_bar() +
  scale_fill_manual(values = ecoregionColorsLight) +
  scale_color_manual(values = ecoregionColors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "none") +
  labs(x = "Level III Ecoregions",
       y = "Count of unique species") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

samplesPerEcoregion

# Panel C, Collection effort per ecoregion ####
# Get the area of each ecoregion and the proportion of the state it covers:
caEcoregionsDataframe <- as.data.frame(caEcoregions)
caEcoregionsDataframe$area <- st_area(caEcoregions) %>%
  str_split_i(pattern = " ",
              i = 1) %>%
  as.numeric()

caEcoregionsDataframe <- caEcoregionsDataframe %>%
  group_by(US_L3NAME) %>%
  summarise(totalArea = sum(area))

caEcoregionsDataframe$proportionOfState <- caEcoregionsDataframe$totalArea / sum(caEcoregionsDataframe$totalArea)

# Get the total number of collecting events in each ecoregion:
collectingEventsPerEcoregion <- samplingLocations %>%
  as.data.frame() %>%
  select(-c(BIN)) %>%
  distinct() %>%
  group_by(ecoregion) %>%
  summarise(totalEvents = n())

collectingEventsPerEcoregion$proportionOfEvents <- collectingEventsPerEcoregion$totalEvents / sum(collectingEventsPerEcoregion$totalEvents)

# Join the data on coverage of state by ecoregion with collecting events by ecoregion:
samplingEffort <- full_join(collectingEventsPerEcoregion,
                            caEcoregionsDataframe,
                            by = c("ecoregion" = "US_L3NAME"))

samplingEffort$samplingEffort <- samplingEffort$proportionOfEvents / samplingEffort$proportionOfState

samplingEffort$bias <- samplingEffort$samplingEffort - 1

samplingEffort <- arrange(samplingEffort,
                          -desc(bias))

# Add a column with ecoregion abbreviations for plotting:
samplingEffort$abbreviation <- sapply(stringr::str_extract_all(samplingEffort$ecoregion,
                                                               '[A-Z]'),
                                      paste0, 
                                      collapse = '')

# Add in the total number of samples per ecoregion, so that the bars can be ordered for plotting:
samplesPerEcoregionDataframe <- sampleEcoregions %>%
  group_by(ecoregion) %>%
  summarise(totalSamples = n())

samplingEffort <- full_join(samplingEffort,
                            samplesPerEcoregionDataframe,
                            by = c("ecoregion" = "ecoregion"))

samplingEffort$ecoregion <- fct_reorder(samplingEffort$ecoregion, 
                                        -samplingEffort$totalSamples,
                                        .desc = FALSE) 

# Plot:
biasInSamplingEffort <- ggplot(data = filter(samplingEffort,
                                             !is.na(ecoregion)), 
                               mapping = aes(x = fct_reorder(abbreviation,
                                                             -totalSamples),
                                             y = bias,
                                             fill = ecoregion,
                                             color = ecoregion)) +
  geom_col() +
  scale_fill_manual(values = ecoregionColorsLight) +
  scale_color_manual(values = ecoregionColors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "none") + 
  guides(colour = guide_legend(ncol = 3),
         fill = guide_legend(ncol = 3)) +
  labs(x = "Level III Ecoregions",
       y = "Bias in sampling effort") 

biasInSamplingEffort

#### Combine and save plots ####
figureTwoPlots <- (allSitesMap | (samplesPerEcoregion / biasInSamplingEffort)) +
  plot_layout(heights = c(2, 1),
              width = c(1.5, 1))

figureTwo <- figureTwoPlots / ggpubr::get_legend(biasInSamplingEffort, position = "bottom")

figureTwo <- figureTwo + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(3, 1))

figureTwo

ggsave("./cibiFigureTwo.png",
       width = 12,
       height = 8,
       units = "in",
       dpi = 900)

ggsave("./cibiFigureTwo.pdf",
       width = 12,
       height = 8,
       units = "in",
       dpi = 900)
