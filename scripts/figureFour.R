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

#### Figure 4, community dissimilarity across sites and ecoregions ####
# Panel A, violin plot of average dissimilarity by same, adjacent, nonadjacent ecoregions ####
# Filter to just BINs from Observatory Sites, and select only columns that we will use:
observatorySiteBINS <- dataJan04 %>%
  filter(projectType == "Observatory") %>%
  select(correctedSiteName,
         `BIN`,
         Lat,
         Lon,
         `Project Code`,
         `Sampling Protocol`,
         `Collection Date.y`,
         `Collection Date Accuracy`)

# Group observations by correctedSiteName, so we can get a single lat/long for each site (fixing the lat/long of varying precisions):
observatorySiteBINS <- observatorySiteBINS %>%
  group_by(correctedSiteName) %>%
  mutate(preciseLat = max(as.character(Lat)),
         preciseLon = max(as.character(Lon))) %>%
  ungroup()

# Create a column that combines site name with lat/long, so that we can check the number of unique sites at this point:
observatorySiteBINS$Location <- paste(observatorySiteBINS$correctedSiteName,
                                      observatorySiteBINS$preciseLat,
                                      observatorySiteBINS$preciseLon)

# Attempt to subset the data to only BINs collected during months that had collecting occur at all observatory sites:
# Convert dates to actual date values:
observatorySiteBINS$collectionDate <- as.Date(observatorySiteBINS$`Collection Date.y`,
                                              "%d-%b-%Y")
# Create a column for month
observatorySiteBINS$collectionMonth <- month(observatorySiteBINS$collectionDate)

# There are no months with all sites having a collecting event, so we will keep collecting month and use that as a control in the statistical testing:
observatorySiteBINS$LocationMonth <- paste(observatorySiteBINS$Location,
                                           observatorySiteBINS$collectionMonth,
                                           sep = "_")

#### Calculate Bray-Curtis dissimilarities for site-month pairs ####
# Turn the observatorySiteBINS into a matrix of species counts per site for calculating dissimilarities:
countsMatrix <- observatorySiteBINS %>%
  select(LocationMonth,
         BIN) 

countsMatrix <- table(countsMatrix) %>%
  as.data.frame.matrix()

# Compute dissimilarity indices:
bray <- vegdist(countsMatrix,
                "bray") 

brayTriangle <- bray
brayTriangle[upper.tri(brayTriangle)] <- NA

# Get a dataframe of the dissimilarity values for each pair of site-month values:
brayDataframe <- reshape2::melt(as.matrix(bray), 
                                varnames = c("row", "col"),
                                na.rm = FALSE)

# Create a column of alphabetically-ordered site-month pairs, so we can remove duplicate rows (e.g. if we have A-B and B-A, keep only A-B):
brayDataframe <- brayDataframe %>%
  rowwise() %>%     
  mutate(sitePair = paste(sort(c(row, 
                                 col)), 
                          collapse = "_")) %>%
  ungroup()  

# Get a column to tell us the sampling month of the first site-month in the pair:
brayDataframe$monthOne <- str_split_i(string = brayDataframe$sitePair,
                                      pattern = "_",
                                      i = 2)
# Now get just the site value:
brayDataframe$row <- str_split_i(string = brayDataframe$sitePair,
                                 pattern = "_",
                                 i = 1)
# Repeat for the second site-month in the pair:
brayDataframe$monthTwo <- str_split_i(string = brayDataframe$sitePair,
                                      pattern = "_",
                                      i = 4)
brayDataframe$col <- str_split_i(string = brayDataframe$sitePair,
                                 pattern = "_",
                                 i = 3)

# Remove the duplicates:
brayDataframe <- brayDataframe %>%
  distinct() %>%
  select(-c(sitePair))

#### Associate each sampling site with an ecoregion: ####
# Generate an sf object for the sampling locations:
locations <- observatorySiteBINS %>%
  distinct()

samplingLocations <- st_as_sf(locations, 
                              coords = c("preciseLon","preciseLat"))
samplingLocations <- st_set_crs(samplingLocations, 
                                "+proj=longlat +datum=WGS84") 
samplingLocations <- st_transform(samplingLocations, 
                                  src = st_crs(samplingLocations),
                                  crs = 3857)

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
rm(ecoregions)
rm(california)

matches <- st_intersects(samplingLocations, 
                         caEcoregions)
matches <- matches %>% 
  lapply(function(x) ifelse(is_empty(x), NA, x)) %>% 
  unlist() 

matchingNames <- caEcoregions$US_L3NAME[matches]
rm(matches)

# Add those ecoregion names back in to the sample dataframe:
samplingLocations$ecoregion <- matchingNames
rm(matchingNames)

# Join the ecoregion data to the sites:
sitesEcoregions <- samplingLocations %>%
  as.data.frame() %>%
  select(c(Location,
           ecoregion)) %>%
  distinct()

#### Add in info about whether sites are identical, in the same ecoregion, in adjacent ecoregions, or in nonadjacent ecoregions ####
# Connect the sampling locations in the Bray-Curtis matrix to their ecoregion:
brayDataframe <- left_join(brayDataframe,
                           sitesEcoregions,
                           by = c("row" = "Location"))
colnames(brayDataframe) <- c("row",
                             "col",
                             "value",
                             "monthOne",
                             "monthTwo",
                             "ecoregion_row")
brayDataframe <- left_join(brayDataframe,
                           sitesEcoregions,
                           by = c("col" = "Location"))
colnames(brayDataframe) <- c("row",
                             "col",
                             "value",
                             "monthOne",
                             "monthTwo",
                             "ecoregion_row",
                             "ecoregion_col")

# Get a column that tells us if the two sites being compared are identical, in the same ecoregion, or in different ecoregions:
brayDataframe <- brayDataframe %>%
  mutate(sameOrDifferent = case_when(row == col ~ "Identical",
                                     ecoregion_row == ecoregion_col ~ "Same",
                                     ecoregion_row != ecoregion_col ~ "Different"))

# For pairs that are in different ecoregions, are those adjacent ecoregions or nonadjacent?
# Merge all the components of the level III ecoregions and get them as a single feature:
levelIIIecoregions <- caEcoregions %>% 
  group_by(US_L3NAME) %>% 
  summarise()

singleFeature <- levelIIIecoregions[1,]

# Write a function to test whether two ecoregions touch or not:
testAdjacent <- function(rowNumberOne,
                         rowNumberTwo) {
  touch <- st_intersects(levelIIIecoregions[rowNumberOne,], 
                         levelIIIecoregions[rowNumberTwo,])
  touch <- touch[[1]] %>%
    as.numeric()
  
  if (rlang::is_empty(touch) == TRUE) {
    touch <- 0
  }
  result <- c(levelIIIecoregions[rowNumberOne,],
              levelIIIecoregions[rowNumberTwo,],
              touch)
  return(result)
}

# Make safe with possibly:
possiblyTestAdjacent <- possibly(testAdjacent)

# Generate a dataframe with all unique permutations of the row numbers in levelIIIecoregions:
rowNumbers <- 1:length(levelIIIecoregions$US_L3NAME)
allCombos <- gtools::permutations(n = length(rowNumbers),
                                  r = 2, 
                                  rowNumbers,
                                  repeats.allowed = TRUE) %>%
  as.data.frame()

# Iterate with map2:
adjacentResults <- purrr::map2(allCombos$V1,
                               allCombos$V2,
                               possiblyTestAdjacent)
adjacentResults <- as.data.frame(do.call(rbind, adjacentResults)) %>%
  as.data.frame()
colnames(adjacentResults) <- c("levelIIIregionOne", 
                               "geometryOne",
                               "levelIIIregionTwo",
                               "geometryTwo",
                               "adjacent")
adjacentResults <- adjacentResults %>%
  rowwise() %>%     
  mutate(ecoregionPair = paste(levelIIIregionOne, 
                               levelIIIregionTwo,
                               sep = "_")) 

adjacentResults$levelIIIregionOne <- str_split_i(adjacentResults$ecoregionPair,
                                                 pattern = "_",
                                                 i = 1)
adjacentResults$levelIIIregionTwo <- str_split_i(adjacentResults$ecoregionPair,
                                                 pattern = "_",
                                                 i = 2)

adjacentResults$adjacent <- as.numeric(adjacentResults$adjacent)

rm(allCombos)

# Now join the dictionary of ecoregion pair adjacency values to brayDataframe:
brayDataframeAdjacent <- full_join(brayDataframe,
                                   adjacentResults,
                                   by = c("ecoregion_row" = "levelIIIregionOne",
                                          "ecoregion_col" = "levelIIIregionTwo")) %>%
  mutate(category = case_when(adjacent == 1 ~ "Adjacent",
                              adjacent == 0 ~ "Nonadjacent")) %>%
  mutate(category = case_when(sameOrDifferent == "Identical" ~ "Identical",
                              sameOrDifferent == "Same" ~ "Same",
                              TRUE ~ category))

brayDataframeAdjacent$category <- factor(brayDataframeAdjacent$category,
                                         levels = c("Same",
                                                    "Adjacent",
                                                    "Nonadjacent",
                                                    "identical"))

# Save the model input data object:
brayDataframeAdjacent %>%
  select(-c(geometryOne,
            geometryTwo)) %>%
  readr::write_csv(file = "./objects/brayDataframeAdjacent.csv",
                   col_names = TRUE)

# Restore the model input data object:
brayDataframeAdjacent <- read_delim(file = "./objects/brayDataframeAdjacent.csv") %>%
  filter(!is.na(category))
brayDataframeAdjacent$category <- factor(brayDataframeAdjacent$category,
                                         levels = c("Same",
                                                    "Adjacent",
                                                    "Nonadjacent"))

#### Test whether being in the same/adjacent/nonadjacent ecoregions impacts community similarity, controlling for sampling month ####
# Run the bayesian regression with the sampling months as a multimembership grouping term:
formula <- bf(value ~ category + (1|mm(monthOne, 
                                       monthTwo, 
                                       row, 
                                       col)),
              phi ~ 1,
              zoi ~ 1,
              coi ~ 1 )

brmBray <- brm(formula,
               family = brms::zero_one_inflated_beta(),
               data = brayDataframeAdjacent,
               cores = 7, 
               init = 0,
               iter = 20000,
               control = list(max_treedepth = 20))

# Save the brm object:
saveRDS(brmBray, 
        file = "./objects/brmBray.rds")

# Restore the brm object:
brmBray <- readRDS(file = "./objects/brmBray.rds")

# If the 95% credible intervals do not overlap zero and do not overlap each other, we'd say we have a significant effect:
brmSummary <- summary(brmBray)
brmSummary

overtCredibleIntervals <- conditional_effects(brmBray, 
                                              method = "fitted", 
                                              effects = "category")
overtCredibleIntervals[["category"]]

brayDataframeAdjacent <- left_join(brayDataframeAdjacent,
                                   overtCredibleIntervals[["category"]],
                                   by = c("category" = "category"))

overtCredibleIntervals[["category"]] <- overtCredibleIntervals[["category"]] %>%
  mutate(xValue = case_when(category == "Same" ~ 1,
                            category == "Adjacent" ~ 2,
                            category == "Nonadjacent" ~ 3))

#### Plot results ####
siteViolinPlot <- ggplot() +
  geom_violin(data = filter(brayDataframeAdjacent,
                            category != "Identical"),
              mapping = aes(x = category,
                            y = value.x),
              alpha = 0.5,
              fill = "#d4d4d4",
              color = "#949494") +
  #  geom_point(alpha = 0.01,
  #            position = position_jitter(seed = 1, 
  #                                       width = 0.3)) + 
  geom_errorbar(data = overtCredibleIntervals[["category"]],
                mapping = aes(x = category,
                              ymin = `lower__`,
                              ymax = `upper__`)) +
  geom_segment(data = overtCredibleIntervals[["category"]],
               mapping = aes(x = xValue - 0.4,
                             xend = xValue + 0.4,
                             y = `estimate__`,
                             yend = `estimate__`),
               color = "red") +
  theme(legend.position = "none") +
  theme_bw() +
  labs(x = "Sites in the same, adjacent, or nonadjacent ecoregions",
       y = "Bray-Curtis dissimilarity index")

siteViolinPlot

#### Panel B, plot a heatmap of ecoregions to see which are most/least dissimilar: ####
# Which ecoregion pairs are most/least dissimilar?
brayDataframeAdjacent <- read_delim(file = "./objects/brayDataframeAdjacent.csv") %>%
  filter(!is.na(category))
brayDataframeAdjacent$category <- factor(brayDataframeAdjacent$category,
                                         levels = c("Same",
                                                    "Adjacent",
                                                    "Nonadjacent"))
brayForHeatmap <- brayDataframeAdjacent

# Create a column of alphabetically-ordered ecoregion pairs, so we can remove duplicate rows (e.g. if we have A-B and B-A, keep only A-B):
brayForHeatmap <- brayForHeatmap %>%
  rowwise() %>%     
  mutate(ecoregionPair = paste(sort(c(ecoregion_row, 
                                      ecoregion_col)), 
                               collapse = "_")) %>%
  ungroup()  

length(unique(brayForHeatmap$ecoregionPair))

# Get a column to tell us the first ecoregion in the pair:
brayForHeatmap$ecoregion_row <- str_split_i(string = brayForHeatmap$ecoregionPair,
                                            pattern = "_",
                                            i = 1)

# Repeat for the second ecoregion in the pair:
brayForHeatmap$ecoregion_col <- str_split_i(string = brayForHeatmap$ecoregionPair,
                                            pattern = "_",
                                            i = 2)

pairwiseSimilarity <- brayForHeatmap %>%
  filter(category != "Identical") %>%
  group_by(ecoregionPair) %>%
  summarise(mean = mean(value))

pairwiseSimilarity$ecoregionOne <- str_split_i(pairwiseSimilarity$ecoregionPair, 
                                               pattern = "_",
                                               i = 1)

pairwiseSimilarity$ecoregionTwo <- str_split_i(pairwiseSimilarity$ecoregionPair, 
                                               pattern = "_",
                                               i = 2)

# Make a second dataframe to fill the bottom diagonal in the heatmap:
pairwiseSimilarityTwo <- pairwiseSimilarity

pairwiseSimilarityTwo$ecoregionOne <- str_split_i(pairwiseSimilarityTwo$ecoregionPair, 
                                                  pattern = "_",
                                                  i = 2)

pairwiseSimilarityTwo$ecoregionTwo <- str_split_i(pairwiseSimilarityTwo$ecoregionPair, 
                                                  pattern = "_",
                                                  i = 1)

pairwiseSimilarity <- rbind(pairwiseSimilarity,
                            pairwiseSimilarityTwo) %>%
  filter(!is.na(mean))

# Add rows for site combos that are missing in the current data:
pairwiseSimilarity <- rbind(pairwiseSimilarity, 
                            c("Cascades_Cascades", NA, "Cascades", "Cascades"))
pairwiseSimilarity$mean <- as.numeric(pairwiseSimilarity$mean)

# Make a dendrogram to cluster ecoregions by similarity:
matrixOfSimilarity <- pairwiseSimilarity %>%
  group_by(ecoregionPair) %>%
  slice(1) %>%
  ungroup() %>%
  select(c("ecoregionOne",
           "ecoregionTwo",
           "mean")) %>%
  pivot_wider(names_from = ecoregionTwo, 
              values_from = mean) %>%
  as.data.frame()

rownames(matrixOfSimilarity) <- matrixOfSimilarity$ecoregionOne
matrixOfSimilarity <- matrixOfSimilarity %>%
  select(-c(ecoregionOne))

dendrogram <- as.dendrogram(hclust(d = dist(x = matrixOfSimilarity)))
dendrogramPlot <- ggdendro::ggdendrogram(data = dendrogram, 
                                         rotate = TRUE)

# Get the order of ecoregions so that we can order the heatmap:
dendrogramOrder <- order.dendrogram(dendrogram)
pairwiseSimilarity$ecoregionOne <- factor(x = pairwiseSimilarity$ecoregionOne,
                                          levels = rownames(matrixOfSimilarity)[dendrogramOrder], 
                                          ordered = TRUE)
pairwiseSimilarity$ecoregionTwo <- factor(x = pairwiseSimilarity$ecoregionTwo,
                                          levels = rownames(matrixOfSimilarity)[dendrogramOrder], 
                                          ordered = TRUE)

# Now plot:
ecoregionHeatMap <- ggplot(pairwiseSimilarity) +
  geom_tile(mapping = aes(x = ecoregionOne, 
                          y = ecoregionTwo, 
                          fill = mean),
            stat = "identity") +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, 
                                                         "YlGnBu"),
                       na.value = "#adadad") +
  coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25,
                                   hjust = 1)) +
  labs(fill = "Mean Bray–Curtis\ndissimilarity",
       x = "",
       y = "" )

ecoregionHeatMap

# Combine the panels ####
library(patchwork)

brayCurtisFigure <- ecoregionHeatMap | siteViolinPlot

brayCurtisFigure + 
  plot_annotation(tag_levels = 'A') 

ggsave("./cibiFigureFour.png",
       width = 12,
       height = 8,
       units = "in",
       dpi = 900)

ggsave("./cibiFigureFour.pdf",
       width = 12,
       height = 8,
       units = "in",
       dpi = 900)

