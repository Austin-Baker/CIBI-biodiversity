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

#### Figure 3: Estimated species richness and sample completeness for Observatory Sites ####
# Read in the completeness estimator results:
siteCompleteness <- readxl::read_xlsx("./objects/chaoEstimatesUpdated.xlsx") %>%
  head(n = 19) %>%
  filter(!is.na(`Site name`))

colnames(siteCompleteness) <- c("site",
                                "bulk samples pre-filter",
                                "bulk samples post-filer",
                                "sample size",
                                "# observed species",
                                "coverage estimate for entire dataset",
                                "CV for entire dataset",
                                "cut-off point",
                                "# observed individuals for rare group",
                                "# observed species for rare group",
                                "estimate of sample coverage for rare group",
                                "estimate of CV for rare group in ACE",
                                "estimate of CV1 for rare group in ACE-1",
                                "number of observed individuals for abundant group",
                                "number of observed species for abundant group",
                                "estimate",
                                "standard error",
                                "95%Lower",
                                "95%Upper",
                                "lower difference",
                                "upper difference",
                                "% completeness",
                                "lat",
                                "long",
                                "null")

siteCompleteness$`% completeness` <- as.numeric(siteCompleteness$`% completeness`)
siteCompleteness$estimate <- as.numeric(siteCompleteness$estimate)
siteCompleteness$`# observed species` <- as.numeric(siteCompleteness$`# observed species`)
siteCompleteness$`95%Lower` <- as.numeric(siteCompleteness$`95%Lower`)
siteCompleteness$`95%Upper` <- as.numeric(siteCompleteness$`95%Upper`)


#### Panel A, Map the sites measured for completeness, with moon plots of completeness ####
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

# Generate an sf object for the samples:
samplingLocations <- st_as_sf(siteCompleteness, 
                              coords = c("long","lat"))
samplingLocations <- st_set_crs(samplingLocations, 
                                "+proj=longlat +datum=WGS84") 
samplingLocations <- st_transform(samplingLocations, 
                                  src = st_crs(samplingLocations),
                                  crs = 3857)
samplingLocations$x <- as.character(samplingLocations$geometry) %>%
  str_split_i( pattern = "\\(",
               i = 2) %>%
  str_split_i( pattern = ",",
               i = 1) %>%
  as.numeric()

samplingLocations$y <- as.character(samplingLocations$geometry) %>%
  str_split_i( pattern = "\\(",
               i = 2) %>%
  str_split_i( pattern = ",",
               i = 2) %>%
  str_split_i( pattern = "\\)",
               i = 1) %>%
  as.numeric()

# Plot:
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

# Make the gibbous map plot:
gibbousMap <- ggplot() + 
  geom_sf(data = caHighRes, 
          mapping = aes(),
          fill = "#e3e3e3",
          linewidth = 0) + 
  gggibbous::geom_moon(data = samplingLocations,
                       mapping = aes(x = x,
                                     y = y,
                                     ratio = `% completeness`), 
                       right = FALSE,
                       size = 5, 
                       fill = "#8E4101") + 
  gggibbous::geom_moon(data = samplingLocations,
                       mapping = aes(x = x,
                                     y = y,
                                     ratio = 1 - `% completeness`), 
                       right = TRUE,
                       size = 5, 
                       fill = "white") +
  ggrepel::geom_text_repel(data = samplingLocations,
                           mapping = aes(x = x,
                                         y = y,
                                         label = paste(site,
                                                       scales::percent(`% completeness`,
                                                                       accuracy = 0.01),
                                                       sep = ", ")),
                           size = 3,
                           point.padding = 10,
                           min.segment.length = unit(0, 
                                                     'lines'),
                           max.overlaps = 30) +
  theme_void() 
gibbousMap

#### Panel B, bar plot of sampled BINs, estimated BINs, and sampling completeness per site ####
# Create a column with percent out of 100, rather than as a decimal:
siteCompleteness$percentComplete <- siteCompleteness$`% completeness` * 100
# Scale that number by the estimated number of species, to be able to plot them all together:
scalingFactor <- 100/(100/max(siteCompleteness$estimate))
siteCompleteness$percentCompleteScaled <- siteCompleteness$`% completeness` * scalingFactor

# Fix the site names so that they plot better:
siteCompleteness <- siteCompleteness %>% 
  mutate(site = case_when(site == "Reinhardt Redwood Regional Park (Oak Woodland)" ~ "Reinhardt Redwood\nRegional Park (Oak Woodland)",
                          site == "Reinhardt Redwood Regional Park (Redwood Grove)" ~ "Reinhardt Redwood\nRegional Park (Redwood Grove)",
                          site == "Mojave National Preserve, Desert Studies Center" ~ "Mojave National Preserve,\nDesert Studies Center",
                          site == "Burns Pinon Ridge Reserve" ~ "Burns Pinon\nRidge Reserve",
                          site == "Granite Mountains Reserve" ~ "Granite Mountains\nReserve",
                          site == "Blue Oak Ranch Reserve" ~ "Blue Oak\nRanch Reserve",
                          site == "Tierra Del Sol SDAA" ~ "Tierra Del\nSol SDAA",
                          site == "Columbine Rd SW of Weed" ~ "Columbine Rd\nSW of Weed",
                          site == "Anza Borrego\nUC Reserve" ~ "Anza Borrego UC Reserve",
                          site == "Santa Cruz Island, Trap B" ~ "Santa Cruz Island,\nTrap B",
                          site == "Lopez Ridge Vernal Pools" ~ "Lopez Ridge\nVernal Pools",
                          TRUE ~ site))

# Order the sites by estimated number of species:
siteCompleteness$site <- fct_reorder(siteCompleteness$site,
                                     siteCompleteness$`# observed species`) 

# Pivot the data longer to be able to plot multiple factors together:
longCompleteness <- siteCompleteness %>%
  select(c("site",
           "# observed species",
           "estimate",
           "percentCompleteScaled",
           "95%Lower",
           "95%Upper")) 

colnames(longCompleteness) <- c("Site",
                                "Observed species",
                                "Estimated species present",
                                "Percent sampling completeness",
                                "95%Lower",
                                "95%Upper")

longCompleteness <- longCompleteness %>%
  pivot_longer(c("Observed species",
                 "Estimated species present",
                 "Percent sampling completeness"),
               names_to = "measure",
               values_to = "value")

# Fix the confidence intervals, since they only apply to the Chao estimates:
longCompleteness <- longCompleteness %>%
  mutate(`95%Lower` = case_when(measure == "Estimated species present" ~ `95%Lower`,
                                TRUE ~ NA),
         `95%Upper` = case_when(measure == "Estimated species present" ~ `95%Upper`,
                                TRUE ~ NA))

# Re-order measures so that observed, then estimated, then percent in terms of bar order
longCompleteness$measure <- factor(longCompleteness$measure,
                                   c("Observed species",
                                     "Estimated species present",
                                     "Percent sampling completeness"))

# Generate a color scheme:
measureColors <- c("#7B802B",
                   "#D7D9B1",
                   "#8E4101")
names(measureColors) <- c("Observed species",
                          "Estimated species present",
                          "Percent sampling completeness")

# Plot number sampled, number estimated, and percent complete all together:
samplingCompletenessBars <- ggplot(data = longCompleteness,
                                   mapping = aes(x = Site, 
                                                 y = value,
                                                 fill = measure,
                                                 ymin = `95%Lower`, 
                                                 ymax = `95%Upper`)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(position = position_dodge(width = 0.9), 
                width = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(colour = "#7B802B"),
        axis.text.y.right = element_text(colour = "#8E4101"),
        axis.title.y = element_text(colour = "#7B802B"),
        axis.title.y.right = element_text(colour = "#8E4101"),
        legend.position = "bottom") + 
  scale_y_continuous(name = "Number of species",
                     limits = c(0, max(longCompleteness$value) * 1.1),
                     sec.axis = sec_axis(~.*1/scalingFactor, 
                                         name = "Estimated sampling\ncompleteness",
                                         labels = scales::percent),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = measureColors) +
  labs(fill = "")

samplingCompletenessBars

#### Panel C, relationship between number of individuals sampled and completeness estimates ####
# Calculate the chao1 estimator for subsets of specimens:
calculateChaoForSubset <- function(subsetSize,
                                   dataFilter) {
  print(paste("Calculating estimator for a subset of",
              subsetSize,
              "samples"))
  
  if (dataFilter == "Observatory") {
    print("Running on observatory data only")
    BOLD_data <- dataJan04 %>%
      filter(projectType == "Observatory")
  } else {
    print("Running on all data")
    BOLD_data <- dataJan04
  }
  
  subsetSamples <- slice_sample(BOLD_data,
                                n = subsetSize, 
                                replace = FALSE)
  
  binCount <- subsetSamples %>% 
    group_by(BIN) %>%
    count()
  
  # Calculate the estimator:
  observatorySiteRichness <- SpadeR::ChaoSpecies(binCount$n,
                                                 "abundance",
                                                 k = 10,
                                                 conf = 0.95)
  
  # Extract the key results:
  dataInfo <- observatorySiteRichness$Basic_data_information %>%
    select(-c(Variable)) %>%
    t() %>%
    as.data.frame()
  
  speciesTable <- as.data.frame(observatorySiteRichness$Species_table) %>%
    as.data.frame()
  
  # Number of observed species:
  observedSpecies <- dataInfo["Value",
                              "    Number of observed species"] %>%
    as.numeric()
  # Chao1 estimator value:
  estimatedCompleteness <- speciesTable["    iChao1 (Chiu et al. 2014)",
                                        "Estimate"] %>%
    as.numeric()
  # Lower bound of 95% confidence interval:
  lower95Interval <- speciesTable["    iChao1 (Chiu et al. 2014)",
                                  "95%Lower"] %>%
    as.numeric()
  # Upper bound of 95% confidence interval:
  upper95Interval <- speciesTable["    iChao1 (Chiu et al. 2014)",
                                  "95%Upper"] %>%
    as.numeric()
  
  # Return those results:
  results <- c(subsetSize,
               observedSpecies,
               estimatedCompleteness,
               lower95Interval,
               upper95Interval,
               dataFilter) 
  
  resultsDataFrame <- as.matrix(t(results))
  write.table(resultsDataFrame, 
              file = "chaoWholeDatasetIterations.csv", 
              sep = ",", 
              col.names = FALSE, 
              row.names = FALSE,
              append = TRUE)
  return(results)
}

possiblyCalculateChaoForSubset <- purrr::possibly(calculateChaoForSubset,
                                                  otherwise = "Error")

# Create a list of subsets for which to calculate the Chao estimator:
subsetValues <- seq(from = 100, 
                    to = 900100, 
                    by = 100000)

# Calculate for all of the subsets with only observatory data:
purrr::map(subsetValues, 
           ~ possiblyCalculateChaoForSubset(.x,
                                            dataFilter = "Observatory"))

# Calculate for all of the subsets with all  data:
purrr::map(subsetValues, 
           ~ possiblyCalculateChaoForSubset(.x,
                                            dataFilter = "all"))

# Read in those results
iterations <- read_csv("./objects/chaoWholeDatasetIterations.csv",
                       col_names = FALSE) 

colnames(iterations) <- c("subsetSize",
                          "observedSpecies",
                          "estimatedCompleteness",
                          "lower95Interval",
                          "upper95Interval",
                          "dataFilter")
iterations$subsetSize <- as.numeric(iterations$subsetSize)
iterations$estimatedCompleteness <- as.numeric(iterations$estimatedCompleteness)

# Also run a rarefaction analysis with iNext, which uses the Chao (rather than iChao) estimator
binCount <- dataJan04 %>% 
  group_by(BIN) %>%
  count()

binCount <- as.data.frame(binCount)
rownames(binCount) <- binCount$BIN
binCount <- binCount %>%
  select(-c(BIN))

rarefaction <- iNEXT::iNEXT(binCount, 
                            q = 0, 
                            datatype = "abundance",
                            endpoint = 1030728*2,
                            knots = 25,
                            se = TRUE)

# Save the rarefaction object:
saveRDS(rarefaction, 
        file = "./objects/globalRarefaction.rds")

# Restore the rarefaction object:
rarefaction <- readRDS(file = "./objects/globalRarefaction.rds")

# Get the sample estimates:
subsampleEstimates <- rarefaction[["iNextEst"]]$size_based

# Generate a simple species accumulation curve for our data:
generateSpeciesAccumulation <- function(subsetSize,
                                        dataFilter) {
  print(paste("Calculating total number of species for a subset of",
              subsetSize,
              "samples"))
  
  if (dataFilter == "Observatory") {
    print("Running on observatory data only")
    BOLD_data <- dataJan04 %>%
      filter(projectType == "Observatory")
  } else {
    print("Running on all data")
    BOLD_data <- dataJan04
  }
  
  subsetSamples <- slice_sample(BOLD_data,
                                n = subsetSize, 
                                replace = FALSE)
  
  binCount <- subsetSamples %>% 
    group_by(BIN) %>%
    count()
  
  numberOfBins <- length(binCount$BIN)
  
  # Return the subset size and number of ensuing BINs:
  results <- c(subsetSize,
               numberOfBins,
               dataFilter) 
  
  resultsDataFrame <- as.matrix(t(results))
  write.table(resultsDataFrame, 
              file = "./objects/speciesAccumulationCurve.csv", 
              sep = ",", 
              col.names = FALSE, 
              row.names = FALSE,
              append = TRUE)
  return(results)
}

possiblyGenerateSpeciesAccumulation <- purrr::possibly(generateSpeciesAccumulation,
                                                       otherwise = "Error")

# Create a list of subsets for which to calculate the Chao estimator:
subsetValues <- subsampleEstimates$m[subsampleEstimates$m <= 927654]

# Calculate for all of the subsets with all  data:
purrr::map(subsetValues, 
           ~ possiblyGenerateSpeciesAccumulation(.x,
                                                 dataFilter = "all"))

# Create a list of subsets for which to calculate the Chao estimator:
subsetValues <- subsampleEstimates$m[subsampleEstimates$m <= 729734]

# Calculate for all of the subsets with only observatory data:
purrr::map(subsetValues, 
           ~ possiblyGenerateSpeciesAccumulation(.x,
                                                 dataFilter = "Observatory"))

# Read in those results
speciesAccumulation <- read_csv("./objects/speciesAccumulationCurve.csv",
                                col_names = FALSE) 

colnames(speciesAccumulation) <- c("subsetSize",
                                   "estimatedCompleteness",
                                   "Analysis")
speciesAccumulation$subsetSize <- as.numeric(speciesAccumulation$subsetSize)
speciesAccumulation$estimatedCompleteness <- as.numeric(speciesAccumulation$estimatedCompleteness)
speciesAccumulation$Method <- "Species accumulation"

# Combine iNext, iChao, and species accumulations:
iNextResults <- subsampleEstimates %>%
  select(c(m,
           qD,
           qD.LCL,
           qD.UCL))
iNextResults$dataFilter <- "iNext"
iNextResults$method <- subsampleEstimates$Method

ourResults <- iterations %>%
  select(subsetSize,
         estimatedCompleteness,
         lower95Interval,
         upper95Interval,
         dataFilter) %>%
  filter(!(dataFilter == "Observatory" &
             subsetSize > 700100))
ourResults$method <- "Rarefaction"

allEstimationResults <- rbind(ourResults,
                              setNames(iNextResults, 
                                       names(ourResults)))

colnames(allEstimationResults) <- c("subsetSize",
                                    "estimatedCompleteness",
                                    "lower95Interval",
                                    "upper95Interval",
                                    "Analysis",
                                    "Method")

allEstimationResults <- plyr::rbind.fill(allEstimationResults,
                                         speciesAccumulation)
allEstimationResults$`Analysis, data type` <- paste(allEstimationResults$Method,
                                                    allEstimationResults$Analysis,
                                                    sep = ", ")
allEstimationResults$`Analysis, data type` <- factor(allEstimationResults$`Analysis, data type`, 
                                                     levels = c("Rarefaction, iNext",
                                                                "Observed, iNext",
                                                                "Extrapolation, iNext",
                                                                "Species accumulation, all",
                                                                "Rarefaction, all",
                                                                "Species accumulation, Observatory",
                                                                "Rarefaction, Observatory"))

# Export as a CSV:
write_csv(allEstimationResults,
          file = "allSamplingCompletenessResults.csv")

# Look at the asymptotic estimator from iNext:
asymptote <- rarefaction[["AsyEst"]]
asymptote

# Plot:
rarefaction <- ggplot(data = allEstimationResults) +
  geom_line(mapping = aes(x = subsetSize,
                          y = estimatedCompleteness,
                          color = `Analysis, data type`)) +
  geom_point(mapping = aes(x = subsetSize,
                           y = estimatedCompleteness,
                           color = `Analysis, data type`,
                           shape = `Analysis, data type`),
             size = 2,
             fill = "white") +
  geom_errorbar(mapping = aes(x = subsetSize,
                              y = estimatedCompleteness,
                              ymin = lower95Interval, 
                              ymax = upper95Interval), 
                width = 0.2,
                position = position_dodge(0.05)) +
  geom_hline(yintercept = 60944.93181,
             linetype = "dotted") +
  scale_shape_manual(values = c("Rarefaction, all" = 21,
                                "Rarefaction, Observatory" = 21,
                                "Rarefaction, iNext" = 21,
                                "Observed, iNext" = 17,
                                "Extrapolation, iNext" = 4,
                                "Species accumulation, Observatory" = 16,
                                "Species accumulation, all" = 16)) +
  scale_color_manual(values = c("Rarefaction, all" = "#2caad4",
                                "Rarefaction, Observatory" = "#72d6c9",
                                "Rarefaction, iNext" = "#0f056e",
                                "Observed, iNext" = "#675fb3",
                                "Extrapolation, iNext" = "#b0abe0",
                                "Species accumulation, Observatory" = "#72d6c9",
                                "Species accumulation, all" = "#2caad4")) +
  theme_bw() +
  theme(axis.text.x = element_text()) +
  labs(x = "Number of sampled individuals",
       y = "Number of species estimated to be present") + 
  scale_x_continuous(labels = scales::label_comma())

rarefaction

#### Combine panels ####
figure3 <- gibbousMap | (samplingCompletenessBars / rarefaction) 

figure3 + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(1.5, 1))

ggsave("./cibiFigureThree.png",
       width = 14,
       height = 10,
       units = "in",
       dpi = 900)
ggsave("./cibiFigureThree.pdf",
       width = 14,
       height = 10,
       units = "in",
       dpi = 900)

