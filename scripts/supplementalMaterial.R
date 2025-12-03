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

#### Create the suppmental table that summarizes counts for BIN/ID: ####
# Get the ID specificity for each sample:
dataJan04 <- dataJan04 %>%
  mutate(idSpecificity = case_when((Class == "NA" & 
                                      Order == "NA" & 
                                      Family == "NA" &
                                      Subfamily == "NA" &
                                      Genus == "NA" &
                                      Species == "NA") ~ "Phylum",
                                   (Order == "NA" & 
                                      Family == "NA" &
                                      Subfamily == "NA" &
                                      Genus == "NA" &
                                      Species == "NA") ~ "Class",
                                   (Family == "NA" &
                                      Subfamily == "NA" &
                                      Genus == "NA" &
                                      Species == "NA") ~ "Order",
                                   (Subfamily == "NA" &
                                      Genus == "NA" &
                                      Species == "NA") ~ "Family",
                                   (Genus == "NA" &
                                      Species == "NA") ~ "Subfamily",
                                   Species == "NA" ~ "Genus",
                                   Species != "NA" ~ "Species",
                                   TRUE ~ "other")) %>%
  mutate(idSpecificity = fct_relevel(idSpecificity, 
                                     "Species", 
                                     "Genus", 
                                     "Subfamily",
                                     "Family",
                                     "Order")) %>%
  mutate(idSpecificityNumeric = case_when(idSpecificity == "Phylum" ~ 0,
                                          idSpecificity == "Class" ~ 1,
                                          idSpecificity == "Order" ~ 2,
                                          idSpecificity == "Family" ~ 3,
                                          idSpecificity == "Subfamily" ~ 4,
                                          idSpecificity == "Genus" ~ 5,
                                          idSpecificity == "Species" ~ 6)) %>% # Some BINs have multiple IDs to multiple levels of specificity; the following code removes less-specific IDs.
  group_by(BIN) %>%
  mutate(idSpecificityNumericMax = max(idSpecificityNumeric)) %>%
  select(-c(idSpecificity)) %>%
  mutate(idSpecificity = case_when(idSpecificityNumericMax == 0 ~ "Phylum",
                                   idSpecificityNumericMax == 1 ~ "Class",
                                   idSpecificityNumericMax == 2 ~ "Order",
                                   idSpecificityNumericMax == 3 ~ "Family",
                                   idSpecificityNumericMax == 4 ~ "Subfamily",
                                   idSpecificityNumericMax == 5 ~ "Genus",
                                   idSpecificityNumericMax == 6 ~ "Species")) %>%
  ungroup()

dataJan04$idSpecificity <- fct_relevel(dataJan04$idSpecificity, 
                                       "Species", 
                                       "Genus", 
                                       "Subfamily",
                                       "Family",
                                       "Order",
                                       "Class")

binCounts <- dataJan04 %>%
  mutate(`Most Specific Identification` = case_when(idSpecificityNumericMax == 0 ~ paste(Phylum, "sp."),
                                                    idSpecificityNumericMax == 1 ~ paste(Class, "sp."),
                                                    idSpecificityNumericMax == 2 ~ paste(Order, "sp."),
                                                    idSpecificityNumericMax == 3 ~ paste(Family, "sp."),
                                                    idSpecificityNumericMax == 4 ~ paste(Subfamily, "sp."),
                                                    idSpecificityNumericMax == 5 ~ paste(Genus, "sp."),
                                                    idSpecificityNumericMax == 6 ~ Species)) %>%
  select(c(BIN,
           Order,
           `Most Specific Identification`)) %>%
  group_by(BIN) %>%
  dplyr::summarise(Order = first(Order),
                   `Most Specific Identification` = first(`Most Specific Identification`),
                   count = n()) %>%
  distinct() %>%
  dplyr::arrange(desc(count))

write_csv(binCounts,
          "supplementaryFile_binCounts.csv")

#### Supplemental figure for uniqueness of species across collecting methods ####
#### Panel A, upset plot ####
numberMalaiseSpecimens <- length(which(dataJan04$`Sampling Protocol` == "Malaise Trap"))
numberNonmalaiseSpecimens <- length(which(dataJan04$`Sampling Protocol` != "Malaise Trap"))
totalSpecimens <- numberMalaiseSpecimens + numberNonmalaiseSpecimens

# Get all combinations of BINs and collecting methods present in the real data:
samplingMethodDataPresent <- dataJan04 %>%
  mutate(`Sampling Protocol` = case_when(`Sampling Protocol` == "Malaise Trap" ~ "Malaise Trap",
                                         TRUE ~ "non-Malaise Trap")) %>%
  select(c(BIN,
           `Sampling Protocol`)) %>%
  distinct()
samplingMethodDataPresent$present <- TRUE

# Get all theoretically possible combinations of BIN and collecting method:
samplingMethodData <- dataJan04 %>%
  mutate(`Sampling Protocol` = case_when(`Sampling Protocol` == "Malaise Trap" ~ "Malaise Trap",
                                         TRUE ~ "non-Malaise Trap")) %>%
  select(c(BIN,
           `Sampling Protocol`)) %>%
  distinct() %>%
  complete(BIN,
           `Sampling Protocol`) %>%
  left_join(samplingMethodDataPresent,
            by = c("BIN" = "BIN",
                   "Sampling Protocol" = "Sampling Protocol")) %>%
  mutate(present = case_when(is.na(present) ~ "FALSE",
                             !is.na(present) ~ "TRUE"))
samplingMethodData$present <- as.logical(samplingMethodData$present)

# Convert that to wide format for ComplexUpset:
samplingMethodDataWide <- samplingMethodData %>%
  pivot_wider(names_from = `Sampling Protocol`,
              values_from = present) %>%
  as.data.frame()

rownames(samplingMethodDataWide) <- samplingMethodDataWide$BIN

samplingMethodDataWide <- samplingMethodDataWide %>%
  dplyr::select(-c("BIN"))

plotData <- samplingMethodDataWide

# Plot:
upsetPlot <- ComplexUpset::upset(plotData,
                                 unique(colnames(samplingMethodDataWide)),
                                 queries=list(upset_query(intersect = c("Malaise Trap", "non-Malaise Trap"),
                                                          color = '#8bcc95',
                                                          fill = '#8bcc95',
                                                          only_components = c('Intersection size')),
                                              upset_query(intersect = c("Malaise Trap"),
                                                          color = '#c5e3bc',
                                                          fill = '#c5e3bc',
                                                          only_components = c('Intersection size',
                                                                              'set_sizes')),
                                              upset_query(intersect = c("non-Malaise Trap"),
                                                          color = '#027822',
                                                          fill = '#027822',
                                                          only_components = c('Intersection size',
                                                                              'set_sizes')),
                                              upset_query(set = "Malaise Trap", 
                                                          fill = '#c5e3bc'),
                                              upset_query(set = "non-Malaise Trap", 
                                                          fill = '#027822'))) +
  ggtitle("Sampling overlap across different\ncollecting methods")

upsetPlot

# Specific numbers:
# How many total bins?
totalBins <- length(samplingMethodDataWide$`Malaise Trap`)

# How many BINs were detected by the malaise traps (could also be detected by other methods)?
malaiseBins <- samplingMethodDataWide %>%
  filter(`Malaise Trap` == TRUE)

# How many BINs were detected ***only*** by the malaise traps?
malaiseOnlyBins <- samplingMethodDataWide %>%
  filter(`Malaise Trap` == TRUE) %>%
  filter(if_all(c(-`Malaise Trap`), 
                ~ . == FALSE))

# How many BINs were detected ***only*** via methods other than malaise traps?
nonMalaiseOnlyBins <- samplingMethodDataWide %>%
  filter(`Malaise Trap` == FALSE) 

# How many BINs were detected by all of the non-malaise trap methods (could have also been detected by the malaise traps)?
nonMalaiseBins <- samplingMethodDataWide %>%
  filter(if_any(c(-`Malaise Trap`), 
                ~ . == TRUE))

#### Panel B, what is the taxonomic breakdown of each type of collecting? ####
test <- dataJan04 %>%
  mutate(`Sampling Protocol` = case_when(`Sampling Protocol` == "Malaise Trap" ~ "Malaise Trap",
                                         TRUE ~ "non-Malaise Trap"))

getProportionsByMethod <- function(specificOrder) {
  nonMalaise <- test %>%
    filter(`Sampling Protocol` == "non-Malaise Trap")
  proportionInNonMalaise <- length(which(nonMalaise$Order == specificOrder)) / length(nonMalaise$Order)
  
  Malaise <- test %>%
    filter(`Sampling Protocol` == "Malaise Trap")
  proportionInMalaise <- length(which(Malaise$Order == specificOrder)) / length(Malaise$Order)
  
  difference <- (proportionInNonMalaise - proportionInMalaise) 
  proportions <- c(specificOrder,
                   proportionInNonMalaise,
                   proportionInMalaise,
                   difference)
  
  return(proportions)
}

proportionsOrders <- purrr::map(unique(dataJan04$Order),
                                getProportionsByMethod)
proportionsOrders <- as.data.frame(do.call(rbind, proportionsOrders)) %>%
  as.data.frame()
colnames(proportionsOrders) <- c("specificOrder",
                                 "proportionInNonMalaise",
                                 "proportionInMalaise",
                                 "difference (non-Malaise minus Malaise)")

proportionsOrders$proportionInNonMalaise <- as.numeric(proportionsOrders$proportionInNonMalaise)
proportionsOrders$proportionInMalaise <- as.numeric(proportionsOrders$proportionInMalaise)
proportionsOrders$`difference (non-Malaise minus Malaise)` <- as.numeric(proportionsOrders$`difference (non-Malaise minus Malaise)`)

proportionsOrders <- proportionsOrders %>%
  arrange(-`difference (non-Malaise minus Malaise)`)
proportionsOrders$specificOrder <- base::factor(proportionsOrders$specificOrder,
                                                levels = proportionsOrders$specificOrder)

# Plot those breakdowns:
proportionsOrdersLong <- proportionsOrders %>%
  select(-c(`difference (non-Malaise minus Malaise)`)) %>%
  pivot_longer(c("proportionInNonMalaise",
                 "proportionInMalaise"),
               names_to = "method",
               values_to = "proportion") %>%
  mutate(majorOrder = case_when(specificOrder == "Diptera" ~ "Major orders",
                                specificOrder == "Hemiptera" ~ "Major orders",
                                specificOrder == "Coleoptera" ~ "Major orders",
                                specificOrder == "Hymenoptera" ~ "Major orders",
                                specificOrder == "Lepidoptera" ~ "Major orders",
                                TRUE ~ "Minor orders"))

ordersByMethod <- ggplot(data = proportionsOrdersLong,
                         mapping = aes(x = specificOrder, 
                                       y = proportion,
                                       fill = method)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = c(proportionInNonMalaise = "#027822",
                               proportionInMalaise = "#c5e3bc")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  facet_wrap(~ majorOrder,
             scales = "free") +
  labs(x = "Order",
       y = "Proportion of specimens",
       title = "Proportion of specimens across collecting methods")

ordersByMethod

#### Combine panels ####
combinedPlots <- (upsetPlot + theme(plot.margin = unit(c(0, 0, 60, 0), "pt"))) / ordersByMethod 

combinedPlots

ggsave("./supplementalFigureOnCollectingMethods.png",
       width = 12,
       height = 12,
       units = "in",
       dpi = 900)

ggsave("./supplementalFigureOnCollectingMethods.pdf",
       width = 12,
       height = 12,
       units = "in",
       dpi = 900)





