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

# Fix some Observatory Site names:
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

# Info about observatory/snapshot is in dataJan04$`Project Code`; 
# If it starts with "CIO" or CIBIP, and in Field ID, if it starts with GMP:
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

# Create a dataframe summarizing counts for BIN/ID:
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

#### Figure 1: Taxonomic breakdown and abundance accumulation curve ####
# Get numbers on percent ID specificity: ####
percentClass <- length(which(dataJan04$idSpecificity == "Class")) / length(dataJan04$idSpecificity)
percentOrder <- length(which(dataJan04$idSpecificity == "Order")) / length(dataJan04$idSpecificity)
percentFamily <- length(which(dataJan04$idSpecificity == "Family")) / length(dataJan04$idSpecificity)
percentSubfamily <- length(which(dataJan04$idSpecificity == "Subfamily")) / length(dataJan04$idSpecificity)
percentGenus <- length(which(dataJan04$idSpecificity == "Genus")) / length(dataJan04$idSpecificity)
percentSpecies <- length(which(dataJan04$idSpecificity == "Species")) / length(dataJan04$idSpecificity)

# Panel A, Plot the specificity of BIN assignment on a per-specimen basis: #### 
specificityOfID <- ggplot(data = dataJan04, 
                          mapping = aes(x = idSpecificity)) +
  geom_bar(aes(y = (..count..)/sum(..count..))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  labs(x = "Identification specificity",
       y = "Percent of specimens") +
  scale_y_continuous(labels = percent,
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(limits = rev)

specificityOfID

# Panel B, abundance of each order: #### 
# For the manuscript, we are aggregating orders with fewer than ~2400 specimens:

# Get a dataframe of unique BINs:
# Process data to get unique BINs:
uniqueBINs <- dataJan04 %>%
  select(-c(`Sample ID`,
            `Project Code`,
            `Field ID`,
            "Identifier",
            "Collectors",
            "Collection Date.y",
            "Country/Ocean",
            "State/Province",
            "Region",
            "Sector",
            "Exact Site",
            "Lat",
            "Lon",
            "Elev",
            "Collection Date Accuracy",
            "Habitat",
            "Sampling Protocol",
            "Tribe",
            "Subspecies",
            "projectType",
            "US_L3NAME",
            "US_L4NAME",
            "EVT_NAME",
            "correctedSiteName",
            "preciseLat",
            "preciseLon",
            idSpecificityNumericMax)) %>%
  group_by(BIN) %>% 
  top_n(1, 
        idSpecificityNumeric) %>%
  distinct() 

uniqueBINsOrders <- uniqueBINs %>%
  select(c(BIN,
           Order)) %>%
  distinct()

ordersAggregated <- uniqueBINsOrders %>%
  group_by(Order) %>%
  dplyr::summarise(count = n()) %>%
  mutate(belowThreshold = case_when(count < 2400 ~ "Minor orders",
                                    TRUE ~ Order)) %>%
  group_by(belowThreshold) %>%
  dplyr::summarise(count = sum(count))

colnames(ordersAggregated) <- c("Order",
                                "Number of unique BINs")

ordersAggregated$Source <- "California Insect\nBarcoding Initiative"

# Do the same thing for GBIF data so that they can be plotted together:
gbifData <- readxl::read_xlsx("GBIF_all_Hexapoda_250821.xlsx") %>%
  filter(taxonRank == "SPECIES" &
           taxonomicStatus == "ACCEPTED") %>%
  select(c(taxonKey,
           scientificName,
           order)) %>%
  distinct() %>%
  group_by(order) %>%
  dplyr::summarise(count = n()) %>%
  mutate(belowThreshold = case_when(order %in% ordersAggregated$Order ~ order,
                                    TRUE ~ "Minor orders")) %>%
  group_by(belowThreshold) %>%
  dplyr::summarise(count = sum(count))

colnames(gbifData) <- c("Order",
                        "Number of unique BINs")

gbifData$Source <- "GBIF"

# Join CIBI and GBIF data:
allOrderAbundanceData <- rbind(ordersAggregated,
                               gbifData)

# Generate a color scheme:
abundanceColors <- c("#595959",
                     "#b3b3b3")
names(abundanceColors) <- c("California Insect\nBarcoding Initiative",
                            "GBIF")

# Plot:
aggregateOrders <- ggplot(data = allOrderAbundanceData, 
                          mapping = aes(x = fct_infreq(Order,
                                                       `Number of unique BINs`),
                                        y = `Number of unique BINs`,
                                        fill = Source)) +
  geom_col(position = position_dodge()) +
  theme_bw() +
  labs(x = "Order",
       y = "Count of unique species") +
  scale_y_continuous(limits = c(0,
                                17000),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = abundanceColors) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "bottom")

aggregateOrders

# Panel C, rank-abundance curve: #### 
# Plot rank abundances:
abundance <- dataJan04 %>%
  count(BIN) %>%
  left_join(uniqueBINsOrders,
            by = c("BIN" = "BIN")) %>%
  distinct() %>%
  arrange(desc(n))

abundance$rank <- as.numeric(rownames(abundance)) - 1

rankAbundance <- ggplot(data = abundance,
                        mapping = aes(x = rank,
                                      y = n)) + 
  geom_line() + 
  geom_point(size = 1,
             shape = 1) + 
  theme_bw() +
  scale_y_log10() +
  labs(x = "Species rank",
       y = "Abundance (log10 scale)")  

rankAbundance

# Combine the three panels for Figure 1 ####
figureOne <- (specificityOfID / aggregateOrders) | rankAbundance

figureOne <- figureOne + plot_annotation(tag_levels = 'A')

figureOne

ggsave("./cibiFigureOne.png",
       width = 9,
       height = 6,
       units = "in",
       dpi = 900)

ggsave("./cibiFigureOne.pdf",
       width = 9,
       height = 6,
       units = "in",
       dpi = 900)
