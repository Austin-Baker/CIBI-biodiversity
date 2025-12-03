library(readxl)
library(dplyr)
library(ggplot2)
library(SpadeR)

# Read in BOLD datasheet
## UPDATE ##
BOLD_data <- read_excel("/Users/austinbaker/Desktop/CIBI_Biodiversity_Manuscript/250626_BOLD_data_clean.xlsx")

# Create list of CIBI codes for a given Observatory Site
## UPDATE ##
OS_list <- c("GMP#55427",	"GMP#55431",	"GMP#55435",	"GMP#55439",	"GMP#55444",	"GMP#55448",	"GMP#55452",	"GMP#55457",	"GMP#55461",	"GMP#55466",	"GMP#55470",	"GMP#55474")

# Create new dataframe with only records matching codes above
OS_data <- BOLD_data[BOLD_data$`Field ID` %in% OS_list, ]

# Ensure that 'Collection Date.y' is in Date format
OS_data$`Collection Date.y` <- as.Date(OS_data$`Collection Date.y`, format = "%Y-%m-%d")

# Sort dataframe by date collected
OS_data_sorted <- OS_data[order(OS_data$`Collection Date.y`), ]


#############################
#### Summary Table ####
# Get unique 'Field ID' values from OS_data_sorted
unique_field_ids <- unique(OS_data_sorted$`Field ID`)

# Create new dataframe to record statistics: # sampling days, # total specimens, # NA, # unique BINs
summary_df <- data.frame(
  Field.ID = unique_field_ids,  # First column is the list of codes
  Collection.Date = as.Date(rep(NA, length(unique_field_ids))),  # Initialize with NAs of correct length
  Total.Specimens = rep(NA, length(unique_field_ids)),  # Placeholder for total specimen count
  Missing.Data = rep(NA, length(unique_field_ids)),  # Placeholder for missing BIN count
  Percent.Missing.Data = rep(NA, length(unique_field_ids)),  # Placeholder for % Missing Data
  Unique.Species = rep(NA, length(unique_field_ids)),  # Placeholder for unique species count (unique BINs)
  Accumulated.BINs = I(vector("list", length(unique_field_ids))),  # Initialize list of lists
  BINs.Accumulated = rep(NA, length(unique_field_ids)),  # Placeholder for count of accumulated BINs
  stringsAsFactors = FALSE  # Ensure strings are not converted to factors
)

#### Fill in Table ####
# Initialize an empty set to accumulate unique BINs over time
accumulated_bins <- character()

# Populate the 'Collection Date', 'Total Specimens', and 'Missing Data' columns
for (i in 1:length(unique_field_ids)) {
  code <- OS_list[i]
  
  # Filter rows for the current Field ID
  code_data <- OS_data_sorted[OS_data_sorted$`Field ID` == code, ]
  
  # If there are any records for this Field ID, copy the 'Collection Date.y', calculate total specimens, and count NA and "NA" in 'BIN'
  if (nrow(code_data) > 0) {
    summary_df$Collection.Date[i] <- code_data$`Collection Date.y`[1]
    summary_df$Total.Specimens[i] <- nrow(code_data)
    summary_df$Missing.Data[i] <- sum(code_data$BIN == "NA")
    summary_df$Percent.Missing.Data[i] <- (summary_df$Missing.Data[i] / summary_df$Total.Specimens[i]) * 100
    
    unique_bins <- unique(code_data$BIN[!is.na(code_data$BIN) & code_data$BIN != "NA"])
    
    accumulated_bins <- unique(c(accumulated_bins, unique_bins))
    
    summary_df$Unique.Species[i] <- length(unique_bins)
    summary_df$Accumulated.BINs[[i]] <- accumulated_bins
    summary_df$BINs.Accumulated[i] <- length(accumulated_bins)
  }
}

# Ensure the dates are displayed correctly in the summary dataframe
summary_df$`Collection.Date` <- as.Date(summary_df$`Collection.Date`, origin = "1970-01-01")

# Make an updated list of Field ID's with <20% missing data
# Filter the summary_df for rows with < 20% missing data
filtered_summary <- summary_df[summary_df$Percent.Missing.Data < 20, ]
print(filtered_summary)

# Create a new OS_list with only the Field IDs that meet the criteria
OS_list_filtered <- filtered_summary$Field.ID

OS_data_filtered <- BOLD_data[BOLD_data$`Field ID` %in% OS_list_filtered, ]

###############################
#### Chao estimator ####

BIN_count_OS <- OS_data_filtered %>% 
  group_by(BIN) %>%
  count()

richness_OS <- ChaoSpecies(BIN_count_OS$n,"abundance",k=10,conf=0.95)

print(richness_OS)