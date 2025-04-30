## Purpose of script: grab BIEN data
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Dr. Diyang Cui & Dr. Lei Song 

# ------------- Setting up --------------
# Install and load packages
library(googledrive)
library(googlesheets4)
library(BIEN)
library(stringr) # for "fuzzy matching" and filtering our target species 
library(ggplot2)
library(sf)
library(dplyr)

# Set a directory for data
here() # first check path
occ_dir <- here("data/occurrences") # create a data folder with relative path
file.path(occ_dir) # double check its absolute path
if (!dir.exists(occ_dir)) dir.create(occ_dir) # if it's not there already, create it


# Authoritization to Google Drive
drive_auth()
gs4_auth(token = drive_token())

#------------- 1. Get species list  -------------
# Get data
ss <- drive_get("Special Status Species")

## E.g. sheet names
ss_meta <- gs4_get(ss) 
sheet_names <- ss_meta$sheets$name

dat <- read_sheet(ss, sheet = sheet_names[1])
head(dat)

# Get species' names
scientific_names_df <- dat[, c("Name (common)", "Name (Latin)")]

# Use Latin name for searching
scientific_names <- scientific_names_df$`Name (Latin)`

#-------------2. Get data from BIEN -------------
# Get boundaries for areas of interest
country_vector<-c("United States","United States","United States")
state_vector<-c("California","California","California")
county_vector<-c("Santa Barbara","Ventura","San Luis Obispo")

# Get all records in the three counties
BIEN_occ<-BIEN_occurrence_county(country=country_vector, state = state_vector, county = county_vector)

# Filter to our target species
filtered_occ <- BIEN_occ %>%
  filter(
    sapply(scrubbed_species_binomial, function(x) {
      # Split each species name and get the first two words
      scrubbed_first_two <- str_split(x, " ", simplify = TRUE)[1:2]
      # Create a regular expression to match the first two words
      regex <- paste(scrubbed_first_two, collapse = " ")
      # Check if any name in scientific_names contains the same first two words
      any(str_detect(scientific_names, regex))
    }) &
      !is.na(latitude) & !is.na(longitude)  # Remove rows with NA in latitude or longitude
  )

# Save out the downloaded data
write.csv(filtered_occ, 
          file.path(occ_dir, "BIEN_3counties_occurrences.csv"),
          row.names = FALSE)

# Upload data to google drive
# Get the target folder first to ensure it exists
speciesObs_folder <- drive_find(pattern = "speciesObs", type = "folder")
if(nrow(speciesObs_folder) > 0) {
  # Upload to Google Drive if folder found
  drive_upload(
    file.path(occ_dir, "BIEN_3counties_occurrences.csv"),
    path = as_id(speciesObs_folder$id[1]),
    name = "BIEN_3counties_occurrences.csv"
  )
} else {
  warning("Could not find 'specieObs' folder in Google Drive. File saved locally only.")
}

#-------------3. Exclude records from GBIF -------------
BIEN_occ <- filtered_occ

# Or read the saved BIEN data directly
# BIEN_occ<-read.csv("replace to your path/BIEN_3counties_occurrences.csv")

# Get datasource and summarize the information
unique(BIEN_occ$datasource)

BIEN_occ %>% group_by(datasource) %>%
  summarize(n=n())

# Exclude records from GBIF
BIEN_occ_subset<-BIEN_occ %>% filter(datasource!="GBIF")

#-------------4. Merge records from GBIF and BIEN -------------
# Read the saved GBIF data 
GBIF_occ_clean<-read.csv("replace to your path/Gbif-plant-clean.csv")

# Check column names from the two datasets 
# and rename them so that they match each other
colnames(GBIF_occ_clean)
colnames(BIEN_occ_subset)
colnames(BIEN_occ_subset)[1:4]=c("species","decimalLatitude","decimalLongitude","eventDate")

# Merge the two datasets
df<-bind_rows(BIEN_occ_subset,GBIF_occ_clean)

# View duplicates
df %>%
  group_by(species, decimalLatitude, decimalLongitude) %>%
  filter(n() > 1)

# Remove duplicates
df_unique <- df %>%
  distinct(species, decimalLatitude, decimalLongitude, .keep_all = TRUE)

# Check the records number for species
df_unique_stat <- df_unique %>%
  group_by(species) %>%
  summarize(n=n())

# Plot the results
ggplot(df_unique_stat, aes(x = reorder(species, -n), y = n)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Species Occurrence Counts",
    x = "Species",
    y = "Number of Records"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = element_text(size = 14)
  )

#-------------5. Keep only one records for a species in a given buffer -------------
# Convert the merged dataset to sf format (spatial data)
df_sf <- df_unique %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Projecting the data
df_sf <- st_transform(df_sf, crs = 32611)

# Calculate 30 m buffer for each record (can change to other resolution)
buffered <- st_buffer(df_sf, dist = 30)

# Merge overlapping buffered points 
# Records with the same group ID are within 30 meters of each other
# Retain only one records for per group
df_sf$group <- as.integer(st_within(df_sf, buffered, sparse = FALSE) %>% 
                            apply(1, function(x) which(x)[1]))

# Convert back to Geographic Coordinate System (GCS)
df_sf_latlon <- st_transform(df_sf, crs = 4326)

df_unique_buffer <- df_sf_latlon %>%
  group_by(species, group) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    longitude = st_coordinates(geometry)[, 1],
    latitude = st_coordinates(geometry)[, 2]
  ) %>%
  st_drop_geometry()

# Save out the downloaded data
write.csv(df_unique_buffer, 
          file.path(occ_dir, "Plant-Merge_Gbif_BIEN_3counties_occurrences.csv"),
          row.names = FALSE)

# Check final records number for species
df_unique_stat <- df_unique_buffer %>%
  group_by(species) %>%
  summarize(n=n())
write.csv(df_unique_stat, 
          file.path(occ_dir, "Plant_occurrences_stat.csv"),
          row.names = FALSE)

# Upload data to google drive
# Get the target folder first to ensure it exists
speciesObs_folder <- drive_find(pattern = "speciesObs", type = "folder")
if(nrow(speciesObs_folder) > 0) {
  # Upload to Google Drive if folder found
  drive_upload(
    file.path(occ_dir, "Plant_occurrences_stat.csv"),
    path = as_id(speciesObs_folder$id[1]),
    name = "Plant_occurrences_stat.csv"
  )
} else {
  warning("Could not find 'specieObs' folder in Google Drive. File saved locally only.")
}