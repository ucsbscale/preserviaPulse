## Purpose of script: grab data from BIEN database
## Authors: GEOG 274
## Date: Spring, 2025
## Credits to: Xue Yan, Diyang Cui

################       1. Get species list   ########################
library(here)
library(googledrive)
library(googlesheets4)
# install.packages("RBIEN")
library(BIEN)
library(dplyr)
library(stringr) # for "fuzzy matching" and filtering our target species

# Set a directory for data
occ_dir <- here("data/occurrences") # create a data folder with relative path
file.path(occ_dir) # double check its absolute path
if (!dir.exists(occ_dir)) dir.create(occ_dir) # if it's not there already, create it

# Authorization
drive_auth()
gs4_auth(token = drive_token())

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

################       2. Get data from BIEN   ########################
# `BIEN_occurrence_county` Returns all occurrences records within a given state/province

country_vector <- rep("United States", 3)
state_vector <- rep("California", 3)
county_vector <- c("Santa Barbara", "Ventura", "San Luis Obispo")

# Get all records in the three counties
BIEN_occ <- BIEN_occurrence_county(
  country = country_vector,
  state = state_vector,
  county = county_vector)

# Filter to our target species
# NOTE: This is a loose match that does not account for specific cases such as synonyms.
# Logic: if the first two words match, then assume two names are the same.
filtered_occ <- BIEN_occ %>% 
  # Remove rows with NA in latitude or longitude
  filter(!is.na(latitude) & !is.na(longitude)) %>% 
  # Split each species name and get the first two words
  # Check if any name in scientific_names contains the same first two words
  filter(tolower(word(scrubbed_species_binomial, 1, 2)) %in% 
           tolower(word(scientific_names, 1, 2)))

raw_fname <- file.path(occ_dir, "BIEN_filtered_occurrences.csv")
write.csv(filtered_occ, raw_fname, row.names = FALSE)

# Get the target folder first to ensure it exists
speciesObs_folder <- drive_find(pattern = "speciesObs", type = "folder")
if(nrow(speciesObs_folder) > 0) {
  # Upload to Google Drive if folder found
  drive_upload(
    raw_fname,
    path = as_id(speciesObs_folder$id[1]),
    name = basename(raw_fname)
  )
} else {
  warning("Could not find 'specieObs' folder in Google Drive. File saved locally only.")
}
