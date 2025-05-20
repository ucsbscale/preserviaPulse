## Purpose of script: make sure final plants and animals occ points have the same format for SSDM
## Authors: GEOG 274
## Date: Spring, 2025

library(here)
library(dplyr)
library(googledrive)
library(googlesheets4)
# --------------------- 1. Post process all animals ----------------------
# note: herps = herps + inverts
## -------------------- 1.1 read in all files -------------------
anim_final_files <- list.files(path = here('data/occurrences/animals'), pattern = '*-cleaned-0519.csv')

anim_df <- data.frame(matrix(ncol=3, nrow=0))
colnames(anim_df) <- c('species', 'geometry', 'taxon')

for(f in anim_final_files){
  thisdf <- read.csv(file.path('data/occurrences/animals', f), sep=';') %>%
    select('species', 'geometry')
  thistaxon = strsplit(f, '-')[[1]][1]
  thisdf$taxon = thistaxon
  anim_df <- rbind(anim_df, thisdf)
}

## -------------------- 1.2 convert geometry to x and y ------------------

options(digits = 16)
anim_df <- anim_df %>%
  mutate(
    coords = gsub("c\\(|\\)", "", geometry),  # remove c() wrapper
    x = as.double(sapply(strsplit(coords, ","), `[`, 1)),  # first coord
    y = as.double(sapply(strsplit(coords, ","), `[`, 2))   # second coord
  ) %>% select(-c(geometry, coords))

name_anim_df <- 'animals-data-ready-occ-0519.csv'
write.table(anim_df, file.path('data/occurrences/animals', name_anim_df), sep=',', row.names = FALSE)

## -------------------- 1.3 upload to gdrive -----------------------
drive_auth()
gs4_auth(token = drive_token())

speciesObs_folder <- drive_find(pattern = "speciesObs", type = "folder")
if(nrow(speciesObs_folder) > 0) {
  # Upload to Google Drive if folder found
  drive_upload(
    file.path('data/occurrences/animals', name_anim_df),
    path = as_id(speciesObs_folder$id[1]),
    name = name_anim_df,
    overwrite = TRUE
  )
} else {
  warning("Could not find 'specieObs' folder in Google Drive. File saved locally only.")
}
