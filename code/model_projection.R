library(SSDM)
library(terra)
library(dplyr)
library(raster)
select <- dplyr::select

Env <- stack("/Users/leisong/Downloads/Stack_Env/final_env_1980_2010_stack.tif")

occ <- read.csv("/Users/leisong/Downloads/Anim_Plant_merge.csv") %>% 
  rename("SPECIES" = species, 
         "LONGITUDE" = x, "LATITUDE" = y) %>% select(-taxon) %>% 
  filter(SPECIES == "Agelaius tricolor")

mods <- ensemble_modelling(c('RF', "GAM", "MAXENT"), occ, Env, 
  Xcol = 'LONGITUDE', Ycol = 'LATITUDE', verbose = TRUE)

