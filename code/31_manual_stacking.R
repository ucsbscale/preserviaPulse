## Purpose of script: manual stacking of binary layers
## Authors: GEOG 274
## Date: Spring, 2025

library(here)
library(sf)
library(terra) # read raster
library(dplyr)
library(rgbif)
library(googledrive)
library(googlesheets4)
library(tidyr)
library(ggplot2)

here() # first check path
mod_dir <- here("data/ModelOutputs")
list.dirs(mod_dir, recursive = FALSE, full.names = FALSE)
for(taxa in c('plants', 'birds', 'mammals', 'inv_herps')){
  cat('========= working on ', taxa, '=============\n')
  taxa_dir <- file.path(mod_dir, taxa, 'current')
  species_folders <- list.dirs(taxa_dir, recursive = FALSE, full.names = FALSE)
  
  for(i in 1:length(species_folders)){
    spp <- species_folders[i]
    cat('reading in ', spp, '\n')
    spp_folder <- file.path(taxa_dir, spp)
    if(taxa == 'birds'){
      bin_raster_path <- file.path(spp_folder, 'Specie', 'Rasters', 'Binary.tif')
      prob_raster_path <- file.path(spp_folder, 'Specie', 'Rasters', 'Probability.tif')
    }
    else{
      bin_raster_path <- file.path(spp_folder, 'Rasters', 'Binary.tif')
      prob_raster_path <- file.path(spp_folder, 'Rasters', 'Probability.tif')
    }
    bin_raster <- rast(bin_raster_path)
    prob_raster <- rast(prob_raster_path)
    if(i==1){
      stacked_bin <- bin_raster
      stacked_prob <- prob_raster
    } else {
      stacked_bin <- stacked_bin + bin_raster
      stacked_prob <- stacked_prob + prob_raster
    }
  }
  
  #plot(stacked_bin)
  #plot(stacked_prob)
  writeRaster(stacked_bin, file.path(taxa_dir, paste0(taxa, '_bin_stacked.tif')), overwrite=TRUE)
  writeRaster(stacked_prob, file.path(taxa_dir, paste0(taxa, '_prob_stacked.tif')), overwrite=TRUE)
  
}
