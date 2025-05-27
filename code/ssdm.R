library(SSDM)
library(raster)
library(terra)
library(here)
library(spatialEco)
library(googledrive)
library(googlesheets4)
library(dplyr)
# ----------------- prep ---------------------

setwd('Desktop/wenxinyang/github/SDM2025')

## -------------- read in and normalize env. data --------------
Env = rast("final_env_1980_2010_stack.tif")

# z-score normalization
means <- cellStats(Env, "mean", na.rm=TRUE)
sds   <- cellStats(Env, "sd",   na.rm=TRUE)

Env_z <- brick(lapply(seq_len(nlayers(Env)), function(i) {
  (Env[[i]] - means[i]) / sds[i]
}))

names(Env_z) <- names(Env)

## ---------------- prep occ data ---------------
occ_pre <- read.csv('Anim_Plant_merge.csv')
occ_pre <- na.omit(occ_pre)
occ_pre$Presence <- 1

## ---------------- prep abs data ----------------
occ_abs <- read.csv("sampled_background_points.csv")
occ_abs <- na.omit(occ_abs)
occ_abs <- occ_abs[,2:3]
occ_abs$species <- target_sp
occ_abs$Presence <- 0

## ---------------- get species name list ------------
drive_auth()
gs4_auth(token = drive_token())
ss <- drive_get("Special Status Species")
dat <- read_sheet(ss, sheet="Birds") %>% 
  filter(`Where Listed?` == 'IRMP') %>% 
  select(`Name (Latin)`, `Name (common)`)

li_species <- unique(dat$`Name (Latin)`)
li_species <- li_species[li_species %in% unique(occ_pre$species)]

# --------------- SDM building ------------------
fitSDM <- function(target_sp){
  cat("====== start working on: ", target_sp, "======\n")
  occ_pre_target <- subset(occ_pre, occ_pre$species == target_sp)
  # merge presence and occurrence data
  occ_target <- rbind(occ_pre_target[,-2], occ_abs)
  # occ_target$Presence <- as.factor(occ_target$Presence)
  
  ##### Set sampling weight
  prNum <- sum(occ_target$Presence == 1) # number of presence records
  bgNum <- sum(occ_target$Presence == 0)
  wt <- ifelse(occ_target$Presence == 1, 1, prNum / bgNum)
  spsize <- c("0" = bgNum, "1" = prNum)
  
  sp_dir <- file.path('mod/test', target_sp)
  dir.create(sp_dir)
  
  set.seed(1)
  SDM_ens <- ensemble_modelling(c('GAM','MAXENT','RF'),
                                occ_target, Env_z, Xcol = 'x', Ycol = 'y', Pcol="Presence",
                                cv = "holdout", cv.param = c(0.7, 10),
                                rep = 1, cores = 10,
                                parmode = "algorithms",
                                save = TRUE, path = sp_dir,
                                #uncertainty = TRUE,
                                SDM.projections = TRUE, # save individual SDM models
                                ensemble.metric = c("AUC"),
                                ensemble.thresh = c(0.5),
                                # axes.metric = "AUC", # access variable importance
                                #final.fit.data = "all", # controls which data go into the final fits, "all" uses every record
                                #bin.thresh = "SES", # thresholding rule, "SES": sensitivity–specificity equality
                                #weight = TRUE, # models are weighted by their performance scores; otherwise each selected model contributes equally
                                gam.args = list(family = binomial(link = "logit"),
                                                weights = wt,
                                                method = "REML"),
                                cta.args = list(ntree = 500,
                                                sampsize = spsize,
                                                replace = TRUE))
  
  saveRDS(SDM_ens, file.path('mod', paste0(target_sp, '.rds')))
  cat('====== finished for species: ', target_sp, '=======\n')
  
}

lapply(li_species, fitSDM)