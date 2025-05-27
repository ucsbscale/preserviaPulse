library(SSDM)
library(raster)
library(here)
library(spatialEco)
library(terra)
library(dismo)
library(rJava)


occ_pre <- read.csv('data/occurrences/animals/animals-data-ready-occ-0519.csv')
colnames(occ_pre)
occ_pre <- na.omit(occ_pre)
occ_pre$Presence <- 1

# select species
i <- 1
target_sp <- unique(occ_pre$species)[i]
print(target_sp)

occ_pre_target <- subset(occ_pre, occ_pre$species == target_sp)
head(occ_pre_target)


##### Abs data
occ_abs <- read.csv("data/occurrences/sampled_background_points.csv")
occ_abs <- na.omit(occ_abs)
occ_abs <- occ_abs[,2:3]
occ_abs$species <- target_sp
occ_abs$Presence <- 0

# merge
occ_target <- rbind(occ_pre_target[,-2], occ_abs)


# --------------------- read in env. vars --------------
Env = rast("data/env/final_env_1980_2010_stack.tif")

# z-score normalization
Env_z = raster.Zscore(Env)

Env_z = stack(Env_z)
Env_z

##### Sampling weight
prNum <- sum(occ_target$Presence == 1) # number of presence records
bgNum <- sum(occ_target$Presence == 0)
wt <- ifelse(occ_target$Presence == 1, 1, prNum / bgNum)
spsize <- c("0" = prNum, "1" = prNum)

set.seed(1)

SDM_ens <- ensemble_modelling(c('GAM','MAXENT'),
                              occ_target, Env_z, Xcol = 'x', Ycol = 'y', Pcol="Presence",
                              cv = "holdout", cv.param = c(0.7, 10),
                              SDM.projections = TRUE,
                              weight = TRUE,
                              save = TRUE,
                              path = file.path('mod/test'), # a folder for each species
                              rep = 1, 
                              # cores = 4,
                              axes.metric = "AUC",
                              gam.args = list(family = binomial(link = "logit"),
                                              weights = wt,
                                              method = "REML"),
                              cta.args = list(ntree = 500,
                                              sampsize = spsize,
                                              replace = TRUE))

length(wt)
#bin.thresh
#SDM.projections


SDM
plot(SDM@projection, main = 'SDM\nwith GLM algorithm')

library(SSDM)
a <- readRDS('data/Agelaiustricolor_ntree500.rds')
b <- readRDS('data/Agelaiustricolor_ntree2.RDS')
plot(a@uncertainty)
a@evaluation
b@evaluation
