library(mgcv)
library(dismo)
library(randomForest)
library(raster)
library(pROC)
library(dplyr)
library(caret)
library(doParallel)

rm(list=ls())
setwd("C:/Users/yzhan/Desktop/courses/geog_274")



# --- Inputs ---
##### Occ data
occ_pre <- read.csv('Anim_Plant_merge_final0616.csv')
occ_pre <- na.omit(occ_pre)
occ_pre$Presence <- 1


##### env_stack: RasterStack of predictors
Env <- stack("full_env_1980_2010_stack.tif")
#Env <- Env[[c(1:3,5:12,19)]] 

Env <- Env[[c(6,10,11,12,17,7,18,9)]] 

names(Env)
plot(Env)

# z score
means <- cellStats(Env, "mean", na.rm=TRUE)
sds   <- cellStats(Env, "sd",   na.rm=TRUE)

env_stack <- brick(lapply(seq_len(nlayers(Env)), function(i) {
  (Env[[i]] - means[i]) / sds[i]
}))

names(env_stack) <- names(Env)

names(env_stack)
plot(env_stack)


##### species list
# plant list
taxon <- "plant"
# shrub
species_lst <- c("Arctostaphylos purissima", 
                 "Cirsium rhothophilum",
                 "Deinandra increscens villosa",
                 "Horkelia cuneata sericea",
                 "Juglans californica", 
                 "Malacothrix saxatilis saxatilis",
                 "Mucronea californica", 
                 "Ribes amarum hoffmannii",
                 "Eriodictyon capitatum",
                 "Astragalus nuttallii nuttallii",
                 "Senecio blochmaniae")

# herb
species_lst <- c("Arctostaphylos purissima",
                 "Astragalus nuttallii nuttallii",
                 "Deinandra increscens villosa",
                 "Horkelia cuneata sericea",
                 "Phacelia hubbyi",
                 "Abronia maritima", 
                 "Cirsium rhothophilum",
                 "Scrophularia atrata", 
                 "Mucronea californica")

# check stats
species_lst %in% unique(occ_pre$species)

occ_pre %>% filter(species %in% species_lst) %>% group_by(species) %>% summarize(obs=n())


##### paralel
#cl = makeCluster(10)
#registerDoParallel(cl)

#foreach( j=1:length(species_lst),.packages=c("raster","mgcv","dismo","randomForest","pROC","dplyr","caret","rsample") ) %dopar% {
for(j in 1:length(species_lst)){
  
  # extract target species
  #target_sp <- "Juglans californica"
  target_sp <- species_lst[j]
  print(target_sp)
  
  occ_pre_target <- subset(occ_pre, occ_pre$species == target_sp)
  head(occ_pre_target)
  
  
  ##### Abs data
  occ_abs <- read.csv("sampled_background_points.csv")
  occ_abs <- na.omit(occ_abs)
  occ_abs <- occ_abs[,2:3]
  occ_abs$species <- target_sp
  occ_abs$Presence <- 0
  
  
  # ---  GAM and RF formula ---
  var_names <- names(env_stack)
  gam_formula <- as.formula(paste("Presence ~", paste0("s(", var_names, ")", collapse = " + ")))
  rf_formula  <- as.formula(paste("as.factor(Presence) ~", paste(var_names, collapse = " + ")))
  
  
  ##### all data model
  occ_data <- rbind(occ_pre_target[,-2], occ_abs)
  
  env_all <- extract(env_stack, occ_data[, c("x", "y")])
  
  all_data <- na.omit(cbind(occ_data, env_all))
  
  # --- GAM ---
  prNum_all <- sum(all_data$Presence == 1) # number of presence records
  bgNum_all <- sum(all_data$Presence == 0)
  wt_all <- ifelse(all_data$Presence == 1, 1, prNum_all / bgNum_all)
  
  gam_model_all <- gam(gam_formula, data = all_data, 
                       family = binomial(link = "logit"), weights = wt_all, method = "REML")
  
  # prob output
  gam_pred_all <- predict(env_stack, gam_model_all, type = "response")
  plot(gam_pred_all)
  
  # stats
  pred_all_gam <- predict(gam_model_all, all_data, type = "response")
  roc_all_gam <- roc(all_data$Presence, pred_all_gam)
  auc_all_gam <- auc(roc_all_gam)
  thresh_all_gam <- coords(roc_all_gam, "best")[1]
  
  # binary output
  binary_all_gam <- gam_pred_all >= as.numeric(thresh_all_gam)
  plot(binary_all_gam)
  
  
  # --- MaxEnt ---
  #maxent_model_all0 <- maxent(env_stack, occ_data[occ_data$Presence == 1, c("x", "y")])
  #maxent_pred_all0 <- predict(env_stack, maxent_model_all0)
  #plot(maxent_pred_all0)
  
  maxent_model_all <- maxent(x=all_data[,-(1:4)], p=all_data$Presence)
  
  # prob output
  maxent_pred_all <- predict(env_stack, maxent_model_all)
  plot(maxent_pred_all)
  
  # stats
  pred_all_maxent <- predict(maxent_model_all, all_data)
  roc_all_maxent <- roc(all_data$Presence, pred_all_maxent)
  auc_all_maxent <- auc(roc_all_maxent)
  thresh_all_maxent <- coords(roc_all_maxent, "best")[1]
  
  # binary output
  binary_all_maxent <- maxent_pred_all >= as.numeric(thresh_all_maxent)
  plot(binary_all_maxent)
  
  
  
  # --- RF ---
  spsize_all <- c("0" = prNum_all, "1" = prNum_all)
  
  set.seed(274)
  rf_model_all <- randomForest(rf_formula, data = all_data, 
                               importance=TRUE, ntree = 1000, sampsize = spsize_all, replace = TRUE)
  
  # prob output
  rf_pred_all <- predict(env_stack, rf_model_all, type = "prob", index = 2)
  plot(rf_pred_all)
  
  # stats
  pred_all_rf <- predict(rf_model_all, all_data, type = "prob")[,2]
  roc_all_rf <- roc(all_data$Presence, pred_all_rf)
  auc_all_rf <- auc(roc_all_rf)
  thresh_all_rf <- coords(roc_all_rf, "best")[1]
  
  # binary output
  binary_all_rf <- rf_pred_all >= as.numeric(thresh_all_rf)
  plot(binary_all_rf)
  
  
  
  
  # --- Variable Importance ---
  ### GAM
  summary(gam_model_all)
  
  # Extract chi-squared statistics
  gam_imp <- summary(gam_model_all)$s.table[, "Chi.sq"]
  
  # Normalize to percent
  gam_imp_percent <- 100 * gam_imp / sum(gam_imp)
  
  # Make a data frame
  gam_varimp <- data.frame(
    Variable = names(gam_imp),
    Importance = gam_imp_percent
  )
  
  # sort  
  gam_varimp <- gam_varimp[order(-gam_varimp$Importance), ]
  print(gam_varimp)
  
  ### MaxEnt
  maxent_varimp <- maxent_model_all@results
  
  # Extract only the variable contributions
  maxent_contrib <- maxent_varimp[grep("contribution", rownames(maxent_varimp)), ,drop = FALSE]
  
  # Clean up names
  var_names <- gsub(".contribution", "", rownames(maxent_contrib))
  maxent_varimp <- data.frame(
    Variable = var_names,
    Importance = maxent_contrib[, 1]
  )
  
  # sort
  maxent_varimp <- maxent_varimp[order(-maxent_varimp$Importance), ]
  print(maxent_varimp)
  
  ### RF
  rf_imp <- importance(rf_model_all)
  
  # Turn into data frame
  rf_varimp <- data.frame(
    Variable = rownames(rf_imp),
    Importance = rf_imp[, 1]  # Usually %IncMSE or IncNodePurity
  )
  
  # sort
  rf_varimp <- rf_varimp[order(-rf_varimp$Importance), ]
  print(rf_varimp)
  
  
  
  
  # --- Output ---
  model_list <- list(
    gam = gam_model_all,
    gam_auc = auc_all_gam,
    gam_thresh = thresh_all_gam,
    gam_varimp = gam_varimp,
    gam_prob = gam_pred_all,
    gam_binary = binary_all_gam,
    
    maxent = maxent_model_all,
    maxent_auc = auc_all_maxent,
    maxent_thresh = thresh_all_maxent,
    maxent_varimp = maxent_varimp,
    maxent_prob = maxent_pred_all,
    maxent_binary = binary_all_maxent,
    
    rf = rf_model_all,
    rf_auc = auc_all_rf,
    rf_thresh = thresh_all_rf,
    rf_varimp = rf_varimp,
    rf_prob = rf_pred_all,
    rf_binary = binary_all_rf
    
  )
  
  saveRDS(model_list, file = paste0("new/",taxon,"/",target_sp, ".RDS"))
  
  
  
}






##### check
p=1
species_lst[p]
test=readRDS(paste0("new/",taxon,"/",species_lst[p], ".RDS"))
#test2=readRDS(paste0(species_lst[p], ".RDS"))
#test$gam
#test$maxent
#test$rf



test$gam_thresh
test$rf_thresh
test$maxent_thresh

plot(test$rf_prob)
plot(test$maxent_prob)
plot(test$gam_prob)

plot(test$rf_binary)
plot(test$maxent_binary)
plot(test$gam_binary)


plot(test2@sdms[[2]]@projection) #maxent
plot(test2@sdms[[2]]@binary) #maxent





#### save tif (avg prob + binary on majority vote)

rds_file=list.files(path="./new",pattern="*.RDS",recursive=T,full.names=T)
rds_file

for(q in 1:length(rds_file)){
  
  print(q)
  
  rds_full=readRDS(rds_file[q])
  
  # ensemble (prob avg)
  all_prob = (rds_full$gam_prob + rds_full$maxent_prob + rds_full$rf_prob)/3
  plot(rast(all_prob))
  writeRaster(all_prob, 
              paste0(substring(rds_file[q],1,nchar(rds_file[q])-4),"_ensprob.tif"),
              overwrite=T
              )
  
  # ensemble (majority vote)
  all_bin_avg = (rds_full$gam_binary + rds_full$maxent_binary + rds_full$rf_binary)/3
  
  all_bin <- ifel(rast(all_bin_avg) > 0.5, 1, 0)
  plot(all_bin)
  writeRaster(all_bin, 
              paste0(substring(rds_file[q],1,nchar(rds_file[q])-4),"_ensbin.tif"),
              overwrite=T
  )
  
}





#####
# stack
bin_tif=list.files(path="./new/herp_invert",pattern="*._ensbin.tif",recursive=T,full.names=T)

bin_stack <- rast(bin_tif)

# calculate pixel-wise sum (i.e., diversity)
diversity_map <- sum(bin_stack, na.rm = TRUE)

# Plot the result
plot(diversity_map, main = "Species Richness (Diversity)")

writeRaster(diversity_map,filename="./new/herp_invert/herp_invert_stack.tif",overwrite=T)




#####
# SSDM
#
j=1
species_lst[j]
a=readRDS(paste0("results_plants/models/",species_lst[j],".RDS"))
plot(a@projection)
plot(a@uncertainty) #different colors

all(values(a@sdms[[2]]@projection)==values(maxent_pred_all))
identical(a@sdms[[2]]@projection,maxent_pred_all)
plot(a@sdms[[2]]@projection)

plot(a@sdms[[1]]@projection)

