library(raster)
library(ENMeval)
library(dplyr)
library(ecospat)
library(caret)


# LOAD ENVIRONMENTAL PREDICTOR VARIABLES derived in script 1 
envs <- stack(list.files("Data/Derived/envs", full.names = TRUE))

# LOAD OCCURENCES RECORD AND FILTER 
load("Data/FINAL_RECS.Rdata")

# Without considering the islands, filter out the island using the geological substrate map 
GEO_ext <-raster::extract(envs[[1]], full[,c('LONGDEC','LATDEC')])
oc <- dplyr::filter(full, !is.na(GEO_ext)) #these are the cleaned occurrences

# Remove FIA data from occs used to fit models (since we compare w/abund at FIA points)
oc <- oc[oc$SOURCE!="FIA",]

# Express the geological substrate map as a factorial variable
envs[[1]] <- as.factor(envs[[1]])

# RUNNING THE MODELS using ENMeval (maxnet)
for (sp in 1:length(unique(oc$CODE))){
  
  # This line to deal with some species that don't work with hinge features...
  problem_spp <- c(275, 513, 539, 541) # (INGLAU, SURMAR, TETCRO, TETURB)

  # These lines are to to select the background points for each species which in our case are the all the other occurrence points 
  focsp <- unique(oc$CODE)[sp]
  occs <- oc[oc$CODE==focsp,2:1]
  bg <- oc[oc$CODE!=focsp,2:1]
  
  #to see on which species the model is working on! 
  message(paste('working on', focsp))
  print(Sys.time())
  
  if(nrow(occs)>10){
    
    if(nrow(occs)>=15){
      partition <- "checkerboard2" #if species have more than 15 occurence points use "checkerboard2" data partition
    } else {
      partition <- "jackknife" #otherwise use jackknife
    }
    
    if(!sp %in% problem_spp) {
      ps <- list(fc = c("L","LQ","LQH","H"), rm = 1:5)
    } else {
      ps <- list(fc = c("L","LQ"), rm = 1:5)
    }
    # HOW ARE 'DUPLICATE' RECORDS BEING TREATED? (AS IN, OCCS IN THE SAME GRID CELL?)
    mod <- ENMevaluate(occs=occs, envs=envs, bg=bg, 
                       algorithm='maxnet', 
                       partitions=partition,
                       categoricals="GEO",
                       partition.settings=list(aggregation.factor=c(5,5)),
                       tune.args=ps,
                       parallel=T)
    # tune.args=list(fc = c("L","H"), rm = 1:5))
    
    filename <- paste0("Data/2025-03-17_ENMeval_results-noFIA/", focsp, ".RDA")
    saveRDS(mod, file=filename)
  }
}

# MODEL SELECTION 

# Load mod files #mod = enmeval result object 
# modfiles <- list.files("Data/2022-12-07_ENMeval_results", full.names = TRUE)
modfiles <- list.files("Data/2025-03-17_ENMeval_results-noFIA", full.names = TRUE)

# Create an empty for the results
res <- list(length=length(modfiles))

# Select the models with the the lowest omission rate at 10p and the highest AUC
for(i in seq_along(modfiles)){
  print(i)
  
  #red in the model files mod
  mod <- readRDS(modfiles[i]) 

  # Select the models with the minimum omission rate and the maximum AUC
  tmpres <- mod@results %>%
    filter(or.10p.avg == min(or.10p.avg, na.rm=T)) %>%
    filter(auc.val.avg == max(auc.val.avg, na.rm=T))

  # Select the models with the continuous Boyce index
  # tmpres <- mod@results[which(mod@results$cbi.train == max(mod@results$cbi.train, na.rm=T)),][1,]
  
  # Take the full names (file path) and extract the just the code of the species 
  tmpres$species <- gsub(".RDA", "", strsplit(modfiles[i],"/")[[1]][3])
  
  res[[i]] <- tmpres[1,]
  
  # Take the prediction map of the selected model
  pred <- mod@predictions[[which(names(mod@predictions)==res[[i]]$tune.args)]]
  names(pred) <- tmpres$species[1]
  
  # Make a vector of 1/0 for occ/bg
  occ_bg <- c(rep(1, nrow(mod@occs)), rep(0, nrow(mod@bg)))
  
  # Make a vector of predicted values at occ and bg points
  pred_vals <- c(raster::extract(pred, mod@occs[,1:2]),
                 raster::extract(pred, mod@bg[,1:2]))
  
  # Find the threshold at max TSS
  tss_thresh <- ecospat.max.tss(pred_vals, occ_bg)$max.threshold
  
  # Impose the max TSS threshold
  pred_thresh <- pred > tss_thresh
  names(pred_thresh) <- tmpres$species[1]
  
  # Save the continuous and thresholded predictions as raster layers
  raster::writeRaster(pred, paste0("Data/SDM_predictions-noFIA/", names(pred)), format='GTiff')
  raster::writeRaster(pred_thresh, paste0("Data/SDM_threshold-noFIA/", names(pred_thresh)), format='GTiff')
}

# Table of all of the results 
resall <- do.call(rbind, res)

# Save model outputs 
write.csv(resall, "Data/Derived/SDM-mod_output-noFIA.csv")


