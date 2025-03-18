## TRY TO SEE IF FOREST AGE MAKES BETTER SDMs...

library(raster)
library(ENMeval)
library(dplyr)
library(ecospat)
library(caret)
library(sf)


# LOAD ENVIRONMENTAL PREDICTOR VARIABLES derived in script 1 
envs <- stack(list.files("Data/Derived/envs", full.names = TRUE))

# Load PR outline shapefile
pr <- st_read("Data/PR_shapes/outline/PR_outline_Project.shp")
pr <- st_transform(pr, CRS("+proj=lcc +lat_0=17.8333333333333 +lon_0=-66.4333333333333 +lat_1=18.0333333333333 +lat_2=18.4333333333333 +x_0=152400.3048 +y_0=0 +datum=NAD27 +units=m +no_defs"))

# Load forest maps
r51 <- raster("Data/Maps_1951-2000/puerto51_sub1_100905_landcov_final.img")
r77 <- raster("Data/Maps_1951-2000/puerto77_sub1_100905_landcov_urbveg_final.img")
r91 <- raster("Data/Maps_1951-2000/pr91_100805_final_quarry_recode2landcov_subset.img")
r00 <- raster("Data/Maps_1951-2000/pr2000_100805_final_quarry_recode2landcov_subset.img")

# Reproject some of r91 and r77 and r00 (with different extent and coordinates) and save them 
# rp77 <- projectRaster(r77, r51, method='ngb')
# writeRaster(rp77, "Data/Maps_1951-2000/Rpm_77.img")
rp77 <- raster("Data/Maps_1951-2000/Rpm_77.img")

# rp91<-projectRaster(r91, r51, method='ngb')
# writeRaster(rp91, "Data/Maps_1951-2000/Rpm_91.img")
rp91 <- raster("Data/Maps_1951-2000/Rpm_91.img")

# rp00 <-projectRaster(r00, r51, method='ngb')
# writeRaster(rp00, "Data/Maps_1951-2000/Rpm_00.img")
rp00 <- raster("Data/Maps_1951-2000/Rpm_00.img")

### RECLASSIFY FOREST AREAS (pixels that are covered in forest)
f51 <- r51 %in% 5
f77 <- rp77 %in% c(5,7)
f91 <- rp91 %in% c(5,7)
f00 <- rp00 %in% c(5,7)

# Stack and mask the forest classes
f <- mask(stack(f51, f77, f91, f00), pr)

age <- sum(f)
age_rp <- projectRaster(age, envs)
envs <- stack(envs, age_rp)
names(envs)[5] <- "Age"

envs <- crop(envs, st_read("Data/PR_shapes/outline/PR_outline_Project.shp"))




# LOAD OCCURENCES RECORD AND FILTER 
load("Data/FINAL_RECS.Rdata")

# Without considering the islands, filter out the island using the geological substrate map 
GEO_ext <-raster::extract(envs[[1]], full[,c('LONGDEC','LATDEC')])
oc <- dplyr::filter(full, !is.na(GEO_ext)) #these are the cleaned occurrences

# Express the geological substrate map as a factorial variable
envs[[1]] <- as.factor(envs[[1]])

# RUNNING THE MODELS using ENMeval (maxnet)
for (sp in 1:length(unique(oc$CODE))){
  
  # This line to skip some species that don't work for some reason...
  if(!sp %in% c(253,387,472,515,541)){
    
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
      
      # HOW ARE 'DUPLICATE' RECORDS BEING TREATED? (AS IN, OCCS IN THE SAME GRID CELL?)
      mod <- ENMevaluate(occs=occs, envs=envs, bg=bg, 
                         algorithm='maxnet', 
                         partitions=partition,
                         categoricals="GEO",
                         partition.settings=list(aggregation.factor=c(5,5)),
                         tune.args=list(fc = c("L","LQ","LQH","H"), rm = 1:5))
      # tune.args=list(fc = c("L","H"), rm = 1:5))
      
      filename <- paste0("Results/2022-12-07_ENMeval_results/", focsp, ".RDA")
      saveRDS(mod, file=filename)
    }
  }
}


mod2 <- ENMevaluate(occs=occs, envs=envs, bg=bg, 
                   algorithm='maxent.jar', 
                   partitions=partition,
                   categoricals="GEO",
                   partition.settings=list(aggregation.factor=c(5,5)),
                   tune.args=list(fc = "LQH", rm = 2))

# MODEL SELECTION 

# Load mod files #mod = enmeval result object 
modfiles <- list.files("Data/2022-12-07_ENMeval_results", full.names = TRUE)

# Create an empty for the results
res <- list()

# Select the models with the the lowest omission rate at 10p and the highest AUC
for(i in seq_along(modfiles)){
  print(i)
  
  #red in the model files mod
  mod <- readRDS(modfiles[i]) 
  
  # Select the models with the minimum omission rate and the maximum AUC
  tmpres <- mod@results %>%
    filter(or.10p.avg == min(or.10p.avg)) %>%
    filter(auc.val.avg == max(auc.val.avg))
  
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
  # raster::writeRaster(pred, paste0("Data/SDM_predictions/", names(pred)), format='GTiff')
  # raster::writeRaster(pred_thresh, paste0("Data/SDM_threshold-boyce/", names(pred_thresh)), format='GTiff')
}

# Table of all of the results 
resall <- do.call(rbind, res)

# Save model outputs 
write.csv(resall, "Derived/SDM-mod_output.csv")


