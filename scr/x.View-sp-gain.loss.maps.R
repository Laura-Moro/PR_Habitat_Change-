
library(terra)

### Full resolution forest cover maps
f <- rast("Data/Derived/forest_cover_stack_native.tif")
f <- rast("Data/Derived/Forest_change_resampled.tif")


### Full resolution SDM predictions
# Import the continuous maps 
Pred_stack_raw <- rast(list.files(path = "Data/SDM_predictions", 
                                   pattern='.tif', all.files=TRUE, full.names=T))

# Import the treshholded maps 
Pred_stack <- rast(list.files(path = "Data/SDM_threshold", 
                               pattern='.tif', all.files=TRUE, full.names=T))

# Assign the names of the species
names(Pred_stack_raw) <- gsub(".tif", "", 
                              list.files(path = "Data/SDM_predictions",
                                         pattern='.tif', all.files=TRUE, full.names=F))
names(Pred_stack) <- gsub(".tif", "", 
                          list.files(path = "Data/SDM_threshold",
                                     pattern='.tif', all.files=TRUE, full.names=F))

for(sp in 1:dim(Pred_stack_raw)[[3]]){
  
  message(paste("working on species", sp, "out of", dim(Pred_stack_raw)[[3]]))
  # Reproject the threshold maps 
  Pred_stack_raw_rp <- terra::project(Pred_stack_raw[[sp]], f[[1]], method='near')
  Pred_stack_rp <- terra::project(Pred_stack[[sp]], f[[1]], method='near')
  
  # tmp_pred <- terra::resample(rast(Pred_stack_raw[[sp]]), rast(f[[1]]), method='ngb')
  
  tmp_pred_stack_raw_forest <- Pred_stack_raw_rp * (f[[4]]-f[[1]])
  tmp_pred_stack_forest <- Pred_stack_rp * (f[[4]]-f[[1]])
  
  sum(values(test_pred_stack_forest)[values(test_pred_stack_forest)<0], na.rm=T)
  sum(values(test_pred_stack_forest)[values(test_pred_stack_forest)>0], na.rm=T)
  
  plot(test_pred_stack_forest<0)
  plot(test_pred_stack_forest>0)
  plot(pr$geometry, add=T)
}






# ===========================

library(terra)

### Full resolution forest cover maps
f <- rast("Data/Derived/Forest_change_resampled.tif")
newforest <- (f_rs[[4]] - f_rs[[1]])

plot(newforest, main="Forest change 1951-2000 (loss, remain, gain)")

### Full resolution SDM predictions
# Import the continuous maps 
Pred_stack_raw <- terra::rast("Data/Derived/raw_stack.tif")

# Import the continuous maps 
Pred_stack_thres <- terra::rast("Data/Derived/thresholded_stack.tif")

# Get species names on the stacks
names(Pred_stack_thres) <-  names(Pred_stack_raw) <- gsub(".tif", "", 
                                                          list.files(path = "Data/SDM_predictions",
                                                                     pattern='.tif', 
                                                                     all.files=TRUE, full.names=F))

rwb <- colorRampPalette(colors = c("red", "white", "blue"))
wb <- colorRampPalette(colors = c("white", "blue"))
rw <- colorRampPalette(colors = c("white", "red"))


# LOAD ENVIRONMENTAL PREDICTOR VARIABLES derived in script 1 
envs <- stack(list.files("Data/Derived/envs", full.names = TRUE))

# Set water geological type to NA
values(envs[[1]])[values(envs[[1]]) %in% 0] <- NA

# Scale envs
envs_scaled <- rast(scale(envs))

# Get scaled envs into same units as the forest map
envs_scaled_rp <- terra::project(envs_scaled, newforest, method='near')

# Make mask of the new forest gains
newforest_mask <- newforest
values(newforest_mask)[values(newforest_mask) <= 0] <- NA

# Look at it...
# par(mfrow=c(2,2), mar=c(4,4,1,1))
# for(i in 1:4){
#   hist(envs_scaled_rp[[i]])
#   hist(mask(envs_scaled_rp[[i]], newforest_mask), add=T, col='red')
# }  

# plot(sum(Pred_stack_thres))

richness_in_new_forests <- sum(Pred_stack_thres * newforest > 0)
values(richness_in_new_forests)[values(richness_in_new_forests)==0] <- NA

plot(pr$geometry, col=rgb(0,0,0,0.05), axes=T)
plot(richness_in_new_forests, col=viridis::inferno(100), add=T)

plot(sum(Pred_stack_thres * newforest < 0))
plot(pr$geometry, add=T)


res_thresh <- res_raw <- matrix(nrow=dim(Pred_stack_raw)[3], ncol=2)

# pdf("/Users/au529793/Desktop/forest-gain-loss-maps.pdf")

for(sp in 1:dim(Pred_stack_raw)[3]){
  
  message(paste("working on species", sp, "out of", dim(Pred_stack_raw)[[3]]))
  
  tmp_pred_raw_forest <- Pred_stack_raw[[sp]] * newforest
  tmp_pred_thresh_forest <- Pred_stack_thres[[sp]] * newforest
  
  values(tmp_pred_raw_forest)[values(tmp_pred_raw_forest) == 0] <- NA
  values(tmp_pred_thresh_forest)[values(tmp_pred_thresh_forest) == 0] <- NA
  
  # plot(tmp_pred_raw_forest, col=rwb(100), main=names(tmp_pred_raw_forest))
  # plot(pr$geometry, add=T)

  res_raw[sp,1] <- sum(values(tmp_pred_raw_forest)[values(tmp_pred_raw_forest) < 0], na.rm=T)
  res_raw[sp,2] <- sum(values(tmp_pred_raw_forest)[values(tmp_pred_raw_forest) > 0], na.rm=T)

  res_thresh[sp,1] <- sum(values(tmp_pred_thresh_forest)[values(tmp_pred_thresh_forest) < 0], na.rm=T)
  res_thresh[sp,2] <- sum(values(tmp_pred_thresh_forest)[values(tmp_pred_thresh_forest) > 0], na.rm=T)

  
  plot(mask(envs_scaled_rp, tmp_pred_raw_forest))
  
  plot(envs_scaled_rp[[1]] * tmp_pred_thresh_forest)

  
}

# dev.off()

# What propo
hist(abs(res_thresh[,1] / rowSums(abs(res_thresh))))



hist(res[,2], breaks=10, col='blue')
hist(abs(res[,1]), add=T, col='red', breaks=10)

plot(res[,1], res[,2])
