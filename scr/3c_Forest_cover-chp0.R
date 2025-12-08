library(raster)
library(terra)
library(sp)
library(sf)

#1 CALCULATE TOTAL FOREST COVER 

# Load PR outline shapefile
pr <- st_read("Data/PR_shapes/outline/PR_outline_Project.shp")
pr <- st_transform(pr, CRS("+proj=lcc +lat_0=17.8333333333333 +lon_0=-66.4333333333333 +lat_1=18.0333333333333 +lat_2=18.4333333333333 +x_0=152400.3048 +y_0=0 +datum=NAD27 +units=m +no_defs"))

# Load forest maps
r51 <- rast("Data/Maps_1951-2000/puerto51_sub1_100905_landcov_final.img")
# r00 <- rast("Data/Maps_1951-2000/pr2000_100805_final_quarry_recode2landcov_subset.img")

### Reproject r00 (different extent and crs) and save
# rp00 <-projectRaster(r00, r51, method='ngb')
# writeRaster(rp00, "Data/Maps_1951-2000/Rpm_00.img")
rp00 <- rast("Data/Maps_1951-2000/Rpm_00.img")

### RECLASSIFY FOREST AREAS (pixels that are covered in forest)
f51 <- r51 %in% 5
f00 <- rp00 %in% c(5,7)

# Stack and mask the forest classes
f <- crop(mask(c(f51, f00), pr), pr)

# writeRaster(f, "Data/Derived/forest_cover_stack_native.tif", overwrite=T)

# plot(sum(f), col=rev(viridis::viridis(5)))
# par(mfrow=c(2,2))
# for(i in 1:4){
#   plot(f[[i]], legend=FALSE)
#   plot(pr$geometry, add=T)
# }

### READ IN SDM LAYERS AND PROCESS--> habitats gains and losses for each species 
#1 Re-project binary maps  

### Import the continuous maps 
# Pred_stack_raw <- rast(list.files(path="Data/SDM_predictions-noFIA", 
#                                   pattern='.tif', all.files=TRUE, full.names=T))
Pred_stack_raw <- rast(list.files(path="Data/SDM_predictions", 
                                  pattern='.tif', all.files=TRUE, full.names=T))

### Import the tresh-holded maps 
# Pred_stack <- rast(list.files(path="Data/SDM_threshold-noFIA", 
#                               pattern='.tif', all.files=TRUE, full.names=T))
Pred_stack <- rast(list.files(path="Data/SDM_threshold", 
                              pattern='.tif', all.files=TRUE, full.names=T))

# Assign the names of the species
# names(Pred_stack_raw) <- gsub(".tif", "", 
#                               list.files(path = "Data/SDM_predictions-noFIA",
#                                          pattern='.tif', all.files=TRUE, 
#                                          full.names=F))
names(Pred_stack_raw) <- gsub(".tif", "", 
                              list.files(path = "Data/SDM_predictions",
                                         pattern='.tif', all.files=TRUE, 
                                         full.names=F))

# names(Pred_stack) <- gsub(".tif", "", 
#                           list.files(path = "Data/SDM_threshold-noFIA",
#                                      pattern='.tif', all.files=TRUE, 
#                                      full.names=F))
names(Pred_stack) <- gsub(".tif", "", 
                          list.files(path = "Data/SDM_threshold",
                                     pattern='.tif', all.files=TRUE, 
                                     full.names=F))

### Reproject the threshold maps to projection system of the Helmer maps
Pred_stack_raw_rp <- terra::project(Pred_stack_raw, crs(r51), method='bilinear')
Pred_stack_rp <- terra::project(Pred_stack, crs(r51), method='mode')

### Mask and crop to smaller extent (no offshore islands)
Pred_stack_raw_rp <- crop(mask(Pred_stack_raw_rp, pr), pr)
Pred_stack_rp <- crop(mask(Pred_stack_rp, pr), pr)

### Save the reprojected raster to use in the Landscape section 
writeRaster(Pred_stack_raw_rp, "Data/Derived/raw_stack-20251022.tif", overwrite=T)
writeRaster(Pred_stack_rp, "Data/Derived/thresholded_stack-20251022.tif", overwrite=T)

##################################################
### Aggregate the forest maps to match SDM output
f_rs <- terra::resample(f, Pred_stack_rp , method="modal")

### Save the reprojected forest change rasters
terra::writeRaster(f_rs, "Data/Derived/Forest_change_resampled_modal-20251022.tif", 
                   overwrite=T)

#########################
### Mask the SDMs with the forest cover maps at the different time points 
Pred_f51_raw <- f_rs[[1]] * Pred_stack_raw_rp 
Pred_f00_raw <- f_rs[[2]] * Pred_stack_raw_rp
Pred_f51 <- f_rs[[1]] * Pred_stack_rp 
Pred_f00 <- f_rs[[2]] * Pred_stack_rp

### Assign names to the predictions Layers
names(Pred_f51_raw) <- names(Pred_stack_raw_rp)
names(Pred_f00_raw) <- names(Pred_stack_raw_rp)
names(Pred_f51) <- names(Pred_stack_rp)
names(Pred_f00) <- names(Pred_stack_rp)

### Save rasters of predicted habitat maps
terra::writeRaster(Pred_f51_raw, "Data/Derived/Pred_f51_raw-modal-20251022.tif")
terra::writeRaster(Pred_f00_raw, "Data/Derived/Pred_f00_raw-modal-20251022.tif")
terra::writeRaster(Pred_f51, "Data/Derived/Pred_f51-modal-20251022.tif")
terra::writeRaster(Pred_f00, "Data/Derived/Pred_f00-modal-20251022.tif")

#########################
### Sum all of the pixels that were forest in 1951, 1977, 1991, 2000
fcover_51_raw <- global(Pred_f51_raw, "sum", na.rm=T)
fcover_00_raw <- global(Pred_f00_raw, "sum", na.rm=T)
habitat_raw <- global(Pred_stack_raw_rp, "sum", na.rm=T)

fcover_51 <- global(Pred_f51, "sum", na.rm=T)
fcover_00 <- global(Pred_f00, "sum", na.rm=T)
habitat <- global(Pred_stack_rp, "sum", na.rm=T)

### Transform into the data frame 
F_cover <- data.frame(sp=names(Pred_stack_rp),
                      
                      total_hab_raw=habitat_raw[,1],
                      fcover_51_raw=fcover_51_raw[,1], 
                      fcover_00_raw=fcover_00_raw[,1],
                      tot_change_raw=(fcover_00_raw - fcover_51_raw)[,1],
                      
                      total_hab=habitat[,1],
                      fcover_51=fcover_51[,1], 
                      fcover_00=fcover_00[,1],
                      tot_change=(fcover_00 - fcover_51)[,1])

write.csv(F_cover, "Data/Derived/3b-output-20251022.csv", row.names=F)


