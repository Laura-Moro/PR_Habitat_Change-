# library(raster)
library(terra)
library(sf)
library(sp)
library(landscapemetrics)
library(dplyr)

### Landscape evaluation for each species 

## Read the species habitat maps 
# Pred_f51 <- rast("Data/Derived/Pred_f51-noFIA.tif")
# Pred_f77 <- rast("Data/Derived/Pred_f77-noFIA.tif")
# Pred_f91 <- rast("Data/Derived/Pred_f91-noFIA.tif")
# Pred_f00 <- rast("Data/Derived/Pred_f00-noFIA.tif")

Pred_f51 <- rast("Data/Derived/Pred_f51-modal.tif")
Pred_f77 <- rast("Data/Derived/Pred_f77-modal.tif")
Pred_f91 <- rast("Data/Derived/Pred_f91-modal.tif")
Pred_f00 <- rast("Data/Derived/Pred_f00-modal.tif")

### Mask out the threshold maps
values(Pred_f51)[values(Pred_f51)==0] <- NA
values(Pred_f77)[values(Pred_f77)==0] <- NA
values(Pred_f91)[values(Pred_f91)==0] <- NA
values(Pred_f00)[values(Pred_f00)==0] <- NA

## Check if the maps are valid
check_landscape(Pred_f51[[1]])
check_landscape(Pred_f77[[1]])
check_landscape(Pred_f91[[1]])
check_landscape(Pred_f00[[1]])

## Compute fragmentation indices
lm51 <- calculate_lsm(Pred_f51, what=c("lsm_c_ai", "lsm_c_enn_mn"), progress=T)
lm77 <- calculate_lsm(Pred_f77, what=c("lsm_c_ai", "lsm_c_enn_mn"), progress=T)
lm91 <- calculate_lsm(Pred_f91, what=c("lsm_c_ai", "lsm_c_enn_mn"), progress=T)
lm00 <- calculate_lsm(Pred_f00, what=c("lsm_c_ai", "lsm_c_enn_mn"), progress=T)

## Merge all landscape metrics across years
lm51$year <- 1951
lm77$year <- 1977
lm91$year <- 1991
lm00$year <- 2000

## Get species names
# modfiles <- list.files("Data/SDM_predictions", full.names = TRUE)
# splist <- unlist(lapply(strsplit(modfiles,"/"), function(x) gsub(".tif", "", x[[3]])))
# names(Pred_f51) <- splist

## Add species names
lm51$sp <- rep(names(Pred_f51), each=length(unique(lm51$metric)))
lm77$sp <- rep(names(Pred_f77), each=length(unique(lm51$metric)))
lm91$sp <- rep(names(Pred_f91), each=length(unique(lm51$metric)))
lm00$sp <- rep(names(Pred_f00), each=length(unique(lm51$metric)))

## Merge
lm <- rbind(lm51, lm77, lm91, lm00)
lm <- lm[,c("sp", "year", "metric", "value")]

## Save full landscape metrics aggregation data
write.csv(lm, file="Data/Derived/Landscape-agg-metrics-20251022.csv", row.names=F)




