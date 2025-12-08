library(terra)
library(sf)
library(landscapemetrics)

### Landscape evaluation for each species 

### Read the species habitat maps 
Pred_f51 <- rast("Data/Derived/Pred_f51-modal-20251022.tif")
Pred_f00 <- rast("Data/Derived/Pred_f00-modal-20251022.tif")

### Mask out the threshold maps
values(Pred_f51)[values(Pred_f51)==0] <- NA
values(Pred_f00)[values(Pred_f00)==0] <- NA

### Check if the maps are valid
check_landscape(Pred_f51[[1]])
check_landscape(Pred_f00[[1]])

### Compute fragmentation indices
lm51 <- calculate_lsm(Pred_f51, what=c("lsm_c_ai", "lsm_c_enn_mn"), progress=T)
lm00 <- calculate_lsm(Pred_f00, what=c("lsm_c_ai", "lsm_c_enn_mn"), progress=T)

### Merge all landscape metrics across years
lm51$year <- 1951
lm00$year <- 2000

### Add species names
lm51$sp <- rep(names(Pred_f51), each=length(unique(lm51$metric)))
lm00$sp <- rep(names(Pred_f00), each=length(unique(lm00$metric)))

### Merge
lm <- rbind(lm51, lm77, lm91, lm00)
lm <- lm[,c("sp", "year", "metric", "value")]

## Save full landscape metrics aggregation data
write.csv(lm, file="Data/Derived/Landscape-agg-metrics-20251022.csv", row.names=F)


# ######################################################
# ### ISLAND-LEVEL FRAGMENTATION STATS
# ######################################################
# 
# ### Read the reprojected forest change rasters
# f <- rast("Data/Derived/Forest_change_resampled_modal-20251022.tif")
# 
# ### Compute fragmentation indices
# lm51 <- calculate_lsm(Pred_f51, what=c("lsm_c_ai", "lsm_c_enn_mn"), progress=T)
# lm00 <- calculate_lsm(Pred_f00, what=c("lsm_c_ai", "lsm_c_enn_mn"), progress=T)
# 
# ### Merge all landscape metrics across years
# lm51$year <- 1951
# lm00$year <- 2000
# 
# ### Add species names
# lm51$sp <- rep(names(Pred_f51), each=length(unique(lm51$metric)))
# lm00$sp <- rep(names(Pred_f00), each=length(unique(lm00$metric)))
# 
# ### Merge
# lm <- rbind(lm51, lm77, lm91, lm00)
# lm <- lm[,c("sp", "year", "metric", "value")]
# 






