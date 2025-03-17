library(raster)
library(sf)
library(sp)
library(landscapemetrics)
library(dplyr)

### Landscape evaluation for each species 

## Read the species habitat maps 
Pred_f51 <- stack("Data/Derived/Pred_f51.tif")
Pred_f77 <- stack("Data/Derived/Pred_f77.tif")
Pred_f91 <- stack("Data/Derived/Pred_f91.tif")
Pred_f00 <- stack("Data/Derived/Pred_f00.tif")

## Check if the maps are valid
check_landscape(Pred_f51)
check_landscape(Pred_f77)
check_landscape(Pred_f91)
check_landscape(Pred_f00)

## Compute all aggregation indices
lm51 <- calculate_lsm(Pred_f51, level="landscape", type="aggregation metric", progress=T)
lm77 <- calculate_lsm(Pred_f77, level="landscape", type="aggregation metric", progress=T)
lm91 <- calculate_lsm(Pred_f91, level="landscape", type="aggregation metric", progress=T)
lm00 <- calculate_lsm(Pred_f00, level="landscape", type="aggregation metric", progress=T)

## Merge all landscape metrics across years
lm51$year <- 1951
lm77$year <- 1977
lm91$year <- 1991
lm00$year <- 2000

## Add species names
lm51$sp <- rep(names(Pred_f51), each=length(unique(lm51$metric)))
lm77$sp <- rep(names(Pred_f51), each=length(unique(lm51$metric)))
lm91$sp <- rep(names(Pred_f51), each=length(unique(lm51$metric)))
lm00$sp <- rep(names(Pred_f51), each=length(unique(lm51$metric)))

## Merge
lm <- rbind(lm51, lm77, lm91, lm00)
lm <- lm[lm$metric != 'iji',]
lm <- lm[,c("sp", "year", "metric", "value")]

## Example way to subset
# lm[lm$year==1951 & lm$metric=='ai',]

## Save full landscape metrics aggregation data
write.csv(lm, file="Data/Derived/Landscape-agg-metrics-20240828.csv", row.names=F)




