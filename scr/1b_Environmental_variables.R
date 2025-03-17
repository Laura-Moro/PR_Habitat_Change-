library(sp)
library(raster)
library(dplyr)
library(sf)
library(terra)

## DERIVING ENVIRONMNETAL VARIBLES 
#1 Precipitation seasonality 
#2 Precipitation in the driest month
#3 Temperature in the hottest month
#4 Geological substrate 
#"/Users/laumo791/Documents/PR/C1/Data/Derived/envs"#

#1 PRECIPITATION SEASONALIY 
ppt <- raster::stack(list.files("Data/ppt", full.names = TRUE)[1:12])/100

#total annual rainfall (R)
R <- raster::raster(list.files("Data/ppt", full.names = TRUE)[13])/100

# Monthly probability distribution for rainfall (pm)
pm <- (ppt)/R

# Get monthly entropy (D)
D <- sum(pm * log2(pm/(1/12)))

# Compute seasonality index by multiplying by standardized R 
# we decide not you standardize the index 
# S_global <- D * R/cellStats(R, max)#

# Save raster
# raster::crs(D) <- "EPSG:4326"
# writeRaster(D,  filename="Data/Derived/envs/PSI.asc", format="ascii")


#2 PRECIPITATION IN THE DRIEST QUARTER 
# Precipitations in the 4 quarters 
PQ1 <- sum(ppt[[1:3]])
PQ2 <- sum(ppt[[4:6]])
PQ3 <- sum(ppt[[7:9]])
PQ4 <- sum(ppt[[10:12]])

# Precipitation of the driest quarter 
Pm <- min(stack(PQ1, PQ2, PQ3, PQ4))

# writeRaster(Pm, filename="Data/Derived/envs/Pm.asc", format="ascii")


#3 MAXIMUM TEMPERATURE IN THE HOTTEST QUARTER 
tmax <- raster::stack(list.files("Data/tmax", full.names = TRUE)[1:12])/100

# Find the warmest quarter 
TM1 <- sum(tmax[[1:3]])/3
TM2 <- sum(tmax[[4:6]])/3
TM3 <- sum(tmax[[7:9]])/3
TM4 <- sum(tmax[[10:12]])/3

# Maximum temperature of the warmest quarter 
TM <- max(stack(TM1, TM2, TM3, TM4))
plot(TM)

# writeRaster(TM, filename="Data/Derived/envs/TM.asc", format="acsii")


#4 GEOLOGICAL SUBSTARETE 

# READ GEOLOGY
geo <- st_read("Data/GEO/geopr_krushensky/geopr_krushensky.shp")

# TRANSFORM GEOLOGY COORDINATE SYSTEM
geo <- st_transform(geo, crs = crs(ppt))

# RASTERIZE BASED ON PPT DATA
geo_r <- terra::rasterize(geo, field="RECLASS1", ppt, fun=min)

# writeRaster(geo_r, filename="Data/Derived/envs/GEO.asc", format="acsii")

env <- stack(D, Pm, TM, geo_r)
names(env) <- c("psi", "pmin", "tmin", "geo")
plot(env)


age <- raster::raster("Data/Maps_1951-2000/iitf_jgr113_puertorico_forestage_zone_reprojectedWGS84.tif")
plot(age)

env <- crop(env, age)
env2 <- resample(env, age)
env3 <- mask(env2, age)
plot(env3)





