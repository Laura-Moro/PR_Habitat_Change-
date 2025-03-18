library(raster)
library(sp)
library(sf)

#1 CALCULATE TOTAL FOREST COVER 
# Load Helmer age map
age <- raster("Data/Maps_1951-2000/iitf_jgr113_puertorico_forestage_zone_reprojectedWGS84.tif")

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

# plot(sum(f), col=rev(viridis::viridis(5)))

# Plot all together
f_plot <- crop(mask(stack(f51, f77, f91, f00), pr), pr)
par(mfrow=c(2,2))
for(i in 1:4){
  plot(f_plot[[i]], legend=FALSE)
  plot(pr$geometry, add=T)
}

### READ IN SDM LAYERS AND PROCESS--> habitats gains and losses for each species 
#1 Re-project binary maps  

# Import the continuous maps 
Pred_stack_raw <- stack(list.files(path = "Data/SDM_predictions-noFIA", 
                               pattern='.tif', all.files=TRUE, full.names=T))

# Import the treshholded maps 
Pred_stack <- stack(list.files(path = "Data/SDM_threshold-noFIA", 
                   pattern='.tif', all.files=TRUE, full.names=T))

# Assign the names of the species
names(Pred_stack_raw) <- gsub(".tif", "", 
                              list.files(path = "Data/SDM_predictions-noFIA",
                                         pattern='.tif', all.files=TRUE, full.names=F))
names(Pred_stack) <- gsub(".tif", "", 
                          list.files(path = "Data/SDM_threshold-noFIA",
                                     pattern='.tif', all.files=TRUE, full.names=F))

# Projection system of the Helmer maps 
newproj <- "+proj=lcc +lat_0=17.8333333333333 +lon_0=-66.4333333333333 
          +lat_1=18.0333333333333 +lat_2=18.4333333333333 +x_0=152400.3048
          +y_0=0 +datum=NAD27 +units=m +no_defs"

# Reproject the threshold maps 
Pred_stack_raw_rp <- projectRaster(Pred_stack_raw, crs=crs(pr), method='ngb')
Pred_stack_rp <- projectRaster(Pred_stack, crs=crs(pr), method='ngb')

# Mask and crop to smaller extent (no offshore islands)
Pred_stack_raw_rp <- crop(mask(Pred_stack_raw_rp, pr), pr)
Pred_stack_rp <- crop(mask(Pred_stack_rp, pr), pr)

# Save the reprojected raster to use in the Landscape section 
writeRaster(Pred_stack_raw_rp, "Data/Derived/raw_stack-noFIA.tif", format="GTiff")
writeRaster(Pred_stack_rp, "Data/Derived/thresholded_stack-noFIA.tif", format="GTiff")


# Resample the forest maps 
f_rs <- raster::resample(f_plot, Pred_stack_rp , method="ngb")

# Mask the SDMs with the forest cover maps at the different time points 
Pred_f51_raw <- f_rs[[1]] * Pred_stack_raw_rp 
Pred_f77_raw <- f_rs[[2]] * Pred_stack_raw_rp
Pred_f91_raw <- f_rs[[3]] * Pred_stack_raw_rp
Pred_f00_raw <- f_rs[[4]] * Pred_stack_raw_rp

Pred_f51 <- f_rs[[1]] * Pred_stack_rp 
Pred_f77 <- f_rs[[2]] * Pred_stack_rp
Pred_f91 <- f_rs[[3]] * Pred_stack_rp
Pred_f00 <- f_rs[[4]] * Pred_stack_rp


# Assign names to the predictions Layers
names(Pred_f51_raw) <- names(Pred_stack_raw_rp)
names(Pred_f77_raw) <- names(Pred_stack_raw_rp)
names(Pred_f91_raw) <- names(Pred_stack_raw_rp)
names(Pred_f00_raw) <- names(Pred_stack_raw_rp)

names(Pred_f51) <- names(Pred_stack_rp)
names(Pred_f77) <- names(Pred_stack_rp)
names(Pred_f91) <- names(Pred_stack_rp)
names(Pred_f00) <- names(Pred_stack_rp)

# Save rasters of predicted habitat maps
writeRaster(Pred_f51_raw, "Data/Derived/Pred_f51_raw-noFIA.tif", format="GTiff")
writeRaster(Pred_f77_raw, "Data/Derived/Pred_f77_raw-noFIA.tif", format="GTiff")
writeRaster(Pred_f91_raw, "Data/Derived/Pred_f91_raw-noFIA.tif", format="GTiff")
writeRaster(Pred_f00_raw, "Data/Derived/Pred_f00_raw-noFIA.tif", format="GTiff")

writeRaster(Pred_f51, "Data/Derived/Pred_f51-noFIA.tif", format="GTiff")
writeRaster(Pred_f77, "Data/Derived/Pred_f77-noFIA.tif", format="GTiff")
writeRaster(Pred_f91, "Data/Derived/Pred_f91-noFIA.tif", format="GTiff")
writeRaster(Pred_f00, "Data/Derived/Pred_f00-noFIA.tif", format="GTiff")

# Sum all of the pixels that were forest in 1951, 1977, 1991, 2000
fcover_51_raw <- cellStats(Pred_f51_raw, 'sum')
fcover_77_raw <- cellStats(Pred_f77_raw, 'sum')
fcover_91_raw <- cellStats(Pred_f91_raw, 'sum')
fcover_00_raw <- cellStats(Pred_f00_raw, 'sum')
habitat_raw <- cellStats(Pred_stack_raw_rp, 'sum')

fcover_51 <- cellStats(Pred_f51, 'sum')
fcover_77 <- cellStats(Pred_f77, 'sum')
fcover_91 <- cellStats(Pred_f91, 'sum')
fcover_00 <- cellStats(Pred_f00, 'sum')
habitat <- cellStats(Pred_stack_rp, 'sum')


# Transform into the data frame 
F_cover <- data.frame(sp=names(Pred_stack_rp),
                      
                      total_hab_raw=habitat_raw,
                      fcover_51_raw=fcover_51_raw, 
                      fcover_77_raw=fcover_77_raw,
                      fcover_91_raw=fcover_91_raw,
                      fcover_00_raw=fcover_00_raw,
                      tot_change_raw=fcover_00_raw - fcover_51_raw,
                      
                      total_hab=habitat,
                      fcover_51=fcover_51, 
                      fcover_77=fcover_77,
                      fcover_91=fcover_91,
                      fcover_00=fcover_00,
                      tot_change=fcover_00 - fcover_51)

rownames(F_cover) <- NULL

write.csv(F_cover, "Data/Derived/3b-output-20250318.csv", row.names = F)



#------------------------------------------------------------------------------#

# 3  FOREST COVER IN THE 6 LIFE ZONES 
# Area in the life zones layer (poligon)
lz <- readOGR("/Users/lauramoro/Documents/PUERTO_RICO/Land_Cover_GAP/Data/Lifezones/lifezones_Project.shp")

# Transform the lifezone layer into a raster
lzr <- rasterize(lz['ECOZONE'], age)
writeRaster(lzr, "/Users/lauramoro/Documents/PUERTO_RICO/Forest_Age_Helmer/Lzr.tiff/", format="GTiff", overwrite=TRUE) #save the raster 

# Classify the 4 lifezone into numerical IDs
lzr@data@attributes[[1]]$lz <- as.numeric(as.factor(lzr@data@attributes[[1]]$ECOZONE))

# Reclassify the the life zones only with id numbers 
lzr2 <- reclassify(lzr, lzr@data@attributes[[1]][,c("ID","lz")])

# Get results (forest cover for each lifezone) using the polygons to mask different areas
result_raster <- matrix(nrow=length(unique(lzr@data@attributes[[1]]$lz)), ncol=4)
rownames(result_raster) <- sort(unique(lzr@data@attributes[[1]]$ECOZONE))
colnames(result_raster) <- c(1951, 1977, 1991, 2000)

## Lifezones and forest cover 
#1 - Total number of pixel of forest at different time points 
for(i in 1:length(unique(lzr@data@attributes[[1]]$lz))){
  print(i)
  tmp <- f *(lzr2 == i)
  result_raster[i,] <- cellStats(tmp, sum)
}
#
plot(c(1951, 1977, 1991, 2000), result_raster[1,], 
     type='l', ylim=c(0, max(result_raster)), lwd=3)

#2 - Proportion of pixel covered in forest (total percentage of each given lifezone coverd in forest)
for(i in 1:length(unique(lzr@data@attributes[[1]]$lz))){
  print(i)
  tmp <- f *(lzr2 == i)
  result_raster[i,] <- cellStats (tmp, sum) / cellStats((lzr2 == i), sum)*100
}


# 2->  PLOT that shows % forest cover change in the different life zones

# choose a color palette
my.palette <- (brewer.pal(n=6, name = 'Set2'))

# plot the the % of forest cover change for one lifezone
plot(c(1951, 1977, 1991, 2000), result_raster[1,], 
     xlab= "Years", ylab="% forest cover",
     type='l', ylim=c(0, max(result_raster)), lwd=3,
     col=my.palette[6])
# complete the plots for the other lifezones 
for(i in 2:6){
  lines(c(1951, 1977, 1991, 2000), result_raster[i,], 
        col=rev(my.palette)[i], lwd=3)
}
# legend 
legend("left", legend=rownames(result_raster),
       col=rev(my.palette), lty=1, bty='n', cex=.75, lwd=3)

# create a map of life zones in the same color palette 
my.palette <- (brewer.pal(n=6, name = 'Set2'))
plot(lzr2, col=rev(my.palette) ,axes=FALSE, box=FALSE, legend=NULL)

legend (legend=rownames(result_raster))
plot(pr, add=T)

plot(T, axes=FALSE, box=FALSE, legend=F)
plot(pr, add=TRUE)


# box plot of forest gain in the different
# load full occurrence data
head(full)
lz_sp <- vector()
for(i in 1:length(Pred)){
  print(i)
  focsp <- substring(Pred[i], 1, 6)
  focsp_xy <- full[full$CODE %in% focsp, c("LONGDEC","LATDEC")]
  lz_sp[i] <- as.numeric(names(sort(table(extract(lzr2, focsp_xy)),
                                    decreasing = T))[1])
}

# generate random data and fill output just for testing!
outmat_abs <- matrix(nrow=length(splayers), ncol=4, 
                     data=rnorm(4*length(splayers)))

# compute change in habitat available from time 4 to time 1
# absolute values 
delta_habitat <- outmat_abs[,4] - outmat_abs[,1]
a<-(delta_habitat+4000) #add a constant for the log trasformation 
# box plot absolute values 
boxplot((a~lz_sp), log="y", ylab= "Habitat Change",
        xlab="Lifezones", names= c('df-S','mf-S','rf-S','wf-LM','wf-S'),
        col=rev(my.palette[-4]))
abline(h=4000, lty=2)

# relative values 
delta_habitat <- outmat_rel[,4] - outmat_rel[,1]
# boxplot
boxplot((delta_habitat~lz_sp), ylab= "Habitat Change",
        xlab="Lifezones", names= c('df-S','mf-S','rf-S','wf-LM','wf-S'),
        col=rev(my.palette[-4]))
abline(h=0, lty=2)

boxplot((B~lz_sp),col=rev(my.palette[-4]))
abline(h=0, lty=2)






