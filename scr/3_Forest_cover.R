library(raster)
library(rgdal)
library(sp)


# Load Helmer age map
age <- raster("Data/Maps_1951-2000/iitf_jgr113_puertorico_forestage_zone_reprojectedWGS84.tif")

#load pr outline shapefiel
pr<- readOGR("Data/PR_shapes/outline/PR_outline_Project.shp")
pr <- spTransform(pr, CRS("+proj=lcc +lat_0=17.8333333333333 +lon_0=-66.4333333333333 +lat_1=18.0333333333333 +lat_2=18.4333333333333 +x_0=152400.3048 +y_0=0 +datum=NAD27 +units=m +no_defs"))

#load forest maps
r51 <- raster("Data/Maps_1951-2000/puerto51_sub1_100905_landcov_final.img")
r77 <- raster("Data/Maps_1951-2000/puerto77_sub1_100905_landcov_urbveg_final.img")
r91 <- raster("Data/Maps_1951-2000/pr91_100805_final_quarry_recode2landcov_subset.img")
r00 <- raster("Data/Maps_1951-2000/pr2000_100805_final_quarry_recode2landcov_subset.img")

#reproject some of r91 and r77 (with different extent and cordinates)
rp77 <- projectRaster(r77, r51, method='ngb')
writeRaster(rp77, "Data/Maps_1951-2000/Rp_77.img", format= 'raster', overwrite=TRUE)

rp91<-projectRaster(r91, r51, method='ngb')
writeRaster(rp91, "Data/Maps_1951-2000/Rp_91.img", format= 'raster')

rp00 <-projectRaster(r00, r51, method='ngb')
writeRaster(rp00, "Data/Maps_1951-2000/Rp_00.img")

### RECLASSIFY FOREST AREAS (pixels that are covered in forest)

f51 <- r51 %in% 5
f77 <- r77 %in% c(5,7)
f91 <- r91 %in% c(5,7)
f00 <- r00 %in% c(5,7)

#stack forest cover maps for the 4 different time point 
f <- mask(stack(f51, f77, f91, f00), pr)

#Plot d all together in one silngle quagrant
f_plot<- crop(mask(stack(f51, f77, f91, f00), pr), buffer(pr))
plot(f_plot, axes=FALSE,box=FALSE,legend=FALSE)
plot(pr, add=T)

#Plot each year separately 
f51<- crop(mask(f51, pr), buffer(pr))
plot(f51, axes=FALSE,box=FALSE,legend=FALSE)
plot(pr, add=T)

f77<- crop(mask(f77, pr), buffer(pr))
plot(f77, axes=FALSE,box=FALSE,legend=FALSE)
plot(pr, add=T)

f91<- crop(mask(f91, pr), buffer(pr))
plot(f91, axes=FALSE,box=FALSE,legend=FALSE)
plot(pr, add=T)

f00<- crop(mask(f00, pr), buffer(pr))
plot(f00, axes=FALSE,box=FALSE,legend=FALSE)
plot(pr, add=T)


### FOREST COVER IN THE 6 LIFE ZONES 
#rea in the life zones layer (poligon)
lz <- readOGR("/Users/lauramoro/Documents/PUERTO_RICO/Land_Cover_GAP/Data/Lifezones/lifezones_Project.shp")

#transform the lifezone layer into a raster
lzr <- rasterize(lz['ECOZONE'], age)
writeRaster(lzr, "/Users/lauramoro/Documents/PUERTO_RICO/Forest_Age_Helmer/Lzr.tiff/", format="GTiff", overwrite=TRUE) #save the raster 

#classify the 4 lifezone into numerical IDs
lzr@data@attributes[[1]]$lz <- as.numeric(as.factor(lzr@data@attributes[[1]]$ECOZONE))

#reclassify the the life zones only with id numbers 
lzr2 <- reclassify(lzr, lzr@data@attributes[[1]][,c("ID","lz")])

# Get results (forest cover for each lifezone) using the polygons to mask different areas
result_raster <- matrix(nrow=length(unique(lzr@data@attributes[[1]]$lz)), ncol=4)
rownames(result_raster) <- sort(unique(lzr@data@attributes[[1]]$ECOZONE))
colnames(result_raster) <- c(1951, 1977, 1991, 2000)

##lifezones and forest cover 
#1-total number of pixel of forest at diffrent time points 
for(i in 1:length(unique(lzr@data@attributes[[1]]$lz))){
  print(i)
  tmp <- f *(lzr2 == i)
  result_raster[i,] <- cellStats(tmp, sum)
}
#
plot(c(1951, 1977, 1991, 2000), result_raster[1,], 
     type='l', ylim=c(0, max(result_raster)), lwd=3)

#2-propotion of pixel covered in forest (total percentage of each given lifezone coverd in forest)
for(i in 1:length(unique(lzr@data@attributes[[1]]$lz))){
  print(i)
  tmp <- f *(lzr2 == i)
  result_raster[i,] <- cellStats (tmp, sum) / cellStats((lzr2 == i), sum)*100
}

# 2->  PLOT that shows % forest cover chnege in the different life zones

# choose a color palette
my.palette <- (brewer.pal(n=6, name = 'Set2'))

#plot the the % of forest cover change for one lifezone
plot(c(1951, 1977, 1991, 2000), result_raster[1,], 
     xlab= "Years", ylab="% forest cover",
     type='l', ylim=c(0, max(result_raster)), lwd=3,
     col=my.palette[6])
#complete the plots for the other lifezones 
for(i in 2:6){
  lines(c(1951, 1977, 1991, 2000), result_raster[i,], 
        col=rev(my.palette)[i], lwd=3)
}
#leggend 
legend("left", legend=rownames(result_raster),
       col=rev(my.palette), lty=1, bty='n', cex=.75, lwd=3)

#creat a map of life zones in the sam ecolor palette 
my.palette <- (brewer.pal(n=6, name = 'Set2'))
plot(lzr2, col=rev(my.palette) ,axes=FALSE, box=FALSE, legend=NULL)


legend (legend=rownames(result_raster))
plot(pr, add=T)

plot(T, axes=FALSE, box=FALSE, legend=F)
plot(pr, add=TRUE)


#box plot of forest gain in the different
# load full occurrence data
head(full)
lz_sp <- vector()
for(i in 1:length(splayers)){
  print(i)
  focsp <- substring(splayers[i], 1, 6)
  focsp_xy <- full[full$CODE %in% focsp, c("LONGDEC","LATDEC")]
  lz_sp[i] <- as.numeric(names(sort(table(extract(lzr2, focsp_xy)),
                                    decreasing = T))[1])
}

# generate random data and fill output just for testing!
outmat_abs <- matrix(nrow=length(splayers), ncol=4, 
                     data=rnorm(4*length(splayers)))

# compute change in habitat available from time 4 to time 1
#absolute values 
delta_habitat <- outmat_abs[,4] - outmat_abs[,1]
a<-(delta_habitat+4000) #add a constant for the log trasformation 
#box plot absolute values 
boxplot((a~lz_sp), log="y", ylab= "Habitat Change",
        xlab="Lifezones", names= c('df-S','mf-S','rf-S','wf-LM','wf-S'),
        col=rev(my.palette[-4]))
abline(h=4000, lty=2)

#relative values 
delta_habitat <- outmat_rel[,4] - outmat_rel[,1]
#boxplot
boxplot((delta_habitat~lz_sp), ylab= "Habitat Change",
        xlab="Lifezones", names= c('df-S','mf-S','rf-S','wf-LM','wf-S'),
        col=rev(my.palette[-4]))
abline(h=0, lty=2)

boxplot((B~lz_sp),col=rev(my.palette[-4]))
abline(h=0, lty=2)



### READ IN SDM LAYERS AND PROCESS--> habitata gains and losses for each species 

sdm_path <- "/Users/laumo791/Documents/PR/C1/SDMs/results/Models/Predictions"
splayers <- list.files(sdm_path)[grepl(".gri", splayers)] #all the species 

raster("/Users/lauramoro/Documents/PUERTO_RICO/Tree_SDMs/Bob_Tree_SDMs/ZANMON.grd")

outmat_abs <- outmat_rel <- matrix(nrow=length(splayers), ncol=4)

thresh <- 0.7

#trial Run 
for(i in sample(1:length(splayers), 3)){   #length(splayers)){
  print(i)
  sp <- crop(raster(paste0(sdm_path, splayers[[i]])), f)
  spp <- raster::resample(sp, f[[1]])
  fsp <- f * (spp > thresh)
  outmat_abs[i,] <- cellStats(fsp, sum)
  outmat_rel[i,] <- cellStats(fsp, sum)/cellStats(spp > thresh, sum)
}

#for the all species layers (267)
for(i in sample(1:length(splayers))){ #length(splayers)){
  print(i)
  sp <- crop(raster(paste0(sdm_path, splayers[[i]])), f)
  spp <- raster::resample(sp, f[[1]])
  fsp <- f * (spp > 0.7)
  outmat_abs[i,] <- cellStats(fsp, sum)
  outmat_rel[i,] <- cellStats(fsp, sum)/cellStats(spp > 0.7, sum)
}

matplot(t(outmat_abs), type='l', lty=1, col = my.palette, )
matplot(t(outmat_rel), type='l', lty=1, log = "y")

plot(t(outmat_abs)[4,]-t(outmat_abs)[1,], 
     t(outmat_rel)[4,]-t(outmat_rel)[1,], bg=rainbow(15),
     pch=21, log='xy')

plot(a[41:267], bg=rainbow(15),type='l',lty=1)

#plot histrogram of species ABSOLUTE habitat change
Abs<-(t(outmat_abs)[4,]-t(outmat_abs)[1,])
#habitat gain
a<-sort(Abs)
hist(log(a[41:267]), col="yello")
# habitata loss
n<-a[1:42]
loss<-(n+4000)
hist(log(loss))

hist((t(outmat_rel)[4,]-t(outmat_rel)[1,]), xlab="Habitat change 1951-2000")

#make an hystorgram with the soecies and the gain/loss of habitat 


### USING SHAPEFILES
# Get results using the polygons to mask different areas

lz <- readOGR("/Users/lauramoro/Documents/PUERTO_RICO/Land_Cover_GAP/Data/Lifezones/lifezones_Project.shp")

result <- matrix(nrow=length(unique(lz$ECOZONE)), ncol=4)
rownames(result) <- unique(lz$ECOZONE)
colnames(result) <- c(1951, 1977, 1991, 2000)

for(i in 1:length(unique(lz$ECOZONE))){
  focal_lz <- unique(lz$ECOZONE)[i]
  tmp <- mask(f, lz[lz$ECOZONE == focal_lz,])
  result[i,] <- cellStats(tmp, sum)
}

plot(c(1951, 1977, 1991, 2000), result[1,], 
     type='l', ylim=c(0,max(result)), lwd=3)

for(i in 2:6){
  lines(c(1951, 1977, 1991, 2000), result[i,], col=i, lwd=3)
}

legend("topleft", legend=rownames(result), col=1:6, lty=1, bty='n', cex=.75, lwd=3)




