library(raster)
library(rgdal)
library(sp)

# 1 CALCULATE TOTAL FOREST COVER 
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

#reproject some of r91 and r77 and r00 (with different extent and cordinates) and save them 

#rp77 <- projectRaster(r77, r51, method='ngb')
#writeRaster(rp77, "Data/Maps_1951-2000/Rpm_77.img")
rp77 <- raster("Data/Maps_1951-2000/Rpm_77.img")

#rp91<-projectRaster(r91, r51, method='ngb')
#writeRaster(rp91, "Data/Maps_1951-2000/Rpm_91.img")
rp91<- raster("Data/Maps_1951-2000/Rpm_91.img")

#rp00 <-projectRaster(r00, r51, method='ngb')
#writeRaster(rp00, "Data/Maps_1951-2000/Rpm_00.img")
rp00<- raster("Data/Maps_1951-2000/Rpm_00.img")

### RECLASSIFY FOREST AREAS (pixels that are covered in forest)
f51 <- r51 %in% 5
f77 <- rp77 %in% c(5,7)
f91 <- rp91 %in% c(5,7)
f00 <- rp00 %in% c(5,7)

#stack and mask the forest classes
f <- mask(stack(f51, f77, f91, f00), pr)

#Plot d all together in one silngle quagrant
f_plot<- crop(mask(stack(f51, f77, f91, f00), pr), buffer(pr))
plot(f_plot, axes=FALSE,box=FALSE,legend=FALSE)
plot(pr, add=T)

### READ IN SDM LAYERS AND PROCESS--> habitats gains and losses for each species 
#1 re-project binary maps  

#import the treshholded maps 
Pred <- list.files(path = "/Users/laumo791/Documents/PR/C1/Results/SDM_threshold", 
                   pattern='.tif', all.files=TRUE, full.names=T)

#make a stack of all the tresholded predicitons
Pred_stack <- stack(Pred)

#take a the names of the species 
names_pred <- gsub(".tif", "", list.files(path = "/Users/laumo791/Documents/PR/C1/Results/SDM_threshold", 
                                          pattern='.tif', all.files=TRUE, full.names=F))

#assign names to the layers of the stack 
names(Pred_stack) <- names_pred

#Projection system of the hemlmer maps 
newproj <- "+proj=lcc +lat_0=17.8333333333333 +lon_0=-66.4333333333333 
          +lat_1=18.0333333333333 +lat_2=18.4333333333333 +x_0=152400.3048
          +y_0=0 +datum=NAD27 +units=m +no_defs"

# reproject the model the threshold maps 
Pred_stack_rp <- projectRaster(Pred_stack, crs=newproj, method='ngb')

#save the re projected raster to use in the Landscape section 
writeRaster(Pred_stack_rp, "/Users/laumo791/Documents/PR/C1/Results/F_stack/t_stack.tif", format= 'GTiff')

names(Pred_stack_rp) <- names_pred

#forest map and here we resample the forest maps 
f_rs <- raster::resample( f, Pred_stack_rp , method="ngb")

#mask the SDMs with the forest cover maps at the different time points 
Pred_f51 <- f_rs[[1]]* Pred_stack_rp 
Pred_f77 <- f_rs[[2]]* Pred_stack_rp
Pred_f91 <- f_rs[[3]]* Pred_stack_rp
Pred_f00 <- f_rs[[4]]* Pred_stack_rp

#assign names to the predictions Layers
names(Pred_f51) <- names_pred
names(Pred_f77) <- names_pred
names(Pred_f91) <- names_pred
names(Pred_f00) <- names_pred

#summ all of the pixels that were forest in 1951, 1977, 1991, 2000

fcover_51<-cellStats(Pred_f51, 'sum')
fcover_77<-cellStats(Pred_f77, 'sum')
fcover_91<-cellStats(Pred_f91, 'sum')
fcover_00<-cellStats(Pred_f00, 'sum')

#total change 
tot_change <- fcover_00 - fcover_51

#transform into the data frame 
df_fcover_51 <- as.data.frame(fcover_51)
df_fcover_77 <- as.data.frame(fcover_77)
df_fcover_91 <- as.data.frame(fcover_91)
df_fcover_00 <-as.data.frame(fcover_00)
df_tot_change <-as.data.frame(tot_change)

#make a data frame of forets cover at the different time frames 
F_cover <- cbind(df_fcover_51, df_fcover_77, df_fcover_91, df_fcover_00, df_tot_change)
write.csv(F_cover, "Data/Derived/Forest_Cover.csv")
F_cover <- read.csv("Data/Derived/Forest_Cover.csv") 


#transform the dataframe into a matrix 
F_mat <- as.matrix(F_cover[2:5]*0.2025)
#change row and column names 
rownames(F_mat) <- F_cover$X
colnames(F_mat) <- c(1951, 1977, 1991, 2000)

#plot total forest cover change 
matplot(t(F_mat), type='l', lty=1, xaxt='n', xlab = "Time (years)", ylab = "climaticaly suitable Habitat (km2)", cex.lab=2, cex.axis=1.8)
axis(side=1, at= c(1, 2, 3, 4), labels = c(1951, 1977, 1991, 2000), cex.axis=1.8)


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


#------------------------------------------------------------------------------#

# 3  FOREST COVER IN THE 6 LIFE ZONES 
#Area in the life zones layer (poligon)
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






