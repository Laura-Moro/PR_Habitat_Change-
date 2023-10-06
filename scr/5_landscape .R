library(raster)
library(rgdal)
library(sp)
library(landscapemetrics)

#lanscape evaluation for each species 
#1 select the forest cove

# Load Helmer age map
age <- raster("Data/Maps_1951-2000/iitf_jgr113_puertorico_forestage_zone_reprojectedWGS84.tif")

#load pr outline shapefiel
pr<- readOGR("Data/PR_shapes/outline/PR_outline_Project.shp")
pr <- spTransform(pr, CRS("+proj=lcc +lat_0=17.8333333333333 +lon_0=-66.4333333333333 +lat_1=18.0333333333333 +lat_2=18.4333333333333 +x_0=152400.3048 +y_0=0 +datum=NAD27 +units=m +no_defs"))

#load forest maps
r51 <- raster("Data/Maps_1951-2000/puerto51_sub1_100905_landcov_final.img")
rp77 <- raster("Data/Maps_1951-2000/Rpm_77.img")
rp91<- raster("Data/Maps_1951-2000/Rpm_91.img")
rp00<- raster("Data/Maps_1951-2000/Rpm_00.img")

#slect only the forest classes 
f51 <- r51 %in% 5
f77 <- rp77 %in% c(5,7)
f91 <- rp91 %in% c(5,7)
f00 <- rp00 %in% c(5,7)

#stack and mask the forest classes
f <- mask(stack(f51, f77, f91, f00), pr)

#check if the maps are 
check_landscape(f)

#aggregation index for the forest cover at each time point 
lsm_l_ai(f51)
lsm_l_ai(f77)
lsm_l_ai(f91)
lsm_l_ai(f00)

#SDMs
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

#asssing names to the raster stack layers 
names(Pred_stack_rp) <- names_pred

#here we resample the forest maps 
f_rs <- raster::resample( f, Pred_stack_rp , method="ngb")

#overaly by multipy the forest cover with the species models one forest map at the time 
Pred_f51 <- f_rs[[1]]* Pred_stack_rp 
Pred_f77 <- f_rs[[2]]* Pred_stack_rp
Pred_f91 <- f_rs[[3]]* Pred_stack_rp
Pred_f00 <- f_rs[[4]]* Pred_stack_rp

#assign names to the 
names(Pred_f51) <- names_pred
names(Pred_f77) <- names_pred
names(Pred_f91) <- names_pred
names(Pred_f00) <- names_pred

#agregation index using Landscape matrix for the thresh holded maps overlayed with the forest cover 
AI_Pred_f51 <-lsm_l_ai(Pred_f51)
AI_Pred_f77 <-lsm_l_ai(Pred_f77)
AI_Pred_f91 <-lsm_l_ai(Pred_f91)
AI_Pred_f00 <-lsm_l_ai(Pred_f00)

#nearest neiburgh distance 
D_Pred_f51 <-lsm_l_enn_mn(Pred_f51)
D_Pred_f77 <-lsm_l_enn_mn(Pred_f77)
D_Pred_f91 <-lsm_l_enn_mn(Pred_f91)
D_Pred_f00 <-lsm_l_enn_mn(Pred_f00)



#make a data frame and add names to the species Ai 
AI_total <- bind_cols(F_cover$X,AI_Pred_f51$value, AI_Pred_f77$value, AI_Pred_f91$value, AI_Pred_f00$value)
as.data.frame(AI_total)

#make a dataframe of the neerest neighboir mean distance 
D_total <- bind_cols(F_cover$X, D_Pred_f51$value, D_Pred_f77$value, D_Pred_f91$value, D_Pred_f00$value)
colnames(D_total) <- c("Code", "D_51", "D_71", "D_91", "D_00")
as.data.frame(D_total)

#Save the the AI 
write.csv(AI_total, "Data/Derived/AI_total.csv")
AI_total <- read.csv("Data/Derived/AI_total.csv", sep=";")

#plot ovverall change in ha
AI_mat <- as.matrix(AI_total[2:5])
#change row and column names 
rownames(AI_mat) <- AI_total$CODE
colnames(AI_mat) <- c(1951, 1977, 1991, 2000)

#plot total forest cover change 
matplot(t(AI_mat), type='l', lty=1, xaxt='n', xlab = "Time (years)", ylab = "Aggregation index", cex.lab=2, cex.axis=1.8)
axis(side=1, at= c(1, 2, 3, 4), labels = c(1951, 1977, 1991, 2000), cex.axis=1.8)

D_mat <- as.matrix(D_total[2:5])
matplot(t(D_mat), type='l', lty=1, xaxt='n', xlab = "Time (years)", ylab = "distance", , cex.lab=2, cex.axis=1.8)
axis(side=1, at= c(1, 2, 3, 4), labels = c(1951, 1977, 1991, 2000), cex.axis=1.8)


delta <- D_total$D_00 - D_total$D_51

plot(delta)







