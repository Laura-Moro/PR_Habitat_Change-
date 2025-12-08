library(raster)
library(terra)
library(pracma) # for area of polygons (niche breadth calculation)

# ========================================================
# ========================================================
# 
# ### Import data 
# df <- read.csv("Data/Derived/6b-output-20250829.csv")
# # The amount of habitat change favored low wood density species
# # Species with higher wood density had, on average, lower gains of available habitat 
# plot(df$WD, df$tot_change_raw)
# abline(lm(df$tot_change_raw ~ df$WD))
# summary(lm(df$tot_change_raw ~ df$WD))
# 
# # Higher values of PC1 - phylogenetic imputation were more acquisitive traits
# plot(df$pca1_pi, df$wd.z)
# 
# # Whereas higher values of PC1 raw where more conservative traits
# plot(df$pca1, df$wd.z)
# 
# # Species with more aquisitive strategies gained more total habitat
# plot(df$pca1_pi, df$tot_change_raw)
# abline(lm(df$tot_change_raw ~ df$pca1_pi))
# summary(lm(df$tot_change_raw ~ df$pca1_pi))
# plot(df$pca1, df$tot_change_raw)
# abline(lm(df$tot_change_raw ~ df$pca1))
# summary(lm(df$tot_change_raw ~ df$pca1))
# 
# # Species with more aquisitive traits had higher gains of habitat even proportional to their already generally larger habitat in 1951
# 
# col <- rev(RColorBrewer::brewer.pal(11, "Spectral"))[cut(df$wd.z_pi, 11)]
# plot(I(df$tot_change_raw/df$fcover_51_raw) ~ df$pca1_pi, pch=21, bg=col)
# summary(lm(I(df$tot_change_raw/df$fcover_51_raw) ~ df$pca1_pi))
# abline(lm(I(df$tot_change_raw/df$fcover_51_raw) ~ df$pca1_pi))

# plot(I(df$tot_change_raw/df$fcover_51_raw) ~ df$pca1, pch=21, bg=col)
# summary(lm(I(df$tot_change_raw/df$fcover_51_raw) ~ df$pca1))
# abline(lm(I(df$tot_change_raw/df$fcover_51_raw) ~ df$pca1))

# Was there a bias in terms of proportional gains of habitat for species that have different niche positions (i.e., moist forest species had proportainally higher gains of habitat than dry forest species)??
# ========================================================
# ========================================================

# Load PR outline shapefile
pr_wgs <- st_read("Data/PR_shapes/outline/PR_outline_Project.shp")
pr <- st_transform(pr_wgs, CRS("+proj=lcc +lat_0=17.8333333333333 +lon_0=-66.4333333333333 +lat_1=18.0333333333333 +lat_2=18.4333333333333 +x_0=152400.3048 +y_0=0 +datum=NAD27 +units=m +no_defs"))

# Load the predictions 
# pred <- stack(list.files(path="Data/SDM_predictions", pattern='.tif', all.files=T, full.names=T))

# Load occurrence records for extraction of values at points
load("Data/FINAL_RECS.Rdata")

# Load SDM predictions and forest change maps that are in same extent / projection
# pred_rs <- stack("Data/Derived/raw_stack.tif")

### Full resolution forest cover maps
f_rs <- rast("Data/Derived/Forest_change_resampled.tif")
newforest <- (f_rs[[4]] - f_rs[[1]])

# Make mask of the new forest gains
newforest_mask <- lostforest_mask <- newforest
values(newforest_mask)[values(newforest_mask) <= 0] <- NA
values(lostforest_mask)[values(lostforest_mask) >= 0] <- NA
values(lostforest_mask) <- abs(values(lostforest_mask))

### Plot and compute richness potential in gain/lost forest
# plot(newforest_mask)
# plot(lostforest_mask)
# richness_in_new_forests <- sum(Pred_stack_thres * newforest_mask)
# richness_in_lost_forests <- sum(Pred_stack_thres * lostforest_mask)
### Plot potential richness in new forests
# plot(pr$geometry, col=rgb(0,0,0,0.05), axes=T)
# plot(richness_in_new_forests, col=viridis::inferno(100), add=T)
### Plot potential richness in lost forests
# plot(pr$geometry, col=rgb(0,0,0,0.05), axes=T)
# plot(richness_in_lost_forests, col=viridis::inferno(100), add=T)

### Plot forest change map
# plot(newforest, main="Forest change 1951-2000 (loss, remain, gain)")

### Full resolution SDM predictions
# Import the continuous maps 
Pred_stack_raw <- terra::rast("Data/Derived/raw_stack.tif")

# Import the continuous maps 
Pred_stack_thres <- terra::rast("Data/Derived/thresholded_stack.tif")

# Get species names vector
modfiles <- list.files("Data/SDM_predictions", full.names = TRUE)
splist <- unlist(lapply(strsplit(modfiles,"/"), function(x) gsub(".tif", "", x[[3]])))
names(Pred_stack_raw) <- names(Pred_stack_thres) <- splist

# LOAD ENVIRONMENTAL PREDICTOR VARIABLES derived in script 1b
envs <- rast(list.files("Data/Derived/envs", full.names = TRUE))
envs <- crop(envs, pr_wgs)

# Set water geological type to NA
values(envs[[1]])[values(envs[[1]]) %in% 0] <- NA

# Get scaled envs into same units as the forest map
envs_scaled <- scale(terra::project(envs, newforest, method='near'))
envs_scaled_newforest <- scale(mask(terra::project(envs, newforest, method='near'), newforest_mask))

# res_thresh <- res_raw <- matrix(nrow=dim(Pred_stack_raw)[3], ncol=2)

### ***------------------------------*****
### ***----- SCRIPT IN LOOP NEEDS TO BE ADJUSTED AFTERWARDS FOR PROPER NICHE CALCULATIONS
### ***----- ALSO NEEDS TO BE CHECKED - WHY DOES NICHE POSITION GO TO THE SAME MAX ON EITHER SIDE OF THE MOST COMMON VALUE, EVEN IF MAXSUIT VALUES ARE MUCH HIGHER ?
### ***------------------------------*****

### Initialize results holders
nafill <- rep(NA, length(splist))

maxsuit <- data.frame(row.names=splist, maxsuit1=nafill, 
                      maxsuit2=nafill, maxsuit3=nafill, maxsuit4=nafill)
maxsuit_occs <- data.frame(row.names=splist, maxsuit1=nafill, 
                           maxsuit2=nafill, maxsuit3=nafill, maxsuit4=nafill)

np_sdm <- data.frame(row.names=splist, np1=nafill, np2=nafill, np3=nafill, np4=nafill)
np_sdm_island <- data.frame(row.names=splist, np1=nafill, np2=nafill, np3=nafill, np4=nafill)

np_occs <- data.frame(row.names=splist, np1=nafill, np2=nafill, np3=nafill, np4=nafill)
np_occs_island <- data.frame(row.names=splist, np1=nafill, np2=nafill, np3=nafill, np4=nafill)

nb_sdm <- data.frame(row.names=splist, nb1=nafill, nb2=nafill, nb3=nafill, nb4=nafill)
nb_occs <- data.frame(row.names=splist, nb1=nafill, nb2=nafill, nb3=nafill, nb4=nafill)


for (sp in 1:nlayers(pred)) {
  
  message(paste("=== Working on species", sp, "of", nlayers(pred), "==="))
  
  p_unscaled <- p <- Pred_stack_raw[[sp]]
  v <- values(p)
  
  # Skip the species if the prediction is all NA
  if(!all(is.nan(v))){
  
    ### ***------------------------------*****
    # FOR NICHE BREADTH, DOES THIS MAKE SENSE?  TO RELATIVIZE AT THE SPECIES LEVEL OR ACCOUNT FOR LOWER PREVALENCE ACROSS SPECIES?
    # Re-scale habitat suitability values between 0-1
    values(p) <- (v - min(v, na.rm=T))/(max(v, na.rm=T) - min(v, na.rm=T))

    # Define thresholds for niche breadth calculation
    thresholds <- seq(0, 1, 0.01)
    
    # Initialize holder for range of env values in each bin of suitability
    out <- out_unscaled <- matrix(ncol=8, nrow=length(thresholds))
    
    # Loop through suitability bins to get range of env values in each
    for (j in 1:nrow(out)){
      message(j)
      tmp <- ((p >= thresholds[j]) * 1) * envs_scaled_rp_newforest 
      tmp_unscaled <- ((p_unscaled >= thresholds[j]) * 1) * envs_scaled_rp_newforest 
      
      # for(l in 2:dim(tmp)[3]){
      #   values(tmp[[l]]) <- ifelse(values(tmp[[l]]) %in% 0, NA, values(tmp[[l]]))
      # }
      
      out[j,1:2] <- sum(!unique(values(tmp[[1]])) %in% c(0, NA))
      out[j,1:2] <- sum(!unique(values(tmp_unscaled[[1]])) %in% c(0, NA))
      
      out[j,3:8] <- unlist(global(tmp[[2:4]], range, na.rm=T))
      out[j,3:8] <- unlist(global(tmp_unscaled[[2:4]], range, na.rm=T))
    }
    
    ### Compute niche position as median of env values in grid cells with max suitability
    maxvals <- which(values(p) == 1)
    if(length(maxvals)>1){
      maxsuit[sp,2:4] <- apply(values(envs_scaled[[2:4]])[maxvals,], 2, median)
    } else {
      maxsuit[sp,2:4] <- values(envs_scaled[[2:4]])[maxvals,]
    }
    
    # Get the soil category with most cells in the max suitability prediction
    x <- as.numeric(names(sort(table(values(envs_scaled[[1]])[maxvals]), decreasing=T))[1])
    if(length(x)==1){ 
      
      maxsuit[sp, 1] <- x 
      
      # Proportion of grid cells in each soil category in max suit cells
      samp <- table(values(envs_scaled[[1]])[maxvals])/sum(table(values(envs_scaled[[1]])[maxvals]))
      
      # Overall proportion of grid cells in each soil category
      all <- table(values(envs_scaled[[1]])) / sum(table(values(envs_scaled[[1]])))
      
      # Bind the two rows above, take squared sum of absolute diffs of proportions
      np[sp, 1] <- sum((apply(rbind(samp[match(names(all), names(samp))], all), 2, diff))^2, na.rm=T)
      
    }
    
    ### Compute niche position as difference between median env values in grid cells with max suitability and the median of env values in the landscape. Uses the empirical cumulative density function.
    ecdf_list <- apply(values(envs_scaled_newforest[[2:4]]), 2, ecdf)
    np[sp,2:4] <- (unlist(Map(function(f, v) f(v), ecdf_list, maxsuit[sp,2:4])) - 0.5)^2
    
    
    ### Niche breadth computations
    for(r in 2:4){
      out2 <- out[,1:2]
      out2[,1] <- as.numeric(global(envs_scaled[[r]], min, na.rm=T))
      out2[,2] <- as.numeric(global(envs_scaled[[r]], max, na.rm=T))
      
      lo <- rbind(1:2, 3:4, 5:6, 7:8)[r,][1]
      hi <- rbind(1:2, 3:4, 5:6, 7:8)[r,][2]
      
      a1 <- polyarea(y=c(out[,lo], rev(out[,hi])), x=c(1:101, 101:1)) 
      a2 <- polyarea(y=c(out2[,1], rev(out2[,2])), x=c(1:101, 101:1))
      nb[sp, r] <- a1/a2
      
      # Soil niche breadth is the average 
      nb[sp, 1] <- sum(out[,1])/sum(rep(6, nrow(out)))
      
    }
  }
}


###  NEED TO PASS THIS ONWARDS>>>

# Visualize niche position vs. environment at most suitable habitat
plot(maxsuit[,1], np[,1])
plot(maxsuit[,2], np[,2])
plot(maxsuit[,3], np[,3])
plot(maxsuit[,4], np[,4])

plot(maxsuit[,1], nb[,1])
plot(maxsuit[,2], nb[,2])
plot(maxsuit[,3], nb[,3])
plot(maxsuit[,4], nb[,4])


# Negative correlation between overall niche position and niche breadth
plot(rowMeans(np), rowMeans(nb), pch=21, bg=pcacols)
cor.test(rowSums(np), rowSums(nb), use='p')

# Merge overall NP / NB metrics with the main dataframe
df$np <- rowMeans(np)[match(df$sp, rownames(np))]
df$nb <- rowMeans(nb)[match(df$sp, rownames(nb))]



### Statistical model
m1 <- lm(tot_change_raw.z ~ wd.z_pi + np + nb, data=df)
summary(m1)
car::vif(m1)








wdcols <- rev(RColorBrewer::brewer.pal(11, "Spectral"))[cut(df$wd.z_pi,11)]
wdcols <- rev(viridis::plasma(20))[cut(df$wd.z_pi,20)]
wdcols <- rev(viridis::plasma(20))[cut(df$sla.log.z_pi,20)]
pcacols <- (RColorBrewer::brewer.pal(11, "Spectral"))[cut(df$pca1_pi,11)]
pcacols <- rev(viridis::plasma(20))[cut(df$pca1_pi,20)]


par(mfrow=c(2,2), mar=c(4,4,1,1))

# plot(df$np, df$nb, pch=21, bg=pcacols, lwd=0.5, log='x')
# plot(df$np, df$wd.z, pch=21, bg=pcacols, lwd=0.5)
# plot(df$nb, df$wd.z, pch=21, bg=pcacols, lwd=0.5)

plot(df$np, df$tot_change_raw.z, pch=21, bg=wdcols, lwd=0.5)
cor.test(df$np, df$tot_change_raw, use='p')

plot(df$nb, df$tot_change_raw.z, pch=21, bg=wdcols, lwd=0.5)
cor.test(df$NB, df$tot_change_raw, use='p')

plot(df$wd.z, df$tot_change_raw.z, pch=21, bg=wdcols, lwd=0.5)
cor.test(df$tot_change_raw, df$WD, use='p')


df$fcover_00_raw.z <- as.vector(scale(df$fcover_00_raw))

par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(df$fcover_51_raw.z, df$tot_change_raw.z, pch=21, bg=wdcols, lwd=0.5, 
     xlab="Available habitat in 1951", 
     ylab="Increase in available habitat 1951-2000")

hist(df$tot_change_raw)


df$np.z <- as.vector(scale(df$np))
df$nb.z <- as.vector(scale(df$nb))

m1 <- lm(tot_change_raw.z ~ np.z + nb.z + pca1_pi, data=df)
summary(m1)
car::vif(m1)

par(mfrow=c(2,2), mar=c(4,4,1,1))

plot(df$np.z, df$tot_change_raw.z, pch=21, bg=wdcols, lwd=0.5, 
     xlab="Niche position (marginality)", ylab="Scaled habitat change")
lines(seq(-5,5,0.1), lwd=2, lty=2,
      predict(m1, newdata=data.frame(np.z=seq(-5,5,0.1), nb.z=0, pca1_pi=0)))

plot(df$nb.z, df$tot_change_raw.z, pch=21, bg=wdcols, lwd=0.5,
     xlab="Niche breadth", ylab="Scaled habitat change")
lines(seq(-5,5,0.1), lwd=2,
      predict(m1, newdata=data.frame(np.z=0, nb.z=seq(-5,5,0.1), pca1_pi=0)))

plot(df$pca1_pi, df$tot_change_raw.z, pch=21, bg=wdcols, lwd=0.5,
     xlab="Trait PCA 1", ylab="Scaled habitat change")
lines(seq(-5,5,0.1), lwd=2,
      predict(m1, newdata=data.frame(np.z=0, nb.z=0, pca1_pi=seq(-5,5,0.1))))



plot(df$tot_change_raw, df$tpa_diam2_2014-df$tpa_diam2_2004, 
     pch=21, bg=wdcols, lwd=0.5, log='')



plot(df$NB.z, df$tpa2004, pch=21, bg=wdcols, lwd=0.5, log='y')
plot(df$NP.z, df$tpa2014, pch=21, bg=wdcols, lwd=0.5, log='y')



plot(c(out[,1], rev(out[,2])), 
     c(1:101, 101:1), type="l")
plot(pred[[sp]])
plot(envs[[2]])    
  

out2 <- out
out2[,1] <- cellStats(envs[[2]], min)
out2[,2] <- cellStats(envs[[2]], max)

plot(c(out2[,1], rev(out2[,2])), 
     c(1:101, 101:1), type="l", col='blue')
lines(c(out[,1], rev(out[,2])), 
     c(1:101, 101:1), type="l", col='red')

polygon(x=c(out2[,1], rev(out2[,2])),
        y=c(1:101, 101:1), col='blue')
polygon(x=c(out[,1], rev(out[,2])),
        y=c(1:101, 101:1), col='red')
abline(v=pptatmaxsuit[sp], lwd=3)



polyarea(y=c(out[,1], rev(out[,2])), x=c(1:101, 101:1))/polyarea(y=c(out2[,1], rev(out2[,2])), x=c(1:101, 101:1))

### MODEL SELECTION 

# Load mod files #mod = enmeval result object 
# modfiles <- list.files("Data/2022-12-07_ENMeval_results", full.names = TRUE)
modfiles <- list.files("Data/2025-03-17_ENMeval_results-noFIA", full.names = TRUE)

# Select the models with the the lowest omission rate at 10p and the highest AUC
for(sp in seq_along(modfiles)){
  message(sp)
  
  #red in the model files mod
  mod <- readRDS(modfiles[sp]) 
  
  # Select the models with the minimum omission rate and the maximum AUC
  tmpres <- mod@results %>%
    filter(or.10p.avg == min(or.10p.avg, na.rm=T)) %>%
    filter(auc.val.avg == max(auc.val.avg, na.rm=T))
  
  # Take the full names (file path) and extract the just the code of the species 
  tmpres$species <- gsub(".RDA", "", strsplit(modfiles[i],"/")[[1]][3])
  
  focmod <- mod@models[[which(names(mod@predictions)==tmpres[1,]$tune.args)]]
  
  plot(focmod, vars = names(focmod$samplemeans), 
       common.scale = T, 
       type = "logistic")
  

}

       
c("link", "exponential", "cloglog", "logistic"))

plot(focmod, vars="Pm", type = "logistic")

m <- data.frame(matrix(focmod$samplemeans, 100, length(focmod$samplemeans), byrow = T))
colnames(m) <- names(mm)

predict(focmod)

focmod$varmin["Pm"]












x <- resample(rast(rall[[1]]), ppt)
plot(x)

plot(rast(rall[[1]]))
plot(ppt)




# Are those disproportate gains reflected in the abundance across the island, corrected for habitat area?



plot(ppt)
plot(rast(f00>0))

