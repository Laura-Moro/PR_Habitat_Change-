################################################
### NICHE CHARACTERISTICS
################################################

library(terra)
library(sf)

### Load occurrence records for extraction of values at points
pr_wgs <- st_read("Data/PR_shapes/outline/PR_outline_Project.shp")
load("Data/FINAL_RECS.Rdata")
occs <- st_as_sf(full, coords=c("LONGDEC","LATDEC"), crs=crs(pr_wgs))

### Load predicted maps
Pred_stack_raw <- rast("Data/Derived/raw_stack-20251022.tif")
Pred_stack_thres <- rast("Data/Derived/thresholded_stack-20251022.tif")

### Forest cover maps
f <- rast("Data/Derived/Forest_change_resampled_near.tif")

### Reproject points and lines
pr_lcc <- st_transform(pr_wgs, crs(f))
occs_lcc <- st_transform(occs, crs(f))

### Initialize results holders
nb <- matrix(nrow=dim(Pred_stack_raw)[3], ncol=4, 
             dimnames=list(names(Pred_stack_raw), 1:4))
nb_scaled <- matrix(nrow=dim(Pred_stack_raw)[3], ncol=4, 
                    dimnames=list(names(Pred_stack_raw), 1:4))
nb_thresh <- matrix(nrow=dim(Pred_stack_raw)[3], ncol=4, 
                    dimnames=list(names(Pred_stack_raw), 1:4))
nb_occs <- matrix(nrow=dim(Pred_stack_raw)[3], ncol=4, 
                  dimnames=list(names(Pred_stack_raw), 1:4))
maxvals <- matrix(nrow=dim(Pred_stack_raw)[3], ncol=4, 
                  dimnames=list(names(Pred_stack_raw), 1:4))
maxvals_occs <- matrix(nrow=dim(Pred_stack_raw)[3], ncol=4, 
                       dimnames=list(names(Pred_stack_raw), 1:4))
np <- matrix(nrow=dim(Pred_stack_raw)[3], ncol=4, 
             dimnames=list(names(Pred_stack_raw), 1:4))
np_occs <- matrix(nrow=dim(Pred_stack_raw)[3], ncol=4, 
                  dimnames=list(names(Pred_stack_raw), 1:4))

### Initialize temporary results holders
breaks <- list()
bin_id <- list()
breaks_nf <- list()
bin_id_nf <- list()

### Prep computations
n_bins <- 100 # No. of bins to use for nb calculations

### LOAD ENVIRONMENTAL PREDICTOR VARIABLES (derived in script 1b)
envs <- rast(list.files("Data/Derived/envs", full.names = TRUE))
envs_rp <- crop(project(envs, crs(f)), f)
envs_rp <- mask(envs_rp, pr_lcc)

### Set water geological type to NA
values(envs_rp[[1]])[values(envs_rp[[1]]) %in% 0] <- NA

### Get scaled envs into same units as the forest map
# envs_rp <- terra::project(envs, newforest, method='bilinear')
envs_scaled <- scale(envs_rp)
# envs_scaled_wgs <- scale(envs)
envs_scaled_newforest <- scale(mask(envs_rp, f, maskvalues=c(1,2,4)))

### Get median env values of new forest cells for np calculations below
median_envs_nf <- global(envs_scaled_newforest, median, na.rm=T)[,1]

### Make niche evenness function for soil categories
niche_evenness <- function(soil_rast, occ_pts) {
  # extract soil categories at occurrence points
  # occ_pts may be SpatVector or a matrix/data.frame of coords; terra::extract handles both
  occ_vals <- terra::extract(soil_rast, occ_pts)[,2]  # extract returns a data.frame: ID, value
  
  # remove NA occurrences (if any fall outside raster/masked)
  occ_vals <- occ_vals[!is.na(occ_vals)]
  if(length(occ_vals) == 0) return(NA_real_)
  
  # get all possible categories in the raster (exclude NA)
  all_vals <- unique(values(soil_rast))
  all_vals <- all_vals[!is.na(all_vals)]
  
  # ensure factor with levels for all categories so zeros are counted
  occ_fac <- factor(occ_vals, levels = all_vals)
  counts <- as.numeric(table(occ_fac))    # includes zeros for categories not used by occurrences
  total <- sum(counts)
  p <- counts / total
  
  # Shannon entropy (natural log)
  p_pos <- p[p > 0]
  H <- -sum(p_pos * log(p_pos))
  
  # Maximum possible entropy = log(K) where K = number of categories considered
  K <- length(all_vals)
  if(K <= 1) return(0)   # if only one category in landscape, evenness is degenerate -> 0
  
  Hmax <- log(K)
  J <- H / Hmax          # Pielou's evenness in [0,1]
  return(J)
}

### Prepare bins
for(r in 1:4){
  r1 <- envs_scaled[[r]]
  r1_nf <- envs_scaled_newforest[[r]]
  breaks[[r]] <- seq(minmax(r1)[1], minmax(r1)[2], length.out=n_bins + 1)
  breaks_nf[[r]] <- seq(minmax(r1_nf)[1], minmax(r1_nf)[2], length.out=n_bins + 1)
  bin_id[[r]] <- cut(values(r1), breaks=breaks[[r]],
                     include.lowest=T, labels=F)
  bin_id_nf[[r]] <- cut(values(r1_nf), breaks=breaks_nf[[r]], 
                        include.lowest=T, labels=F)
}

### Make function that computes the AUC proportion for a given binning
get_auc <- function(bin, vals) {
  max_per_bin <- tapply(vals, bin, median, na.rm = TRUE)
  max_per_bin[is.infinite(max_per_bin)] <- NA
  mean(max_per_bin, na.rm = TRUE)
}

### Run loop across species
for(sp in 1:nrow(np)){
  message(paste("=== Working on species", sp, "of", nlyr(Pred_stack_raw), "==="))
  
  # Process occurrence records values
  # occs <- full[full$CODE %in% names(Pred_stack_raw[[sp]]),]
  # ext <- terra::extract(envs_scaled_wgs, occs[,2:1], ID=F)
  occs <- occs_lcc[occs_lcc$CODE %in% names(Pred_stack_raw[[sp]]),]
  ext <- terra::extract(envs_scaled, occs, ID=F)
  nb_occs[sp,] <- apply(ext, 2, sd, na.rm=T)
  
  # Calculate soil nb using niche evenness function at top of script
  # nb_occs[sp,1] <- niche_evenness(envs$GEO, occs[,2:1])
  nb_occs[sp,1] <- niche_evenness(envs_rp$GEO, occs_lcc)
  
  maxvals_occs[sp,] <- apply(ext, 2, median, na.rm=T)
  
  # Re-scale prediction for scaled niche breadth
  p_unscaled <- p <- Pred_stack_raw[[sp]]
  v <- values(p)
  values(p) <- (v - min(v, na.rm=T))/(max(v, na.rm=T) - min(v, na.rm=T))
  
  # Get thresholded prediction
  p_thresh <- Pred_stack_thres[[sp]]
  
  # Store AUC for max suitability in each env bin, in each env layer
  nb[sp,] <- vapply(bin_id, get_auc, numeric(1), vals=values(p_unscaled))
  nb_scaled[sp,] <- vapply(bin_id, get_auc, numeric(1), vals=values(p))
  # nb_thresh[sp,] <- vapply(bin_id, get_auc, numeric(1), vals=values(p_thresh))
  
  # Niche position calculations: ID cells with max suitability
  p_nf <- terra::mask(p, f, maskvalues=c(1,2,4))
  max_cells <- which(values(p_nf)==max(values(p_nf), na.rm=T))
  
  # Extract environment values from max suitability cells
  vals_other <- values(envs_scaled_newforest)[max_cells,]
    
  # Get median value of cells with max suit, and calc diff from median of new forest cells
  if(!all(is.na(vals_other))){
    if(is.matrix(vals_other)){
      maxvals[sp,] <- apply(vals_other, 2, median, na.rm=T)
      np[sp,] <- apply(rbind(maxvals[sp,], median_envs_nf), 2, diff)
      np_occs[sp,] <- apply(rbind(maxvals_occs[sp,], median_envs_nf), 2, diff)
    }
    if(!is.matrix(vals_other)){
      maxvals[sp,] <- vals_other
      np[sp,] <- apply(rbind(maxvals[sp,], median_envs_nf), 2, diff)
      np_occs[sp,] <- apply(rbind(maxvals_occs[sp,], median_envs_nf), 2, diff)
    }
  }
}


### Run loop across species to recompute total weighted area change
d <- data.frame(sp=names(Pred_stack_raw), hab_change_raw=NA)
Pred_stack_raw51 <- global(mask(Pred_stack_raw, f, maskvalues=c(1,3)), sum, na.rm=T)
Pred_stack_raw00 <- global(mask(Pred_stack_raw, f, maskvalues=c(1,2)), sum, na.rm=T)
d$hab_change_raw <- (Pred_stack_raw00 - Pred_stack_raw51)[,1]


###############################################################
###############################################################
###############################################################
### Import data 
df <- read.csv("Data/Derived/6b-output-20251022.csv")

### Load fragmentation data
lm <- readRDS("Data/Derived/Species_Landscape-metrics-near-20251024.RDA")

### Set color gradient
pc1cols <- (RColorBrewer::brewer.pal(11, "Spectral"))[cut(df$pca1_pi, 11)]

### Look at average absolute np vs nb
plot(rowMeans(abs(np)), rowMeans(nb), pch=21, bg=pc1cols)
cor.test(rowMeans(abs(np)), rowMeans(nb))

## ADD NP AND NB VALUES TO DF
df$np <- as.vector(scale(rowMeans(abs(np), na.rm=T)[match(df$sp, rownames(np))]))
df$nb <- as.vector(scale(rowMeans(nb, na.rm=T)[match(df$sp, rownames(nb))]))
df$nb_scaled <- as.vector(scale(rowMeans(nb_scaled, na.rm=T)[match(df$sp, rownames(nb_scaled))]))
df$pca1_pi.z <- as.vector(scale(df$pca1_pi))
df$pca2_pi.z <- as.vector(scale(df$pca2_pi))

## Select fragmentation metric(s) to add to analysis data
df$enn_51 <- lm$value[lm$year==1951 & lm$metric=='enn_mn'][match(df$sp, lm$sp[lm$year==1951 & lm$metric=='enn_mn'])]
df$enn_00 <- lm$value[lm$year==2000 & lm$metric=='enn_mn'][match(df$sp, lm$sp[lm$year==2000 & lm$metric=='enn_mn'])]
df$enn_change <- (df$enn_00 - df$enn_51)
df$enn_change.z <- scale(df$enn_change)

df$ai_51 <- lm$value[lm$year==1951 & lm$metric=='ai'][match(df$sp, lm$sp[lm$year==1951 & lm$metric=='ai'])]
df$ai_00 <- lm$value[lm$year==2000 & lm$metric=='ai'][match(df$sp, lm$sp[lm$year==2000 & lm$metric=='ai'])]
df$ai_change <- (df$ai_00 - df$ai_51)
df$ai_change.z <- scale(df$ai_change)

df$mesh_51 <- lm$value[lm$year==1951 & lm$metric=='mesh'][match(df$sp, lm$sp[lm$year==1951 & lm$metric=='mesh'])]
df$mesh_00 <- lm$value[lm$year==2000 & lm$metric=='mesh'][match(df$sp, lm$sp[lm$year==2000 & lm$metric=='mesh'])]
df$mesh_change <- (df$mesh_00 - df$mesh_51)
df$mesh_change.z <- scale(df$mesh_change)


df$tot_change_raw.z <- as.numeric(scale(df$tot_change_raw))

par(mfcol=c(2,3))
plot(df$nb, df$enn_change, pch=21, bg=pc1cols)
cor.test(df$nb, df$enn_change, use='p')
plot(df$np, df$enn_change, pch=21, bg=pc1cols)
cor.test(df$np, df$enn_change, use='p')

plot(df$nb, df$ai_change, pch=21, bg=pc1cols)
cor.test(df$nb, df$ai_change, use='p')
plot(df$np, df$ai_change, pch=21, bg=pc1cols)
cor.test(df$np, df$ai_change, use='p')

plot(df$nb, df$mesh_change, pch=21, bg=pc1cols)
cor.test(df$nb, df$mesh_change, use='p')
plot(df$np, df$mesh_change, pch=21, bg=pc1cols)
cor.test(df$np, df$mesh_change, use='p')

### SAVE OUTPUT FROM ABOVE
# saveRDS(df, file="Data/Derived/post-nb-np_calcs-near.RDA")
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################


########################################################################
# Load data
df <- readRDS("Data/Derived/post-nb-np_calcs-near.RDA")

# Set color gradient
pc1cols <- rev(RColorBrewer::brewer.pal(11, "Spectral"))[cut(df$pca1_pi.z, 11)]

### Load and reproject life zone zone map
lz <- st_read("Data/PR_shapes/Lifezones/lifezones_Project.shp")
lz_lcc <- st_transform(lz, crs(f))


########################################################################
### FIGURE 1A
### Map of land cover change

fchange <- rast("Data/Derived/Forest_change_resampled_near.tif")
fchange_wgs <- terra::project(fchange, crs(pr_wgs), method="near")

plot(ext(-67.35, -65.5, 17.8, 18.8), border=NA)
plot(pr_wgs$geometry, axes=T, main="Forest cover change (1951-2000)", add=T)
plot(fchange_wgs, 
     col=c('black', viridis::viridis(4)[c(4,2,3)]),
     legend=F, main="Forest cover change (1951-2000)", add=T)
add_legend("bottomright", 
           legend=c("Nonforest", 
                    "Loss", 
                    "Gain", 
                    "Forest"), bty="n",
           # pt.bg=viridis::viridis(4)[c(1,4,2,3)], 
           pt.bg=c('black', viridis::viridis(4)[c(4,2,3)]),
           pch=22, pt.cex=2, cex=0.85)

plot(lz$geometry, border='white', add=T)
plot(pr_wgs$geometry, add=T)

########################################################################
# Compute percentage of total island in each category
(round(cbind(global(fchange==(1), sum, na.rm=T),
             global(fchange==(2), sum, na.rm=T),
             global(fchange==(3), sum, na.rm=T),
             global(fchange==(4), sum, na.rm=T)) / sum(!is.na(values(fchange))), 2))

values(fchange) <- as.factor(values(fchange))

fchangemap_sdf <- mask(fchange, lz_lcc[lz_lcc$ECOZONE=="df-S",])
fchangemap_smf <- mask(fchange, lz_lcc[lz_lcc$ECOZONE=="mf-S",])
fchangemap_swf <- mask(fchange, lz_lcc[lz_lcc$ECOZONE=="wf-S",])
fchangemap_lmwf <- mask(fchange, lz_lcc[lz_lcc$ECOZONE=="wf-LM",])
fchangemap_srf <- mask(fchange, lz_lcc[lz_lcc$ECOZONE=="rf-S",])
fchangemap_lmrf <- mask(fchange, lz_lcc[lz_lcc$ECOZONE=="rf-LM",])

change_table <- rbind(table(values(fchangemap_sdf))/sum(!is.na(values(fchangemap_sdf))),
                      table(values(fchangemap_smf))/sum(!is.na(values(fchangemap_smf))),
                      table(values(fchangemap_swf))/sum(!is.na(values(fchangemap_swf))),
                      table(values(fchangemap_lmwf))/sum(!is.na(values(fchangemap_lmwf))),
                      table(values(fchangemap_srf))/sum(!is.na(values(fchangemap_srf))),
                      table(values(fchangemap_lmrf))/sum(!is.na(values(fchangemap_lmrf))))

(change_table <- round(change_table, 2))

# Table 1. LULCC by life zone
round(100*(tapply(st_area(lz), lz$ECOZONE, sum) / sum(st_area(lz))), 1)

# Compute pixels gained / lost per species (*SLOW*)
gl.df <- data.frame(sp=rownames(np), g=NA, l=NA)
for(sp in 1:nrow(np)){
  message(paste("=== Working on species", sp, "of", nlyr(Pred_stack_raw), "==="))
  pred_thresh <- terra::resample(Pred_stack_thres[[sp]], f_rs, method="near")
  gl.df$l[sp] <- global(lostforest_mask * pred_thresh, sum, na.rm=T)
  gl.df$g[sp] <- global(newforest_mask * pred_thresh, sum, na.rm=T)
}

# Summary of gain / loss total at species level
range(unlist(gl.df$g)/unlist(gl.df$l))
mean(unlist(gl.df$g)/unlist(gl.df$l))
sd(unlist(gl.df$g)/unlist(gl.df$l))


########################################################################
### FIGURE 3
### Map of species richness

richnessmap <- sum(Pred_stack_thres)
# richnessmap_mask <- (richnessmap * newforest_mask)
# richnessmap_rp <- terra::project(richnessmap, crs(pr_wgs))

# richnessmap_highres <- terra::resample(richnessmap, f, method="bilinear")
# richnessmap_highres_mask_rp <- terra::project(richnessmap_highres * newforest_mask, pr_wgs)
# richnessmap_highres_rp <- terra::project(richnessmap_highres, pr_wgs)

richnessmap_mask <- mask(richnessmap, fchange, maskvalues=c(1,2,4))
richnessmap_mask_rp <- terra::project(richnessmap_mask, crs(pr_wgs), method='near')

# COMPUTE RICHNESS IN GAIN / LOST FORESTS BY LIFE ZONE
lz_richness_df <- matrix(nrow=6, ncol=3)
for(i in 1:length(unique(lz$ECOZONE))){
  message(i)
  foclz <- sort(unique(lz_lcc$ECOZONE))[i]
  lz_richness_df[i,] <- unlist(c(foclz,
                                 global(mask(richnessmap_mask, 
                                             lz_lcc[lz_lcc$ECOZONE==foclz,]), min, na.rm=T),
                                 global(mask(richnessmap_mask, 
                                             lz_lcc[lz_lcc$ECOZONE==foclz,]), max, na.rm=T)))
}
lz_richness_df


### FIGURE 3
plot(pr_wgs$geometry, col=rgb(0,0,0,1), axes=T, xlim=c(-67.25, -65.4))
# plot(richnessmap_rp, col=(viridis::inferno(100)[10:100]), add=T)
plot(richnessmap_mask_rp, col=(viridis::inferno(100)[10:100]), add=T)
plot(lz$geometry, border='white', add=T)
plot(pr_wgs$geometry, add=T)

# plot(pr$geometry, col=rgb(0,0,0,1), axes=T)#, xlim=c(-67.25, -65.4))
# plot(richnessmap_mask, col=(viridis::inferno(100)[10:100]), add=T)
# # plot(richnessmap_highres_mask_rp , col=(colorRamps::matlab.like(100)[5:95]), add=T)
# plot(lz_lcc$geometry, border='white', add=T)
# plot(pr$geometry, add=T)




########################################################################
### FIGURE 1 C, D
### Compute and plot amount and fragmentation metrics by life zone


par(mfrow=c(2,2), mar=c(4,6,2,1))

lzcols <- rep(c("grey","palegoldenrod","palegreen","seagreen2",
                "skyblue","royalblue1","purple3"), each=2)

dens <- rep(c(40,100), times=7)

b <- barplot(rbind(landscape_lz$enn_mn_51, landscape_lz$enn_mn_00)/1000, 
             beside=T, col=lzcols, names.arg=landscape_lz$ls, 
             las=2, density=dens, ylab="Mean euclidean NN dist. (km)",
             xlab="Life zone")
abline(v=3.5, lty=2)
axis(1, at=b, labels=NA)
legend("topright", legend=c(1951, 2000), 
       bty="n",  pt.cex=2, density=c(40,100), fill='grey', cex=1.2)


barplot(rbind(landscape_lz$ai_51,
              landscape_lz$ai_00), beside=T, col=lzcols,
        names.arg=landscape_lz$ls, las=2, ylim=c(0,100), dens=dens,
        ylab="Aggregation index", xlab="Life zone")
abline(v=3.5, lty=2)
axis(1, at=b, labels=NA)



### Plot the number and area of patches
barplot(rbind(landscape_lz$npatch_51,
              landscape_lz$npatch_00), beside=T, col=lzcols,
        names.arg=landscape_lz$ls, las=2,  dens=dens,
        ylab="Number of patches", xlab="Life zone")
abline(v=3.5, lty=2)
axis(1, at=b, labels=NA)

barplot(rbind(landscape_lz$patch_area_51,
              landscape_lz$patch_area_00), beside=T, col=lzcols,
        names.arg=landscape_lz$ls, las=2,  dens=dens,
        ylab="Mean patch area (ha)", xlab="Life zone")
abline(v=3.5, lty=2)
axis(1, at=b, labels=NA)







# plot(x, col=c('gold', 'black', 'blue'),
#      legend=F, main="Forest cover change (1951-2000)")
# add_legend("bottomright", legend=c("Loss", "Stable", "Gain"), bty="n",
#            pt.bg=c('gold', 'black', 'blue'), pch=22, pt.cex=2)


########################################################################
### FIGURE 1.B
### Density plot of environmental conditions in different LC categories

#total annual rainfall (R)
R <- rast(list.files("Data/ppt", full.names = TRUE)[13])/100

envtmp <- terra::resample(R, richnessmap_highres_mask_rp)

values(stable_forest)[values(stable_forest)==0] <- NA

ppt_newforest <- terra::mask(terra::project(R, newforest_mask), newforest_mask)
ppt_lostforest <- terra::mask(terra::project(R, lostforest_mask), lostforest_mask)
ppt_stable <- terra::mask(terra::project(R, stable_forest), stable_forest)

# Set up the plot area
x <- values(envtmp)[!is.na(values(envtmp))]
d <- stats::density(x, adjust=3)
d$y <- d$y * length(x) * (d$x[2] - d$x[1])

# Plot the frequency-scaled curve
# par(mar=c(6,6,2,2))
plot(d, type = "n", ylab = "Frequency", log='', main=NA, 
     xlab="Mean annual precipitation (mm)", axes=F, xlim=c(100,5000))
axis(1); axis(2)
polygon(d, col="grey50", border=NA)
# abline(v=median(x), lwd=2, lty=2)

legend("topright", legend=c("All island", 
                            "Loss", 
                            "Gain", 
                            "Stable forest"), bty="n",
       pt.bg=c('black', viridis::viridis(4)[c(4,2,3)]),
       pch=22, pt.cex=2, cex=0.85)

# Redraw the density curve on top
xgain <- values(ppt_newforest)[!is.na(values(ppt_newforest))]
d_gain <- stats::density(xgain, adjust=3)
d_gain$y <- d_gain$y * length(xgain) * (d_gain$x[2] - d_gain$x[1])
polygon(d_gain, col=scales::alpha(viridis::viridis(4)[2], 1),
        #border=viridis::viridis(4)[2], lwd=2
        border=NA)
# abline(v=median(xgain), lwd=2, lty=2, col=viridis::viridis(4)[2])

# Redraw the density curve on top
xstable <- values(ppt_stable)[!is.na(values(ppt_stable))]
d_stable <- stats::density(xstable, adjust=3)
d_stable$y <- d_stable$y * length(xstable) * (d_stable$x[2] - d_stable$x[1])
polygon(d_stable, col=scales::alpha(viridis::viridis(4)[3], 1),
        # border=viridis::viridis(4)[3], lwd=2
        border=NA
)
# abline(v=median(xstable), lwd=2, lty=2, col=viridis::viridis(4)[3])

# Draw the loss density
xloss <- values(ppt_lostforest)[!is.na(values(ppt_lostforest))]
d_loss <- stats::density(values(ppt_lostforest)[!is.na(values(ppt_lostforest))])
d_loss$y <- d_loss$y * length(xloss) * (d_loss$x[2] - d_loss$x[1])
polygon(d_loss, col=scales::alpha(viridis::viridis(4)[4], 1),
        # border=viridis::viridis(4)[4], lwd=2
        border=NA
)
# abline(v=median(xloss), lwd=2, lty=2, col=viridis::viridis(4)[4])




########################################################################
### Compile and summarize eval stats from SDMs

library(ENMeval)
modfiles <- list.files("Data/2022-12-07_ENMeval_results", full.names = TRUE)
splist <- unlist(lapply(strsplit(modfiles,"/"), function(x) gsub(".RDA", "", x[[3]])))
res <- list()

for(i in seq_along(modfiles)){
  message(paste("=== Working on species", sp, "of", nlyr(Pred_stack_raw), "==="))
  # read in the model files mod
  mod <- readRDS(modfiles[i])
  # Select the models with the minimum omission rate and the maximum AUC
  tmpres <- mod@results %>%
    filter(or.10p.avg == min(or.10p.avg, na.rm=T)) %>%
    filter(auc.val.avg == max(auc.val.avg, na.rm=T))
  # Save to results
  noccs <- nrow(mod@occs)
  sp <- splist[i]
  res[[i]] <- cbind(tmpres[1,], noccs, sp)
}

res <- do.call(rbind, res)

mean(res$auc.train)
sd(res$auc.train)
mean(res$or.10p.avg)
sd(res$or.10p.avg)


########################################################################
### FIGURE 2

library(factoextra)

### DO TRAIT PCA with phyloimputation
res.pca_pi <- prcomp(~ df$ss.log.z_pi +
                       wd.z_pi + 
                       thk.log.z_pi +
                       la.log.z_pi +
                       sla.log.z_pi +
                       maxht.z_pi, 
                     data=df)

# Add the first two principal components as multivariate traits
res.ind_pi <- get_pca_ind(res.pca_pi)$coord[,1:2]
rownames(res.ind_pi) <- df$sp[as.numeric(rownames(res.ind_pi))]
df$pca1_pi <- res.ind_pi[match(df$sp, rownames(res.ind_pi)),1]
df$pca2_pi <- res.ind_pi[match(df$sp, rownames(res.ind_pi)),2]
rownames(res.pca_pi$rotation) <- c("Seed", "WD", "Thk", "LA", "SLA", "Maxht")

### Habitat in 1951 vs. Gain
pccols <- (RColorBrewer::brewer.pal(11, "Spectral"))[cut(df$pca1_pi, 11)]

par(mfrow=c(2,2), mar=c(5,5,2,2))

plot(df$pca1_pi, df$pca2_pi, pch=16, 
     col=scales::alpha("grey", 0.5),
     # col=scales::alpha(pccols, 0.5),
     xlim=c(-5,5), ylim=c(-4,4),
     xlab="PC 1 (42%)", ylab="PC 2 (22%)")
arrows(0, 0, 
       5*res.pca_pi$rotation[,1],
       5*res.pca_pi$rotation[,2],
       len=0.075, lwd=1.5, col="blue")
text(4*res.pca_pi$rotation[,1], 
     4*res.pca_pi$rotation[,2], 
     rownames(res.pca_pi$rotation), 
     col='black', font=2, pos=c(1,2,1,4,1,3), adj=2)
mtext("A", line=-1.75, adj=0.025, cex=1.5)

df$fcover_51_raw.z <- scale(df$fcover_51_raw)[,1]
plot(df$fcover_51_raw.z, df$tot_change_raw.z, 
     pch=21, bg=pccols,
     xlab="Potential habitat 1951", ylab="Gain in potential habitat 1951-2000")
legend("bottomright", legend=c("High PC1", "", "","", "Low PC1"),
       pt.bg=rev(RColorBrewer::brewer.pal(11, "Spectral"))[c(1,3,5,8,11)], 
       pch=21, bty='n', pt.cex = 1.5, cex=0.75)
mtext("B", line=-1.75, adj=0.025, cex=1.5)



# par(mfrow=c(2,2), mar=c(4,4,1,1))
# plot(df$np, df$nb, pch=21, bg=cols,
#      xlab="Niche position (marginality)", ylab="Niche breadth")
# cor.test(df$np, df$nb)
# 
# plot(df$pca1_pi.z, df$np, pch=21, bg=cols,
#      xlab="PC1", ylab="Niche position (marginality)")
# cor.test(df$pca1_pi.z, df$np)
# 
# plot(df$pca1_pi.z, df$nb, pch=21, bg=cols,
#      xlab="PC1", ylab="Niche breadth")
# cor.test(df$pca1_pi.z, df$nb)


########################################################################
### RESULTS - CORRELATION OF NP, NB, TRAITS
cor.test(df$np, df$pca1_pi.z)
cor.test(df$nb, df$pca1_pi.z)

# Figure S1

par(mfrow=c(2,2), mar=c(5,5,2,1))
plot(df$np, df$nb, pch=21, bg=pc1cols,
     xlab="Niche position", ylab="Niche breadth")
mtext("A", adj=0.05, cex=1.5, line=0.1)
cor.test(df$np, df$nb, use='p')
plot(df$pca1_pi.z, df$nb, pch=21, bg=pc1cols,
     xlab="PCA axis 1", ylab="Niche breadth")
mtext("B", adj=0.05, cex=1.5, line=0.1)
cor.test(df$pca1_pi.z, df$nb, use='p')
plot(df$np, df$pca1_pi.z, pch=21, bg=pc1cols,
     xlab="Niche position", ylab="PCA axis 1")
mtext("C", adj=0.05, cex=1.5, line=0.1)
cor.test(df$np, df$pca1_pi.z, use='p')



########################################################################
### FIGURE 4

pccols <- (RColorBrewer::brewer.pal(11, "Spectral"))[cut(df$pca1_pi, 11)]


par(mfrow=c(2,3), mar=c(4,4,1,1), oma=c(4,2,1,1))
plot(df$np, df$tot_change_raw, pch=21, bg=pccols, 
     xlab="Niche position", ylab="Potential habitat change", 
     # xlim=c(-2.5,3.3), ylim=c(-2.3,1.7)
)
abline(h=0, lty=2)
mtext("A", adj=0.025)

plot(df$nb, df$tot_change_raw, pch=21, bg=pccols,
     xlab="Niche breadth", ylab="Potential habitat change")
abline(h=0, lty=2)
mtext("B", adj=0.025)

plot(df$pca1_pi.z, df$tot_change_raw, pch=21, bg=pccols,
     xlab="PCA axis 1", ylab="Potential habitat change", 
     # xlim=c(-3.3,3.5), ylim=c(-2.4,1.7)
)
abline(h=0, lty=2)
mtext("C", adj=0.025)

### Panels for enn
plot(df$np, df$enn_change/1000, pch=21, bg=pccols, 
     xlab="Niche position", ylab="Patch distance change (km)",
     # ylim=c(-2,2)
)
abline(h=0, lty=2)
mtext("D", adj=0.025)

plot(df$nb, df$enn_change/1000, pch=21, bg=pccols,
     xlab="Niche breadth", ylab="Patch distance change (km)", 
     # ylim=c(-2,2)
)
abline(h=0, lty=2)
mtext("E", adj=0.025)

plot(df$pca1_pi.z, df$enn_change/1000, pch=21, bg=pccols,
     xlab="PCA axis 1", ylab="Patch distance change (km)", 
     # ylim=c(-2,2)
)
abline(h=0, lty=2)
mtext("F", adj=0.025)

# plot(df$np, df$ai_change, pch=21, bg=pccols, 
#      xlab="Niche position", ylab="Aggregation change")
# abline(h=0, lty=2)
# mtext("G", adj=0.025)
# 
# plot(df$nb, df$ai_change, pch=21, bg=pccols,
#      xlab="Niche breadth", ylab="Aggregation change")
# abline(h=0, lty=2)
# mtext("H", adj=0.025)
# 
# plot(df$pca1_pi.z, df$ai_change, pch=21, bg=pccols,
#      xlab="PCA axis 1", ylab="Aggregation change")
# abline(h=0, lty=2)
# mtext("I", adj=0.025)



df$tot_change_raw.z

df$pca1_pi.z

# For model, remove rows where NA for traits
df2 <- df[!is.na(df$pca1_pi.z),]
nrow(df2) # model is fit for 413 species 

m1 <- lm(tot_change_raw.z ~ np + nb + pca1_pi.z, data=df2)
m2 <- lm(tot_change_raw.z ~ np + pca1_pi.z, data=df2)
m3 <- lm(tot_change_raw.z ~ nb + pca1_pi.z, data=df2)
m4 <- lm(tot_change_raw.z ~ np + nb, data=df2)
m5 <- lm(tot_change_raw.z ~ np, data=df2)
m6 <- lm(tot_change_raw.z ~ nb, data=df2)
m7 <- lm(tot_change_raw.z ~ pca1_pi.z, data=df2)

AIC(m1,m2,m3,m4,m5,m6,m7)

summary(m1)



# For model, remove rows where NA for traits
df2 <- df[!is.na(df$pca1_pi.z),]
nrow(df2) # model is fit for 413 species 

n1 <- lm(enn_change.z ~ np + nb + pca1_pi.z, data=df2)
n2 <- lm(enn_change.z ~ np + pca1_pi.z, data=df2)
n3 <- lm(enn_change.z ~ nb + pca1_pi.z, data=df2)
n4 <- lm(enn_change.z ~ np + nb, data=df2)
n5 <- lm(enn_change.z ~ np, data=df2)
n6 <- lm(enn_change.z ~ nb, data=df2)
n7 <- lm(enn_change.z ~ pca1_pi.z, data=df2)

AIC(n1,n2,n3,n4,n5,n6,n7)

summary(n1)







sp <- 1
plot(nb[sp,], type='l', ylim=c(0,1))
lines(nb_scaled[sp,], col='red')
lines(nb_thresh[sp,], col='blue')
lines(nb_occs[sp,]/2, col='green')









plot(pr$geometry, col='grey')
plot(sum(Pred_stack_thres * newforest_mask), 
     col=viridis::plasma(100), add=T, axes=T)


plot(Pred_stack_raw[[1]])
plot(Pred_stack_raw[[2]])


v <- Pred_stack_thres * df$WD[match(names(Pred_stack_thres), df$sp)]
v <- terra::mask(v, v, maskvalues=0)

par(mfrow=c(2,2), mar=c(4,6,1,1))
plot(pr$geometry, col='lightgrey', axes=T)
plot(mean(v, na.rm=T) * newforest_mask, add=T)

hist(mean(v, na.rm=T) * newforest_mask, 
     xlab="wood density (g / cm^3)", main="Mean wood density")




v1 <- values(envs_scaled_newforest)
v2 <- values(sum(Pred_stack_thres * newforest_mask))

par(mfrow=c(2,2), mar=c(4,6,1,1))

plot(jitter(v1[,1]), v2, col=viridis::viridis(100)[cut(v2, 100)], pch=16, cex=0.25,
     xlab="Soil type", ylab="Predicted richness")
plot(v1[,2], v2, col=viridis::viridis(100)[cut(v2, 100)], pch=16, cex=0.25,
     xlab="Pm", ylab="Predicted richness")
plot(v1[,3], v2, col=viridis::viridis(100)[cut(v2, 100)], pch=16, cex=0.25,
     xlab="PSI", ylab="Predicted richness")
plot(v1[,4], v2, col=viridis::viridis(100)[cut(v2, 100)], pch=16, cex=0.25,
     xlab="TM", ylab="Predicted richness")















