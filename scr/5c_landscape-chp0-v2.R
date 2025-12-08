################################################
### FULL RESOLUTION METRICS (AGGREGATED BELOW)
################################################

library(terra)
library(sf)
library(landscapemetrics)

### Load forest maps
r51 <- rast("Data/Maps_1951-2000/puerto51_sub1_100905_landcov_final.img")
r00 <- rast("Data/Maps_1951-2000/pr2000_100805_final_quarry_recode2landcov_subset.img")
r00 <- terra::project(r00, r51, method="near")

### Classify forest / non-forest pixels
f51 <- r51 %in% 5
f00 <- r00 %in% c(5,7)

### Load PR outline shapefile
pr_wgs <- st_read("Data/PR_shapes/outline/PR_outline_Project.shp")
pr_lcc <- st_transform(pr_wgs, crs(f51))

### Stack and mask the forest classes
f <- crop(mask(c(f51, f00), pr_lcc), pr_lcc)

### Get a forest change map
newforest <- (f[[2]] - f[[1]])
stable_forest <- (f[[2]]==1 & f[[1]]==1)
stable_nonforest <- (f[[2]]==0 & f[[1]]==0)
fchange <- ((stable_nonforest*-1) + stable_forest + f[[2]])+2

plot(fchange, col=c('black', viridis::viridis(4)[c(4,2,3)]))
plot(lz_lcc$geometry, border='white', add=T)

### Load and reproject life zone zone map
lz <- st_read("Data/PR_shapes/Lifezones/lifezones_Project.shp")
lz_lcc <- st_transform(lz, crs(fchange))

### Assume the polygon attribute with zone ID is called "zone_id"
zone_rast <- rasterize(vect(lz_lcc), fchange, field="ECOZONE")

### Combine forest change map and life zones; unique ID for each combination
combined <- as.factor(fchange + zone_rast*4)

### Make a reclassification table 
rcl <- data.frame(id=1:24,
                  change=rep(1:4, times=6),
                  change_cat=rep(c("nonforest","loss","gain","forest"), times=6),
                  zone=rep(levels(as.factor(lz_lcc$ECOZONE)), each=4),
                  zoneorder=c(rep(1,4), rep(2,4), rep(5,4), rep(6,4), rep(4,4), rep(3,4)),
                  for51=rep(c(F,T,F,T), times=6),
                  for00=rep(c(F,F,T,T), times=6))
rcl$zone_for51 <- as.numeric(as.factor(paste(rcl$zone, rcl$for51)))
rcl$zone_for00 <- as.numeric(as.factor(paste(rcl$zone, rcl$for00)))

### Subset change map to forest-only in 1951 and 2000
f51_mask <- mask(combined, combined, maskvalues=rcl$id[!rcl$for51])
f51_mask_noRF <- mask(combined, combined, maskvalues=rcl$id[!rcl$for51 | grepl("rf", rcl$zone)])
f00_mask <- mask(combined, combined, maskvalues=rcl$id[!rcl$for00])
f00_mask_noRF <- mask(combined, combined, maskvalues=rcl$id[!rcl$for00 | grepl("rf", rcl$zone)])

### Reclassify
f51_reclass <- classify(f51_mask, rcl[,c("id","zone_for51")])
f51_reclass_noRF <- classify(f51_mask_noRF, rcl[,c("id","zone_for51")])
f00_reclass <- classify(f00_mask, rcl[,c("id","zone_for00")])
f00_reclass_noRF <- classify(f00_mask_noRF, rcl[,c("id","zone_for00")])

### Metrics to calculate on all classes and full landscape (all classes) 
lsm1 <- calculate_lsm(combined, what=c("lsm_c_ca", "lsm_l_ta"))
lsm1 <- cbind(lsm1, rcl[match(lsm1$class, rcl$id),])
lsm1 <- lsm1[order(lsm1$zoneorder),]

### Metrics to calculate on all forest in 51 or all forest in 00
m51 <- calculate_lsm(f51_reclass, what=c("lsm_c_ca", "lsm_c_np", "lsm_c_area_mn", 
                                    "lsm_c_enn_mn", "lsm_c_ai", "lsm_c_mesh",
                                    "lsm_l_area_mn", "lsm_l_np", 
                                    "lsm_l_ai", "lsm_l_mesh"))
lsm2 <- lsm_l_enn_mn(f51_reclass_noRF)
m51 <- rbind(m51, lsm2)

m00 <- calculate_lsm(f00_reclass, what=c("lsm_c_ca", "lsm_c_np", "lsm_c_area_mn", 
                                         "lsm_c_enn_mn", "lsm_c_ai", "lsm_c_mesh", 
                                         "lsm_l_area_mn", "lsm_l_np", 
                                         "lsm_l_ai", "lsm_l_mesh"))
lsm3 <- lsm_l_enn_mn(f00_reclass_noRF)
m00 <- rbind(m00, lsm3)

m51$year <- 1951
m00$year <- 2000

lm <- rbind(m51, m00)

lm$zone <- rcl$zone[match(lm$class, rcl$zone_for51)]
lm$zoneorder <- rcl$zoneorder[match(lm$class, rcl$zone_for51)]

lm <- lm[order(lm$zoneorder),]

saveRDS(lm, file="Data/Derived/Landscape_metrics_full_resolution.RDA")


############################################
############# DIAGNOSTIC PLOTS #############
### PLOT IT!
par(mfrow=c(2,2), mar=c(4,6,2,1))
lzcols <- rep(c("grey","palegoldenrod","palegreen","seagreen2",
                "skyblue","royalblue1","purple3"), each=2)
dens <- rep(c(40,100), times=7)


get_metric <- function(metric){
  tmp <- rbind(c(lm$value[lm$level=="landscape" & lm$metric==metric & lm$year==1951],
               lm$value[lm$level=="class" & lm$metric==metric & lm$year==1951]),
             c(lm$value[lm$level=="landscape" & lm$metric==metric & lm$year==2000],
               lm$value[lm$level=="class" & lm$metric==metric & lm$year==2000]))
  return(tmp)
  }

### Plot the area of patches
b <- barplot(log10(get_metric("area_mn")), beside=T, col=lzcols,
        names.arg=c("All", unique(lm$zone[!is.na(lm$zone)])), las=2,  dens=dens,
        ylab="Mean patch area (ha)", xlab="Life zone")
abline(v=3.5, lty=2)
axis(1, at=b, labels=NA)

### Plot the number of patches
barplot(log10(get_metric("np")), beside=T, col=lzcols,
        names.arg=c("All", unique(lm$zone[!is.na(lm$zone)])), las=2,  dens=dens,
        ylab="Number of patches", xlab="Life zone")
abline(v=3.5, lty=2)
axis(1, at=b, labels=NA)

### Plot the nearest neighbor distance
barplot(log10(get_metric("enn_mn")), beside=T, col=lzcols, 
        names.arg=c("All", unique(rcl$zone)), las=2, density=dens, 
        ylab="Mean euclidean NN dist. (km)", xlab="Life zone")
abline(v=3.5, lty=2)
axis(1, at=b, labels=NA)
legend("topright", legend=c(1951, 2000), 
       bty="n",  pt.cex=2, density=c(40,100), fill='grey', cex=1.2)

### Plot the aggregation index
barplot(get_metric("ai"), beside=T, col=lzcols, 
        names.arg=c("All", unique(rcl$zone)), las=2, density=dens, 
        ylab="Aggregation index", xlab="Life zone", ylim=c(0,100))
abline(v=3.5, lty=2)
axis(1, at=b, labels=NA)



### Plot the area of patches
plot(get_metric("np"), get_metric("area_mn"), pch=21, bg=lzcols,
     xlab="Mean patch area (ha)", ylab="Number of patches", log='xy')
arrows(get_metric("np")[1,], get_metric("area_mn")[1,], 
       get_metric("np")[2,], get_metric("area_mn")[2,], 
       col=lzcols[seq(1,14,2)], len=0.1, lwd=3)
points(get_metric("np"), get_metric("area_mn"), pch=21, bg=lzcols)
arrows(get_metric("np")[1,], get_metric("area_mn")[1,], 
       get_metric("np")[2,], get_metric("area_mn")[2,], 
       len=0.1)
legend("bottomleft", legend=c("All", unique(rcl$zone)),
       bty='n', pch=21, pt.bg=lzcols[seq(1,14,2)])



############################################
############# VALIDATION PLOTS #############

### Make a plot of total area in each land cover class during 1951 and 2000, at the full island and life zone levels
### For the aggregated version, this will be a % of the unaggregated version in each category.

change_cols <- c('black', viridis::viridis(4)[c(4,2,3)])

lz_changes <- function(d=NULL, foczone="df-S", prop=T){
  tmp <- d[d$zone %in% foczone,]
  tmp2 <- tmp$value[match(1:4, tmp$change)]
  tmp2[is.na(tmp2)] <- 0
  if(prop){return(100*(tmp2/sum(tmp2)))}
  if(!prop){return(tmp2)}
}

par(mfrow=c(2,2), mar=c(5,5,2,1))

b <- barplot(
  cbind(tapply(lsm1$value, lsm1$change, sum, na.rm=T),
        lz_changes(lsm1, foczone="df-S", prop=F),
        lz_changes(lsm1, foczone="mf-S", prop=F),
        lz_changes(lsm1, foczone="wf-S", prop=F),
        lz_changes(lsm1, foczone="wf-LM", prop=F),
        lz_changes(lsm1, foczone="rf-S", prop=F),
        lz_changes(lsm1, foczone="rf-LM", prop=F)), 
  beside=T, col=change_cols, xlab="Lifezone", ylab="Area (ha)")
axis(1, at=b, labels=NA)
axis(1, at=colMeans(b[2:3,]), labels=c("All", unique(lsm1$zone))[1:7], tick=F)
abline(v=5.5, lty=2)

b <- barplot(cbind(100*tapply(lsm1$value, lsm1$change, sum, na.rm=T)/
                     lsm1$value[lsm1$level=="landscape"],
        lz_changes(lsm1, foczone="df-S"),
        lz_changes(lsm1, foczone="mf-S"),
        lz_changes(lsm1, foczone="wf-S"),
        lz_changes(lsm1, foczone="wf-LM"),
        lz_changes(lsm1, foczone="rf-S"),
        lz_changes(lsm1, foczone="rf-LM")), 
  beside=T, col=change_cols, xlab="Lifezone", ylim=c(0,100), ylab="Area (%)")
axis(1, at=b, labels=NA)
axis(1, at=colMeans(b[2:3,]), labels=c("All", unique(lsm1$zone))[1:7], tick=F)
abline(v=5.5, lty=2)

plot.new()
plot.new()
legend("topleft", legend=c("Stable non-forest", 
                           "Loss", 
                           "Gain", 
                           "Stable forest"), bty="n",
       pt.bg=change_cols, pch=22, pt.cex=3, cex=1.5)


#############################################################################################
#############################################################################################
#############################################################################################

################################################
### AGGREGATED RESOLUTION METRICS
################################################

### Load predicted maps
Pred_stack_raw <- rast("Data/Derived/raw_stack-20251022.tif")

### Aggregate forest map to SDM resolution
f_agg <- terra::resample(f, Pred_stack_raw[[1]], method="near")
# f_agg <- terra::resample(f, Pred_stack_raw[[1]], method="mode")
# f_agg <- terra::resample(f, Pred_stack_raw[[1]], method="mean")>0.3
# f_agg <- terra::resample(f, Pred_stack_raw[[1]], method="mean")>0.7

### Get a forest change map
newforest_agg <- (f_agg[[2]] - f_agg[[1]])
stable_forest_agg <- (f_agg[[2]]==1 & f_agg[[1]]==1)
stable_nonforest_agg <- (f_agg[[2]]==0 & f_agg[[1]]==0)
fchange_agg <- ((stable_nonforest_agg*-1) + stable_forest_agg + f_agg[[2]])+2

writeRaster(fchange_agg, file="Data/Derived/Forest_change_resampled_near.tif")

plot(fchange_agg, col=c('black', viridis::viridis(4)[c(4,2,3)]))
plot(lz_lcc$geometry, border='white', add=T)


### Rasterize the life zone map
zone_rast <- rasterize(vect(lz_lcc), fchange_agg, field="ECOZONE")

### Combine forest change map and life zones; unique ID for each combination
combined_agg <- as.factor(fchange_agg + zone_rast*4)

### Subset change map to forest-only in 1951 and 2000
f51_mask <- mask(combined_agg, combined_agg, 
                 maskvalues=rcl$id[!rcl$for51])
f51_mask_noRF <- mask(combined_agg, combined_agg, 
                      maskvalues=rcl$id[!rcl$for51 | grepl("rf", rcl$zone)])
f00_mask <- mask(combined_agg, combined_agg, 
                 maskvalues=rcl$id[!rcl$for00])
f00_mask_noRF <- mask(combined_agg, combined_agg, 
                      maskvalues=rcl$id[!rcl$for00 | grepl("rf", rcl$zone)])

### Reclassify
f51_reclass <- classify(f51_mask, rcl[,c("id","zone_for51")])
f51_reclass_noRF <- classify(f51_mask_noRF, rcl[,c("id","zone_for51")])
f00_reclass <- classify(f00_mask, rcl[,c("id","zone_for00")])
f00_reclass_noRF <- classify(f00_mask_noRF, rcl[,c("id","zone_for00")])

### Metrics to calculate on all classes and full landscape (all classes) 
lsm_agg <- calculate_lsm(combined_agg, what=c("lsm_c_ca", "lsm_l_ta"))
lsm_agg <- cbind(lsm_agg, rcl[match(lsm_agg$class, rcl$id),])
lsm_agg <- lsm_agg[order(lsm_agg$zoneorder),]

### Metrics to calculate on all forest in 51 or all forest in 00
m51_agg <- calculate_lsm(f51_reclass, what=c("lsm_c_ca", "lsm_c_np", "lsm_c_area_mn", 
                                         "lsm_c_enn_mn", "lsm_c_ai", 
                                         "lsm_l_area_mn", "lsm_l_np", 
                                         "lsm_l_ai","lsm_l_para_mn","lsm_c_para_mn"))
lsm2_agg <- lsm_l_enn_mn(f51_reclass_noRF)
m51_agg <- rbind(m51_agg, lsm2_agg)

m00_agg <- calculate_lsm(f00_reclass, what=c("lsm_c_ca", "lsm_c_np", "lsm_c_area_mn", 
                                         "lsm_c_enn_mn", "lsm_c_ai", 
                                         "lsm_l_area_mn", "lsm_l_np", 
                                         "lsm_l_ai","lsm_l_para_mn","lsm_c_para_mn"))
lsm3_agg <- lsm_l_enn_mn(f00_reclass_noRF)
m00_agg <- rbind(m00_agg, lsm3_agg)

m51_agg$year <- 1951
m00_agg$year <- 2000

lm_agg <- rbind(m51_agg, m00_agg)

lm_agg$zone <- rcl$zone[match(lm_agg$class, rcl$zone_for51)]
lm_agg$zoneorder <- rcl$zoneorder[match(lm_agg$class, rcl$zone_for51)]

lm_agg <- lm_agg[order(lm_agg$zoneorder),]

############################################
############# DIAGNOSTIC PLOTS #############
### PLOT IT!
par(mfrow=c(2,2), mar=c(4,6,2,1))
lzcols <- rep(c("grey","palegoldenrod","palegreen","seagreen2",
                "skyblue","royalblue1","purple3"), each=2)
dens <- rep(c(40,100), times=7)


get_metric <- function(d, metric){
  tmp <- rbind(c(d$value[d$level=="landscape" & d$metric==metric & d$year==1951],
                 d$value[d$level=="class" & d$metric==metric & d$year==1951]),
               c(d$value[d$level=="landscape" & d$metric==metric & d$year==2000],
                 d$value[d$level=="class" & d$metric==metric & d$year==2000]))
  return(tmp)
}


### Density plot of environmental conditions in different LC categories
#total annual rainfall (R)
R <- rast(list.files("Data/ppt", full.names = TRUE)[13])/100

envtmp <- terra::project(R, fchange_agg)

ppt_newforest <- mask(mask(terra::project(R, fchange_agg), fchange_agg,
                           maskvalues=c(1,2,4)), !is.na(zone_rast), maskvalues=0)
ppt_lostforest <- mask(mask(terra::project(R, fchange_agg), fchange_agg, 
                              maskvalues=c(1,3,4)), !is.na(zone_rast), maskvalues=0)
ppt_stable <-  mask(mask(terra::project(R, fchange_agg), fchange_agg, 
                           maskvalues=c(1,2,3)), !is.na(zone_rast), maskvalues=0)

# Set up the plot area
x <- values(envtmp)[!is.na(values(envtmp))]
d <- stats::density(x, adjust=3)
d$y <- d$y * length(x) * (d$x[2] - d$x[1])

# Plot the frequency-scaled curve
# par(mar=c(6,6,2,2))
plot(d, type = "n", ylab = "Frequency", log='', main=NA, 
     xlab="Mean annual precipitation (mm)", axes=F, xlim=c(100,5000))
axis(1); axis(2)
polygon(d, col="grey", border=NA)
# abline(v=median(x), lwd=2, lty=2)

legend("topright", legend=c("All island", 
                            "Loss", 
                            "Gain", 
                            "Stable forest"), bty="n",
       pt.bg=c('grey', viridis::viridis(4)[c(4,2,3)]),
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

### Plot the area of patches
b <- barplot(get_metric(lm_agg, "area_mn"), beside=T, col=lzcols,
             names.arg=c("All", unique(lm_agg$zone[!is.na(lm_agg$zone)])), las=2,  dens=dens,
             ylab="Mean patch area (ha)", xlab="Life zone")
abline(v=3.5, lty=2)
axis(1, at=b, labels=NA)

### Plot the number of patches
barplot(get_metric(lm_agg, "np"), beside=T, col=lzcols,
        names.arg=c("All", unique(lm_agg$zone[!is.na(lm_agg$zone)])), las=2,  dens=dens,
        ylab="Number of patches", xlab="Life zone")
abline(v=3.5, lty=2)
axis(1, at=b, labels=NA)
legend("topright", legend=c(1951, 2000), 
       bty="n",  pt.cex=2, density=c(40,100), fill='grey', cex=1.2)

### Plot the nearest neighbor distance
barplot(get_metric(lm_agg, "enn_mn"), beside=T, col=lzcols, 
        names.arg=c("All", unique(lm_agg$zone[!is.na(lm_agg$zone)])), las=2,  dens=dens,
        ylab="Mean euclidean NN dist. (km)", xlab="Life zone")
abline(v=3.5, lty=2)
axis(1, at=b, labels=NA)

### Plot the aggregation index
# barplot(get_metric(lm_agg, "ai"), beside=T, col=lzcols, 
#         names.arg=c("All", unique(lm_agg$zone[!is.na(lm_agg$zone)])), las=2,  dens=dens,
#         ylab="Aggregation index", xlab="Life zone")
# abline(v=3.5, lty=2)
# axis(1, at=b, labels=NA)

### Plot the mean perimeter-area ratio
# barplot(get_metric(lm_agg, "para_mn"), beside=T, col=lzcols,
#         names.arg=c("All", unique(lm_agg$zone[!is.na(lm_agg$zone)])), las=2,  dens=dens,
#         ylab="Perimeter-area ratio", xlab="Life zone")
# abline(v=3.5, lty=2)
# axis(1, at=b, labels=NA)




### Plot the area of patches
# plot(get_metric(lm_agg, "np"), get_metric(lm_agg, "area_mn"), pch=21, bg=lzcols,
#      ylab="Mean patch area (ha)", xlab="Number of patches", log='xy')
# arrows(get_metric(lm_agg, "np")[1,], get_metric(lm_agg, "area_mn")[1,], 
#        get_metric(lm_agg, "np")[2,], get_metric(lm_agg, "area_mn")[2,], 
#        col=lzcols[seq(1,14,2)], len=0.1, lwd=3)
# points(get_metric(lm_agg, "np"), get_metric(lm_agg, "area_mn"), pch=21, bg=lzcols)
# arrows(get_metric(lm_agg, "np")[1,], get_metric(lm_agg, "area_mn")[1,], 
#        get_metric(lm_agg, "np")[2,], get_metric(lm_agg, "area_mn")[2,], 
#        len=0.1)
# legend("bottomleft", legend=c("All", unique(rcl$zone)),
#        bty='n', pch=21, pt.bg=lzcols[seq(1,14,2)])
# 


############################################
############# VALIDATION PLOTS #############

### Make a plot of total area in each land cover class during 1951 and 2000, at the full island and life zone levels
### For the aggregated version, this will be a % of the unaggregated version in each category.

change_cols <- c('black', viridis::viridis(4)[c(4,2,3)])

lz_changes <- function(d=lsm_agg, foczone="df-S", prop=T){
  tmp <- d[d$zone %in% foczone,]
  tmp2 <- tmp$value[match(1:4, tmp$change)]
  tmp2[is.na(tmp2)] <- 0
  if(prop){return(100*(tmp2/sum(tmp2)))}
  if(!prop){return(tmp2)}
}


par(mfrow=c(2,2), mar=c(5,5,2,1))

######### ABSOLUTE AREA OF LIFEZONE IN EACH CHANGE CLASS #########
b <- barplot(
  cbind(tapply(lsm_agg$value, lsm_agg$change, sum, na.rm=T),
        lz_changes(lsm_agg, foczone="df-S", prop=F),
        lz_changes(lsm_agg, foczone="mf-S", prop=F),
        lz_changes(lsm_agg, foczone="wf-S", prop=F),
        lz_changes(lsm_agg, foczone="wf-LM", prop=F),
        lz_changes(lsm_agg, foczone="rf-S", prop=F),
        lz_changes(lsm_agg, foczone="rf-LM", prop=F)), 
  beside=T, col=change_cols, xlab="Lifezone", ylab="Area (ha)")
axis(1, at=b, labels=NA)
axis(1, at=colMeans(b[2:3,]), labels=c("All", unique(lsm1$zone))[1:7], tick=F)
abline(v=5.5, lty=2)

legend("topright", legend=c("Stable non-forest", 
                           "Loss", 
                           "Gain", 
                           "Stable forest"), bty="n",
       pt.bg=change_cols, pch=22, pt.cex=2)

######### PROPORTION OF LIFEZONE IN EACH CHANGE CLASS #########
b <- barplot(cbind(100*tapply(lsm_agg$value, lsm_agg$change, sum, na.rm=T)/
                     lsm_agg$value[lsm_agg$level=="landscape"],
                   lz_changes(lsm_agg, foczone="df-S"),
                   lz_changes(lsm_agg, foczone="mf-S"),
                   lz_changes(lsm_agg, foczone="wf-S"),
                   lz_changes(lsm_agg, foczone="wf-LM"),
                   lz_changes(lsm_agg, foczone="rf-S"),
                   lz_changes(lsm_agg, foczone="rf-LM")), 
             beside=T, col=change_cols, xlab="Lifezone", ylim=c(0,100), ylab="Area (%)")
axis(1, at=b, labels=NA)
axis(1, at=colMeans(b[2:3,]), labels=c("All", unique(lsm1$zone))[1:7], tick=F)
abline(v=5.5, lty=2)

######### ABSOLUTE DIFFERENCE COMPARED TO FULL RESOLUTION VALUES #########
lz_changes_abs <- function(d=lsm_agg, foczone="df-S", fd=lsm1){
  tmp <- d[d$zone %in% foczone,]
  tmpf <- fd[fd$zone %in% foczone,]
  tmp2 <- tmp$value[match(1:4, tmp$change)]
  tmp2[is.na(tmp2)] <- 0
  tmpf2 <- tmpf$value[match(1:4, tmpf$change)]
  tmpf2[is.na(tmpf2)] <- 0
  return(tmp2-tmpf2)
}

b <- barplot(
  cbind(tapply(lsm_agg$value, lsm_agg$change, sum, na.rm=T) - 
          tapply(lsm1$value, lsm1$change, sum, na.rm=T),
        lz_changes_abs(lsm_agg, foczone="df-S", fd=lsm1),
        lz_changes_abs(lsm_agg, foczone="mf-S", fd=lsm1),
        lz_changes_abs(lsm_agg, foczone="wf-S", fd=lsm1),
        lz_changes_abs(lsm_agg, foczone="wf-LM", fd=lsm1),
        lz_changes_abs(lsm_agg, foczone="rf-S", fd=lsm1),
        lz_changes_abs(lsm_agg, foczone="rf-LM", fd=lsm1)), 
  beside=T, col=change_cols, xlab="Lifezone", ylab="Area (ha)")
axis(1, at=b, labels=NA)
axis(1, at=colMeans(b[2:3,]), labels=c("All", unique(lsm1$zone))[1:7], tick=F)
abline(v=5.5, lty=2)


######### PROPORTIONAL DIFFERENCE COMPARED TO FULL RESOLUTION VALUES #########
lz_changes_prop <- function(d=lsm_agg, foczone="df-S", fd=lsm1){
  tmp <- d[d$zone %in% foczone,]
  tmpf <- fd[fd$zone %in% foczone,]
  tmp2 <- tmp$value[match(1:4, tmp$change)]
  tmp2[is.na(tmp2)] <- 0
  tmpf2 <- tmpf$value[match(1:4, tmpf$change)]
  tmpf2[is.na(tmpf2)] <- 0
  return(100*(tmp2-tmpf2)/tmpf2)
}

b <- barplot(
  cbind((tapply(lsm_agg$value, lsm_agg$change, sum, na.rm=T)-tapply(lsm1$value, lsm1$change, sum, na.rm=T))/tapply(lsm1$value, lsm1$change, sum, na.rm=T),
        lz_changes_prop(lsm_agg, foczone="df-S", fd=lsm1),
        lz_changes_prop(lsm_agg, foczone="mf-S", fd=lsm1),
        lz_changes_prop(lsm_agg, foczone="wf-S", fd=lsm1),
        lz_changes_prop(lsm_agg, foczone="wf-LM", fd=lsm1),
        lz_changes_prop(lsm_agg, foczone="rf-S", fd=lsm1),
        lz_changes_prop(lsm_agg, foczone="rf-LM", fd=lsm1)), 
  beside=T, col=change_cols, xlab="Lifezone", ylab="% Diff from full resolution")
axis(1, at=b, labels=NA)
axis(1, at=colMeans(b[2:3,]), labels=c("All", unique(lsm1$zone))[1:7], tick=F)
abline(v=5.5, lty=2)



################################################
### AGGREGATED SPECIES LSM METRICS
################################################

### Load predicted maps
Pred_stack_thres <- rast("Data/Derived/thresholded_stack-20251022.tif")
fchange_agg <- rast("Data/Derived/Forest_change_resampled_near.tif")

### Combine forest change map and life zones; unique ID for each combination
combined_agg <- as.factor(fchange_agg + zone_rast*4)

### Subset change map to forest-only in 1951 and 2000
f51_mask <- mask(combined_agg, combined_agg, 
                 maskvalues=rcl$id[!rcl$for51])
f51_mask_noRF <- mask(combined_agg, combined_agg, 
                      maskvalues=rcl$id[!rcl$for51 | grepl("rf", rcl$zone)])
f00_mask <- mask(combined_agg, combined_agg, 
                 maskvalues=rcl$id[!rcl$for00])
f00_mask_noRF <- mask(combined_agg, combined_agg, 
                      maskvalues=rcl$id[!rcl$for00 | grepl("rf", rcl$zone)])

### Reclassify
f51_reclass <- classify(f51_mask, rcl[,c("id","zone_for51")])
f00_reclass <- classify(f00_mask, rcl[,c("id","zone_for00")])


predthresh <- mask(Pred_stack_thres, Pred_stack_thres, maskvalues=0)
predthresh_51 <- mask(predthresh, f51_reclass)
predthresh_00 <- mask(predthresh, f00_reclass)

### Metrics to calculate on all classes and full landscape (all classes) 
sp_lsm_51 <- calculate_lsm(predthresh_51, what=c("lsm_l_area_mn", "lsm_l_np", 
                                              "lsm_l_ai", "lsm_l_enn_mn", "lsm_l_mesh"))
sp_lsm_00 <- calculate_lsm(predthresh_00, what=c("lsm_l_area_mn", "lsm_l_np", 
                                                 "lsm_l_ai", "lsm_l_enn_mn", "lsm_l_mesh"))



sp_lsm_51$year <- 1951
sp_lsm_00$year <- 2000
sp_lm <- rbind(sp_lsm_51, sp_lsm_00)
sp_lm$sp <- names(predthresh_00)[sp_lm$layer]
sp_lm <- sp_lm[order(sp_lm$sp, sp_lm$metric, sp_lm$year),]

# Drop values for NN for species with fewer than 4 patches
sp_lm$value[sp_lm$metric=="enn_mn"][sp_lm$value[sp_lm$metric=="np"]<4] <- NA

saveRDS(sp_lm, "Data/Derived/Species_Landscape-metrics-near-20251024.RDA")





plot(Pred_stack_thres[["MICPAC"]])




