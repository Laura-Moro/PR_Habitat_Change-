library(rgdal)
library(raster)
library(ENMeval)
library(dplyr)
library(ecospat)
library(caret)


# LOAD ENVIRONMENTAL PREDICTOR VARIABLES derived in script 1 
envs <- stack(list.files("Data/Derived/envs", full.names = TRUE))

#LOAD OCCURENCES RECORD AND FILTER 
load("Data/FINAL_RECS.Rdata")

#without considering the islands -> filter out the island using the geological substrate map 
GEO_ext <-raster::extract(envs[[1]], full[,c('LONGDEC','LATDEC')])
oc <- dplyr::filter(full, !is.na(GEO_ext)) #these are the cleaned occurencese

# espress the geological substrate map as a factorial variable
envs[[1]] <- as.factor(envs[[1]])

#RUNNING THE MODELS univing ENMeval (maxnet)
for (sp in 1:length(unique(oc$CODE))){
  
  # This line to skip some species that dont work
  if(!sp %in% c(253,387,472,515,541)){
  
    #these lines are to to select the backgorund points for each species whic in our case are the all the othe occurence points 
    focsp <- unique(oc$CODE)[sp]
    occs <- oc[oc$CODE==focsp,2:1]
    bg <- oc[oc$CODE!=focsp,2:1]
    
    #to see on which species the model is working on! 
    message(paste('working on', focsp))
    print(Sys.time())
    
    if(nrow(occs)>10){
      
      if(nrow(occs)>=15){
        partition <- "checkerboard2" #if species have more than 15 occurence points use "checkerboard2" data partition
      } else {
        partition <- "jackknife" #otherwise use jackknife
      }
      
# HOW ARE 'DUPLICATE' RECORDS BEING TREATED? (AS IN, OCCS IN THE SAME GRID CELL?)
      mod <- ENMevaluate(occs=occs, envs=envs, bg=bg, 
                         algorithm='maxnet', 
                         partitions=partition,
                         categoricals="GEO",
                         partition.settings=list(aggregation.factor=c(5,5)),
                         tune.args=list(fc = c("L","LQ","LQH","H"), rm = 1:5))
      # tune.args=list(fc = c("L","H"), rm = 1:5))
      
      filename <- paste0("Results/2022-12-07_ENMeval_results/", focsp, ".RDA")
      saveRDS(mod, file=filename)
    }
  }
}

#MODEL SELECTION 

load 
#load mod files #mod = enmeval result object 
modfiles <- list.files("/Users/laumo791/Documents/PR/C1/Results/2022-12-07_ENMeval_results", full.names = TRUE)
#create empty lists 
res <- list() #this is for the results 

#select the models with the the lowest omission rate at 10p and the highest AUC
for(i in seq_along(modfiles)){
  print(i)
  
  #red in the model files mod
  mod <- readRDS(modfiles[i]) 

  # select the models ith the minimum omission rate and the maximum AUC   
  tmpres <- mod@results %>% 
    filter(or.10p.avg == min(or.10p.avg)) %>% 
    filter(auc.val.avg == max(auc.val.avg))
 
  #take the full names (file path) and extract the just the code of the species 
  species <- gsub(".RDA", "", strsplit(modfiles[i],"/")[[1]][9])
  
  tmpres$species <- species
  
  res[[i]] <- tmpres[1,]
  
  #take the predictions maps of the selected models
  pred <- mod@predictions[[which(names(mod@predictions)==res[[i]]$tune.args)]]
  names(pred) <- species
  
  # vector of 1/0 for occ/bg
  occ_bg <- c(rep(1, nrow(mod@occs)), rep(0, nrow(mod@bg)))
  
  # vector of predicted values at occ and bg points
  pred_vals <- c(raster::extract(pred, mod@occs[,1:2]),
                 raster::extract(pred, mod@bg[,1:2]))
  
  # find the threshold at max TSS
  tss_thresh <- ecospat.max.tss(pred_vals, occ_bg)$max.threshold
  
  # impose the max TSS threshold
  pred_thresh <- pred > tss_thresh
  names(pred_thresh) <- species 
  
  #Save the continuous and thresholded predictions (raster layers)! 
  raster::writeRaster(pred, paste0("/Users/laumo791/Documents/PR/C1/Results/SDM_predictions/", names(pred)), format='GTiff')
  raster::writeRaster(pred_thresh, paste0("/Users/laumo791/Documents/PR/C1/Results/SDM_threshold/", names(pred_thresh)), format='GTiff')
}


#tables of all of the results 
resall <- do.call(rbind, res)
#Save model outputs 
write.csv(resall, "/Users/laumo791/Documents/PR/C1/SDMs/results/mod_output.csv")








