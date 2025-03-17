# Explore SDM results

library(ENMeval)

x <- readRDS(list.files("Data/2022-12-07_ENMeval_results", full.names=T)[3])

x
x@occs

plot(x@predictions[[20]])
points(x@occs)

evalplot.stats(x, stats=c("cbi.train"), x.var="rm", color.var="fc")
evalplot.stats(x, stats=c("or.10p"), x.var="rm", color.var="fc")
evalplot.stats(x, stats=c("or.10p"), x.var="rm", color.var="fc")


# select the models ith the minimum omission rate and the maximum AUC   
tmpres <- x@results %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))

x@results %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))

x@results %>% filter(cbi.train == max(cbi.train, na.rm=T))

which(x@results$cbi.train == min(x@results$cbi.train, na.rm=T))

which(x@results[x@results$or.10p.avg == min(x@results$or.10p.avg, na.rm=T),]
      [x@results$auc.val.avg == min(x@results$auc.val.avg, na.rm=T),])

par(mfrow=c(2,1))
plot(x@predictions[[1]])
plot(x@predictions[[4]])

par(mfrow=c(1,1))
plot(raster::extract(x@predictions[[1]], x@occs[,c(1,2)]),
     raster::extract(x@predictions[[4]], x@occs[,c(1,2)]))
abline(0,1)


pred <- x@predictions[[4]]

# vector of 1/0 for occ/bg
occ_bg <- c(rep(1, nrow(x@occs)), rep(0, nrow(x@bg)))

# vector of predicted values at occ and bg points
pred_vals <- c(raster::extract(pred, x@occs[,1:2]),
               raster::extract(pred, x@bg[,1:2]))

# find the threshold at max TSS
tss_thresh <- ecospat.max.tss(pred_vals, occ_bg)$max.threshold

# impose the max TSS threshold
pred_thresh <- pred > tss_thresh

par(mfrow=c(2,1))
plot(pred)
plot(pred_thresh)
points(x@occs)
