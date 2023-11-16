#exploratory analysis FAI Abundance data
library(quantreg)
library(bayesQR)
library(bayesplot)
library(rstan)
library(Brq)

#import data 
df <- read.csv("Data/Derived/master_final.csv")

#tranfdorm forest cover abbundance in km2 each grid cell of the raster was 450 m resolution 
#df$fcover_51 <- df$fcover_51*0.2025
#df$fcover_00 <- df$fcover_00*0.2025
#df$tot_change <-df$tot_change*0.2025

#transform Data scale
df$fcover_51.std<- as.vector(scale(df$fcover_51))
df$tot_change.std <- as.vector(scale(df$tot_change))
df$D_51.std <- as.vector(scale(df$D_51))
df$PC1.std <- as.vector(scale(df$PC1))

#look at varibles correlation
cor(df$fcover_51.std, df$tot_change.std)
cor(df$tot_change.std, df$D_51.std)
cor(df$fcover_51.std, df$D_51.std)

#plot variables distribution
hist(df$tpa_2014)
hist(df$fcover_51)
hist(df$fcover_51.std)
hist(df$tot_change)
hist(df$tot_change.std)
hist(df$D_51)
hist(df$D_51.std)
hist(df$PC1)
hist(df$PC1.std)

#Run first model considering the effect of fragmentation and habitat amount 
#bayesQR
#deinfine
prior <- prior(tpa_2014 ~ fcover_51.std + tot_change.std + D_51.std,data=df)
#fit the model 
mod1.1 <- bayesQR(tpa_2014 ~ fcover_51.std + tot_change.std + D_51.std, data=df, 
                   quantile= c(0.95), normal.approx=F, prior=prior ,ndraw = 100000)
res_mod1.1 <-summary(mod1.1, keep=200, burnin = 10000)

#check for model convergence 
#keep the estimated values into a matrix 
conv <- (mod1.1[[1]][["betadraw"]])
#remove the first 10000
conv<-conv[-c(1:10000),]
#check for the Rhat 
rstan::Rhat(conv) # --> 2.5 it seems is not converging well eventhough I have increased the iterations

#tryin the runnig the model with different Bayesis package Brq
mod1.2 <-Brq (tpa_2014 ~ fcover_51.std + tot_change.std + D_51.std, data=df, 
           tau=0.95, runs =100000, burn =10000, thin=5)
res_mod1.2 <- summary(mod1.2) #the results are bit different

#try to use the predict function
nd <- data.frame("tot_change.std"=seq(-2,3,length.out=272),
                 "fcover_51.std"=mean(df$fcover_51.std))

nd_51 <- data.frame("fcover_51.std"=seq(-2,3,length.out=272),
                    "tot_change.std"=mean(df$tot_change.std))
#ideally this how it wousls be but not working with the model output 
pred <- predict(mod1.2, newdata = nd, type= percentile)
pred <- predict(mod1.1, newdata = nd, type= percentile)

#plot-> I'm still haveing trubles to plot the trend line 

#run second model including traits
#remove data that don't have trait values 
df_t <- na.omit(df)
#define priors
prior_2 <- prior(tpa_2014 ~ fcover_51.std + tot_change.std + D_51.std + PC1.std  + tot_change.std * PC1.std , data=df_t)
#run the model 
mod2 <- bayesQR(tpa_2014 ~ fcover_51.std +
                              tot_change.std + 
                              PC1.std +
                              D_51.std +
                              tot_change.std * PC1.std , 
                              data=df_t, 
                              quantile= c(0.95) ,
                              normal.approx=F,
                              prior = prior_2,
                              ndraw = 100000)
#look at the results 
res_mod2 <- summary(mod2, burnin=5000)
#look at iterations convergence 
conv <- (mod2[[1]][["betadraw"]])
conv<-conv[-c(1:10000),]
#check for the Rhat 
rstan::Rhat(conv)

########## you can stop here 

#Run a quantile regression Abundace ~ forest cover 1951
mod1 <- rq(tpa_2014 ~ fcover_51.std + tot_change.std, data=df_f, tau = 0.95)
res.mod1 <- summary(mod1)

par(mfrow=c(1,2))

nd <- data.frame("tot_change.std"=seq(-2,3,length.out=155),
                 "fcover_51.std"=mean(df_f$fcover_51.std))
nd_51 <- data.frame("fcover_51.std"=seq(-2,3,length.out=155),
                 "tot_change.std"=mean(df_f$tot_change.std))
                 
pred <-predic(mod1, newdata = nd, type= percentile)
pred_51 <- predict.rq(mod1, newdata = nd_51, type= percentile)

plot(df_f$tot_change.std  , df_f$tpa_2014 , pch = 16,
     main = " Abundace ~ Suitable habitat change",
     xlab = " Suitable habitat change",
     ylab = " log Abundance (TPA)")
lines(nd$tot_change.std, pred)

plot(df_f$fcover_51.std, df_f$tpa_2014 , pch = 16,
    main = " Abundace ~ Suitable habitat 1951",
    xlab = " Suitable habitat 1951",
    ylab = " Abundance (TPA)")
lines(nd_51$fcover_51.std, pred_51)


# model 2 
mod2<- rq(tpa_2014~ 
            fcover_51.std + 
            tot_change.std + 
            tot_change.std * PC1.std , data=df_f, tau = 0.95)
res.mod2 <- summary(mod2)

nd_l <- data.frame("fcover_51.std"=mean(df_f$fcover_51.std),
                   "tot_change.std"=seq(-2,3,length.out=155),
                   "PC1.std"=-2.223796)

nd_h <- data.frame("fcover_51.std"=mean(df_f$fcover_51.std),
                  "tot_change.std"=seq(-2,3,length.out=155), 
                  "PC1.std"=7.783721)

pred_l <- predict.rq(mod2, newdata = nd_l)
pred_h <- predict.rq(mod2, newdata = nd_h)

data$Colour[data$col_name2>=3]="red"
data$Colour[data$col_name2<=1]="blue"

rbPal <- colorRampPalette(c('red','blue'))
dat$Col <- rbPal(10)[df_f$PC1.std,]

plot(df_f$tot_change.std, df_f$tpa_2014, pch = 16,
     main = " Abundace ~ Suitable habitat change ",
     xlab = " Suitable habitat change",
     ylab = " Species Abundance (TPA)", 
     col=df_f$PC1.std) 
lines(nd_l$tot_change.std, pred_l, col="chartreuse4")
lines(nd_h$tot_change.std, pred_h, col="chocolate1")
legend(-1.6, 82, bg="transparent", legend=c("min PC1", "max PC1"),
       col=c("chartreuse4", "chocolate1"), lty = 1, box.lty=0)

#########################################3
#trying Bayesian quantile regression 

#make a new data frame without NA values 
df <- data.frame(df_f$tpa_2014,df_f$fcover_51.std, df_f$tot_change.std)
colnames(df)[1] <- "tpa_2014"
colnames(df)[2] <- "fcover_51.std"
colnames(df)[3] <- "tot_change.std"
df<- na.omit(df)

#rund second model 
prior_2 <- prior(tpa_2014 ~ fcover_51.std + tot_change.std + PC1.std + tot_change.std * PC1.std , data=df_f)
mod2_bp <- bayesQR(tpa_2014 ~ fcover_51.std + tot_change.std + PC1.std + tot_change.std * PC1.std , data=df_f, 
                   quantile= c(0.95) ,prior=prior_2 ,ndraw = 100000, keep=5)
summary(mod2_bp, burnin=10000)
plot(mod2_bp, plottype = "trace", burnin=10000)

#adding distance between patches in 1951
prior_2 <- prior(tpa_2014 ~ fcover_51.std + tot_change.std + D_51.std, data=df_f)
mod3_bp <- bayesQR(tpa_2014 ~ fcover_51.std + D_51.std + tot_change.std , data=df_f, 
                   quantile= c(0.95) ,prior=prior ,ndraw = 20000)
summary(mod3_bp)




#Adding the species# 
# will all of the trait imputed data# 

#transform Data
df_uf$fcover_51.std<- as.vector(scale(df_uf$fcover_51))
df_uf$tot_change.std <- as.vector(scale(df_uf$tot_change))
df_uf$PC1.std <- as.vector(scale(df_uf$PC1))

cor(df_uf$fcover_51.std, df_uf$tot_change.std)
cor(df_uf$fcover_51.std, df_uf$PC1.std)
cor(df_uf$tot_change.std, df_uf$PC1.std)

#Run a quantile regression Abundace ~ forest cover 1951
prior_1 <- prior(tpa_2014 ~ fcover_51.std + tot_change.std, data=df_f)
mod1_BP <- bayesQR(tpa_2014 ~ fcover_51.std + tot_change.std, data=df_uf, 
                   quantile= c(0.95) ,prior=prior ,ndraw = 20000)
res.mod1 <- summary(mod1_BP)

prior_2 <- prior(tpa_2014 ~ fcover_51.std + tot_change.std + PC1.std + tot_change.std * PC1.std , data=df_f)
mod2_BP <- bayesQR(tpa_2014 ~ fcover_51.std + tot_change.std + PC1.std + tot_change.std * PC1.std , data=df_t, 
                   quantile= c(0.95),ndraw = 20000)

summary(mod2_BP)




















