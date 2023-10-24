#exploratory analysis FAI Abundance data
library(quantreg)
library(ggplot2)


#import data 
df_f <-read.csv("Data/Derived/master_f.csv")
df_uf <- read.csv("Data/Derived/master_uf.csv")

#tranfdorm forest cover abbundance in km2 each grid cell of the raster was 450 m resolution 
#df$fcover_51 <- df$fcover_51*0.2025
#df$fcover_00 <- df$fcover_00*0.2025
#df$tot_change <-df$tot_change*0.2025

#transform Data
df_f$tpa_2014.z <- log(df_f$tpa_2014)
df_f$fcover_51.std<- as.vector(scale(df_f$fcover_51))
df_f$tot_change.std <- as.vector(scale(df_f$tot_change))
df_f$D_51.std <- as.vector(scale(df_f$D_51))
df_f$PC1.std <- as.vector(scale(df_f$PC1))

cor(df_f$fcover_51.std, df_f$tot_change.std)
cor(df_f$tot_change.std, df_f$D_51.std)
cor(df_f$fcover_51.std, df_f$D_51.std)

hist(df_f$fcover_51)
hist(df_f$fcover_51.std)
hist(df_f$tot_change)
hist(df_f$tot_change.std)
hist(df_f$D_51)
hist(df_f$D_51.std)
hist(df_f$PC1)
hist(df_f$PC1.std)

#Run a quantile regression Abundace ~ forest cover 1951
mod1 <- rq(tpa_2014.z ~ fcover_51.std + tot_change.std, data=df_f, tau = 0.95)
res.mod1 <- summary(mod1, se= 'boot')

#par(mfrow=c(1,3))

nd <- data.frame("tot_change.std"=seq(-2,3,length.out=155),
                 #BOB: You can simply pass a single value instead of repeating like this
                 #"fcover_51.std"=rep(mean(df_f$fcover_51.std),155))
                 "fcover_51.std"=mean(df_f$fcover_51.std))
                 
pred <-predict.rq(mod1, newdata = nd, type= percentile)

plot(df_f$tot_change.std , df_f$tpa_2014.z , pch = 16,
    main = " Abundace ~ suitable habitat change",
    xlab = " Suitable habitat change",
    ylab = " log Abundance (TPA)")

# BOB: You cannot use 'abline' like this - you need to use 'lines'
# abline(df_f$tot_change.std, pred)
lines(nd$tot_change.std, pred)


# model 2 
# BOB: you had df_f$... in the variable names. So, then the variable names in the model are saved as "df_f$..." but you were naming the variables in your new data without the "df_f$" prefix.  Since you specify the data in the rq() call, you do not need to include this prefix.
mod2<- rq(tpa_2014.z~ 
            fcover_51.std + 
            tot_change.std + 
            tot_change.std * PC1 , data=df_f, tau = 0.95)
res.mod2 <- summary(mod2, se= 'boot')

nd_l <- data.frame("fcover_51.std"=mean(df_f$fcover_51.std),
                   "tot_change.std"=seq(-2,3,length.out=155),
                   "PC1"=-2.223796)

nd_h <- data.frame("fcover_51.std"=mean(df_f$fcover_51.std),
                  "tot_change.std"=seq(-2,3,length.out=155), 
                  "PC1"=7.783721)

pred_l <- predict.rq(mod2, newdata = nd_l)
pred_h <- predict.rq(mod2, newdata = nd_h)

plot(df_f$tot_change.std, df_f$tpa_2014.z, pch = 16,
     main = " Abundace ~ forest cover change",
     xlab = " Suitanle habitat change",
     ylab = " log Species Abundance (TPA)")

# BOB: Again here, you were using 'abline()' whereas you should do it like this
lines(nd_l$tot_change.std, pred_l)
lines(nd_h$tot_change.std, pred_h)


# BOB: I have not gone past this point...

mod5 <-  rq(df_f$tpa_2014.z ~  df_f$fcover_51.std + df_f$tot_change.std + df_f$D_change.std, data=df_f, tau = 0.95)
res.mod5 <- summary(mod, se= pred_l)

mod6 <-  rq(df_f$tpa_2014.z ~  
              df_f$fcover_51.std + 
              df_f$tot_change.std + 
              df_f$D_change.std +
              df_f$PC1+
              df_f$PC1 * df_f$tot_change.std +
              df_f$D_change.std * df_f$tot_change.std,
              data=df_f, tau = 0.95)
res.mod6 <- summary(mod6, se= 'boot')


# will all of the trait imputed data 

#transform Data
df_uf$tpa_2014.z <- log(df_uf$tpa_2014)
df_uf$tpa_2014.std <- scale(df_uf$tpa_2014.z)
df_uf$fcover_51.std<- scale(df_uf$fcover_51)
df_uf$tot_change.std <- scale(df_uf$tot_change)
df_uf$D_change <- (df_uf$D_00 - df_uf$D_51)
df_uf$D_change.std <- scale(df_uf$D_change)
df_uf$PC1.std <- scale(df_uf$PC1)

#Run a quantile regression Abundace ~ forest cover 1951
mod3 <- rq(df_uf$tpa_2014.z ~  df_uf$fcover_51.std + df_uf$tot_change.std , data=df_uf, tau = 0.95)
res.mod3 <- summary(mod3, se= 'boot')

mod <- rq(df_uf$tpa_2014.z ~  df_uf$fcover_51.std + df_uf$tot_change.std + df_uf$D_change.std, data=df_uf, tau = 0.95)


plot(df_uf$tpa_2014.z ~ df_uf$tot_change.std + df_uf$fcover_51.std, data=df_uf, pch = 16,
     main = " Abundace ~ suitable habitat change",
     xlab = " Suitable habitat change",
     ylab = " log Abundance (TPA)")
par(xpd=FALSE)
abline(rq(df_uf$tpa_2014.z ~ df_uf$tot_change.std, tau = 0.95, data=df_uf),col = "Red")
abline(rq(df_uf$tpa_2014.z  ~ df_uf$fcover_51.std,  tau = 0.95, data=df_uf), col = "blue")
legend("bottomright", legend = c("suitable habitat change","habitat in 1951"), col = c("red", "blue"), lty = 2)

mod4<- rq(df_uf$tpa_2014.z~ 
            df_uf$fcover_51.std + 
            df_uf$tot_change.std + 
            df_uf$PC1 +
            df_uf$tot_change.std*df_uf$PC1 , data=df_uf, tau = 0.95)
res.mod4 <- summary(mod4, se= 'boot')

plot(df_uf$tpa_2014 ~ df_uf$tot_change, data=df_uf, pch = 16,
     main = " Abundace ~ forest cover change",
     xlab = " Suitanle habitat change",
     ylab = " log Species Abundance (TPA)" , log='y')
par(xpd = "")
abline(rq(df_uf$tpa_2014~ df_uf$tot_change, tau = 0.95, data=df_uf, log='y'), col = "Red")
abline(rq(df_uf$tpa_2014  ~ df_uf$fcover_51,  tau = 0.95, data=df_uf), col = "blue")
abline(rq(df_uf$tpa_2014  ~ df_uf$PC1,  tau = 0.95, data=df_uf), col = "blue")

mod5 <- rq(df_uf$tpa_2014.z ~  df_uf$fcover_51.std + df_uf$tot_change.std + df_uf$D_change.std, data=df_uf, tau = 0.95)

mod6 <_ 















#assing names 
A <- df$tpa_2014.z
dF <- df$tot_change
F51 <- df$fcover_51
PC1 <- df$PC1
PC2 <- df$PC2

#only forest cover 
mod2 <- lm(A ~ dF * F51)
summary(mod2)

#interaction with the traits (PC1)
mod3 <- lm(A ~ F51+ dF+ PC1+ dF*PC1 + F51*PC1)
summary(mod3)
plot(dF,A) #,cex= PC1)
abline(lm(A ~ dF, data=df),col = "blue")
abline(lm(A ~ F51,data=df), col=" Red")
abline(lm(A ~ PC1,data=df))
legend("topright", legend = c("F51","df","PC1"), col = c("red", "blue", "black"), lty = 2)



#try with untrasformed data 
mod5 <- lm(df$tpa_2014~   df$tot_change + 
             df$fcover_51+ 
             df$PC1+
             df$tot_change * df$PC1 +
             df$fcover_51 * df$PC1 +
             df$tot_change * df$fcover_51)
summary(mod5)








            








           
          
           
           
             
             
             


























