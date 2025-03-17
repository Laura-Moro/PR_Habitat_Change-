library(quantreg)
library(bayesQR)
# library(bayesplot)
# library(rstan)
# library(Brq)

## Import data 
df <- read.csv("Data/Derived/6b-output-20250314.csv")

## Load fragmentation data
lm <- read.csv("Data/Derived/Landscape-agg-metrics-20240828.csv")

## Select fragmentation metric(s) to add to analysis data
unique(lm$metric)
df$F <- lm$value[lm$year==1951 & lm$metric=='enn_mn']

## Transform forest cover in km2 each grid cell of the raster (was 450 m resolution)
df$fcover_51 <- df$fcover_51 * 0.2025
df$fcover_77 <- df$fcover_77 * 0.2025
df$fcover_91 <- df$fcover_91 * 0.2025
df$fcover_00 <- df$fcover_00 * 0.2025
df$tot_change <-df$tot_change * 0.2025
df$total_hab <-df$total_hab * 0.2025

# Change of forest from 1951-1977
df$old_change <- df$fcover_77 - df$fcover_51
df$old_change_raw <- df$fcover_77_raw - df$fcover_51_raw

## Transform data scale
df$tpa2014.z <- as.vector(scale(df$tpa2014))
df$tpa_gs_2014.z <- as.vector(scale(df$tpa_gs_2014))
df$tpa_diam2_2014.z <- as.vector(scale(df$tpa_diam2_2014))
df$tpa_seedlings_2014.z <- as.vector(scale(df$tpa_seedlings_2014))

df$fcover_51.z <- as.vector(scale(df$fcover_51))
df$fcover_51_raw.z <- as.vector(scale(df$fcover_51_raw))
df$tot_change.z <- as.vector(scale(df$tot_change))
df$tot_change_raw.z <- as.vector(scale(df$tot_change_raw))
df$total_hab.z <- as.vector(scale(df$total_hab))
df$total_hab_raw.z <- as.vector(scale(df$total_hab_raw))
df$F.z <- as.vector(scale(df$F))
df$old_change.z <- as.vector(scale(df$old_change))
df$old_change_raw.z <- as.vector(scale(df$old_change_raw))

## Drop outliers of fragmentation
df <- df[df$F.z < 5,]

## Look at correlations
cor(df$fcover_51.z, df$tot_change.z)
cor(df$fcover_51.z, df$tot_change.z)
cor(df$tot_change.z, df$F.z)
cor(df$fcover_51.z, df$F.z)
cor(df$fcover_51.z, df$total_hab.z)

## Plot variable distributions
hist(df$tpa2014)
hist(scale(log10(df$tpa2014.z)))
hist(df$fcover_51)
hist(df$fcover_51.z)
hist(df$tot_change)
hist(df$tot_change.z)
hist(df$F)
hist(df$F.z)
hist(df$pca1)


plot(df$old_change_raw.z, df$tpa_diam2_2014)

plot(tpa_diam2_2014 ~ old_change_raw.z, data=df)
m1 <- rq(tpa_diam2_2014 ~ old_change_raw.z + fcover_51_raw.z + F.z, data=df, tau=0.95)
# m1 <- rq(tpa_diam2_2014 ~ old_change_raw.z, data=df, tau=0.95)
summary(m1, se='boot')
abline(m1)



m1 <- rq(tpa_diam2_2014 ~ fcover_51.z + old_change.z + F.z, data=df, tau=0.95)
summary(m1, se='boot')

df2 <- df[df$nPlots2014 > 5 & df$fcover_51.z < 2,]

plot(tpa_gs_2014 ~ fcover_51.z, data=df2)
m1 <- rq(tpa_gs_2014 ~ fcover_51.z + old_change.z + F.z, data=df2, tau=0.95)
summary(m1, se='boot')
abline(m1)

plot(tpa_gs_2014 ~ fcover_51.z, data=df2)
m1 <- rq(tpa_gs_2014 ~ fcover_51.z, data=df2, tau=0.95)
summary(m1, se='boot')
abline(m1)

plot(tpa_gs_2014 ~ old_change.z, data=df2)
m1 <- rq(tpa_gs_2014 ~ old_change.z, data=df2, tau=0.95)
summary(m1, se='boot')
abline(m1)

plot(tpa_diam2_2014 ~ F.z, data=df2)
m1 <- rq(tpa_diam2_2014 ~ F.z, data=df2, tau=0.95)
summary(m1, se='boot')
abline(m1)



m1 <- rq(nPlots2014 ~ fcover_51.z + old_change.z + F.z, data=df, tau=0.95)
summary(m1, se='boot')

plot(df$fcover_51.z, df$nPlots2014, log='')


head(df$sp[rev(order(df$fcover_51.z))],10)


plot(df$fcover_51.z, df$tpa_diam2_2014, log='')
plot(df$fcover_51.z, df$tpa_diam2_2014, log='y')


plot(df$old_change.z, df$tpa_diam2_2009, log='')

m2 <- lm(log(tpa_gs_2009) ~ fcover_51.z + old_change.z + total_hab.z + F.z, data=df)
summary(m2)

car::vif(m2)


plot(df$tpa2019 ~ df$maxht.z, log='y')
mx <- lm(log(df$tpa2019) ~ df$maxht.z)
summary(mx)
abline(mx)

plot(log(tpa2014) ~ old_change.z, data=df)
plot(log(tpa2014) ~ F.z, data=df)
plot((tpa_diam2_2004-tpa_diam2_2014) ~ fcover_51.z, data=df)
plot((tpa_diam2_2004-tpa_diam2_2019) ~ tot_change.z, data=df)
abline(0.22, 0.39)


library(ENMeval)
x <- readRDS("Data/2022-12-07_ENMeval_results/POIFLO.RDA")
plot(x@predictions[[5]])
points(x@occs)


#Run first model considering the effect of fragmentation and habitat amount 
#bayesQR
#deinfine
prior <- prior(tpa2014 ~ fcover_51_raw.z + tot_change_raw.z + F.z, data=df)

#fit the model 
mod1.1 <- bayesQR(tpa2014 ~ fcover_51_raw.z + tot_change_raw.z + F.z, data=df, 
                  quantile= c(0.95), normal.approx=T, prior=prior, ndraw=50000)
(summary(mod1.1, keep=500, burnin=10000))

plot(tpa2014 ~ fcover_51_raw.z, data=df)
plot(tpa2014 ~ tot_change_raw.z, data=df)
abline(18.70, 6.35)
plot(tpa2014 ~ F.z, data=df)



# Run first model considering the effect of fragmentation and habitat amount 
mod0 <- rq(tpa2014.z ~ fcover_51.z, data=df, tau=0.95)
mod1 <- rq(tpa2014.z ~ fcover_51.z + tot_change.z + F.z, data=df, tau=0.95)
summary(mod0)
summary(mod1)

plot(df$tpa2014-df$tpa2004)
abline(h=0)

y <- df$tpa2019 - df$tpa2014

plot(df$fcover_00, y)


plot(tpa2014.z ~ fcover_51.z, data=df)
abline(mod1)

cor(df[,!colnames(df) %in% c("X","sp","species")], use='p')

mod1 <- lm(tpa2014.z ~ fcover_51.z + tot_change.z + F.z, data=df)
car::vif(mod1)


### bayesQR
# define prior
prior <- prior(tpa2014.z ~ fcover_51.z + tot_change.z + F.z, data=df)

# fit the model 
mod1.1 <- bayesQR(tpa2014.z ~ fcover_51.z + tot_change.z + F.z, data=df, 
                   quantile= c(0.95), prior=prior, ndraw=50000)
res_mod1.1 <- summary(mod1.1, keep=200, burnin=5000)
res_mod1.1

















#### BOB TEST GROUND
# prior <- prior(tpa_2014 ~ fcover_51.std + tot_change.std + D_51.std, data=df)

df$tpa2014.log.z <- scale(log10(df$tpa2014))

plot(df$fcover_51.z, df$tpa2014.z)

m0_2014 <- bayesQR(tpa2014.z ~ fcover_51.z + tot_change.z + F.z, data=df, 
                   quantile= c(0.95), normal.approx=T, ndraw=10000)

summary(m0_2014)
plot(df$fcover_51.z, df$tpa2014.z)

plot(tpa_2014 ~ fcover_51.std, data=df)
abline(12.856, 2.976)
abline(20.25, 7.18)

plot(m0_2014, plottype='hist', burnin=1000, keep=200)



m1 <- bayesQR(tpa_2014 ~ fcover_51.std + tot_change.std + D_51.std, data=df, 
                  quantile= c(0.95), normal.approx=T, ndraw=10000)

summary(m0a, keep=200, burnin=1000)
summary(m0b, keep=200, burnin=1000)


summary(m0, keep=200, burnin=1000)
summary(m1, keep=200, burnin=1000)

plot(m0, plottype='hist', burnin=1000, keep=200)
plot(m1, plottype='hist', burnin=1000, keep=200)

par(mfrow=c(2,2))
plot(m0, plottype='trace', burnin=1000, keep=200)
plot(m1, plottype='trace', burnin=1000, keep=200, col=2)



str(m0)


# Simulate data from heteroskedastic regression
set.seed(66)
n <- 200
X <- runif(n=n,min=0,max=10)
X <- X
y <- 1 + 2*X + rnorm(n=n, mean=0, sd=.6*X)

# Estimate series of quantile regressions with adaptive lasso
# NOTE: to limit execution time of the example, ndraw is set
#       to a very low value. Set value to 5000 for a better
#       approximation of the posterior distirubtion.
out <- bayesQR(y~X, quantile=c(.5), alasso=TRUE, ndraw=5000)

# Initiate plot
## Plot datapoints
plot(X, y, main="", cex=.6, xlab="X")




# Run a quantile regression Abundace ~ forest cover 1951

mod1 <- rq(tpa_2014 ~ fcover_51.std + tot_change.std , data=df, tau=0.95)
summary(mod1)

mod0 <- rq(tpa_2014 ~ 1, data=df, tau = 0.95)
summary(mod0) 

rho <- function(u, tau=.5)u*(tau - (u < 0))

V <- sum(rho(mod1$resid, mod1$tau))
V0 <- sum(rho(mod0$resid, mod0$tau))

R1 <- 1-V/V0
R1





m0 <- bayesQR(tpa_2014 ~ 1, data=df, 
              quantile= c(0.95), normal.approx=F, ndraw=10000)
m1 <- bayesQR(tpa_2014 ~ fcover_51.std + tot_change.std + D_51.std, data=df, 
              quantile= c(0.95), normal.approx=F, ndraw=10000)

summary(m0)
summary(m1)



library(brms)
fit <- brm(bf(tpa_2014 ~ fcover_51.std + tot_change.std + D_51.std, quantile = 0.95), 
           data = df, family = asym_laplace())
summary(fit)

predict(fit)


plot(tpa_2014 ~ fcover_51.std, data=df, log='y')
abline(20.95, 10.14)


fit_0 <- brm(bf(tpa_2014 ~ 1, quantile = 0.95), data = df[!is.na(df$PC1.std),], family = asym_laplace())

df$tpa_2014.z <- scale(df$tpa_2014, scale=F)

fit_1 <- brm(bf(tpa_2014.z ~ fcover_51.std + tot_change.std + D_51.std, quantile = 0.95), 
           data = df, family = asym_laplace())

plot(df$fcover_51, df$tpa_2014.z)

pp_check(fit_1)


hist(do.call(c, lapply(as_draws(fit_1), function(x) x$b_fcover_51.std)))



pred0 <- matrix(ncol=nrow(df), nrow=4*length(as_draws(fit_1)[[1]]$b_Intercept))
for(i in 1:nrow(df)){
  pred0[,i] <- do.call(c, lapply(as_draws(fit_0), function(x) x$Intercept))
}

pred1 <- matrix(ncol=nrow(df), nrow=4*length(as_draws(fit_1)[[1]]$b_Intercept))
for(i in 1:nrow(df)){
  pred1[,i] <- do.call(c, lapply(as_draws(fit_1), function(x) x$Intercept)) + 
  do.call(c, lapply(as_draws(fit_1), function(x) x$b_fcover_51.std)) * df$fcover_51.std[i] + 
  do.call(c, lapply(as_draws(fit_1), function(x) x$b_tot_change.std)) * df$tot_change.std[i] + 
  do.call(c, lapply(as_draws(fit_1), function(x) x$b_D_51.std)) * df$D_51.std[i]
}

plot(df$tpa_2014, apply(pred0, 2, median))
plot(df$tpa_2014, apply(pred1, 2, median))

resid0 <- abs(df$tpa_2014 - apply(pred0, 2, median))
resid1 <- abs(df$tpa_2014 - apply(pred1, 2, median))

rho <- function(u, tau=.5)u*(tau - (u < 0))

V <- sum(rho(resid1[!is.na(resid1)], 0.95))
V0 <- sum(rho(resid0[!is.na(resid0)], 0.95))

R1 <- 1-V/V0
R1






# tryin the runnig the model with different Bayesis package Brq
mod1.2 <-Brq(tpa_2014 ~ fcover_51.std + tot_change.std + D_51.std, data=df, 
              tau=0.95, runs =100000, burn =10000, thin=5)
summary(mod1.2) #the results are bit different







# Simulate data from heteroskedastic regression
set.seed(66)
n <- 200
X <- runif(n=n,min=0,max=10)
X <- X
y <- 1 + 2*X + rnorm(n=n, mean=0, sd=.6*X)

# Estimate series of quantile regressions with adaptive lasso
# NOTE: to limit execution time of the example, ndraw is set
#       to a very low value. Set value to 5000 for a better
#       approximation of the posterior distirubtion.
out <- bayesQR(y~X, quantile=c(.05,.25,.5,.75,.95), alasso=TRUE, ndraw=500)

# Analyze 5 quantiles using default prior
# NOTE: to limit execution time of the example, ndraw is set
#       to a very low value. Set value to 5000 for a better
#       approximation of the posterior distirubtion.
out1 <- bayesQR(y ~ X, quantile=c(.95), ndraw=5000, normal.approx = T)
out2 <- bayesQR(y ~ X, quantile=c(.95), ndraw=5000, normal.approx = F)

# Check traceplot of first variable of .75 quantile regression 
plot(out2, var=1, plottype="trace")

# Check posterior histogram of first variable of .5 quantile regression 
plot(out1, var=2, plottype="hist")
plot(out2, var=2, plottype="hist")

# Create default quantile plot of first variable
plot(out, var=1, plottype="quantile")

# Create quantile plot of second variable with 90% credible interval
plot(out, var="X", credint=c(.05, .95), plottype="quantile", main="This is an example")


plot(X, y)
abline()



#### END BOB TEST GROUND


#check for model convergence 
#keep the estimated values into a matrix 
conv <- (mod1.1[[1]][["betadraw"]])
#remove the first 10000
conv<-conv[-c(1:10000),]
#check for the Rhat 
rstan::Rhat(conv) # --> 2.5 it seems is not converging well eventhough I have increased the iterations




#tryin the runnig the model with different Bayesis package Brq
mod1.2 <-Brq(tpa_2014 ~ fcover_51.std + tot_change.std + D_51.std, data=df, 
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
mod1 <- rq(tpa_2014 ~ fcover_51.std + tot_change.std + D_51.std, data=df, tau=0.95)
summary(mod1)


mod1 <- rq(tpa_2014 ~ fcover_51.std + tot_change.std + D_51.std, data=df, tau=0.95)
summary(mod1)



df$y <- rowSums(cbind(df$tpa_2014, df$tpa_2019))

plot(df$fcover_51.std, df$tpa_2014)


plot(df$fcover_51.std, df$tpa_2014)


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

df_f <- df[!is.na(df$PC1),]
mod2<- rq(tpa_2014 ~ 
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
                   quantile= c(0.95), ndraw = 100000, keep=5)
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





plot(df$pca1, df$nPlots2014)
plot(df$fcover_51, df$nPlots2019, log='xy')

cor.test(df$pca1, df$nPlots2014)





m1 <- rq(nPlots2014 ~ fcover_51.z + old_change.z + F.z, data=df, tau=0.95)
summary(m1, se='boot')

plot(df$fcover_51.z, df$nPlots2014)
plot(df$old_change.z, df$nPlots2014)
plot(df$F.z, df$nPlots2014)
plot(df$wd.z, df$nPlots2014)

lm1 <- lm(nPlots2014 ~ fcover_51.z + old_change.z + F.z, data=df)
summary(lm1)

plot(nPlots2014 ~ fcover_51.z, data=df, log='y')
abline(lm1)

plot(nPlots2014 ~ old_change.z, data=df)
abline(lm1)


jtools::summ(m1, vifs=T, confint=T)
car::vif(lm1)







