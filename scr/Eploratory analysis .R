#exploratory analysis FAI Abundance data
library(quantreg)

#import data 
Abco<- read.csv("Data/Derived/Abundance_Fcover_AI.csv", sep = ";")

##Plot different variables against abbundance
#tranfdorm forest cover abbundance in km2 each grid cell of the raster was 450 m resolution 
Abco$fcover_51 <- Abco$fcover_51*0.2025
Abco$fcover_00 <- Abco$fcover_00*0.2025
Abco$tot_change <-Abco$tot_change *0.2025

#Run a quantile regression Abundace ~ forest cover 1951
fitAF_51<- rq(log10(Abco$tpa_2014) ~ Abco$fcover_51, data=Abco, tau = 0.95)
#plot the quadratic regression 
plot(log10(Abco$tpa_2014) +2 ~ Abco$fcover_51, data=Abco, tau = 0.95, pch = 16, 
     main = " Abundace ~ forest cover 1951",
     xlab = " Habitat amount 1951",
     ylab = "Species Abundance (TPA) 2014 ")
abline(rq(log10(Abco$tpa_2014)+2  ~ Abco$fcover_51,tau = 0.95, data=Abco), col = "red")
abline(lm(log10(Abco$tpa_2014)+2  ~ Abco$fcover_51, tau = 0.5, data=Abco), col = "blue")
abline(rq(log10(Abco$tpa_2014) +2 ~ Abco$fcover_51,tau = 0.05, data=Abco), col = "red")
legend("topright", legend = c("rq","lm","rq"), col = c("red", "blue",  ), lty = 2)

#Run a quantile regression Abundance ~ forest cover 1951
fitAF_00<- rq(log10(Abco$tpa_2014) ~ Abco$fcover_00, data=Abco, tau = 0.95)
#plot the quadratic regression 
plot(log10(Abco$tpa_2014) +2 ~ Abco$fcover_00, data=Abco, tau = 0.95, pch = 16, 
     main = " Abundace ~ forest cover 2000",
     xlab = " Habitat amount 2000",
     ylab = "Species Abundance (TPA) 2014 ")
abline(rq(log10(Abco$tpa_2014)+2  ~ Abco$fcover_00,tau = 0.95, data=Abco), col = "red")
abline(lm(log10(Abco$tpa_2014)+2  ~ Abco$fcover_00, tau = 0.5, data=Abco), col = "blue")
abline(rq(log10(Abco$tpa_2014) +2 ~ Abco$fcover_00,tau = 0.05, data=Abco), col = "red")
legend("topright", legend = c("rq","lm"), col = c("red", "blue"), lty = 2)

#quantile regression of variation abbundnace ~ variation in habbitat 
#variation in abbundnace 
abundance <- Abco$tpa_2014 - Abco$tpa_2004
fitAF<- rq(log10(abundance) ~ Abco$tot_change, data=Abco, tau = 0.95)
#plot the quadratic regression 
plot(log10(abundance) +2 ~ Abco$tot_change, data=Abco, tau = 0.95, pch = 16, 
     main = " Abundace chage 2004-2014 ~  Habitat change -1951-2000 ",
     xlab = "  Habitat change -1951-2000",
     ylab = "Abundace chage 2004-2014 ")
abline(rq(log10(abundance)+2  ~ Abco$tot_change,tau = 0.95, data=Abco), col = "red")
abline(lm(log10(abundance)+2  ~ Abco$tot_change, tau = 0.5, data=Abco), col = "blue")
abline(rq(log10(abundance) +2 ~ Abco$tot_change ,tau = 0.05, data=Abco),col = "red")
legend("topright", legend = c("rq","lm"), col = c("red", "blue"), lty = 2)


#Run a quantile regression Abundance ~ Fragmentation 1951
fitAI_00<- rq(log10(Abco$tpa_2014) ~ Abco$AI_51, data=Abco, tau = 0.95)
#plot the quadratic regression 
plot(log10(Abco$tpa_2014) +2 ~ Abco$AI_51, data=Abco, tau = 0.95, pch = 16, 
     main = " Abundace ~ fragmentation 1951",
     xlab = " Habitat fragmentation 1951",
     ylab = "Species Abundance (TPA) 2014 ")
abline(rq(log10(Abco$tpa_2014)+2  ~ Abco$AI_51,tau = 0.95, data=Abco), col = "red")
abline(lm(log10(Abco$tpa_2014)+2  ~ Abco$AI_51, tau = 0.5, data=Abco), col = "blue")
abline(rq(log10(Abco$tpa_2014) +2 ~ Abco$AI_51,tau = 0.05, data=Abco), col = "red")
legend("topright", legend = c("rq","lm"), col = c("red", "blue"), lty = 2)


#Run a quantile regression Abundance ~ Fragmentation 2000
fitAI_00<- rq(log10(Abco$tpa_2014) ~ Abco$AI_00, data=Abco, tau = 0.95)
#plot the quadratic regression 
plot(log10(Abco$tpa_2014) +2 ~ Abco$AI_00, data=Abco, tau = 0.95, pch = 16, 
     main = " Abundace ~ fragmentation 2000",
     xlab = " Habitat fragmentation 2000",
     ylab = "Species Abundance (TPA) 2014 ")
abline(rq(log10(Abco$tpa_2014)+2  ~ Abco$AI_00,tau = 0.95, data=Abco), col = "red")
abline(lm(log10(Abco$tpa_2014)+2  ~ Abco$AI_00, tau = 0.5, data=Abco), col = "blue")
abline(rq(log10(Abco$tpa_2014) +2 ~ Abco$AI_00,tau = 0.05, data=Abco), col = "red")
legend("topright", legend = c("rq","lm"), col = c("red", "blue"), lty = 2)


####SEE VARIATION IN ABBUNDANCE 

#make a matrix of the to see general trends in species abundance change 
A_mat <- Abco[, c(4,6,8,9)]
matplot(t(A_mat), type='l', lty=1)

#plot abbundance at the start and the and of the time period 
plot(log10(A_mat$tpa_2004)+1.5, log10(A_mat$tpa_2014)+2)
abline(lm(log10(A_mat$tpa_2014)~ log10(A_mat$tpa_2004)))

#simple regression of abbundnace change 
#make a matrix with only the abundnace data
A_mat <- Abco[, c(3,4,6,8,9)]
#change column names 
colnames(A_mat) <- c('CODE','2004','2009','2014','2019')
#transform data into long format 
A_mat <- A_mat %>%
  gather(key = "time", value = "TPA", "2004", "2009", "2014", "2019")
A_mat$CODE <- as.factor(A_mat$CODE)
A_mat$time <- as.numeric(A_mat$time)
str(A_mat)
# fit the linear t model 
lm_A_mat <- lm(A_mat$TPA ~ A_mat$time, data=A_mat)
summary(lm_A_mat)

lm_A_Sp <- lm(A_mat$TPA~ A_mat$CODE, data=A_mat)
summary(lm_A_Sp)
#save the model coefficents and intercept for each species 
res_lm_A_Sp<-res_lm_A_Sp[["coefficients"]]
coeffiecents <- as.data.frame(res_lm_A_Sp)


#repeated measure anova 
# Gather columns of TPA into long format
# Convert id and time into factor variables
#make a spreadsheet with only the TPA values 
A_mat <- Abco[, c(3,4,6,8,9)]
# Gather columns of TPA into long format i have exluded some of the years
#for too few data mabe aggregate the data? 
# Convert id and time into factor variables

#sumamry statistics 
A_mat  %>%
  group_by(time) %>%
  get_summary_stats(TPA, type = "mean_sd")
#visualize data --> there is not realy a fluctiation between years 
boxplot(log(TPA) ~time, A_mat)

#check for outliers 
A_mat %>%
  group_by(time) %>%
  identify_outliers(TPA)
#check normality assumptions 
A_mat %>%
  group_by(time) %>%
  shapiro_test(TPA)
#runn the repeated measure anova 
res.aov <- anova_test(data = A_mat, dv = TPA, wid = CODE, within = time)
get_anova_table(res.aov)

##Plot different variables against abbundance
#tranfdorm forest cover abbundance in km2 each grid cell of the raster was 450 m resolution 
Abco$fcover_51 <- Abco$fcover_51*0.2025
Abco$fcover_00 <- Abco$fcover_00*0.2025
Abco$tot_change <-Abco$tot_change *0.2025

#Run a quantile regression Abundace ~ forest cover 1951
fitAF_51<- rq(log10(Abco$tpa_2014) ~ Abco$fcover_51, data=Abco, tau = 0.95)
#plot the quadratic regression 
plot(log10(Abco$tpa_2014) +2 ~ Abco$fcover_51, data=Abco, tau = 0.95, pch = 16, 
     main = " Abundace ~ forest cover 1951",
     xlab = " Habitat amount 1951",
     ylab = "Species Abundance (TPA) 2014 ")
abline(rq(log10(Abco$tpa_2014)+2  ~ Abco$fcover_51,tau = 0.95, data=Abco), col = "red")
abline(lm(log10(Abco$tpa_2014)+2  ~ Abco$fcover_51, tau = 0.5, data=Abco), col = "blue")
abline(rq(log10(Abco$tpa_2014) +2 ~ Abco$fcover_51,tau = 0.05, data=Abco), col = "red")
legend("topright", legend = c("rq","lm","rq"), col = c("red", "blue", ), lty = 2)

#Run a quantile regression Abundance ~ forest cover 1951
fitAF_00<- rq(log10(Abco$tpa_2014) ~ Abco$fcover_00, data=Abco, tau = 0.95)
#plot the quadratic regression 
plot(log10(Abco$tpa_2014) +2 ~ Abco$fcover_00, data=Abco, tau = 0.95, pch = 16, 
     main = " Abundace ~ forest cover 2000",
     xlab = " Habitat amount 2000",
     ylab = "Species Abundance (TPA) 2014 ")
abline(rq(log10(Abco$tpa_2014)+2  ~ Abco$fcover_00,tau = 0.95, data=Abco), col = "red")
abline(lm(log10(Abco$tpa_2014)+2  ~ Abco$fcover_00, tau = 0.5, data=Abco), col = "blue")
abline(rq(log10(Abco$tpa_2014) +2 ~ Abco$fcover_00,tau = 0.05, data=Abco), col = "red")
legend("topright", legend = c("rq","lm"), col = c("red", "blue"), lty = 2)

#quantile regression of variation abbundnace ~ variation in habbitat 
#variation in abbundnace 
abundance <- Abco$tpa_2014 - Abco$tpa_2004
fitAF<- rq(log10(abundance) ~ Abco$tot_change, data=Abco, tau = 0.95)
#plot the quadratic regression 
plot(log10(abundance) +2 ~ Abco$tot_change, data=Abco, tau = 0.95, pch = 16, 
     main = " Abundace chage 2004-2014 ~  Habitat change -1951-2000 ",
     xlab = "  Habitat change -1951-2000",
     ylab = "Abundace chage 2004-2014 ")
abline(rq(log10(abundance)+2  ~ Abco$tot_change,tau = 0.95, data=Abco), col = "red")
abline(lm(log10(abundance)+2  ~ Abco$tot_change, tau = 0.5, data=Abco), col = "blue")
abline(rq(log10(abundance) +2 ~ Abco$tot_change ,tau = 0.05, data=Abco), col = "red")
legend("topright", legend = c("rq","lm"), col = c("red", "blue"), lty = 2)


#Run a quantile regression Abundance ~ Fragmentation 1951
fitAI_00<- rq(log10(Abco$tpa_2014) ~ Abco$AI_51, data=Abco, tau = 0.95)
#plot the quadratic regression 
plot(log10(Abco$tpa_2014) +2 ~ Abco$AI_51, data=Abco, tau = 0.95, pch = 16, 
     main = " Abundace ~ fragmentation 1951",
     xlab = " Habitat fragmentation 1951",
     ylab = "Species Abundance (TPA) 2014 ")
abline(rq(log10(Abco$tpa_2014)+2  ~ Abco$AI_51,tau = 0.95, data=Abco), col = "red")
abline(lm(log10(Abco$tpa_2014)+2  ~ Abco$AI_51, tau = 0.5, data=Abco), col = "blue")
abline(rq(log10(Abco$tpa_2014) +2 ~ Abco$AI_51,tau = 0.05, data=Abco), col = "red")
legend("topright", legend = c("rq","lm"), col = c("red", "blue"), lty = 2)


#Run a quantile regression Abundance ~ Fragmentation 2000
fitAI_00<- rq(log10(Abco$tpa_2014) ~ Abco$AI_00, data=Abco, tau = 0.95)
#plot the quadratic regression 
plot(log10(Abco$tpa_2014) +2 ~ Abco$AI_00, data=Abco, tau = 0.95, pch = 16, 
     main = " Abundace ~ fragmentation 2000",
     xlab = " Habitat fragmentation 2000",
     ylab = "Species Abundance (TPA) 2014 ")
abline(rq(log10(Abco$tpa_2014)+2  ~ Abco$AI_00,tau = 0.95, data=Abco), col = "red")
abline(lm(log10(Abco$tpa_2014)+2  ~ Abco$AI_00, tau = 0.5, data=Abco), col = "blue")
abline(rq(log10(Abco$tpa_2014) +2 ~ Abco$AI_00,tau = 0.05, data=Abco),col = "red")
legend("topright", legend = c("rq","lm"), col = c("red", "blue"), lty = 2)














