#exploratory analysis FAI Abundance data
library(quantreg)

#import data 
AFF <- read.csv("Data/Derived/Abundance_Fcover_AI.csv", sep = ";")
master_data <-read.csv("Data/Derived/Master_data.csv")

##Plot different variables against abbundance
#tranfdorm forest cover abbundance in km2 each grid cell of the raster was 450 m resolution 
AFF$fcover_51 <- AFF$fcover_51*0.2025
AFF$fcover_00 <- AFF$fcover_00*0.2025
AFF$tot_change <-AFF$tot_change*0.2025

tot_AI <- AFF$AI_00 - AFF$AI_51

#make a correlation matrix of our variables 

#log transform Data
AFF$tpa_2004.z <- log(AFF$tpa_2004)
AFF$tpa_2014.z <- log(AFF$tpa_2014)
abundance <- AFF$tpa_2014.z - AFF$tpa_2004.z

#Run a quantile regression Abundace ~ forest cover 1951
fitAF_51<- rq( AFF$tpa_2014.z ~ AFF$fcover_51,  data=AFF, tau = 0.95)

model.null = rq(AFF$tpa_2014.z ~ 1,
                data = AFF,
                tau = 0.95)
#look at the P-value of the model 
anova(fitAF_51, model.null)

#plot the quantile regression Abundace ~ forest cover 1951
plot(AFF$tpa_2014.z ~ AFF$fcover_51, data=AFF, pch = 16,
     main = " Abundace ~ forest cover 1951",
     xlab = " Habitat amount 1951",
     ylab = " log Species Abundance (TPA)")
abline(rq((AFF$tpa_2014.z) ~ AFF$fcover_51, tau = 0.95, data=AFF), col = "Red")
abline(lm(AFF$tpa_2014.z  ~ AFF$fcover_51, data=AFF), col = "blue")
legend("topright", legend = c("rq","lm"), col = c("red", "blue"), lty = 2)

#Run a quantile regression Abundance ~ forest cover 2000
fitAF_00<- rq(AFF$tpa_2014.z ~ AFF$fcover_00, data=AFF, tau = 0.95)

model.null = rq(AFF$tpa_2014.z ~ 1,
                data = AFF,
                tau = 0.95)
#look at the P-value of the model 
anova(fitAF_00, model.null)


#plot the quadratic regression 
plot(AFF$tpa_2014.z~ AFF$fcover_00, data=AFF, pch = 16, 
     main = " Abundace ~ forest cover 2000",
     xlab = " Habitat amount 2000",
     ylab = "Species Abundance (TPA) 2014 ")
abline(rq(AFF$tpa_2014.z ~ AFF$fcover_00, tau = 0.95, data=AFF), col = "red")
abline(lm(AFF$tpa_2014.z  ~ AFF$fcover_00, data=AFF), col = "blue")
legend("bottomright", legend = c("rq","lm"), col = c("red", "blue"), lty = 2)

#quantile regression of Abundnace ~ variation in habbitat 
fitAF<- rq(AFF$tpa_2014.z ~ AFF$tot_change, data=AFF, tau = 0.95)

model.null = rq(AFF$tpa_2014.z ~ 1,
                data = AFF,
                tau = 0.95)
#look at the P-value of the model 
anova(fitAF, model.null)

#plot the quadratic regression 
plot((AFF$tpa_2014.z) ~ AFF$tot_change, data=AFF, pch = 16, 
     main = " Abundace 2014 ~  Habitat change -1951-2000 ",
     xlab = "  Habitat change -1951-2000",
     ylab = "Abundace 2014 ")
abline(rq(AFF$tpa_2014.z  ~ AFF$tot_change, tau = 0.95, data= AFF), col = "red")
abline(lm(AFF$tpa_2014.z  ~ AFF$tot_change, data= AFF), col = "blue")
legend("topright", legend = c("rq","lm"), col = c("red", "blue"), lty = 2)


#Run a quantile regression Abundance ~ Fragmentation 1951
fitAI_51<- rq(AFF$tpa_2014.z ~ AFF$AI_51, data=AFF, tau = 0.95)
anova(fitAI_51, model.null)

#plot the quadratic regression 
plot(AFF$tpa_2014.z ~ AFF$AI_51, data=AFF, pch = 16, 
     main = " Abundace ~ Connectivity 1951",
     xlab = " Habitat Connectivity 1951",
     ylab = "log Species Abundance (TPA) ")
abline(rq(AFF$tpa_2014.z ~ AFF$AI_51, tau = 0.95, data=AFF), col = "red")
abline(lm(AFF$tpa_2014.z  ~ AFF$AI_51, data= AFF), col = "blue")
legend("topright", legend = c("rq","lm"), col = c("red", "blue"), lty = 2)


#Run a quantile regression Abundance ~ Fragmentation 2000
fitAI_00<- rq(AFF$tpa_2014.z ~ AFF$AI_00, data=AFF, tau = 0.95)
anova(fitAI_00, model.null)

#plot the quadratic regression 
plot(AFF$tpa_2014.z ~ AFF$AI_00, data=AFF, pch = 16, 
     main = " Abundace ~ Connectivity 2000",
     xlab = " Habitat Connectivity 2000",
     ylab = "Species Abundance (TPA) 2014 ")
abline(rq(AFF$tpa_2014.z ~ AFF$AI_00, tau = 0.95, data=AFF), col = "red")
abline(rq(AFF$tpa_2014.z  ~ AFF$AI_00, tau = 0.5 ,data= AFF), col = "blue")
legend("topright", legend = c("rq","lm"), col = c("red", "blue"), lty = 2)

#scale data 
AFF$tpa_2014.z <- scale(AFF$tpa_2014.z)
AFF$AI_51.z <- scale(AFF$AI_51)
AFF$AI_00.z<- scale(AFF$AI_00)
AFF$fcover_51.z <- scale(AFF$fcover_51)
AFF$fcover_00.z <- scale(AFF$fcover_00)
tot_AI <- AFF$AI_00.z - AFF$AI_51.z
tot_F <- AFF$fcover_00.z - AFF$fcover_51.z

#interaction between abbundance and connectivity 
Fit_int_AIF51<- lm(AFF$tpa_2014.z ~ AFF$AI_51.z * AFF$fcover_51.z)
Fit_int_AIF00<- lm(AFF$tpa_2014.z ~ AFF$AI_00.z *AFF$fcover_00.z)

mod1 <- lm(abundance ~ tot_AI* tot_F * master_data$PC1)
summary(mod1)


#log transform and scale 
master_data$tpa_2014.z <- log10(master_data$tpa_2014)
master_data$tpa_2014.z <- scale(master_data$tpa_2014.z)
master_data$tpa_2004.z <- log10(master_data$tpa_2004)
master_data$tpa_2004.z <- scale(master_data$tpa_2004.z)
master_data$fcover_51.z <- scale(master_data$fcover_51)
master_data$fcover_00.z <- scale(master_data$fcover_00)
master_data$AI_51.z <- scale(master_data$AI_51)
master_data$AI_00.z <- scale(master_data$AI_00)


master_data$tot_abundance <- master_data$tpa_2004.z - master_data$tpa_2014.z
master_data$tot_F <- master_data$fcover_00.z - master_data$fcover_51.z
master_datatot_AI <- master_data$AI_00.z - master_data$AI_51.z


mod1 <- lm(master_data$tot_abundance ~ master_data$tot_F * master_datatot_AI * master_data$PC1)
summary(mod)

#abbundace and traits 
abundance <- master_data$tpa_2014.z 
habitat_ammount_1951 <- master_data$fcover_51.z 
habitat_ammount_2000 <- master_data$fcover_00.z 
connectivity_1951 <- master_data$AI_51.z 
connectivity_2000 <- master_data$AI_00.z 
PC1  <- master_data$PC1
PC2 <- master_data$PC2


mod51 <- lm(abundance ~ habitat_ammount_1951 * connectivity_1951 *PC1) 

mod00 <- lm(abundance ~ habitat_ammount_2000 * connectivity_2000 *PC1) 

mod51_2 <- lm(abundance ~ habitat_ammount_1951 * connectivity_1951 *PC2) 

mod51_2<- lm(abundance ~ habitat_ammount_2000 * connectivity_2000 *PC2) 

mod_all <-  lm(abundance ~ habitat_ammount_1951 * 
                      connectivity_1951 *
                      habitat_ammount_2000 * 
                      connectivity_2000 
                      *PC1) 
     





mod00 <- lm(master_data$tpa_2014.z ~ 
             master_data$fcover_00.z *
             master_data$AI_00.z *
             master_data$PC1)

mod51_2 <- lm(master_data$tpa_2014.z ~ 
              master_data$fcover_51.z *
              master_data$AI_51.z *
              master_data$PC2) 

mod00_2 <- lm(master_data$tpa_2014.z ~ 
              master_data$fcover_00.z *
              master_data$AI_00.z *
              master_data$PC2)






             
             
             

#try to look at the interactions 
mod4 <- lm(master_data$tpa_2014.z ~ 
             master_data$fcover_51.z +
             master_data$AI_51.z +
             master_data$fcover_51.z * master_data$PC1 +
             master_data$AI_51.z * master_data$PC1)









           
          
           
           
             
             
             


























