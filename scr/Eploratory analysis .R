#exploratory analysis Forets Abbundance forets cover 

Abco<- read.csv("/Users/laumo791/Documents/PR/C1/Results/Abundace_Cover.csv", sep=";")

mod <- raster("/Users/laumo791/Documents/PR/C1/Results/SDM_threshold/trsh_stack.tif")
names(mod) <-

#tranfdorm forest cover abbundance in km2 each grid cell of the raster was 450 m resolution 
Abco$fcover_51 <- Abco$fcover_51*0.2025
Abco$fcover_00 <- Abco$fcover_00*0.2025
Abco$tot_change <-Abco$tot_change *0.2025

plot((Abco$fcover_51),(Abco$TPA), type="h")
plot((Abco$fcover_00),Abco$TPA, type="h")

#see if variables are correlated 
cor.test(Abco$fcover_00, Abco$TPA)
cor.test(Abco$tot_change, Abco$TPA)

#scale data
Abundance_Cover$fcover_00 <- scale(Abundance_Cover$fcover_00)
Abundance_Cover$fcover_51 <- scale(Abundance_Cover$fcover_51)
Abundance_Cover$TPA <- scale(Abundance_Cover$TPA)

fit_abbundance <- lm(Abco$TPA ~ Abco$fcover_51 + Abco$fcover_00)
plot(fit_abbundance)
summary(fit_abbundance)

#
plot(log10(Abco$fcover_51), log10(Abco$TPA),
     xlab = " Forest cover 1951",
     ylab = " Species abbundance")
abline(lm(fit_abbundance), col = "red")

plot(log(Abco$fcover_00), Abco$TPA,
    xlab = " Forest cover 200",
    ylab = " Species abbundance")
abline(lm(Abco$TPA ~ Abco$fcover_00), col = "red")


fit_abbundance_2<- lm(Abco$TPA ~Abco$tot_change)
plot(fit_abbundance_2)
summary(fit_abbundance_2)

plot(Abco$tot_change, Abco$TPA,
     xlab = " Forest change",
     ylab = " Species abbundance")
abline(lm(Abco$TPA ~ Abco$fcover_00), col = "red")

ggplot(fit_abbundance, aes(x=log(Abco$fcover_00), y=log(Abco$TPA)))+
  geom_point() +
  geom_smooth(method=lm, se=FALSE)

