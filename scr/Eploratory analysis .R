#exploratory analysis Forets Abbundance forets cover 

Abco<- read.csv("Data/Derived/Abundance_Fcover_AI.csv")

#tranfdorm forest cover abbundance in km2 each grid cell of the raster was 450 m resolution 
Abco$fcover_51 <- Abco$fcover_51*0.2025
Abco$fcover_00 <- Abco$fcover_00*0.2025
Abco$tot_change <-Abco$tot_change *0.2025

#plot forets cover at bigging and the end and speces abbundance tpa 
plot(Abco$fcover_51, Abco$fcover_00,cex=(Abco$tpa_2019)/6, 
     xlab="Forest Cover 1951 (km2)", ylab = "Forest Cover 2000 (km2)")
abline(0,1, col="Red")

plot(Abco$AI_51/100, Abco$AI_00/100, cex=(Abco$tpa_2019)/6, 
     xlab="Landscape aggreagtion index (AI) 2000", ylab = "Landscape aggreagtion index (AI) 1951")

plot(Abco$fcover_51, log10(Abco$tpa_2014)+2, pch=21, bg="darkGreen",
     xlab = " Habitat amount 1951",
     ylab = "Species Abbundance (TPA) 2014 ",
     cex.lab=2, cex.axis=1.8)

#variation in abbundnace 
abundance <- Abco$tpa_2014 - Abco$tpa_2004
plot(Abco$tot_change, log10(abundance)+2, pch=21, bg="darkGreen", 
     xlab = " Habitat amount change 1951-2000",
     ylab = "Species Abbundance change (TPA) 2004 - 2014 ", cex.lab=2, cex.axis=1.8)
axis(side=2, at= c(1, 2, 3, 4, 5), labels = c(0, 2, 3, 4 ), cex.axis=1.8)




#plot forest cover in 1951 against current appundnace 
plot(Abco$fcover_51, log10(Abco$tpa_2019)
plot(Abco$fcover_00, Abco$tpa_2019)
    
text(Abco$fcover_51,log10(Abco$tpa_2019), labels=round(Abco$SCIENTIFIC_NAME), cex=(Abco$tpa_2019)/6, pos=1, offset=0.5, cex=0.7)

#plot forest cover in 2000 against current appundnace 
plot(Abco$fcover_00,log10(Abco$tpa_2019), cex=(Abco$tpa_2019)/6)
abline(lm(Abco$tpa_2009~Abco$fcover_00))

plot(Abco$tot_change,Abco$tpa_2019, cex=(Abco$tpa_2019)/6)
text(Abco$tot_change, labes=Abco$CODE, pos=1, offset=0.5,cex=0.7)


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

