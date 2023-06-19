##PREPARING Trait Data 
library(dplyr)
library(readxl)

#1 SPecies Traits 
#compiling the Traint spreas sheets 

traits_E <- read.csv("Data/Traits/Species_Traits.csv")
guanica <- read.csv("Data/Traits/Guanica_seed_mass.csv")
Traits_fix<- read.csv("Data/Traits/Fixed_traits.csv", sep = ";" )
master_data <- read.csv("Data/Derived/Abundance_Fcover_AI.csv")

#take the fixed traits data set and filter the species we have abundance data for 
#it is only 168? 
Traits_fix<- filter(Traits_fix, Traits_fix$Code %in% 
                      master_data$CODE)

#add species names to the traits data set 
Traits_fix$species<- master_data$SCIENTIFIC_NAME[match(Traits_fix$Code, master_data$CODE)]

#take the data on the seed mass form the Helmer data set (extact only infomation on seed mass form the df)
Traits_fix$seed_mass <- traits_E$Seed_wt_avg_Kew_or_other_g_per_1000[match(Traits_fix$species,traits_E$PLANTS_Accepted_Name)]
Traits_fix$seed_mass<- as.numeric(Traits_fix$seed_mass)

#Filter the master data in order to plot some thrits and abbundance 
master_data<- filter(master_data, master_data$CODE %in% 
                       Traits_fix$Code)

# Try out some traits and abbundance 
plot(Traits_fix$THK, (master_data$tpa_2014))
abline(lm(master_data$tpa_2014 ~ Traits_fix$THK))

plot(Traits_fix$SLA.wp,(master_data$tpa_2014))
abline(lm(master_data$tpa_2014 ~ Traits_fix$SLA.wp))

plot(Traits_fix$MAXHT, master_data$tpa_2014)
abline(lm(master_data$tpa_2014 ~ Traits_fix$MAXHT))

plot(Traits_fix$WD, master_data$tpa_2014)
abline(lm(master_data$tpa_2014 ~ Traits_fix$WD))

plot((Traits_fix$seed_mass), log(master_data$tpa_2014))
abline(lm(master_data$tpa_2014 ~ Traits_fix$seed_mass))

# when is not log transformed you can see the tendency! 
plot((Traits_fix$seed_mass), (master_data$tpa_2014))
abline(lm(master_data$tpa_2014 ~ Traits_fix$seed_mass))

plot(Traits_fix$LP.mass, log(master_data$tpa_2014) , cex=(master_data$fcover_51)/1500)
abline(lm(master_data$tpa_2014 ~ Traits_fix$LP.mas))


plot(master_data$fcover_51, Traits_fix$LP.mass, cex=master_data$tpa_2014/10)
abline(lm(Traits_fix$LP.mas ~ master_data$fcover_51))

plot(Traits_fix$LP.mass, master_data$tpa_2014, cex=(master_data$fcover_51)/1500)
abline(lm(master_data$tpa_2014~Traits_fix$LP.mas))












