library(dplyr)
library(rstatix)
library(tidyverse)
library(ggplot2)

### THIS script needs to be reorganzie! some parts are missing 

#load species list
Sp_list <- read.csv("Data/PR_Trees_species_list.csv" , sep = ";")
abundance <- read.csv("Data/FIA/FIA_abundnace.csv", sep = ";" )

#format species names 
# paste togher species name  
Sp_list$SPECIES <- paste(Sp_list$GENUS, Sp_list$SPECIES, sep=" ")
#trnsform everythig to lowercase 
Sp_list$SPECIES <- tolower(Sp_list$SPECIES)
#transform first letter of the genus capital 
Sp_list$SPECIES <- paste0(toupper(substr(Sp_list$SPECIES, 1, 1)), substr(Sp_list$SPECIES, 2, nchar(Sp_list$SPECIES)))

#making a data frame only with the species name and the species code 
sp_code <- Sp_list[4:5]
#making a data frame with the species name and the 
TPA_species <- abundance[5:6]
colnames(TPA_species )[1]  <- "SPECIES" 

#filetr data frame
sp_code<- filter(sp_code, sp_code$SPECIES %in% TPA_species$SPECIES)
sum(duplicated(sp_code$SPECIES))

#filetr data frame
TPA_species <- filter(TPA_species, TPA_species$SPECIES %in% sp_code$SPECIES)
#see if there are duplicates 
sum(duplicated(TPA_species$SPECIES))
#Remove duplicates--- why? I cannot do it? 
TPA<- TPA_species[-185,]

Abco<- read.csv("/Users/laumo791/Documents/PR/C1/Results/Abundace_Cover.csv", sep=";")

# TPA ALL OTHER YEARS 
TPA_04 <- PR_Trees_sp[c(1:283), c(1,4,5)]
TPA_04<- filter(TPA_04 , TPA_04$SCIENTIFIC_NAME %in% 
                  sp_code$SPECIES, .preserve = TRUE)


TPA_08 <- PR_Trees_sp[c(284:313), c(1,4,5)]
TPA_08<- filter(TPA_08, if(TPA_08$SCIENTIFIC_NAME %in% 
                  sp_code$SPECIES),


TPA_09 <- PR_Trees_sp[c(314:601), c(1,4,5)]
TPA_09 <- filter(TPA_09 , TPA_09$SCIENTIFIC_NAME %in% 
                  sp_code$SPECIES)


TPA_13 <- PR_Trees_sp[c(602:625), c(1,4,5)]
TPA_13 <- filter(TPA_13 , TPA_13$SCIENTIFIC_NAME %in% 
                   sp_code$SPECIES)

TPA_14 <- PR_Trees_sp[c(626:917), c(1,4,5)]
TPA_14<- filter(TPA_14 , TPA_09$SCIENTIFIC_NAME %in% 
                   sp_code$SPECIES)












#FOREST COVER 
#take the predictions stacks for the different time perios (scrips 5_landscape)
Pred_f51
Pred_f77 
Pred_f91 
Pred_f00 

#summ all of the pixels that were forets in 51
fcover_51<-cellStats(Pred_f51, 'sum')
#summ all of the pixels that were forets in 00
fcover_00<-cellStats(Pred_f00, 'sum')
#total change 
tot_change <- fcover_00-fcover_51

min(fcover_51)
min(fcover_00)
min(tot_change)

matplot(fcover_51, type='h', lty=1)
matplot(fcover_00, type='h', lty=1)
matplot(a, type='h', lty=1)

#plot abbundance against forest cover 
#load data of TREE ABUNDANCE(FIA)
#load data on forest cover 
Tree_abun <- read.csv("/Users/laumo791/Documents/PR/C1/Results/alberi.csv", sep = ";")
df_fcover_51 <- read.csv("/Users/laumo791/Documents/PR/C1/Results/df_fcover_51.csv", sep = ";")
df_fcover_00 <-read.csv("/Users/laumo791/Documents/PR/C1/Results/df_fcover_00.csv", sep = ";")
df_tot_change <- read.csv("/Users/laumo791/Documents/PR/C1/Results/df_tot_change.csv", sep = ";")

#filter data to consider only the species we Have abbundance data for 
ff51 <- filter(df_fcover_51,  df_fcover_51$CODE %in% 
                       Tree_abun$CODE)
#filter data to consider only the species we have abbundance data for  
ff00 <- filter(df_fcover_00,  df_fcover_00$CODE %in% 
                 Tree_abun$CODE)
#filter data to consider only the species we have abbundance data for 
ffChange <- filter(df_tot_change,  df_tot_change$CODE %in% 
                     Tree_abun$CODE)
#filter data to consider only the species we have abbundance data for 
Tree_abun <- filter(Tree_abun,  Tree_abun$CODE %in% 
                      ff00$CODE)
#make one data frame with abundance and forest cover data 
Abundace_Cover <- cbind (Tree_abun, ff51$fcover_51, ff00$fcover_00, ffChange$tot_change)

#tranfdorm forest cover abbundance in km2
Abundace_Cover$`ff51$fcover_51` <- Abundace_Cover$`ff51$fcover_51`*0.2025
Abundace_Cover$`ff00$fcover_00` <- Abundace_Cover$`ff00$fcover_00`*0.2025
Abundace_Cover$`ffChange$tot_change` <- Abundace_Cover$`ffChange$tot_change`*0.2025

#adjust file in excel 
write.csv(Abundace_Cover, "/Users/laumo791/Documents/PR/C1/Results/Abundace_Cover.csv" )
Abco<- read.csv("/Users/laumo791/Documents/PR/C1/Results/Abundace_Cover.csv", sep = ",")

 
