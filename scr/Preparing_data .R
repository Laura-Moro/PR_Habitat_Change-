library(dplyr)
library(tidyverse)
library(ggplot2)

#Preparing a master_uf file with all of the 

#load species list
Sp_list <- read.csv("Data/PR_Trees_species_list.csv")
FIA <- read.csv("Data/FIA/FIA_abundnace.csv")
F_cover <-read.csv("Data/Derived/Forest_Cover.csv", sep=",")
AI_total <- read.csv("Data/Derived/AI_total.csv", sep=";")
D_total <- read.csv("Data/Derived/D_total.csv")
traits <- read.csv("Data/Derived/Trait_complete.csv")

#format species names 
# paste togher species name  
Sp_list$SPECIES <- paste(Sp_list$GENUS, Sp_list$SPECIES, sep=" ")
#transform everything to lowercase 
Sp_list$SPECIES <- tolower(Sp_list$SPECIES)
#transform first letter of the genus capital 
Sp_list$SPECIES <- paste0(toupper(substr(Sp_list$SPECIES, 1, 1)), substr(Sp_list$SPECIES, 2, nchar(Sp_list$SPECIES)))

FIA$code <- Sp_list$SP.CODE[match(FIA$SCIENTIFIC_NAME,Sp_list$SPECIES)]

#making a data frame only with the species name and the species code 
x <- filter(F_cover,F_cover$X %in% FIA$code)
x$Species<- FIA$SCIENTIFIC_NAME[match(x$X,FIA$code)]
x <- x[, c(1, 7,2,3,4,5,6)] #this I need to ajugst 
colnames(x)[1]<- "CODE"

#add columns with the FIA abundance data to the data fame X
d2004 <- FIA[FIA$YEAR==2004,]
x$tpa_2004<- d2004$TPA[match(x$CODE, d2004$code)]

d2008 <- FIA[FIA$YEAR==2008,]
x$tpa_2008 <- d2008$TPA[match(x$CODE, d2008$code)]

d2009 <- FIA[FIA$YEAR==2009,]
x$tpa_2009 <- d2009$TPA[match(x$CODE, d2009$code)]

d2013<- FIA[FIA$YEAR==2013,]
x$tpa_2013 <- d2013$TPA[match(x$CODE, d2013$code)]

d2014<- FIA[FIA$YEAR==2014,]
x$tpa_2014 <- d2014$TPA[match(x$CODE, d2014$code)]

d2019<- FIA[FIA$YEAR==2019,]
x$tpa_2019<- d2019$TPA[match(x$CODE, d2019$code)]

#add the forest cover 
x$fcover_51<- F_cover$fcover_51[match(x$CODE, F_cover$X)]
x$fcover_77<- F_cover$fcover_77[match(x$CODE, F_cover$X)]
x$fcover_91<- F_cover$fcover_91[match(x$CODE, F_cover$X)]
x$fcover_00<- F_cover$fcover_00[match(x$CODE, F_cover$X)]
x$tot_change<- F_cover$tot_change[match(x$CODE, F_cover$X)]

#add Aggregation index 
x$AI_51<- AI_total$AI_51[match(x$CODE, AI_total$CODE)]
x$AI_77<- AI_total$AI_77[match(x$CODE, AI_total$CODE)]
x$AI_91<- AI_total$AI_91[match(x$CODE, AI_total$CODE)]
x$AI_00<- AI_total$AI_00[match(x$CODE, AI_total$CODE)]

#euclidean mean distance 
x$D_51<- D_total$D_51[match(x$CODE, D_total$Code)]
x$D_77<- D_total$D_71[match(x$CODE, D_total$Code)]
x$D_91<- D_total$D_91[match(x$CODE, D_total$Code)]
x$D_00<- D_total$D_00[match(x$CODE, D_total$Code)]

#save the data frame with Abundance, forest cover, AI , euclidean distance between patches
write.csv(x,"Data/Derived/master.csv")

#Add trait data and make a new data frame
master_f <- read.csv ("Data/Derived/master.csv")
master_uf <- read.csv ("Data/Derived/master.csv")
Trait_f <- read.csv("Data/Derived/Trait_complete_f.csv")
Trait_uf <- read.csv("Data/Derived/Trait_complete_uf.csv")

master_f$WD <- Trait_f$WD[match(master_f$CODE, Trait_f$code)]
master_f$THK<- Trait_f$THK[match(master_f$CODE, Trait_f$code)]
master_f$LA.wp <- Trait_f$LA.wp[match(master_f$CODE, Trait_f$code)]
master_f$SLA.wp <- Trait_f$SLA.wp[match(master_f$CODE, Trait_f$code)]
master_f$MAXHT <- Trait_f$MAXHT[match(master_f$CODE, Trait_f$code)]
master_f$seed_mass <- Trait_f$seed_mass[match(master_f$CODE, Trait_f$code)]
master_f$PC1 <- Trait_f$Dim1[match(master_f$CODE, Trait_f$code)]
master_f$PC2 <- Trait_f$Dim2[match(master_f$CODE, Trait_f$code)]

master_f<- filter(master_f, master_f$CODE %in% Trait_f$code)

#add columns with the species trait data 
master_uf$WD <- Trait_uf$WD[match(master_uf$CODE, Trait_uf$code)]
master_uf$THK<- Trait_uf$THK[match(master_uf$CODE, Trait_uf$code)]
master_uf$LA.wp <- Trait_uf$LA.wp[match(master_uf$CODE, Trait_uf$code)]
master_uf$SLA.wp <- Trait_uf$SLA.wp[match(master_uf$CODE, Trait_uf$code)]
master_uf$LP.mass<- Trait_uf$LP.mass[match(master_uf$CODE, Trait_uf$code)]
master_uf$MAXHT <- Trait_uf$MAXHT[match(master_uf$CODE, Trait_uf$code)]
master_uf$seed_mass <- Trait_uf$seed_mass[match(master_uf$CODE, Trait_uf$code)]
master_uf$PC1 <- Trait_uf$Dim1[match(master_uf$CODE, Trait_uf$code)]
master_uf$PC2 <- Trait_uf$Dim2[match(master_uf$CODE, Trait_uf$code)]

write.csv(master_f, "Data/Derived/master_f.csv")
write.csv(master_uf, "Data/Derived/master_uf.csv")











 
