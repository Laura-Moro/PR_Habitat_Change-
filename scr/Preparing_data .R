library(dplyr)
library(tidyverse)
library(ggplot2)

#Preparing a master file with all of the 

#load species list
Sp_list <- read.csv("Data/PR_Trees_species_list.csv" , sep = ";")
ABCO <-read.csv("Data/Derived/Abundace_Cover.csv", sep = ";" ) # this is the Species we have data of the FAI to be maybe revisioned 
FIA <- read.csv("Data/FIA/FIA_abundnace.csv")
F_cover <-read.csv("Data/Derived/Forest_Cover.csv", sep=",")
AI_total <- read.csv("Data/Derived/AI_total.csv", sep=";")
traits <- read.csv("Data/Derived/Trait_complete.csv")

#format species names 
# paste togher species name  
Sp_list$SPECIES <- paste(Sp_list$GENUS, Sp_list$SPECIES, sep=" ")
#transform everything to lowercase 
Sp_list$SPECIES <- tolower(Sp_list$SPECIES)
#transform first letter of the genus capital 
Sp_list$SPECIES <- paste0(toupper(substr(Sp_list$SPECIES, 1, 1)), substr(Sp_list$SPECIES, 2, nchar(Sp_list$SPECIES)))

#making a data frame only with the species name and the species code 
x <- ABCO[, c(2,4)] #this I need to ajugst 

#add code to the the  Abundance data 
FIA$code <- x$CODE[match(FIA$SCIENTIFIC_NAME, x$SCIENTIFIC_NAME)]

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
x$fcover_51<- F_cover$fcover_51[match(x$CODE, F_cover$SPECIES)]
x$fcover_77<- F_cover$fcover_77[match(x$CODE, F_cover$SPECIES)]
x$fcover_91<- F_cover$fcover_91[match(x$CODE, F_cover$SPECIES)]
x$fcover_00<- F_cover$fcover_00[match(x$CODE, F_cover$SPECIES)]
x$tot_change<- F_cover$tot_change[match(x$CODE, F_cover$SPECIES)]

#add Aggregation index 
x$AI_51<- AI_total$AI_51[match(x$CODE, AI_total$CODE)]
x$AI_77<- AI_total$AI_77[match(x$CODE, AI_total$CODE)]
x$AI_91<- AI_total$AI_91[match(x$CODE, AI_total$CODE)]
x$AI_00<- AI_total$AI_00[match(x$CODE, AI_total$CODE)]

#save the data frame
write.csv(x, "Data/Derived/Abundance_Fcover_AI.csv")

#add Trait data 
Master_data <- read.csv("Data/Derived/Abundance_Fcover_AI.csv",sep=";")
trait <- read.csv("Data/Derived/Trait_complete.csv")

#filter fo 
Master_data<- filter(Master_data, Master_data$CODE %in% trait$code)

#add columns with the species trait data 
Master_data$WD <- trait$WD[match(Master_data$CODE, trait$code)]
Master_data$THK<- trait$THK[match(Master_data$CODE, trait$code)]
Master_data$LA.wp <- trait$LA.wp[match(Master_data$CODE, trait$code)]
Master_data$SLA.wp <- trait$SLA.wp[match(Master_data$CODE, trait$code)]
Master_data$LP.mass<- trait$LP.mass[match(Master_data$CODE, trait$code)]
Master_data$MAXHT <- trait$MAXHT[match(Master_data$CODE, trait$code)]
Master_data$seed_mass <- trait$seed_mass[match(Master_data$CODE, trait$code)]
Master_data$PC1 <- trait$Dim1[match(Master_data$CODE, trait$code)]
Master_data$PC2 <- trait$Dim2[match(Master_data$CODE, trait$code)]

write.csv(Master_data, "Data/Derived/Master_data.csv")
write.csv(Master_data, "Data/Derived/Master_data_uf.csv")











 
