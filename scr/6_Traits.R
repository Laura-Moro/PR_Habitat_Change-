##PREPARING DATA 
library(dplyr)
library(readxl)

#1 SPecies Traits 
#compiling the Traint spreas sheets 

traits_E <- read.csv("Data/Traits/Species_Traits.csv")
abundance <- read.csv("Data/FIA/FIA_abundnace.csv", sep = ";" )
guanica <- read.csv("Data/Traits/Guanica_seed_mass.csv")
Traits_fix<- read.csv("Data/Traits/Fixed_traits.csv", sep = ";" )



#extract species in both lists 
sp_a <- Abundance$SCIENTIFIC_NAME
sp_t <- Traits$PLANTS_Accepted_Name

#species in commmon between the tow 
Species_in_common_df <- intersect(sp_a, sp_t)

#filter the original data 
filtered_df <- filter(Traits, Traits$PLANTS_Accepted_Name %in% 
                        Abundance$SCIENTIFIC_NAME)

filtered_df <- write.csv(filtered_df, "Data/Traits/Traits_abbundnace.csv")




