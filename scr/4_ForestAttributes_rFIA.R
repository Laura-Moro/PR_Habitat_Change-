## Load rFIA package
library(rFIA)

#Download FIA data for Puerto Rico 
PR_Trees <- getFIA(states = 'PR', dir = '/Users/laumo791/Documents/PR/C1/PR_Habitat_Change-/Data/FIA/NEW')

#load data
PR_Trees <- readFIA("Data/FIA/NEW") #fix the download bj
head(PR_Trees)

#Only estimates for the most recent inventory year
MR_PR_Trees <- clipFIA(PR_Trees)#mostRecent = TRUE) #subset of the most recent dara (MR)
MR_PR_Trees_tpa <- tpa(MR_PR_Trees)
head(MR_PR_Trees_tpa)

#All Inventory Years Available
PR_Trees_tpa<-tpa(PR_Trees)
head(PR_Trees_tpa)

# Gorup data by species
PR_Trees_sp <- tpa(PR_Trees , bySpecies = TRUE)
PR_Trees_sp_total <- tpa(MR_PR_Trees , bySpecies = TRUE, totals = TRUE, treeDomain = TRUE, areaDomain = TRUE)
head(PR_Trees_sp)
View(PR_Trees_sp)

#Abundnce Data are saved in Results/FIA

FIA_abbundance <- write.csv(PR_Trees_sp, "Data/FIA/FIA_abundnace.csv")















                      







