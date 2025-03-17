library(rFIA)
library(sf)

# Load the species list
Sp_list <- read.csv("Data/Traits/PR_Trees_species_list.csv")

# Load data output from previous script
indata <- read.csv("Data/Derived/3b-output-20250314.csv") 

# Load PR outline shapefile
pr <- st_read("Data/PR_shapes/outline/PR_outline_Project.shp")

# Download FIA data for Puerto Rico (one time only)
# PR_Trees <- getFIA(states = 'PR', dir = 'Data/FIA/bob')

# Load data (excluding points off of mainland PR (e.g., Vieques, Culebra))
prfia <- clipFIA(readFIA("Data/FIA/bob"), mask=pr, mostRecent=F) 

# Estimate total population size by species
(pr_tpa <- tpa(prfia, bySpecies=TRUE, totals=TRUE))


# Compare abundance estimates from default TPA ('live') with 'growing stock'
pr_tpa2 <- pr_tpa[match(paste(pr_tpa_gs$YEAR, pr_tpa_gs$COMMON_NAME), 
                        paste(pr_tpa$YEAR, pr_tpa$COMMON_NAME)),]
plot(pr_tpa2$TPA, pr_tpa_gs$TPA, log='xy')
abline(0,1,col='blue', lwd=3)

# Compare abundance estimates from default TPA ('live') with 'seedlings'
pr_tpa_seedlings <- seedling(prfia, bySpecies=TRUE, totals=TRUE)
pr_tpa3 <- pr_tpa[match(paste(pr_tpa_seedlings$YEAR, pr_tpa_seedlings$COMMON_NAME), 
                        paste(pr_tpa$YEAR, pr_tpa$COMMON_NAME)),]
plot(pr_tpa3$TPA, pr_tpa_seedlings$TPA, log='xy')
abline(0,1,col='blue', lwd=3)


# Growing stock: live stems > 5 in. DBH which contain at least one 8 ft merchantable log
(pr_tpa_gs <- tpa(prfia, bySpecies=TRUE, totals=TRUE, treeType="gs"))#treeDomain=DIA>2))

# Growing stock: live stems > 5 in. DBH which contain at least one 8 ft merchantable log
(pr_tpa_diam2 <- tpa(prfia, bySpecies=TRUE, totals=TRUE, treeDomain=DIA>2))


pr_tpa$TPA_gs <- pr_tpa_gs$TPA[match(paste(pr_tpa$YEAR, 
                                           pr_tpa$COMMON_NAME),
                                     paste(pr_tpa_gs$YEAR, 
                                           pr_tpa_gs$COMMON_NAME))]

pr_tpa$TPA_diam2 <- pr_tpa_diam2$TPA[match(paste(pr_tpa$YEAR, 
                                                 pr_tpa$COMMON_NAME),
                                     paste(pr_tpa_diam2$YEAR, 
                                           pr_tpa_diam2$COMMON_NAME))]

pr_tpa$TPA_seedlings <- pr_tpa_seedlings$TPA[match(paste(pr_tpa$YEAR, 
                                                         pr_tpa$COMMON_NAME), 
                                                   paste(pr_tpa_seedlings$YEAR, 
                                                         pr_tpa_diam2$COMMON_NAME))]

plot(pr_tpa[,c("TPA","TPA_gs","TPA_diam2","TPA_seedlings")])


### Format species names 
# paste together species names
Sp_list$binom <- tolower(paste(Sp_list$GENUS, Sp_list$SPECIES, sep=" "))
# transform first letter of the genus capital 
Sp_list$binom <- paste0(toupper(substr(Sp_list$binom, 1, 1)), 
                        substr(Sp_list$binom, 2, nchar(Sp_list$binom)))
# Give sp code to FIA data
pr_tpa$code <- Sp_list$SP.CODE[match(pr_tpa$SCIENTIFIC_NAME, Sp_list$binom)]

# Separate by year
pr_tpa_2004 <- pr_tpa[pr_tpa$YEAR==2004,]
pr_tpa_2009 <- pr_tpa[pr_tpa$YEAR==2009,]
pr_tpa_2014 <- pr_tpa[pr_tpa$YEAR==2014,]
pr_tpa_2019 <- pr_tpa[pr_tpa$YEAR==2019,]

# Merge abundance data with habitat data from previous script
indata$tpa2004 <- pr_tpa_2004$TPA[match(indata$sp, pr_tpa_2004$code)]
indata$tpa2009 <- pr_tpa_2009$TPA[match(indata$sp, pr_tpa_2009$code)]
indata$tpa2014 <- pr_tpa_2014$TPA[match(indata$sp, pr_tpa_2014$code)]
indata$tpa2019 <- pr_tpa_2019$TPA[match(indata$sp, pr_tpa_2019$code)]

indata$tpa_gs_2004 <- pr_tpa_2004$TPA_gs[match(indata$sp, pr_tpa_2004$code)]
indata$tpa_gs_2009 <- pr_tpa_2009$TPA_gs[match(indata$sp, pr_tpa_2009$code)]
indata$tpa_gs_2014 <- pr_tpa_2014$TPA_gs[match(indata$sp, pr_tpa_2014$code)]
indata$tpa_gs_2019 <- pr_tpa_2019$TPA_gs[match(indata$sp, pr_tpa_2019$code)]

indata$tpa_diam2_2004 <- pr_tpa_2004$TPA_diam2[match(indata$sp, pr_tpa_2004$code)]
indata$tpa_diam2_2009 <- pr_tpa_2009$TPA_diam2[match(indata$sp, pr_tpa_2009$code)]
indata$tpa_diam2_2014 <- pr_tpa_2014$TPA_diam2[match(indata$sp, pr_tpa_2014$code)]
indata$tpa_diam2_2019 <- pr_tpa_2019$TPA_diam2[match(indata$sp, pr_tpa_2019$code)]

indata$tpa_seedlings_2004 <- pr_tpa_2004$TPA_seedlings[match(indata$sp, 
                                                             pr_tpa_2004$code)]
indata$tpa_seedlings_2009 <- pr_tpa_2009$TPA_seedlings[match(indata$sp, 
                                                             pr_tpa_2009$code)]
indata$tpa_seedlings_2014 <- pr_tpa_2014$TPA_seedlings[match(indata$sp, 
                                                             pr_tpa_2014$code)]
indata$tpa_seedlings_2019 <- pr_tpa_2019$TPA_seedlings[match(indata$sp, 
                                                             pr_tpa_2019$code)]

# Merge FIA plot occupancy data with habitat data from previous script
indata$nPlots2004 <- pr_tpa_2004$nPlots_TREE[match(indata$sp, pr_tpa_2004$code)]
indata$nPlots2009 <- pr_tpa_2009$nPlots_TREE[match(indata$sp, pr_tpa_2009$code)]
indata$nPlots2014 <- pr_tpa_2014$nPlots_TREE[match(indata$sp, pr_tpa_2014$code)]
indata$nPlots2019 <- pr_tpa_2019$nPlots_TREE[match(indata$sp, pr_tpa_2019$code)]

# Save analysis data with abundance data
write.csv(indata, "Data/Derived/4b-output-20250314.csv", row.names = F)
