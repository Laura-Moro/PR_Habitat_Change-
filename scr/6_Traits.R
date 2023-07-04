##PREPARING Trait Data 
library(dplyr)
library(factoextra)

#Species Traits 
#compiling the Traint spreas sheets 

traits_E <- read.csv("Data/Traits/Species_Traits.csv")
guanica <- read.csv("Data/Traits/Guanica_seed_mass.csv")
Traits_fix<- read.csv("Data/Traits/Fixed_traits.csv", sep = ";" )
master_data <- read.csv("Data/Derived/Abundance_Fcover_AI.csv",sep = ";")

#take the fixed traits data set and filter the species we have abundance data for 
#it is only 168? 
Traits_fix<- filter(Traits_fix, Traits_fix$Code %in% master_data$CODE)

#add species names to the traits data set 
Traits_fix$species <- master_data$SCIENTIFIC_NAME [match(Traits_fix$Code, master_data$CODE)]
#add species codes to the helmer data set 
traits_E$CODE <- master_data$CODE[match(traits_E$PLANTS_Accepted_Name, master_data$SCIENTIFIC_NAME)]


#units: seed mass of kew is in for 1000 seeds and in g, seed mass of guanica is in for 1 seed in mg--> they are equicalent  
#take the data on the seed mass form the Helmer data set (extact only infomation on seed mass form the df
Traits_fix$seed_mass <- traits_E$Seed_wt_avg_Kew_or_other_g_per_1000[match(Traits_fix$species,traits_E$PLANTS_Accepted_Name)]
Traits_fix$seed_mass <- as.numeric(Traits_fix$seed_mass)

#!!!!!add seed mass form guanica # not very sure how to do this!!!!
Traits_fix$seed_mass <- guanica$SS.MG[match(Traits_fix$Code, guanica$SPECIES)]

#Filter the master data in order to plot abbundance and and traits 
master_data<- filter(master_data, master_data$CODE %in% Traits_fix$Code)

# Try plotting out some traits and abundance 
# leaf thickness 
plot(Traits_fix$THK, (master_data$tpa_2014))
abline(lm(master_data$tpa_2014 ~ Traits_fix$THK))
#Specific leaf area 
plot(Traits_fix$SLA.wp,(master_data$tpa_2014))
abline(lm(master_data$tpa_2014 ~ Traits_fix$SLA.wp))
#maximum height 
plot(Traits_fix$MAXHT, master_data$tpa_2014)
abline(lm(master_data$tpa_2014 ~ Traits_fix$MAXHT))
#wood density 
plot(Traits_fix$WD, master_data$tpa_2014)
abline(lm(master_data$tpa_2014 ~ Traits_fix$WD))
#sead mass 
plot((Traits_fix$seed_mass), log(master_data$tpa_2014))
abline(lm(master_data$tpa_2014 ~ Traits_fix$seed_mass))

# when is not log transformed you can see the tendency
plot(Traits_fix$seed_mass, master_data$tpa_2014)
abline(lm(master_data$tpa_2014 ~ Traits_fix$seed_mass))

plot(Traits_fix$WD,(master_data$tpa_2014) , cex=(master_data$fcover_51)/1500)
abline(lm(master_data$tpa_2014 ~ Traits_fix$WD))

plot(master_data$fcover_51, Traits_fix$LP.mass, cex=master_data$tpa_2014/10)
abline(lm(Traits_fix$LP.mas ~ master_data$fcover_51))

plot(Traits_fix$LP.mass, master_data$tpa_2014, cex=(master_data$fcover_51)/1500)
abline(lm(master_data$tpa_2014~Traits_fix$LP.mas))

#look at the specie in the multivariate space
#preparedata take only the values for the 6 variables 
data_PCA <- Traits_fix [, c(1:5,7,10)]

#Transform into a matrix 
PCA_mat <- as.data.frame(data_PCA)

#run the pca (i think it is dropping the NA values) 
# maybe try to center the variables 
res.pca <- prcomp(~ data_PCA$WD + 
                    data_PCA$THK + 
                    data_PCA$LA.wp +
                    data_PCA$SLA.wp +
                    data_PCA$MAXHT +
                    data_PCA$seed_mass, 
                    data=data_PCA, scale = TRUE)

#look at the scree plot 
fviz_eig(res.pca)

#plot the pca 
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_biplot(res.pca, repel = TRUE,
                col.ind = "cos2",
                radient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                col.var = "contrib",
)

#look at the results 

#eigenvalues
eig.val <- get_eigenvalue(res.pca)

#results for variables 
res.var <- get_pca_var(res.pca)
res.var$coord #Coordinates for the variables
res.var$cor #Correlations between variables and dimensions
res.var$cos2 #Quality of rapresentation 
res.var$contrib #contribution to the PCs 

#results for individuals 
res.ind <- get_pca_ind(res.pca)
res.ind$coord # Coordinates of distribution of each species in the multivariate space 
res.ind$contrib #contribution to the PCs 
res.ind$cos2 #quality of representation 

#run model of interaction between traits and abundance and forets cover 
Abbundace <- master_data$tpa_2014
forest_1951 <- master_data$fcover_51
forest_2000 <- master_data$fcover_00


#WOOD DENSITY 
fit_Wd_51<- lm(scale(Abbundace) ~ scale(master_data$fcover_51) + 
                 scale(data_PCA$WD) +  
               scale(master_data$fcover_51*data_PCA$WD))
summary(fit_Wd_51)

fit_Wd_00<- lm(Abbundace ~ master_data$fcover_00 + 
                 data_PCA$WD +  
                 master_data$fcover_00*data_PCA$WD)
summary(fit_Wd_00)

#LEAVES THIKNESS
#1951
fit_thk_51<- lm(Abbundace ~ master_data$fcover_51 + 
                 data_PCA$THK +  
                 master_data$fcover_51*data_PCA$THK)
summary(fit_thk_51)
#2000
fit_thk_00<- lm(Abbundace ~ master_data$fcover_00 + 
                 data_PCA$THK +  
                 master_data$fcover_00*data_PCA$THK)
summary(fit_thk_00)

#LEAVES AREA
#1951
fit_la_51<- lm(Abbundace ~ master_data$fcover_51 + 
                  data_PCA$LA.wp +  
                  master_data$fcover_51*data_PCA$THK)
summary(fit_la_51)
#2000
fit_la_00<- lm(Abbundace ~ master_data$fcover_00 + 
                  data_PCA$LA.wp +  
                  master_data$fcover_00*data_PCA$THK)
summary(fit_la_00)


#SPECIFIC LEAVES AREA
#1951
fit_sla_51<- lm(Abbundace ~ master_data$fcover_51 + 
                 data_PCA$SLA.wp +  
                 master_data$fcover_51*data_PCA$SLA.wp)
summary(fit_sla_51)
#2000
fit_la_00<- lm(Abbundace ~ master_data$fcover_00 + 
                 data_PCA$SLA.wp +  
                 master_data$fcover_00*data_PCA$SLA.wp)
summary(fit_sla_00)

#MAXIMUM HIGHT 
#1951
fit_mh_51<- lm(Abbundace ~ master_data$fcover_51 + 
                  data_PCA$MAXHT +  
                  master_data$fcover_51*data_PCA$MAXHT)
summary(fit_mh_51)
#2000
fit_mh_00<- lm(Abbundace ~ master_data$fcover_00 + 
                 data_PCA$MAXHT +  
                 master_data$fcover_00*data_PCA$MAXHT)
summary(fit_mh_00)


#SEED MASS
#1951
fit_sm_51<- lm(Abbundace ~ master_data$fcover_51 + 
                 data_PCA$seed_mass +  
                 master_data$fcover_51*data_PCA$seed_mass)
summary(fit_sm_51)
#2000
fit_sm_00<- lm(Abbundace ~ master_data$fcover_00 + 
                 data_PCA$seed_mass +  
                 master_data$fcover_00*data_PCA$seed_mass)
summary(fit_sm_00)














































