##PREPARING Trait Data 
library(dplyr)
library(factoextra)
library(tibble)

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

data_complete<-bind_cols(master_data, data_PCA)


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

#results for individuals 
res.ind <- get_pca_ind(res.pca)
res.ind$coord # Coordinates of distribution of each species in the multivariate space 

#add the two principal components as multivariate traits 
df <-as.data.frame(res.ind$coord)

#add the code name to the PCA results 
df <- tibble::rowid_to_column(df, var ="id")
data_PCA <-tibble::rowid_to_column(data_PCA,var ="id")
df$CODE<- data_PCA$Code[match(df$id,data_PCA$id)]

#take only the PCA results for the first 2 dimentions 
df_P<- df[,c(2,3,8)]

#just for now I'm filtering the data so that I can use the match function and add the data pcw1 data to the dataset 
data_complete<- filter(data_complete, data_complete$CODE %in% df_P$CODE)
#add PC1(dim1)
data_complete$Dim1 <- df_P$Dim.1 [match(df_P$CODE, data_complete$CODE)]
#add PC2(dim2)
data_complete$Dim2 <- df_P$Dim.2 [match(df_P$CODE, data_complete$CODE)]

write.csv(data_complete, "Data/Derived/data_complete_10_07")


#Try running some models 

#Transform the data 
data_complete$log_abund <- log10(data_complete$tpa_2014)
data_complete$fcover_51.z <- scale(data_complete$fcover_51)
data_complete$fcover_00.z <- scale(data_complete$fcover_00)
data_complete$fcover_77.z <- scale(data_complete$fcover_77)
data_complete$fcover_91.z <- scale(data_complete$fcover_91)
data_complete$AI_51.z <- scale (data_complete$AI_51)
data_complete$AI_77.z <- scale (data_complete$AI_77)
data_complete$AI_91.z <- scale (data_complete$AI_91)
data_complete$AI_00.z <- scale (data_complete$AI_00)
data_complete$wd.z <- scale(data_complete$WD)
data_complete$LA.z <- scale(data_complete$LA.wp)
data_complete$THK.z <- scale(data_complete$THK)
data_complete$SLA.z<- scale(data_complete$SLA.wp)
data_complete$MAXHT.z <- scale(data_complete$MAXHT)
data_complete$seed_mass.z <- scale(data_complete$seed_mass)
data_complete$Dim1.z <- scale(data_complete$Dim1)
data_complete$Dim2.z <- scale(data_complete$Dim2)


#wood density 
fit_Wd_51<- lm(log_abund ~ fcover_51.z + 
                 fcover_00.z + 
                 AI_51.z + 
                 AI_00.z + 
                 wd.z + 
                 fcover_51.z * wd.z, data=data_complete)
summary(fit_Wd_51)

#Leaf area 
fit_LA_51<- lm(log_abund ~ fcover_51.z + 
                 fcover_00.z + 
                 AI_51.z + 
                 AI_00.z + 
                 LA.z + 
                 fcover_51.z * LA.z, data=data_complete)
summary(fit_LA_51)


#leaf thikness
fit_THK_51<- lm(log_abund ~ fcover_51.z + 
                 fcover_00.z + 
                 AI_51.z + 
                 AI_00.z + 
                 THK.z + 
                 fcover_51.z * THK.z, data=data_complete)
summary(fit_THK_51)

#Specific leaf area 
fit_SLA_51<- lm(log_abund ~ fcover_51.z + 
                 fcover_00.z + 
                 AI_51.z + 
                 AI_00.z + 
                  SLA.z + 
                 fcover_51.z * SLA.z, data=data_complete)
summary(fit_SLA_51)

#maximum hight 
fit_MAXHT_51<- lm(log_abund ~ fcover_51.z + 
                  fcover_00.z + 
                  AI_51.z + 
                  AI_00.z + 
                  MAXHT.z + 
                  fcover_51.z * MAXHT.z, data=data_complete)
summary(fit_MAXHT_51)

#seed mass
fit_seed_mass_51<- lm(log_abund ~ fcover_51.z + 
                    fcover_00.z + 
                    AI_51.z + 
                    AI_00.z + 
                    seed_mass.z + 
                    fcover_51.z * seed_mass.z, data=data_complete)
summary(fit_seed_mass_51)


fit_seed_mass_51<- lm(log_abund ~ fcover_51.z + 
                        fcover_00.z + 
                        AI_51.z + 
                        AI_00.z + 
                        seed_mass.z + 
                        fcover_51.z * seed_mass.z, data=data_complete)
summary(fit_seed_mass_51)

fit_Dim1_51<- lm(log_abund ~ fcover_51.z + 
                        fcover_00.z + 
                        AI_51.z + 
                        AI_00.z + 
                        Dim1.z + 
                        fcover_51.z * Dim1.z, data=data_complete)
summary(fit_Dim1_51)


fit_Dim2_51<- lm(log_abund ~ fcover_51.z + 
                   fcover_00.z + 
                   AI_51.z + 
                   AI_00.z + 
                   Dim2.z + 
                   fcover_51.z * Dim2.z, data=data_complete)
summary(fit_Dim2_51)

fit_Dim2_AI51<- lm(log_abund ~ fcover_51.z + 
                   fcover_00.z + 
                   AI_51.z + 
                   AI_00.z + 
                   Dim2.z + 
                   AI_51.z * Dim2.z, data=data_complete)
summary(fit_Dim2_AI51)

fit_Dim2_AI00<- lm(log_abund ~ fcover_51.z + 
                   fcover_00.z + 
                   AI_51.z + 
                   AI_00.z + 
                   Dim2.z + 
                     AI_00.z * Dim2.z, data=data_complete)
summary(fit_Dim2_AI00)


































