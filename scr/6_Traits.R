library(dplyr)
library(tidyverse)
library(factoextra)
library(tibble)
library(ape)
library(Rphylopars)

#1 
#compiling the Trait data sheets 
traits_E <- read.csv("Data/Traits/Species_Traits.csv")
Traits_fix<- read.csv("Data/Traits/Fixed_traits.csv", sep = ";" )
Master <- read.csv("Data/Derived/Master.csv")

#add species names to the traits data set 
Traits_fix$species <- Master$Species [match(Traits_fix$Code, Master$CODE)]
#add species codes to the helmer data set 
traits_E$CODE <- Master$CODE[match(traits_E$PLANTS_Accepted_Name, Master$species)]

#units: seed mass of kew is in for 1000 seeds and in g, seed mass of guanica is in for 1 seed in mg--> they are equicalent  
#take the data on the seed mass form the Helmer data set (extact only infomation on seed mass form the df
Traits_fix$seed_mass <- traits_E$Seed_wt_avg_Kew_or_other_g_per_1000[match(Traits_fix$species,traits_E$PLANTS_Accepted_Name)]
Traits_fix$seed_mass <- as.numeric(Traits_fix$seed_mass)

Trait_f <- filter(Traits_fix, Traits_fix$Code %in% Master$CODE)
Trait_f<- Trait_f[c(1,9,2:5,7,10)]
write_csv(Trait_f, "Data/Derived/Trait_f.csv")

#traits we have at least one data one trait values per species
Trait_uf <- Master[c(2,3)]
Trait_uf$WD <- Trait_f$WD[match(Trait_uf$CODE, Trait_f$Code)]
Trait_uf$LA.wp <- Trait_f$LA.wp[match(Trait_uf$CODE, Trait_f$Code)]
Trait_uf$SLA.wp <- Trait_f$SLA.wp[match(Trait_uf$CODE, Trait_f$Code)]
Trait_uf$MAXHT <- Trait_f$MAXHT[match(Trait_uf$CODE, Trait_f$Code)]
Trait_uf$THK <- Trait_f$THK[match(Trait_uf$CODE, Trait_f$Code)]
Trait_uf$seed_mass <- Trait_f$seed_mass[match(Trait_uf$CODE, Trait_f$Code)]
Trait_uf<- Trait_uf[c(2,1,3,7,4,5,6,8)]
write.csv(Trait_uf, "Data/Derived/Trait_uf.csv")

##2 Perfor trait imputation 

#load treated trait values! i had to remove some specie than were not presnt in the phylogeny 
Trait_f <-read.csv("Data/Derived/Trait_f.csv")
Trait_uf <- read.csv("Data/Derived/Trait_uf.csv")

#reorder columns and elimiante columns
colnames(Trait_uf)[2] <- "species"

T_f <- Trait_f[,2:8]
T_uf <- Trait_uf[,c (-3,-1)]

#remove space 
T_f$species <- sub(" ", "_", T_f$species)
T_uf$species <- sub(" ", "_", T_uf$species)

#removing unmatchig rows with the phylogeny 
f = c("Chrysophyllum_pauciflorum","Coccoloba_diversifolia", "Erythroxylum_rotundifolium", "Henriettea_squamulosum", "Licaria_brittoniana", "Maytenus_elongata", "Persea_americana", "Pithecellobium_unguis-cati", "Psychotria_berteriana", "Rochefortia_acanthophora", "Schefflera_morototonii")
f<-as.data.frame(f)
uf = c ("Alsophila_portoricensis","Cestrum_laurifolium", "Chrysophyllum_pauciflorum", "Citrus_x sinensis" , "Coccoloba_diversifolia", "Erythroxylum_rotundifolium", "Eugenia_cordata", "Henriettea_squamulosum", "Licaria_brittoniana", "Maytenus_elongata", "Maytenus_ponceana", "Myrcia_citrifolia", "Palicourea_crocea", "Persea_americana", "Pithecellobium_unguis-cati", "Prosopis_pallida", "Psidium_guajava", "Psychotria_berteriana","Rochefortia_acanthophora", "Schefflera_morototonii","Trema_micrantha")
uf<-as.data.frame(uf)

T_f <- T_f[!(T_f$species %in% f$f), ]
T_uf <- T_uf[!(T_uf$species %in% uf$uf),]

#save record of imputed data for T_uf  
T_f_imo <- data.frame(T_f$species, is.na(T_f)[,2:7])
colnames(T_f) <- c("species", "WD", "THK", "LA.wp", "SLA.wp", "MAXHT", "seed_mass")

#save record of imputed data for T_uf  
T_uf_imp <- data.frame(T_uf $species, is.na(T_uf )[,2:7])
colnames(T_uf) <- c("species", "WD", "THK", "LA.wp", "SLA.wp", "MAXHT", "seed_mass")

#Scale data 
T_uf$WD <- scale(T_uf $WD)
T_uf$THK <- scale(T_uf $THK)
T_uf$LA.wp <- scale(T_uf $LA.wp)
T_uf$SLA.wp <- scale(T_uf $SLA.wp)
T_uf$MAXHT <- scale(T_uf $MAXHT)
T_uf$seed_mass <- scale(T_uf $seed_mass)

#transform data
T_f$WD <- scale(T_f$WD)
T_f$THK <- scale(T_f$THK)
T_f$LA.wp <- scale(T_f$LA.wp)
T_f$SLA.wp <- scale(T_f$SLA.wp)
T_f$MAXHT <- scale(T_f$MAXHT)
T_f$seed_mass <- scale(T_f$seed_mass)

#phylogeny
myTree <- ape::read.tree(file = "Data/traits/phylogeny.new")
#choose the first phylogeny 
Phy<- myTree[[1]]

#trip the phylogeny with the species we have data
Tree_f <- keep.tip(Phy, T_f$species)
Tree_uf <- keep.tip(Phy , T_uf$species)

#assing 0.0001 values to the branch leght --> otherwise the phylopars does not work
Tree_f$edge.length[which(Tree_f$edge.length==0)]=0.0001
Tree_uf$edge.length[which(Tree_uf$edge.length==0)]=0.0001

#preform phylo imputation 
IMP_f <-phylopars(trait_data = T_f, tree=Tree_f, phylo_correlated = T)

IMP_uf <- phylopars(trait_data = T_uf , tree= Tree_uf , phylo_correlated = T)

#view imputed data 
I_f<-IMP_f$anc_recon[1:157,]
I_f <-as.data.frame(I_f)
I_f <- tibble::rownames_to_column(I_f, "Species")
I_f$Species <- sub("_", " ", I_f$Species)
I_f$code <- Trait_f$Code[match(I_f$Species,Trait_f$species)]
I_f<- I_f[,c(8,1,2,3,4,5,6,7)]
write.csv(I_f, "Data/Derived/Trait_f_I.csv")

#view imputed data 
#the first 157 rows are the thaits imputed means 
I_uf<-IMP_uf$anc_recon[1:260,]
I_uf <-as.data.frame(I_uf)
I_uf<- tibble::rownames_to_column(I_uf, "species")
I_uf$species <- sub("_", " ", I_uf$species)
I_uf$code <- Trait_uf$CODE[match(I_uf$species ,Trait_uf$species)]
I_uf<- I_uf[,c(8,1,2,3,4,5,6,7)]
write.csv(I_uf, "Data/Derived/Trait_uf_I.csv")


##3
#look at the specie in the multivariate space PCA
PCA_f <- read.csv("Data/Derived/Trait_f_I.csv")
PCA_uf <- read.csv("Data/Derived/Trait_uf_I.csv")


#run the pca (i think it is dropping the NA values) 
# maybe try to center the variables 
res.pca_f <- prcomp(~ PCA_f$WD + 
                    PCA_f$THK + 
                    PCA_f$LA.wp +
                    PCA_f$SLA.wp +
                    PCA_f$MAXHT +
                    PCA_f$seed_mass, 
                  data=PCA_f, scale = TRUE)

res.pca_uf <- prcomp(~ PCA_uf$WD + 
                    PCA_uf$THK + 
                    PCA_uf$LA.wp +
                    PCA_uf$SLA.wp +
                    PCA_uf$MAXHT +
                    PCA_uf$seed_mass, 
                  data=PCA_uf, scale = TRUE)

#look at the scree plot 
fviz_eig(res.pca_f)
fviz_eig(res.pca_uf)

#plot the pca 
fviz_pca_ind(res.pca_f,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
 
fviz_pca_var(res.pca_f,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(res.pca_uf,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_biplot(res.pca, repel = TRUE,
                col.ind = "cos2",
                radient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                col.var = "contrib",
)

#results for individuals 
res.ind_f<- get_pca_ind(res.pca_f)
res.ind_f$coord # Coordinates of distribution of each species in the multivariate space 

#results for individuals 
res.ind_uf<- get_pca_ind(res.pca_uf)
res.ind_uf$coord # Coordinates of distribution of each species in the multivariate space 


#add the two principal components as multivariate traits 
df_f <-as.data.frame(res.ind_f$coord)
df_uf <-as.data.frame(res.ind_uf$coord)


#add the code name to the PCA results 
df_f <- tibble::rowid_to_column(df_f, var ="id")
PCA_f <-tibble::rowid_to_column(PCA_f,var ="id")
df_f$code<- PCA_f$code[match(df_f$id,PCA_f$id)]
df_f<- df_f[,c(8,2,3)]

df_uf <- tibble::rowid_to_column(df_uf, var ="id")
PCA_uf <-tibble::rowid_to_column(PCA_uf,var ="id")
df_uf$code<- PCA_uf$code[match(df_uf$id,PCA_uf$id)]
df_uf<- df_uf[,c(8,2,3)]


#add principal components 1 and 2 to the  
PCA_f$Dim1 <- df_f$Dim.1 [match(df_f$code,PCA_f$code )]
PCA_f$Dim2<- df_f$Dim.2 [match(df_f$code,PCA_f$code )]

PCA_uf$Dim1 <- df_uf$Dim.1 [match(df_uf$code,PCA_uf$code)]
PCA_uf$Dim2<- df_uf$Dim.2 [match(df_uf$code,PCA_uf$code)]

#save data 
write.csv(PCA_f,"Data/Derived/Trait_complete_f.csv")
write.csv(PCA_uf,"Data/Derived/Trait_complete_uf.csv")













