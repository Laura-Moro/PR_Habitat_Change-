library(dplyr)
library(ape)
library(Rphylopars)

#load treated trait values! i had to remove some specie than were not presnt in the phylogeny 
trait <-read.csv("Data/traits/Traits_IMP.csv", sep=";")
df <- read.csv("Data/Derived/Master_data_uf.csv", sep=";")

#reorder columns and elimiante columns
trait<-trait[, c(10,3,4,5,6,8,11)]
df <- df [, c(1, 18,20,21,22,23,19)]
colnames(df)[1] <- "species"

#remove space 
trait$species <- sub(" ", "_", trait$species)
df$species <- sub(" ", "_", df$species)

#save record of imputed data for df 
trait_i <- data.frame(trait$species, is.na(trait)[,2:7])
colnames(trait_i) <- c("species", "WD", "THK", "LA.wp", "SLA.wp", "MAXHT", "seed_mass")

#save record of imputed data for df 
df_i <- data.frame(df$species, is.na(df)[,2:7])
colnames(df_i) <- c("species", "WD", "THK", "LA.wp", "SLA.wp", "MAXHT", "seed_mass")

#Scale data 
df$WD <- scale(df$WD)
df$THK <- scale(df$THK)
df$LA.wp <- scale(df$LA.wp)
df$SLA.wp <- scale(df$SLA.wp)
df$MAXHT <- scale(df$MAXHT)
df$seed_mass <- scale(df$seed_mass)

#transform data
trait$WD <- scale(trait$WD)
trait$THK <- scale(trait$THK)
trait$LA.wp <- scale(trait$LA.wp)
trait$SLA.wp <- scale(trait$SLA.wp)
trait$MAXHT <- scale(trait$MAXHT)
trait$seed_mass <- scale(trait$seed_mass)

#phylogeny
myTree <- ape::read.tree(file = "Data/traits/phylogeny.new")
#choose the first phylogeny 
Phy<- myTree[[1]]

#trip the phylogeny with the species we have data
Tree_filter <- keep.tip(Phy, trait$species)
Tree <- keep.tip(Phy , df$species)

#assing 0.0001 values to the branch leght --> otherwise the phylopars does not work
Tree_filter$edge.length[which(Tree_filter$edge.length==0)]=0.0001
Tree$edge.length[which(Tree$edge.length==0)]=0.0001

#preform phylo imputation 
IMP_t<-phylopars(trait_data = trait, tree=Tree_filter, phylo_correlated = T)

IMP_df<- phylopars(trait_data = df, tree= Tree , phylo_correlated = T)

#view imputed data 
IMP_trait<-IMP_t$anc_recon[1:157,]
IMP_trait <-as.data.frame(IMP_trait)
IMP_trait <- tibble::rownames_to_column(IMP_trait, "Species")


#view imputed data 
#the first 157 rows are the thaits imputed means 
IMP_trait2<-IMP_df$anc_recon[1:205,]
IMP_trait2 <-as.data.frame(IMP_trait2)
IMP_trait2 <- tibble::rownames_to_column(IMP_trait2, "Species")








