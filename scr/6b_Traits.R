library(dplyr)
library(tidyverse)
library(factoextra)
library(tibble)
library(ape)
library(Rphylopars)

## 1. Compile the trait data

## Load data
# Species data
Sp_list <- read.csv("Data/Traits/PR_Trees_species_list.csv")
# Helmer trait data (seed mass)
traits_E <- read.csv("Data/Traits/Species_Traits.csv")
# Our species mean trait data
Traits_fix <- read.csv("Data/Traits/Fixed_traits.csv", sep = ";" )
# Output from previous script to add trait data into
data <- read.csv("Data/Derived/4b-output-20250314.csv")

## Get binomial species name
Sp_list$binom <- tolower(paste(Sp_list$GENUS, Sp_list$SPECIES, sep=" "))
Sp_list$binom <- paste0(toupper(substr(Sp_list$binom, 1, 1)), 
                        substr(Sp_list$binom, 2, nchar(Sp_list$binom)))
Sp_list$species <- sub(" ", "_", Sp_list$binom)
data$species <- Sp_list$species[match(data$sp, Sp_list$SP.CODE)]

# Add species codes to the helmer data set 
traits_E$CODE <- Sp_list$SP.CODE[match(traits_E$PLANTS_Accepted_Name, Sp_list$binom)]

# Add seed mass data
# Units: seed mass of kew is in for 1000 seeds and in g, 
# Seed mass of Guanica is in for 1 seed in mg --> they are equivalent
# Swenson seed mass is g per seed so * 1000
data$seedmass <- as.numeric(traits_E$Seed_wt_avg_Kew_or_other_g_per_1000[match(data$sp, traits_E$CODE)])

# ADD GUANICA SEED MASS
guanica <- read.csv("Data/Traits/Guanica_seed_mass.csv")
data$ss.guanica <- guanica$SS.MG[match(data$sp, guanica$SPECIES)]

# ADD LFDP SEED MASS
swenson <- read.csv("Data/Traits/Swenson.UmanaEtAl.dryad.csv")
data$ss.swenson <- (swenson$seed[match(data$sp, swenson$species.code)]) * 1000

# Choose seed mass
data$ss <- dplyr::coalesce(data$ss.swenson, data$ss.guanica, data$seedmass)

# Merge species-mean trait data with input data
data <- cbind(data, Traits_fix[match(data$sp, Traits_fix$Code), c("WD","THK","LA.wp","SLA.wp","MAXHT")])
rownames(data) <- NULL

data$ss.log <- log10(data$ss)
data$thk.log <- log10(data$THK)
data$la.log <- log10(data$LA.wp)
data$sla.log <- log10(data$SLA.wp)

data$wd.z <- scale(data$WD)
data$ss.log.z <- scale(data$ss.log)
data$thk.log.z <- scale(data$thk.log)
data$la.log.z <- scale(data$la.log)
data$sla.log.z <- scale(data$sla.log)
data$maxht.z <- scale(data$MAXHT)

### 2. Trait imputation 

# Subset temporary data for imputation
tdata <- data[,c("species","ss.log.z", "wd.z", "thk.log.z", "la.log.z", "sla.log.z", "maxht.z")]

# Phylogeny
phy <- ape::read.tree(file = "Data/traits/phylogeny.new")[[1]]

# Prune the phylogeny / trait data to the species we have data / in the phylo
phy <- keep.tip(phy, phy$tip.label[phy$tip.label %in% tdata$species])
tdata_phy <- tdata[tdata$species %in% phy$tip.label,]

# Give small (0.0001) values to 0 branch lengths otherwise phylopars does not work
phy$edge.length[which(phy$edge.length==0)] <- 0.0001

# Preform phylo imputation
tdat_pi <- phylopars(trait_data=tdata_phy, tree=phy, phylo_correlated=T)

# Merge imputed trait data back to full dataset
imp_trt <- as.data.frame(tdat_pi$anc_recon[match(tdata_phy$species, rownames(tdat_pi$anc_recon)),])
names(imp_trt) <- paste0(names(imp_trt), "_pi")
data <- cbind(data, imp_trt[match(data$species, rownames(imp_trt)),])


### 3. Trait Ordination
# Run the PCA
res.pca <- prcomp(~ ss.log.z +
                    wd.z + 
                    thk.log.z +
                    la.log.z +
                    sla.log.z +
                    maxht.z, 
                  data=data)

res.pca_pi <- prcomp(~ ss.log.z_pi +
                    wd.z_pi + 
                    thk.log.z_pi +
                    la.log.z_pi +
                    sla.log.z_pi +
                    maxht.z_pi, 
                  data=data)

## Scree plot 
# fviz_eig(res.pca)

## Plot the pca 
fviz_pca_ind(res.pca,
             geom="point",
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Extract results for species
res.ind <- get_pca_ind(res.pca)$coord[,1:2]
res.ind_pi <- get_pca_ind(res.pca_pi)$coord[,1:2]

# Add the first two principal components as multivariate traits
data$pca1 <- res.ind[match(data$species, rownames(res.ind)),1]
data$pca2 <- res.ind[match(data$species, rownames(res.ind)),2]
data$pca1_pi <- res.ind_pi[match(data$species, rownames(res.ind_pi)),1]
data$pca2_pi <- res.ind_pi[match(data$species, rownames(res.ind_pi)),2]

# Save data 
write.csv(data, "Data/Derived/6b-output-20250314.csv")

