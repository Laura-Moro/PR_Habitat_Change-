library(ape)
library(Rphylopars)

#trait data
trait <-read.csv("Data/traits/Traits_IMP.csv", sep=";")
#reorder columns
trait<-trait[, c(10,3,4,5,6,7,8,11)]
#remove spacve 
trait$species <- sub(" ", "_", trait$species)

#phylogeny
myTree <- ape::read.tree(file = "Data/traits/phylogeny.new")


#choose the first phylogeny 
tree<- myTree[[1]]
tree3 <- myTree[[3]]


tree$edge.length[which(tree$edge.length==0)]=0.0001


IMP_data<-phylopars(trait_data=trait, tree=tree, phylo_correlated = FALSE)

#view imputed data 
IMP_data[["anc_recon"]]




