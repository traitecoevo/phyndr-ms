library(ape)
library(plotrix)

calculate.overlap<-function(tree.path,data.species){
  phy<-read.tree(tree.path)
  tree.species<-phy$tip.label
  counts<-data.frame(dataset=tree.path,num.sp.in.tree=length(tree.species),num.sp.in.data=length(data.species),overlap=sum(data.species%in%tree.species))
  return(counts)
}

wood<-read.csv("GlobalWoodinessDatabase.csv")
sr<-calculate.overlap("Vascular_Plants_rooted.dated.tre",wood$gs)


wood<-read.csv("GlobalWoodinessDatabase.csv")
sr<-calculate.overlap("Conifer-timetree.tre",wood$gs)

##
library(phyndr)
library(TaxonLookup)

t <- read.tree("Vascular_Plants_rooted.dated.tre")

plant_lookup[-duplicated(plant_
## build taxonomy from order names
gen <- sapply(t$tip.label, function(x) strsplit(x, split="_")[[1]][1])
## drop tips of genera that don't match
to_drop <- names(gen[-which(gen %in% plant_lookup$genus)])
phy <- drop.tip(t, to_drop)
gen <- gen[-which(names(gen) %in% to_drop)]
tax <- lapply(as.character(gen), function(x) subset(plant_lookup, genus == x))
tax <- do.call(rbind, tax)
rownames(tax) <- phy$tip.label
