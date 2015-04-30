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



