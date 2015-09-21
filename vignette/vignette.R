## # Matching phylogenetic and comparative data with phyndr and taxonlookup

## *In this vignette, we will be using previously published data sets that are publicly available. If you want to try phyndr out using these data sets, we use the R package [remake](https://github.com/richfitz/remake) to facilitate downloading these directly from the canonical source. To use this feature, you will need to download or clone the [GitHub repository](https://github.com/mwpennell/phyndr-ms), and run the R session from the root of the directory.*


## ## Introduction
## We sought to address a simple (and perhaps a little mundane) but common problem in comparative analyses: a mismatch between species for which there is (phylo)genetic information and those that have been measured for some trait of interest. While a number of data imputation approaches have been suggested to get around this problem, our strategy with phyndr is much simpler: use some external information, such as a topological hypothesis or taxonomic resource, to swap out species in the tree that don't have trait data with "phylogenetically equivalent" species that do.

## We have implemented the taxon swapping algorithms in the phyndr R package. To facilitate the use of taxonomic knowledge to make the swaps, at least in land plants, we have built a second R package, taxonlookup, that contains a curated, versioned, and dynamic taxonomic resource. taxonlookup is interoperable with phyndr but can also be used as a stand-alone tool for a wide variety of ecological and evolutionary applications. (As mentioned above, taxonlookup currently only contains a taxonomy of land plants, but the plumbing and infrastructure was designed to be general; if people with taxonomic expertise in other groups would be interested in contributing to this project, we would be thrilled!)


## ## Preliminaries
## Install and load packages
library(ape)
library(phyndr)
library(taxonlookup)

## ## Using a taxononmy with phyndr

## We will first load in a recent phylogenetic tree from [Magallon et al. 2015](http://onlinelibrary.wiley.com/doi/10.1111/nph.13264/abstract), which includes 798 taxa, sampled from across Angiosperms. Following the terminology of our algorithm, this is the chronogram.

angio_phy <- read.tree("../source_data/magallon.tre")
angio_phy

## For the trait data set, we will use a large data base of plant growth form (i.e., woody v. herbaceous) compiled by [Zanne et al. 2014](http://www.nature.com/nature/journal/v506/n7486/abs/nature12872.html)

wood_dat <- read.csv("../source_data/wood.csv", row.names=1)
head(wood_dat)

## ### With taxonlookup

## taxonlookup is, on the surface, a very simple package. It essentially has two functions. First, `plant_lookup` loads the taxonomy into R.
head(plant_lookup())
## The first time this function is run, it will download and compile the latest version of the taxonomy. This may take some time, but subsequent calls (even in different R sessions) will be essentially instantaenous (see below).

## The second function is `lookup_table`, which compiles a taxonomic resource for a set of species. By default this uses the taxonomy produced by `plant_lookup()` but an alternative taxonomy can be supplied to the argument `lookup_table=`.

## We want to build a taxonomic table that includes all species in the phylogenetic tree **AND** all species in the trait data. If there is no taxonomic information for some species in the data set, this is fine; it only means that no swaps involving this taxa will be permitted.

## Create a vector of taxa for the union of names in the tree and the data
angio_spp <- unique(c(angio_phy$tip.label, rownames(wood_dat)))

## Get a taxonomic lookup table for this vector of names
angio_tax <- lookup_table(angio_spp, by_species=TRUE)
head(angio_tax)
## The `by_species=TRUE` flag (which is not the default) provides a taxonomic classification for each species. This is what is needed by `phyndr_taxonomy`. The default `by_species=FALSE` avoids redundancies and subsets the taxonomic table by higher groups that are represented in the species list. For example,
angio_tax_alt <- lookup_table(angio_spp)
head(angio_tax_alt)

## Let's subset the taxonomy so we only consider the genus, family and order columns
angio_tax <- angio_tax[,c("genus", "family", "order")]

## To run the phyndr algorithm, we need to supply the chronogram, the trait data with rownames set to species names, and the taxonomic table. Since we are using the taxonomic version of the algorithm, we use `phyndr_taxonomy`
angio_phyndr <- phyndr_taxonomy(angio_phy, rownames(wood_dat), angio_tax)

## The returned object in `angio_phyndr`
angio_phyndr
## contains a phylogenetic tree in `ape::phylo` format with additional attributes containing the permissable swaps for every taxa. The resulting tree can be of the same size or smaller than the supplied chronogram since no new splits have been induced.
str(angio_phyndr)

## ### With a supplied taxonomy

## As mentioned above, the taxonomic resources in taxonlookup are currently only available for land plants. Therefore, if we want to use the taxonomic version of phyndr with some other group of organisms (or you are working with plants and want to use an alternative taxonomy), we need to supply our own taxonomy from somewhere else. One convenient way to obtain a taxonomy is to query online databases. A number of packages have been developed to facilitate this.

## Here we are going to use the R interface to the [Open Tree of Life](http://opentreeoflife.org/) [API](https://github.com/OpenTreeOfLife/opentree/wiki/Open-Tree-of-Life-APIs) [rotl](https://github.com/ropensci/rotl) to obtain a tree of mammals from a study by [Meredith et al. 2011](http://www.sciencemag.org/content/334/6055/521.short).

## To use rotl we first need to install the following packages
library(rotl)

## To pull down the Meredith tree use the function `rotl::get_study_tree`
mamm_phy <- read.tree("../source_data/meredith.tre")

## And we are going to pull down a data set of basal metabolic rate for mammals from a compilation by [McNab 2008](http://www.sciencedirect.com/science/article/pii/S1095643308007782)
bmr_dat <- read.csv("../source_data/bmr.csv", row.names=1)

## To get a taxonomy for this group, we are going to the [taxize](https://github.com/ropensci/taxize) to query the [NCBI Taxonomy Database](http://www.ncbi.nlm.nih.gov/taxonomy)
library(taxize)

## Get a vector of taxa that occur in either the data or the tree
mamm_spp <- unique(c(mamm_phy$tip.label, rownames(bmr_dat)))

## Get NCBI ID numbers for all taxa
ids <- get_uid(mamm_spp, verbose=FALSE)

## Get classification info for all of these
cls <- classification(ids)

## The warnings indicate that some taxa in our data were not found in the taxonomic database. There are lots of reasons this could be the case but for now, we will be conservative and simply disallow any swaps that involve taxa for which which we cannot find taxonomic information.
no_swap <- which(is.na(names(cls)))
mamm_spp <- mamm_spp[-no_swap]
cls[no_swap] <- NULL

## Put the taxonomy together and select the columns of interest
mamm_tax <- cbind(cls)
rownames(mamm_tax) <- mamm_spp
mamm_tax <- mamm_tax[,c("genus", "family")]
head(mamm_tax)

## **Note: phyndr is dumb. It assumes that every column in your taxonomic table corresponds to a taxonomic rank and that these occur in ascending order (i.e., genus -> family -> order, etc.).**

## Now we can use `phyndr_taxonomy` just as before
mamm_phyndr <- phyndr_taxonomy(mamm_phy, rownames(bmr_dat), mamm_tax)
mamm_phyndr

## There are plenty of other databases to query to obtain taxonomies (see the [taxize tutorial](https://ropensci.org/tutorials/taxize_tutorial.html)) and alternative approaches to querying these. For example, one could get the mammal taxonomy directly from the tip labels using the `gbresolve` function in [geiger](https://github.com/mwpennell/geiger-v2) that is usually used as part of the [congruification approach to dating phylogenies](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12051/abstract).

## ### By genus name

## It's even easier if we just want to restrict the swaps to species occuring in the same genus. Assuming that the tip labels and rownames of the data but are the species names in `genus_species` format, we can just skip the taxonomy curation step and use `phyndr_genus`
mamm_phyndr_genus <- phyndr_genus(mamm_phy, rownames(bmr_dat))
mamm_phyndr_genus

## ## Using a topology with phyndr

## If we assume that taxonomies reflect evolutionary relationships (which they should!), then it is fair to consider a taxonomic table as a way of representing a phylogenetic tree, albeit a very unresolved one, with polytomies at each named level. Indeed, using taxize it is possible to construct a tree from the a classification scheme
mamm_class_tree <- class2tree(cls)
mamm_class_tree

## If we have some more detailed information on evolutionary relationships, we could then leverage this in the same way that we did taxonomies (though the algorithm is slightly different; see manuscript text for details). However, we still want to use the tree we have built for the purposes of our analyses presumably because it is more reliable or robust than previous topological hypotheses.

## ### With a supplied topology

## The phylogeny used in [Zanne et al.](http://www.nature.com/nature/journal/v506/n7486/abs/nature12872.html) was much larger the Magallon phylogeny but was constructed using far fewer genes. In this case, we might want to use the Magallon tree for our analysis and the Zanne tree to inform our swaps.

## Load the Zanne tree into the workspace
zae_phy <- read.tree("../source_data/zanne.tre")

## And use the function `phyndr_topology`
angio_zae_phyndr <-  phyndr_topology(angio_phy, rownames(wood_dat), zae_phy)
angio_zae_phyndr

## ## Additional features of taxonlookup

## ### Getting species counts

head(plant_lookup(include_counts = TRUE))

## This feature will include a column with the number of (non-hybrid) species per genus in the lookup table, which is useful for both diversity weighting in graphics and for diversification analyses.

## ### Versioning

## Because taxonomy is (and always has been) a dynamic field, the best available data will always be changing.  taxonlookup is a dynamic resource built on two other dynamic web resources:  1. [The Plant List v1.1.](http://www.theplantlist.org/) for accepted genera to families and species richness within each genera.  Note that we do not consider hybrids (e.g. Genus X species) as distinct species for this count while the plant list summary statistics do, so the the counts from this package will not line up exactly with the ones on the TPL website. And 2. [APWeb](http://www.mobot.org/MOBOT/research/APweb/) for family-level synonymies and family-to-order for all vascular plant families. Note that there is not currently order-level information available for Bryophytes.  When either of these resource changes we will release a new version of the data underlying taxonlookup.  By default these function use the most recent data.  However, for the purposes of reproducability, taxonlookup also makes it easy to access older versions.

head(plant_lookup())

## downloads the most recent version.  You can see the version number with:

pl<-plant_lookup()
plant_lookup_version_current()

## to load a previous version, simply specify it with the function `plant_lookup`

head(plant_lookup(version="0.2.1"))

## this will download the data from the appropriate github release, and should allow for easy access to both new and old versions and allow for easier reproducibility.
