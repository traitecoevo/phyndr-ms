## # Matching phylogenetic and comparative data with phyndr and taxonlookup

## *In this vignette, we will be using previously published data sets that are publicly available. If you want to try phyndr out using these data sets, we use the R package [remake](https://github.com/richfitz/remake) to facilitate downloading these directly from the canonical source. To use this feature, you will need to download or clone the [GitHub repository](https://github.com/mwpennell/phyndr-ms), and run the R session from the root of the directory.*


## ## Introduction
## We sought to address a simple (and perhaps a little mundane) but common problem in comparative analyses: a mismatch between species for which there is (phylo)genetic information and those that have been measured for some trait of interest. While a number of data imputation approaches have been suggested to get around this problem, our strategy with phyndr is much simpler: use some external information, such as a topological hypothesis or taxonomic resource, to swap out species in the tree that don't have trait data with "phylogenetically equivalent" species that do.

## We have implemented the taxon swapping algorithms in the `phyndr` R package. To facilitate the use of taxonomic knowledge to make the swaps, at least in land plants, we have built a second R package `taxonlookup` that contains a curated, versioned, and dynamic taxonomic resource. `taxonlookup` is interoperable with `phyndr` but can also be used as a stand-alone tool for a wide variety of ecological and evolutionary applications. (As mentioned above, `taxonlookup` currently only contains a taxonomy of land plants, but the plumbing and infrastructure was designed to be general; if people with taxonomic expertise in other groups would be interested in contributing to this project, we would be thrilled!)


## ## Preliminaries
## Install and load packages

## Get `remake`
## install.packages("devtools")
devtools::install_github("richfitz/remake")
library(remake)

## Install `phyndr`
devtools::install_github("richfitz/phyndr")
library(phyndr)

## N.B. `phyndr` imports `ape` but does not load it into the `NAMESPACE`. Reading, writing, and manipulating `phylo` objects requires loading in ape separately.
library(ape)

## Install `taxonlookup`
devtools::install_github("wcornwell/taxonlookup")
library(taxonlookup)


## ## Using a taxononmy with phyndr

## We will first load in a recent phylogenetic tree from [Magallon et al. 2015](http://onlinelibrary.wiley.com/doi/10.1111/nph.13264/abstract), which includes 798 taxa, sampled from across Angiosperms. Following the terminology of our algorithm, this is the chronogram.

mag_phy <- make("magallon_tree_modified")
mag_phy

## For the trait data set, we will use a large data base of plant growth form (i.e., woody v. herbaceous) compiled by [Zanne et al. 2014](http://www.nature.com/nature/journal/v506/n7486/abs/nature12872.html)

wood_dat <- make("woody_data")
head(wood_dat)

## ### With taxonlookup

## `taxonlookup` is, on the surface, a very simple package. It essentially has two functions. First, `plant_lookup` loads the taxonomy into R.
head(plant_lookup())
## The first time this function is run, it will download and compile the latest version of the taxonomy. This may take some time, but subsequent calls (even in different R sessions) will be essentially instantaenous.

## The second function is `lookup_table`, which compiles a taxonomic resource for a set of species. By default this uses the taxonomy produced by `plant_lookup()` but an alternative taxonomy can be supplied to the argument `lookup_table=`.

## We want to build a taxonomic table that includes all species in the phylogenetic tree **AND** all species in the trait data. If there is no taxonomic information for some species in the data set, this is fine; it only means that no swaps involving this taxa will be permitted.

## Create a vector of taxa for the union of names in the tree and the data
spp <- unique(c(mag_phy$tip.label, rownames(wood_dat)))

## Get a taxonomic lookup table for this vector of names
tax <- lookup_table(spp, by_species=TRUE)
head(tax)
## The `by_species=TRUE` flag (which is not the default) provides a taxonomic classification for each species. This is what is needed for `phyndr`. The default `by_species=FALSE` avoids redundancies and subsets the taxonomic table by higher groups that are represented in the species list. For example,
tax_2 <- lookup_table(spp)
head(tax_2)

## Let's subset the taxonomy so we only consider the genus, family and order columns
tax <- tax[,c("genus", "family", "order")] 

## To run the phyndr algorithm, you need to supply the chronogram, the trait data with rownames set to species names, and the taxonomic table. Since we are using the taxonomic version of the algorithm, we use `phyndr_taxonomy`
mag_phyndr <- phyndr_taxonomy(mag_phy, wood_dat, tax)

## The returned object in `mag_phyndr`
mag_phyndr
## contains a phylogenetic tree in `ape::phylo` format with additional attributes containing the permissable swaps for every taxa. The resulting tree can be of the same size or smaller than the supplied chronogram since no new splits have been induced.
str(mag_phyndr)

## ### With a supplied taxonomy

## If you don't work on land plants, or if you don't like our taxonomy (don't worry, we won't be offended),

## ### By genus name

## ## Using a topology with phyndr

## ### With a supplied topology

## ### With the OpenTree of Life API

## ## Additional features of taxonlookup

## 1. include_counts

## 2. versioning

