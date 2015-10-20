# Matching phylogenetic and comparative data with phyndr and taxonlookup

*In this vignette, we will be using previously published data sets that are publicly available. If you want to try phyndr out using these data sets, we use the R package [remake](https://github.com/richfitz/remake) to facilitate downloading these directly from the canonical source. To use this feature, you will need to download or clone the [GitHub repository](https://github.com/mwpennell/phyndr-ms), and run the R session from the root of the directory.*


## Introduction
We sought to address a simple (and perhaps a little mundane) but common problem in comparative analyses: a mismatch between species for which there is (phylo)genetic information and those that have been measured for some trait of interest. While a number of data imputation approaches have been suggested to get around this problem, our strategy with phyndr is much simpler: use some external information, such as a topological hypothesis or taxonomic resource, to swap out species in the tree that don't have trait data with "phylogenetically equivalent" species that do.

We have implemented the taxon swapping algorithms in the phyndr R package. To facilitate the use of taxonomic knowledge to make the swaps, at least in land plants, we have built a second R package, taxonlookup, that contains a curated, versioned, and dynamic taxonomic resource. taxonlookup is interoperable with phyndr but can also be used as a stand-alone tool for a wide variety of ecological and evolutionary applications. (As mentioned above, taxonlookup currently only contains a taxonomy of land plants, but the plumbing and infrastructure was designed to be general; if people with taxonomic expertise in other groups would be interested in contributing to this project, we would be thrilled!)


## Preliminaries
Install and load packages

```r
library(ape)
library(phyndr)
library(taxonlookup)
```

## Using a taxononmy with phyndr

We will first load in a recent phylogenetic tree from [Magallon et al. 2015](http://onlinelibrary.wiley.com/doi/10.1111/nph.13264/abstract), which includes 798 taxa, sampled from across Angiosperms. Following the terminology of our algorithm, this is the chronogram.


```r
angio_phy <- read.tree("../source_data/magallon.tre")
angio_phy
```

```
## 
## Phylogenetic tree with 798 tips and 777 internal nodes.
## 
## Tip labels:
## 	Abatia_spp., Prockia_spp., Dovyalis_spp., Flacourtia_spp., Idesia_spp., Populus_spp., ...
## Node labels:
## 	Sper000, Angi001, , , , Eudi026, ...
## 
## Rooted; includes branch lengths.
```

For the trait data set, we will use a large data base of plant growth form (i.e., woody v. herbaceous) compiled by [Zanne et al. 2014](http://www.nature.com/nature/journal/v506/n7486/abs/nature12872.html)


```r
wood_dat <- read.csv("../source_data/wood.csv", row.names=1)
head(wood_dat)
```

```
##                     woodiness woodiness.count
## Aa_paleacea                 H           1;0;0
## Abarema_adenophora          W           0;0;2
## Abarema_barbouriana         W           0;0;2
## Abarema_campestris          W           0;0;1
## Abarema_curvicarpa          W           0;0;1
## Abarema_jupunba             W          0;0;14
```

### With taxonlookup

taxonlookup is, on the surface, a very simple package. It essentially has two functions. First, `plant_lookup` loads the taxonomy into R.

```r
head(plant_lookup())
```

```
##       genus       family       order       group
## 1    Acorus    Acoraceae    Acorales Angiosperms
## 2 Albidella Alismataceae Alismatales Angiosperms
## 3    Alisma Alismataceae Alismatales Angiosperms
## 4   Astonia Alismataceae Alismatales Angiosperms
## 5 Baldellia Alismataceae Alismatales Angiosperms
## 6  Burnatia Alismataceae Alismatales Angiosperms
```
The first time this function is run, it will download and compile the latest version of the taxonomy. This may take some time, but subsequent calls (even in different R sessions) will be essentially instantaenous (see below).

The second function is `lookup_table`, which compiles a taxonomic resource for a set of species. By default this uses the taxonomy produced by `plant_lookup()` but an alternative taxonomy can be supplied to the argument `lookup_table=`.

We want to build a taxonomic table that includes all species in the phylogenetic tree **AND** all species in the trait data. If there is no taxonomic information for some species in the data set, this is fine; it only means that no swaps involving this taxa will be permitted.

Create a vector of taxa for the union of names in the tree and the data

```r
angio_spp <- unique(c(angio_phy$tip.label, rownames(wood_dat)))
```

Get a taxonomic lookup table for this vector of names

```r
angio_tax <- lookup_table(angio_spp, by_species=TRUE)
head(angio_tax)
```

```
##                      genus     family        order       group
## Abatia_spp.         Abatia Salicaceae Malpighiales Angiosperms
## Prockia_spp.       Prockia Salicaceae Malpighiales Angiosperms
## Dovyalis_spp.     Dovyalis Salicaceae Malpighiales Angiosperms
## Flacourtia_spp. Flacourtia Salicaceae Malpighiales Angiosperms
## Idesia_spp.         Idesia Salicaceae Malpighiales Angiosperms
## Populus_spp.       Populus Salicaceae Malpighiales Angiosperms
```
The `by_species=TRUE` flag (which is not the default) provides a taxonomic classification for each species. This is what is needed by `phyndr_taxonomy`. The default `by_species=FALSE` avoids redundancies and subsets the taxonomic table by higher groups that are represented in the species list. For example,

```r
angio_tax_alt <- lookup_table(angio_spp)
head(angio_tax_alt)
```

```
##        genus     family        order       group
## 1     Abatia Salicaceae Malpighiales Angiosperms
## 2    Prockia Salicaceae Malpighiales Angiosperms
## 3   Dovyalis Salicaceae Malpighiales Angiosperms
## 4 Flacourtia Salicaceae Malpighiales Angiosperms
## 5     Idesia Salicaceae Malpighiales Angiosperms
## 6    Populus Salicaceae Malpighiales Angiosperms
```

Let's subset the taxonomy so we only consider the genus, family and order columns

```r
angio_tax <- angio_tax[,c("genus", "family", "order")]
```

To run the phyndr algorithm, we need to supply the chronogram, the trait data with rownames set to species names, and the taxonomic table. Since we are using the taxonomic version of the algorithm, we use `phyndr_taxonomy`

```r
angio_phyndr <- phyndr_taxonomy(angio_phy, rownames(wood_dat), angio_tax)
```

The returned object in `angio_phyndr`

```r
angio_phyndr
```

```
## 
## Phylogenetic tree with 768 tips and 747 internal nodes.
## 
## Tip labels:
## 	genus::Abatia, genus::Prockia, genus::Dovyalis, genus::Flacourtia, genus::Idesia, genus::Populus, ...
## Node labels:
## 	Sper000, Angi001, , , , Eudi026, ...
## 
## Rooted; includes branch lengths.
```
contains a phylogenetic tree in `ape::phylo` format with additional attributes containing the permissable swaps for every taxa. The resulting tree can be of the same size or smaller than the supplied chronogram since no new splits have been induced.

```r
str(angio_phyndr)
```

```
## List of 6
##  $ edge       : int [1:1514, 1:2] 769 770 771 772 773 774 775 776 777 778 ...
##  $ Nnode      : int 747
##  $ tip.label  : chr [1:768] "genus::Abatia" "genus::Prockia" "genus::Dovyalis" "genus::Flacourtia" ...
##  $ edge.length: num [1:1514] 190.65 0.923 1.012 1.748 1.167 ...
##  $ node.label : chr [1:747] "Sper000" "Angi001" "" "" ...
##  $ clades     :List of 530
##   ..$ genus::Abatia           : chr [1:8] "Abatia_americana" "Abatia_boliviana" "Abatia_canescens" "Abatia_mexicana" ...
##   ..$ genus::Prockia          : chr [1:2] "Prockia_crucis" "Prockia_flava"
##   ..$ genus::Dovyalis         : chr [1:7] "Dovyalis_abyssinica" "Dovyalis_caffra" "Dovyalis_hebecarpa" "Dovyalis_longispina" ...
##   ..$ genus::Flacourtia       : chr [1:10] "Flacourtia_degeneri" "Flacourtia_flavescens" "Flacourtia_indica" "Flacourtia_inermis" ...
##   ..$ genus::Idesia           : chr "Idesia_polycarpa"
##   ..$ genus::Populus          : chr [1:48] "Populus_acuminata" "Populus_adenopoda" "Populus_alba" "Populus_angustifolia" ...
##   ..$ genus::Salix            : chr [1:204] "Salix_acmophylla" "Salix_acutifolia" "Salix_alaxensis" "Salix_alba" ...
##   ..$ genus::Poliothyrsis     : chr "Poliothyrsis_sinensis"
##   ..$ genus::Casearia         : chr [1:101] "Casearia_aculeata" "Casearia_acuminata" "Casearia_altiplanensis" "Casearia_arborea" ...
##   ..$ genus::Lunania          : chr [1:9] "Lunania_cubensis" "Lunania_dentata" "Lunania_dodecandra" "Lunania_ekmanii" ...
##   ..$ genus::Lozania          : chr [1:4] "Lozania_klugii" "Lozania_mediana" "Lozania_mutisiana" "Lozania_pittieri"
##   ..$ genus::Goupia           : chr "Goupia_glabra"
##   ..$ genus::Hybanthus        : chr [1:6] "Hybanthus_concolor" "Hybanthus_floribundus" "Hybanthus_monoepetalus" "Hybanthus_prunifolius" ...
##   ..$ genus::Leonia           : chr [1:3] "Leonia_crassa" "Leonia_cymosa" "Leonia_glycycarpa"
##   ..$ genus::Viola            : chr [1:72] "Viola_acuminata" "Viola_adunca" "Viola_albida" "Viola_arcuata" ...
##   ..$ genus::Rinorea          : chr [1:25] "Rinorea_anguifera" "Rinorea_apiculata" "Rinorea_ardisiiflora" "Rinorea_bengalensis" ...
##   ..$ genus::Paropsia         : chr [1:2] "Paropsia_madagascariensis" "Paropsia_vareciformis"
##   ..$ genus::Passiflora       : chr [1:11] "Passiflora_biflora" "Passiflora_coccinea" "Passiflora_cookii" "Passiflora_edulis" ...
##   ..$ genus::Turnera          : chr [1:4] "Turnera_hindsiana" "Turnera_panamensis" "Turnera_subulata" "Turnera_ulmifolia"
##   ..$ genus::Kiggelaria       : chr "Kiggelaria_africana"
##   ..$ genus::Pangium          : chr "Pangium_edule"
##   ..$ genus::Erythrospermum   : chr [1:3] "Erythrospermum_acuminatissimum" "Erythrospermum_candidum" "Erythrospermum_phytolaccoides"
##   ..$ genus::Hydnocarpus      : chr [1:15] "Hydnocarpus_alpina" "Hydnocarpus_castanea" "Hydnocarpus_curtisii" "Hydnocarpus_gracilis" ...
##   ..$ genus::Archytaea        : chr "Archytaea_vahlii"
##   ..$ genus::Bonnetia         : chr [1:5] "Bonnetia_celiae" "Bonnetia_crassa" "Bonnetia_ptariensis" "Bonnetia_roraimae" ...
##   ..$ genus::Calophyllum      : chr [1:56] "Calophyllum_amblyphyllum" "Calophyllum_angulare" "Calophyllum_australianum" "Calophyllum_austrocoriaceum" ...
##   ..$ genus::Mammea           : chr [1:8] "Mammea_africana" "Mammea_americana" "Mammea_bongo" "Mammea_ebboro" ...
##   ..$ genus::Mesua            : chr [1:7] "Mesua_assamica" "Mesua_ferrea" "Mesua_ferruginea" "Mesua_grandis" ...
##   ..$ genus::Vismia           : chr [1:11] "Vismia_baccifera" "Vismia_billbergiana" "Vismia_brasiliensis" "Vismia_cayennensis" ...
##   ..$ genus::Marathrum        : chr [1:2] "Marathrum_oxycarpum" "Marathrum_rubrum"
##   ..$ genus::Garcinia         : chr [1:78] "Garcinia_acuminata" "Garcinia_adinantha" "Garcinia_afzelii" "Garcinia_archboldiana" ...
##   ..$ genus::Irvingia         : chr [1:4] "Irvingia_gabonensis" "Irvingia_grandifolia" "Irvingia_malayana" "Irvingia_robur"
##   ..$ genus::Klainedoxa       : chr [1:2] "Klainedoxa_gabonensis" "Klainedoxa_latifolia"
##   ..$ genus::Ochthocosmus     : chr [1:4] "Ochthocosmus_barrae" "Ochthocosmus_congolensis" "Ochthocosmus_longipedicellatus" "Ochthocosmus_sessiliflorus"
##   ..$ genus::Spathiostemon    : chr "Spathiostemon_javensis"
##   ..$ genus::Ricinus          : chr "Ricinus_communis"
##   ..$ genus::Dalechampia      : chr [1:2] "Dalechampia_scandens" "Dalechampia_spathulata"
##   ..$ genus::Lasiocroton      : chr [1:3] "Lasiocroton_bahamensis" "Lasiocroton_macrophyllus" "Lasiocroton_trelawniensis"
##   ..$ genus::Conceveiba       : chr [1:4] "Conceveiba_guianensis" "Conceveiba_martiana" "Conceveiba_pleiostemona" "Conceveiba_rhytidocarpa"
##   ..$ genus::Euphorbia        : chr [1:137] "Euphorbia_abyssinica" "Euphorbia_acanthothamnos" "Euphorbia_adiantoides" "Euphorbia_akenocarpa" ...
##   ..$ genus::Homalanthus      : chr [1:6] "Homalanthus_novoguineensis" "Homalanthus_nutans" "Homalanthus_populifolius" "Homalanthus_populneus" ...
##   ..$ genus::Hura             : chr [1:2] "Hura_crepitans" "Hura_polyandra"
##   ..$ genus::Pimelodendron    : chr [1:3] "Pimelodendron_amboinicum" "Pimelodendron_griffithianum" "Pimelodendron_zoanthogyne"
##   ..$ genus::Codiaeum         : chr [1:3] "Codiaeum_ludovicianum" "Codiaeum_peltatum" "Codiaeum_variegatum"
##   ..$ genus::Trigonostemon    : chr "Trigonostemon_verrucosus"
##   ..$ genus::Croton           : chr [1:107] "Croton_actinophyllum" "Croton_adspersus" "Croton_alabamensis" "Croton_antanosiensis" ...
##   ..$ genus::Hevea            : chr [1:7] "Hevea_brasiliensis" "Hevea_guianensis" "Hevea_lutea" "Hevea_microphylla" ...
##   ..$ genus::Moultonianthus   : chr "Moultonianthus_leembruggianus"
##   ..$ genus::Suregada         : chr [1:6] "Suregada_aequoreum" "Suregada_boiviniana" "Suregada_glomerulata" "Suregada_multiflora" ...
##   ..$ genus::Endospermum      : chr [1:9] "Endospermum_banghamii" "Endospermum_chinense" "Endospermum_diadenum" "Endospermum_macrophyllum" ...
##   ..$ genus::Omphalea         : chr [1:4] "Omphalea_diandra" "Omphalea_oleifera" "Omphalea_queenslandiae" "Omphalea_triandra"
##   ..$ genus::Tetrorchidium    : chr [1:5] "Tetrorchidium_didymostemon" "Tetrorchidium_gabonense" "Tetrorchidium_macrophyllum" "Tetrorchidium_rotundatum" ...
##   ..$ genus::Clutia           : chr [1:6] "Clutia_affinis" "Clutia_augustifolia" "Clutia_brassii" "Clutia_paxii" ...
##   ..$ genus::Cespedesia       : chr [1:2] "Cespedesia_macrophylla" "Cespedesia_spathulata"
##   ..$ genus::Sauvagesia       : chr [1:6] "Sauvagesia_africana" "Sauvagesia_angustifolia" "Sauvagesia_calophylla" "Sauvagesia_elata" ...
##   ..$ genus::Ochna            : chr [1:16] "Ochna_arborea" "Ochna_ciliata" "Ochna_gambleoides" "Ochna_holstii" ...
##   ..$ genus::Luxemburgia      : chr "Luxemburgia_ciliosa"
##   ..$ genus::Quiina           : chr [1:11] "Quiina_albiflora" "Quiina_amazonica" "Quiina_florida" "Quiina_guianensis" ...
##   ..$ genus::Touroulia        : chr "Touroulia_guianensis"
##   ..$ genus::Caryocar         : chr [1:12] "Caryocar_amygdaliferum" "Caryocar_amygdaliforme" "Caryocar_barbinerve" "Caryocar_brasiliense" ...
##   ..$ genus::Erythroxylum     : chr [1:72] "Erythroxylum_amazonicum" "Erythroxylum_amplifolium" "Erythroxylum_amplum" "Erythroxylum_ampullaceum" ...
##   ..$ genus::Bruguiera        : chr [1:8] "Bruguiera_caryophylloides" "Bruguiera_conjugata" "Bruguiera_cylindrica" "Bruguiera_exaristata" ...
##   ..$ genus::Carallia         : chr [1:4] "Carallia_borneensis" "Carallia_brachiata" "Carallia_euryoides" "Carallia_integerrima"
##   ..$ genus::Rhizophora       : chr [1:7] "Rhizophora_apiculata" "Rhizophora_candelaria" "Rhizophora_harrisonii" "Rhizophora_mangle" ...
##   ..$ genus::Crossostylis     : chr [1:7] "Crossostylis_biflora" "Crossostylis_harveyi" "Crossostylis_multiflora" "Crossostylis_parksii" ...
##   ..$ genus::Cassipourea      : chr [1:12] "Cassipourea_congensis" "Cassipourea_elliottii" "Cassipourea_elliptica" "Cassipourea_euryoides" ...
##   ..$ genus::Paradrypetes     : chr "Paradrypetes_subintegrifolia"
##   ..$ genus::Ctenolophon      : chr [1:2] "Ctenolophon_englerianus" "Ctenolophon_parvifolius"
##   ..$ genus::Atuna            : chr [1:2] "Atuna_elliptica" "Atuna_racemosa"
##   ..$ genus::Hirtella         : chr [1:25] "Hirtella_adenophora" "Hirtella_americana" "Hirtella_bicornis" "Hirtella_bullata" ...
##   ..$ genus::Licania          : chr [1:96] "Licania_alba" "Licania_albicans" "Licania_apetala" "Licania_arborea" ...
##   ..$ genus::Euphronia        : chr [1:2] "Euphronia_guianensis" "Euphronia_hirtelloides"
##   ..$ genus::Dichapetalum     : chr [1:27] "Dichapetalum_angolense" "Dichapetalum_axillare" "Dichapetalum_barteri" "Dichapetalum_bellum" ...
##   ..$ genus::Tapura           : chr [1:9] "Tapura_acreana" "Tapura_bouquetiana" "Tapura_capitulifera" "Tapura_colombiana" ...
##   ..$ genus::Trigonia         : chr [1:5] "Trigonia_boliviana" "Trigonia_floribunda" "Trigonia_nivea" "Trigonia_rugosa" ...
##   ..$ genus::Balanops         : chr [1:4] "Balanops_australiana" "Balanops_pancheri" "Balanops_pedicellata" "Balanops_vieillardii"
##   ..$ genus::Hugonia          : chr [1:2] "Hugonia_castaneifolia" "Hugonia_platysepala"
##   ..$ genus::Linum            : chr [1:8] "Linum_arboreum" "Linum_catharticum" "Linum_marginale" "Linum_perenne" ...
##   ..$ genus::Galearia         : chr [1:3] "Galearia_aristifera" "Galearia_celebica" "Galearia_filiformis"
##   ..$ genus::Panda            : chr "Panda_oleosa"
##   ..$ genus::Microdesmis      : chr [1:5] "Microdesmis_afrodecandra" "Microdesmis_caseariifolia" "Microdesmis_keayana" "Microdesmis_pierlotiana" ...
##   ..$ genus::Drypetes         : chr [1:60] "Drypetes_alba" "Drypetes_amazonica" "Drypetes_americana" "Drypetes_angustifolia" ...
##   ..$ genus::Dicella          : chr "Dicella_nucifera"
##   ..$ genus::Malpighia        : chr [1:11] "Malpighia_albiflora" "Malpighia_coccigera" "Malpighia_emarginata" "Malpighia_glabra" ...
##   ..$ genus::Thryallis        : chr [1:2] "Thryallis_latifolia" "Thryallis_longifolia"
##   ..$ genus::Tetrapterys      : chr [1:16] "Tetrapterys_acutifolia" "Tetrapterys_chloroptera" "Tetrapterys_crispa" "Tetrapterys_discolor" ...
##   ..$ genus::Bergia           : chr [1:3] "Bergia_ammannioides" "Bergia_glomerata" "Bergia_texana"
##   ..$ genus::Elatine          : chr [1:3] "Elatine_hexandra" "Elatine_hydropiper" "Elatine_triandra"
##   ..$ genus::Tetracoccus      : chr [1:3] "Tetracoccus_dioicus" "Tetracoccus_hallii" "Tetracoccus_ilicifolius"
##   ..$ genus::Austrobuxus      : chr [1:9] "Austrobuxus_carunculatus" "Austrobuxus_clusiaceus" "Austrobuxus_huerlimannii" "Austrobuxus_megacarpus" ...
##   ..$ genus::Micrantheum      : chr [1:3] "Micrantheum_demissum" "Micrantheum_ericoides" "Micrantheum_hexandrum"
##   ..$ genus::Dissiliaria      : chr [1:3] "Dissiliaria_baloghioides" "Dissiliaria_laxinervis" "Dissiliaria_muelleri"
##   ..$ genus::Petalostigma     : chr [1:5] "Petalostigma_banksii" "Petalostigma_pachyphyllum" "Petalostigma_pubescens" "Petalostigma_quadriloculare" ...
##   ..$ genus::Podocalyx        : chr "Podocalyx_loranthoides"
##   ..$ genus::Lachnostylis     : chr "Lachnostylis_bilocularis"
##   ..$ genus::Heywoodia        : chr "Heywoodia_lucens"
##   ..$ genus::Phyllanthus      : chr [1:87] "Phyllanthus_abnormis" "Phyllanthus_acidus" "Phyllanthus_acuminatus" "Phyllanthus_amarus" ...
##   ..$ genus::Humiria          : chr [1:2] "Humiria_balsamifera" "Humiria_procera"
##   ..$ genus::Vantanea         : chr [1:10] "Vantanea_barbourii" "Vantanea_compacta" "Vantanea_cupularis" "Vantanea_depleta" ...
##   .. [list output truncated]
##  - attr(*, "class")= chr [1:2] "phyndr" "phylo"
##  - attr(*, "order")= chr "cladewise"
```

### With a supplied taxonomy

As mentioned above, the taxonomic resources in taxonlookup are currently only available for land plants. Therefore, if we want to use the taxonomic version of phyndr with some other group of organisms (or you are working with plants and want to use an alternative taxonomy), we need to supply our own taxonomy from somewhere else. One convenient way to obtain a taxonomy is to query online databases. A number of packages have been developed to facilitate this.

For this example, we are going to use a tree of mammals from a study by [Meredith et al. 2011](http://www.sciencemag.org/content/334/6055/521.short).

```r
mamm_phy <- read.tree("../source_data/meredith.tre")
```

And we are going to pull down a data set of basal metabolic rate for mammals from a compilation by [McNab 2008](http://www.sciencedirect.com/science/article/pii/S1095643308007782)

```r
bmr_dat <- read.csv("../source_data/bmr.csv", row.names=1)
```

To get a taxonomy for this group, we are going to the [taxize](https://github.com/ropensci/taxize) to query the [NCBI Taxonomy Database](http://www.ncbi.nlm.nih.gov/taxonomy)

```r
library(taxize)
```

Get a vector of taxa that occur in either the data or the tree

```r
mamm_spp <- unique(c(mamm_phy$tip.label, rownames(bmr_dat)))
```

Get NCBI ID numbers for all taxa

```r
ids <- get_uid(mamm_spp, verbose=FALSE)
```

Get classification info for all of these

```r
cls <- classification(ids)
```

The warnings indicate that some taxa in our data were not found in the taxonomic database. There are lots of reasons this could be the case but for now, we will be conservative and simply disallow any swaps that involve taxa for which which we cannot find taxonomic information.

```r
no_swap <- which(is.na(names(cls)))
mamm_spp <- mamm_spp[-no_swap]
cls[no_swap] <- NULL
```

Put the taxonomy together and select the columns of interest

```r
mamm_tax <- cbind(cls)
rownames(mamm_tax) <- mamm_spp
mamm_tax <- mamm_tax[,c("genus", "family")]
head(mamm_tax)
```

```
##                                    genus            family
## Tachyglossus_aculeatus      Tachyglossus    Tachyglossidae
## Ornithorhynchus_anatinus Ornithorhynchus Ornithorhynchidae
## Caenolestes_fuliginosus      Caenolestes     Caenolestidae
## Didelphis_virginiana           Didelphis       Didelphidae
## Dromiciops_gliroides          Dromiciops Microbiotheriidae
## Notoryctes_typhlops           Notoryctes      Notoryctidae
```

**Note: phyndr is dumb. It assumes that every column in your taxonomic table corresponds to a taxonomic rank and that these occur in ascending order (i.e., genus -> family -> order, etc.).**

Now we can use `phyndr_taxonomy` just as before

```r
mamm_phyndr <- phyndr_taxonomy(mamm_phy, rownames(bmr_dat), mamm_tax)
mamm_phyndr
```

```
## 
## Phylogenetic tree with 82 tips and 81 internal nodes.
## 
## Tip labels:
## 	Tachyglossus_aculeatus, Ornithorhynchus_anatinus, Didelphis_virginiana, Dromiciops_gliroides, genus::Notoryctes, genus::Dasyurus, ...
## 
## Rooted; includes branch lengths.
```

There are plenty of other databases to query to obtain taxonomies (see the [taxize tutorial](https://ropensci.org/tutorials/taxize_tutorial.html)) and alternative approaches to querying these. For example, one could get the mammal taxonomy directly from the tip labels using the `gbresolve` function in [geiger](https://github.com/mwpennell/geiger-v2) that is usually used as part of the [congruification approach to dating phylogenies](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12051/abstract).

### By genus name

It's even easier if we just want to restrict the swaps to species occuring in the same genus. Assuming that the tip labels and rownames of the data but are the species names in `genus_species` format, we can just skip the taxonomy curation step and use `phyndr_genus`

```r
mamm_phyndr_genus <- phyndr_genus(mamm_phy, rownames(bmr_dat))
mamm_phyndr_genus
```

```
## 
## Phylogenetic tree with 73 tips and 72 internal nodes.
## 
## Tip labels:
## 	Tachyglossus_aculeatus, Ornithorhynchus_anatinus, Didelphis_virginiana, Dromiciops_gliroides, genus::Notoryctes, genus::Dasyurus, ...
## 
## Rooted; includes branch lengths.
```

## Using a topology with phyndr

If we assume that taxonomies reflect evolutionary relationships (which they should!), then it is fair to consider a taxonomic table as a way of representing a phylogenetic tree, albeit a very unresolved one, with polytomies at each named level. Indeed, using taxize it is possible to construct a tree from the a classification scheme

```r
mamm_class_tree <- class2tree(cls)
mamm_class_tree
```

```
## 
## Phylogenetic tree with 681 tips and 680 internal nodes.
## 
## Tip labels:
## 	Tachyglossus aculeatus, Ornithorhynchus anatinus, Caenolestes fuliginosus, Didelphis virginiana, Dromiciops gliroides, Notoryctes typhlops, ...
## 
## Rooted; includes branch lengths.
```

If we have some more detailed information on evolutionary relationships, we could then leverage this in the same way that we did taxonomies (though the algorithm is slightly different; see manuscript text for details). However, we still want to use the tree we have built for the purposes of our analyses presumably because it is more reliable or robust than previous topological hypotheses.

### With a supplied topology

The phylogeny used in [Zanne et al.](http://www.nature.com/nature/journal/v506/n7486/abs/nature12872.html) was much larger the Magallon phylogeny but was constructed using far fewer genes. In this case, we might want to use the Magallon tree for our analysis and the Zanne tree to inform our swaps.

Load the Zanne tree into the workspace

```r
zae_phy <- read.tree("../source_data/zanne.tre")
```

And use the function `phyndr_topology`

```r
angio_zae_phyndr <-  phyndr_topology(angio_phy, rownames(wood_dat), zae_phy)
angio_zae_phyndr
```

```
## 
## Phylogenetic tree with 244 tips and 235 internal nodes.
## 
## Tip labels:
## 	Lacistema_aggregatum, Rinorea_pubiflora, Malesherbia_linearifolia, Acharia_tragodes, Cratoxylum_cochinchinense, Hypericum_perforatum, ...
## Node labels:
## 	Sper000, Angi001, , , , Eudi026, ...
## 
## Rooted; includes branch lengths.
```

## Additional features of taxonlookup

### Getting species counts


```r
head(plant_lookup(include_counts = TRUE))
```

```
##   number.of.species     genus       family       order       group
## 1                 2    Acorus    Acoraceae    Acorales Angiosperms
## 2                 1 Albidella Alismataceae Alismatales Angiosperms
## 3                 8    Alisma Alismataceae Alismatales Angiosperms
## 4                 1   Astonia Alismataceae Alismatales Angiosperms
## 5                 3 Baldellia Alismataceae Alismatales Angiosperms
## 6                 1  Burnatia Alismataceae Alismatales Angiosperms
```

This feature will include a column with the number of (non-hybrid) species per genus in the lookup table, which is useful for both diversity weighting in graphics and for diversification analyses.

### Versioning

Because taxonomy is (and always has been) a dynamic field, the best available data will always be changing.  taxonlookup is a dynamic resource built on two other dynamic web resources:  1. [The Plant List v1.1.](http://www.theplantlist.org/) for accepted genera to families and species richness within each genera.  Note that we do not consider hybrids (e.g. Genus X species) as distinct species for this count while the plant list summary statistics do, so the the counts from this package will not line up exactly with the ones on the TPL website. And 2. [APWeb](http://www.mobot.org/MOBOT/research/APweb/) for family-level synonymies and family-to-order for all vascular plant families. Note that there is not currently order-level information available for Bryophytes.  When either of these resource changes we will release a new version of the data underlying taxonlookup.  By default these function use the most recent data.  However, for the purposes of reproducability, taxonlookup also makes it easy to access older versions.


```r
head(plant_lookup())
```

```
##       genus       family       order       group
## 1    Acorus    Acoraceae    Acorales Angiosperms
## 2 Albidella Alismataceae Alismatales Angiosperms
## 3    Alisma Alismataceae Alismatales Angiosperms
## 4   Astonia Alismataceae Alismatales Angiosperms
## 5 Baldellia Alismataceae Alismatales Angiosperms
## 6  Burnatia Alismataceae Alismatales Angiosperms
```

downloads the most recent version.  You can see the version number with:


```r
pl<-plant_lookup()
plant_lookup_version_current()
```

```
## [1] "0.3.1"
```

to load a previous version, simply specify it with the function `plant_lookup`


```r
head(plant_lookup(version="0.2.1"))
```

```
##       genus       order       group
## 1    Acorus    Acorales Angiosperms
## 2 Albidella Alismatales Angiosperms
## 3    Alisma Alismatales Angiosperms
## 4   Astonia Alismatales Angiosperms
## 5 Baldellia Alismatales Angiosperms
## 6  Burnatia Alismatales Angiosperms
```

this will download the data from the appropriate github release, and should allow for easy access to both new and old versions and allow for easier reproducibility.
