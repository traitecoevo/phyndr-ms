load_magallon_tree <- function(filename)
    read.nexus(filename)


clean_magallon_tiplabel <- function(tree){
    tmp <- tree$tip.label
    for (i in seq_len(length(tmp))){
        gen <- strsplit(tmp[i], split="_")[[1]]
        if (length(gen) == 1)
            tmp[i] <- paste(gen, "spp.", sep="_")

    }
    tree$tip.label <- tmp
    tree
}

write_mag_tree <- function(tree)
  write.tree(tree, "source_data/magallon.tre")

write_woody_dat <- function(wood){
  write.csv(wood, "source_data/wood.csv")
}

build_zae_tree <- function(filename){
    tmp <- tempdir()
    oo <- options(warn=2)
    on.exit(options(oo))
    file <- 'PhylogeneticResources/Vascular_Plants_rooted.dated.tre'
    unzip(filename, file, junkpaths=TRUE, exdir=tmp)
    t <- read.tree(file.path(tmp, basename(file)))
    ## write to file for use in vignette
    write.tree(t, "source_data/zanne.tre")
    t
}

load_woody_data <- function(filename){
    w <- read.csv(filename)
    rownames(w) <- gsub(" ", "_", w$gs)
    w[,c("woodiness", "woodiness.count")]
}

build_bmr_dat <- function(filename){
  b <- read_excel(filename, skip=6)
  bmr <- b[,c("Mass (g)", "BMR (W)")]
  colnames(bmr) <- c("mass", "bmr")
  rownames(bmr) <- sapply(b[,"Genus Species"], function(x) gsub(" ", "_", x))
  ## write to file for use in vignette
  write.csv(bmr, "source_data/bmr.csv")
  bmr
}

build_meredith_tree <- function(filename){
    t <- read.tree(filename)

    ## Fix up names with rotl version because species epithet missing from file
    drop <- c("Leporidae", "Hydrochaeris", "Agouti",
    "Emballonuridae", "Hippopotamidae", "Geomyidae",
    "Gliridae", "Sciuridae")
    t <- phyndr:::drop_tip(t, drop)
    mer <- rotl::get_study_tree(study_id="pg_1428", tree_id="tree2855")
    tmp <- sapply(t$tip.label, function(x) grep(x, mer$tip.label))
    t$tip.label <- mer$tip.label[unlist(tmp)]
    ## write to file for use in vignette
    write.tree(t, "source_data/meredith.tre")
    t
}
