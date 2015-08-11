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

load_zae_tree <- function(filename){
    tmp <- tempdir()
    oo <- options(warn=2)
    on.exit(options(oo))
    file <- 'PhylogeneticResources/Vascular_Plants_rooted.dated.tre'
    unzip(filename, file, junkpaths=TRUE, exdir=tmp)
    read.tree(file.path(tmp, basename(file)))
}

load_woody_data <- function(filename){
    w <- read.csv(filename)
    rownames(w) <- gsub(" ", "_", w$gs)
    w[,c("woodiness", "woodiness.count")]
}

bmr_dat_xlsx <- function(filename){
  b <- read_excel(filename, skip=6)
  rownames(b) <- sapply(b[,"Genus Species"], function(x) gsub(" ", "_", x))
  b <- b[,c("Mass (g)", "BMR (W)")]
  colnames(b) <- c("mass", "bmr")
  b
}

load_meredith_tree <- function(filename){
    t <- read.tree(filename)

    ## Fix up names with rotl version because species epithet missing from file
    drop <- c("Leporidae", "Hydrochaeris", "Agouti")
    t <- phyndr:::drop_tip(t, drop)
    mer <- rotl::get_study_tree(study_id="pg_1428", tree_id="tree2855")
    tmp <- sapply(t$tip.label, function(x) grep(x, mer$tip.label))
    t$tip.label <- mer$tip.label[unlist(tmp)]
    t
}
