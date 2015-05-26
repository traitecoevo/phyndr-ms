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
    w
}


