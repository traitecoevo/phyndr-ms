phyndr_taxonomy_wrap <- function(tree, dat){
    tax <- lookup_table(unique(c(tree$tip.label, rownames(dat))),
                        by_species=TRUE)
    tax <- tax[,c("genus", "family", "order")]
    new <- phyndr_taxonomy(tree, rownames(dat), tax)
    list(new_tree=new, old_tree=tree)
}


fig_trees <- function(res_1, res_2){
    cols <- tree_cols()
    par(mfrow=c(1,2))
    fig_tree_single(res_1, cols)
    fig_tree_single(res_2, cols)
}

tree_cols <- function(){
    cols_flatui <- flatui_colour_scheme()
    list(swap = c(cols_flatui$blue,cols_flatui$orange))
}

fig_tree_single <- function(res, cols){
    old <- res$old_tree
    new <- res$new_tree
    idx  <- which(new$tip.label %in% old$tip.label) 
    same <- new$tip.label[idx]
    swap <- new$tip.label[-idx]
    labs <- c(rep(0, length(same)), rep(1,length(swap)))
    dat  <- data.frame(swap=labs)
    rownames(dat) <- c(same, swap)
    dat  <- dat[new$tip.label, 1, drop=FALSE]
    diversitree::trait.plot(new, dat=dat, cols=cols,
                            w=1/10, cex.lab=0.001, col.lab="white")
}
