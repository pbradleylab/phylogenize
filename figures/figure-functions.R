last_elem <- function(x) x[length(x)]

make_qvs <- function(tbl,
                     .nestby="phylum",
                     .pvcol="p.value",
                     .qvcol="q.value") {
    .qv <- sym(.qvcol)
    .pv <- sym(.pvcol)
    .nb <- sym(.nestby)
    tbl %>%
        group_by(!!(.nb)) %>%
        nest %>%
        mutate(data=map(data, function(x) {
            mutate(x,
                   !!(.qv) := phylogenize:::qvals(!!(.pv), error_to_file=FALSE))
        })) %>%
        unnest
}

make_equivs <- function(tbl, mfx=0.5) {
    mutate(tbl,
           equiv.pv = pmap_dbl(list(effect.size, std.err, df),
                               phylogenize:::equiv_test,
                               min_fx=mfx))
}

process_data <- function(tbl) {
    tbl %>% make_equivs %>% make_qvs %>%
        make_qvs(., .pvcol="equiv.pv", .qvcol="equiv.qv")
}

get_alpha <- function(phylum, gene, db=pz.db, lab=8) {
                                        #    attr(p, "class") <- "phylo"
    p <- db$trees[[phylum]]
    tips <- intersect(colnames(db$gene.presence[[phylum]]),
                      p$tip.label)
    df <- data.frame(g=db$gene.presence[[phylum]][gene, tips])
    if (var(df$g) == 0) { warning(paste0("zero-variance: ", gene)); return(NA) }
    tryCatch({phyloglm(g ~ 1,
                       data=df,
                       phy=keep.tip(p, tips),
                       log.alpha.bound=lab)$alpha },
             error=function(e) { warning(e); NA })
}

ap_tree <- function(x) {
    attr(x, "class") <- "phylo"
    x
}
