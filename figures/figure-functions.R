last_elem <- function(x) x[length(x)]

make_qvs <- function(tbl, .nestby = "phylum", .pvcol = "p.value", .qvcol = "q.value") {
  .qv <- sym(.qvcol)
  .pv <- sym(.pvcol)
  .nb <- sym(.nestby)

  tbl %>%
    group_by(!!(.nb)) %>%
    mutate(across(all_of(.pvcol), ~phylogenize:::qvals(., error_to_file = TRUE), .names = "{.col}_qval")) %>%
    ungroup()
}

make_equivs <- function(tbl, mfx = 0.5) {
  tbl %>%
    rowwise() %>%
    transmute(equiv.pv = phylogenize:::equiv_test(effect.size, std.err, df, min_fx = mfx))
}

process_data <- function(tbl) {
  tbl %>%
    make_equivs() %>%
    make_qvs() %>%
    make_qvs(.pvcol = "equiv.pv", .qvcol = "equiv.qv")
}

get_alpha <- function(phylum, gene, db = pz.db, lab = 8) {
  tryCatch(
    {
      with(
        data.frame(g = db$gene.presence[[phylum]][gene, intersect(colnames(db$gene.presence[[phylum]]), db$trees[[phylum]]$tip.label)]),
        if (var(g) == 0) stop("zero-variance")
        phyloglm(g ~ 1, data = ., phy = keep.tip(db$trees[[phylum]], .), log.alpha.bound = lab)$alpha
      )
    },
    warning = function(w) warning(w),
    error = function(e) { warning(e); NA }
  )
}

ap_tree <- function(x) {
    class(x) <- "phylo"
    return(x)
}