library("data.table")
library("parallel")
library("phylolm")
library("phyloseq")
library("pbapply")
library("qvalue")
library("magrittr")
library("ggtree")
library("ape")
library("fitdistrplus")
library("ggtree")
library("phytools")
library("knitr")
library("seqinr")
library("kableExtra")
library("scales")
library("xml2")

# helper functions

`%btwn%` <- function(x, y) { (x > min(y)) & (x < max(y)) }
`%btwn.inc%` <- function(x, y) { (x >= min(y)) & (x <= max(y)) }
`%intr%` <- function(x,y) intersect(x,y)
min.nonzero <- function(x) min(x[x > 0])

grepv <- function(x, y, ...) grep(x, y, ..., value = TRUE)
get.from <- function(x, y) { sapply(x, function(z) z[[y]]) }
select.which <- function(x,y,z) x[which(x[[y]]==z), ]
select.return <- function(x,y,z,w) select.which(x,y,z)[,w]
subset.return <- function(x, id, y, w) x[which(x[[id]] %in% intersect(x[[id]], y)), w]
get.runs <- function(id, value, m = metadata, rts = run.to.subject) {
  s <- select.return(m, id, value, "subject_id")
  r <- subset.return(rts, "subject_id", s, "run_accession")
  r
}
grep.runs <- function(id, x, m = metadata, rts = run.to.subject) {
  s <- m[grep(x, m[[id]]), "subject_id"]
  r <- subset.return(rts, "subject_id", s, "run_accession")
  r
}
`%withnames%` <- function(x, y) { names(x) <- y; x }
nw <- function(x) { names(which(x)) }
count.each <- function(x, na.rm = FALSE) {
  u <- unique(x)
  simplify2array(lapply.across.names(u, function(y) sum(x == y, na.rm = na.rm)))
}

fastread <- function(location) {
  # rownames are useful
  master <- data.frame(fread(location, header = T))
  rn <- master[,1]
  rest <- master[,-1]
  rownames(rest) <- rn
  return(data.matrix(rest))
}
lapply.across.names <- function(X, FUN, ...) {
    r <- lapply(X, FUN, ...)
  names(r) <- X
    r
}
pblapply.across.names <- function(X, FUN, ...) {
    r <- pblapply(X, FUN, ...)
  names(r) <- X
    r
}
small.merge <- function(x, y, ...) {
  z <- merge(x, y, by = 0)
  r <- z[, -1]
  rownames(r) <- z[, 1]
  r
}

logit <- function(x) (log(x / (1 - x)))
logistic <- function(x) exp(x) / (1 + exp(x))


qtruncated <- function(x, lim = c(0.1, 0.9)) {
  # truncate a vector by quantiles
  qlim <- quantile(na.omit(x), lim)
  y <- x
  y[which(x < qlim[1])] <- qlim[1]
  y[which(x > qlim[2])] <- qlim[2]
  y
}

truncated <- function(x, lim = c(logit(0.001), logit(0.25))) {
  # truncate a vector by values
  y <- x
  y[which(x < lim[1])] <- lim[1]
  y[which(x > lim[2])] <- lim[2]
  y
}

average.runs.subject <- function(phenos, mtx = species.relabs) {
  accessions.per.sub <- lapply(phenos$subject_id %>% unique, function(x) (phenos %index% (subject_id == x))$run_accession)
  sapply(accessions.per.sub, function(x) rowMeans(mtx[, intersect(colnames(mtx), x)]))
}

`%index%` <- function(mtx, predicate) { index.by(mtx, predicate) }
`%subset%` <- function(mtx, predicate) { substitute(quote(subset(mtx, predicate))) }
index.by <- function(mtx, predicate) { attach(mtx, warn.con = FALSE); mtx[which(with(mtx, eval(quote(predicate)))),] }

# Annotation (FIGfams and taxa)

fig.annot.slow <- function(x) {
  sapply(x, function(y) { 
    r <- figfam.to.fxn[which(figfam.to.fxn$figfam == y), ]$"function"
    if (length(r) == 0) { return(NA) } else { return(r) }
  })
}

gene.annot <- function(x) {
  merge(data.table(x),
        gene.to.fxn,
        by.x = colnames(data.table(x))[1],
        by.y = "gene",
        all = T)[x]$"function" %withnames% x
}

tax.annot.slow <- function(tns) {
  sapply(tns, function(tn) taxonomy$species[which(taxonomy$cluster == tn)[1]])
}

tax.annot <- function(tns, taxonomy) {
  dtt <- data.table(tns)
  merge(dtt, taxonomy, by.x = colnames(dtt)[1], by.y = "cluster", all = T)[tns]$"species" %>% as.character %withnames% tns
}

# Significance testing

make.pos.sig <- function(sigs, signs, cut = "strong") {
  lapply.across.names(names(sigs), function(x) intersect(sigs[[x]][[cut]], nw(signs[[x]] > 0)))
}
qvals <- function(x) qvalue(na.omit(x), pfdr = T, pi0.method = "bootstrap")$qvalues

fdr.bh <- function(x) p.adjust(x, "BH")
fdr.by <- function(x) p.adjust(x, "BY")

make.sigs <- function(results, cuts = c(strong = 0.05, med = 0.1, weak = 0.25), method = qvals, 
  exclude = NULL) {
  lapply.across.names(names(results), function(x) { 
    lapply(cuts, function(cut) {
      if (!is.null(exclude)) {
      valid <- setdiff(colnames(results[[x]]), 
        exclude[[x]])
      } else { valid <- colnames(results[[x]]) }
      tested <- na.omit(results[[x]][2, valid])
      tryCatch(nw(method(tested) <= cut), error = function(e) character(0))
      })
  })
}

make.signs <- function(results) {
  lapply(results, function(r) {
    sign(r[1, ]) %>% na.omit
  })
}


# Misc additional power functions (see also plm-power-functions.R)

calc.alpha.power <- function(pvs, null, alt, alpha = 0.05, filter = NULL) {
  reject = which(pvs <= alpha)
  if (!is.null(filter)) { alt <- intersect(alt, filter); null <- intersect(null, filter); reject <- intersect(reject, filter) }   # allows you to look at e.g. just predicted - fx
  pwr = length(intersect(reject, alt)) / length(alt)
  a = length(intersect(reject, null)) / length(null)
  return(c(r = length(reject) / length(pvs), p = pwr, a = a))
}

# Tree functions

fix.tree <- function(phy) {
  phy <- multi2di(phy)
  phy$edge.length[phy$edge.length==0] <- max(nodeHeights(phy)) * 1e-6 # from Liam Reveil's blog June 23 2015
  phy
}

keep.tips <- function(tree, keeplist) {
  drop.tip(tree, setdiff(tree$tip.label, keeplist))
}

tipToRoot <- function(phy) phy %>% vcv.phylo %>% diag

# Enrichment 

dirxn.enrich <- function(sigs, signs, results, siglevel = "strong", exclude = NULL, dirxn = 1, mapping, ...) {
  lapply(
  c(l1 = "level1", l2 = "level2", l3 = "level3"), function(l) {
    # dirxn = 1 is POSITIVE
    message(l)
    enrichments <- lapply.across.names(names(sigs), function(x) {
      #print(x)
      #print(nw(signs[[x]] > 0))
      #print(intersect(sigs[[x]][[siglevel]], nw(signs[[x]] < 0)))
      enr.wrapper(intersect(sigs[[x]][[siglevel]], nw((dirxn * signs[[x]]) > 0)),
          tested = setdiff(names(na.omit(results[[x]][1,])), exclude[[x]]),
        mapping = mapping[[l]], ...)
    })
  return(enrichments)
  })
}

make.enr.tables <- function(X) { 
  Reduce(rbind, Reduce(c, 
    lapply.across.names(X %>% names, function(x) {
      lapply.across.names(X[[x]] %>% names, function(y) {
        if (!is.null(dim(X[[x]][[y]]$table))) {
          cbind(x, y, X[[x]][[y]]$table)
        } else {
          cbind(x, y, "(none)", "1")
        }
      })}))
    )
}

make.enr.tables.all4 <- function(X) { 
  Reduce(rbind, Reduce(rbind, Reduce(c, 
    lapply.across.names(X %>% names, function(x) {
      lapply.across.names(X[[x]] %>% names, function(y) {
        lapply.across.names(X[[x]][[y]] %>% names, function(z) {
        if (!is.null(dim(X[[x]][[y]][[z]]$table))) {
          cbind(x, y, z, X[[x]][[y]][[z]]$table)
        } else {
          cbind(x, y, z, "(none)", "1")
        }
      })})})))
    )
}



make.plm.enr.tables <- function(enr) {
  Reduce(rbind, zipSapply(enr, function(qcut) {
    cbind(qcut$name, Reduce(rbind, zipSapply(qcut$data, function(level) {
      cbind(level$name, Reduce(rbind, zipSapply(level$data, function(clade) {
        #    print(clade$data$table)
        tb <- data.matrix(clade$data$table)
        if (dim(tb)[1] == 0) { return(matrix(c(clade$name, Enriched = "--", pval = "--"), nr = 1)) }
        if (dim(tb)[2] == 1) { return(cbind(clade$name, t(tb))) }
        cbind(clade$name, tb)
    })
    ))})
  ))})
    )
}
make.sg.enr.tables <- function(enr) {
    Reduce(rbind, zipSapply(enr, function(level) {
      cbind(level$name, Reduce(rbind, zipSapply(level$data, function(clade) {
        #    print(clade$data$table)
        tb <- data.matrix(clade$data$table)
        if (dim(tb)[1] == 0) { return(matrix(c(clade$name, Enriched = "--", pval = "--"), nr = 1)) }
        if (dim(tb)[2] == 1) { return(cbind(clade$name, t(tb))) }
        cbind(clade$name, tb)
    })
    ))})
  )}

perform.enrichments <- function(sigs, signs, results, mapping, dirxn = 1, exclude = NULL, ...) {
  lapply.across.names(c("strong", "med", "weak"), function(x) {
    dirxn.enrich(sigs, signs, results, siglevel = x,
      exclude = exclude, fdr.method = fdr.bh, dirxn = dirxn,
      mapping = mapping, ...)
  })
}


# Parallelization

mcapply <- function(X, MARGIN, FUN, mc.cores = 10, simplify = TRUE, ...) {
  if (MARGIN == 1) {
    mlist <- lapply(seq_len(nrow(X)), function(i) X[i, ])
    names(mlist) <- rownames(X)
  } else if (MARGIN == 2) {
    mlist <- lapply(seq_len(ncol(X)), function(i) X[, i])
    names(mlist) <- colnames(X)
  } else {
    stop("invalid MARGIN value")
  }
  r <- mclapply(mlist, FUN, mc.cores = mc.cores, ...) 
  if (simplify) simplify2array(r) else r
}

pbmcapply <- function(X, MARGIN, FUN, mc.cores = 10, simplify = TRUE, ...) {
  if (MARGIN == 1) {
    mlist <- lapply(seq_len(nrow(X)), function(i) X[i, ])
    names(mlist) <- rownames(X)
  } else if (MARGIN == 2) {
    mlist <- lapply(seq_len(ncol(X)), function(i) X[, i])
    names(mlist) <- colnames(X)
  } else {
    stop("invalid MARGIN value")
  }
  r <- pbmclapply(mlist, FUN, mc.cores = mc.cores, ...) 
  if (simplify) simplify2array(r) else r
}

# NCBI handling

process.ncbitaxa <- function(query) {
  taxa <- lapply(XML::getNodeSet(query, "/TaxaSet/Taxon"), XML::xmlToList) 
  lapply(taxa, function(tax) {
    lin <- tax$LineageEx
    #sn <- lin$Taxon$ScientificName
    sn <- tax$ScientificName
    on <- tax$OtherNames
    on.flat <- sapply(1:length(on), function(nx) {
      xn <- names(on)[nx]
      if (is.null(xn)) { return(NA) }
      if (xn %in% c("Synonym", "EquivalentName")) { return(on[[nx]]) }
      if (xn %in% c("Name")) { return(on[[nx]]$DispName) }
    })
        tmp2 <- Reduce(c, lapply(lin, function(x) {
            a <- c(x$ScientificName)
            names(a) <- x$Rank
            a
        }))
    names(tmp2)[which(names(tmp2) == "superkingdom")] <- "kingdom"
        list(tmp2, c(sn, on.flat))
  })
}


# Actually generate results

result.wrapper <- function(phyla, pheno, tree,
  proteins, clusters, randomize = FALSE, subsample.pheno = FALSE,
  method = matrix.plr, n.cores = 16) {
  lapply.across.names(phyla, function(p) {
    message(p)
    # avoid passing 0.5Gb object to matrix.plr
    if (subsample.pheno) {
      pheno <- pheno.resample(pheno[intersect(names(pheno), clusters[[p]])])
      message(paste0("(resampled to length ", length(pheno), ")"))
    }
    valid <- Reduce(intersect, list(colnames(proteins), clusters[[p]], names(pheno)), tree$tip.label)
    sub.fig <- proteins[, valid]
    if (randomize) {
      sub.fig <- apply(sub.fig, 1, sample) %>% t
      colnames(sub.fig) <- valid
    }
    method(tree, sub.fig, pheno, clusters[[p]] %>% as.character, cores = n.cores)
  })
}

result.wrapper.plm <- function(phyla, pheno, tree,
  proteins = phylum.ff.bin.restricted, clusters = phy.clusters, randomize = FALSE, subsample.pheno = FALSE,
  method = matrix.plm, n.cores = 16, restrict.figfams = NULL, drop.zero.var = FALSE, only.return.names = FALSE, ...) {
  lapply.across.names(phyla, function(p) {
    message(p)
    # avoid passing 0.5Gb object to matrix.plr
    if (class(tree) == "phylo") { tr <- tree } else if (class(tree) == "list") { tr <- tree[[p]] } else { stop("tree must be either an object of class phylo or a list") }
    valid <- Reduce(intersect, list(colnames(proteins[[p]]), clusters[[p]], names(pheno)), tr$tip.label)
    sub.fig <- proteins[[p]][, valid]
    if (!is.null(restrict.figfams)) {
      sub.fig <- sub.fig[intersect(rownames(sub.fig), restrict.figfams), ]
    }
    if (drop.zero.var) {
      fvar <- apply(sub.fig, 1, var)
      message(paste0(sprintf("%.01f", 100 * mean(na.omit(fvar == 0))), "% dropped [no variance]"))
      sub.fig <- sub.fig[which(fvar > 0), ]
    }
    #print(dim(sub.fig))
    if (randomize) {
      sub.fig <- apply(sub.fig, 1, sample) %>% t
      colnames(sub.fig) <- valid
    }
    if (!only.return.names) {
      method(tr, sub.fig, pheno, clusters[[p]] %>% as.character, cores = n.cores, ...)
    } else {
      rownames(sub.fig)
    }
  })
}

matrix.plm <- function(tree, mtx, pheno, restrict, cores = 8) {
  message("phylogenetic linear model")
  restrict <- as.character(restrict)
  valid.cols <- Reduce(intersect, list(colnames(mtx), restrict, names(pheno)))
  restricted.tree <- keep.tips(tree, valid.cols)
  cl <- makeCluster(cores)
  clusterCall(cl, source, file="plm-functions.R")
  clusterExport(cl, c("restricted.tree", "mtx", "pheno", "valid.cols"), 
    envir = environment())
  r <- parApply(cl, mtx, 1, function(m) {
    m1 <- m[valid.cols]
    p1 <- pheno[valid.cols]
    names(m1) <- valid.cols
    names(p1) <- valid.cols
    x <- tryCatch(
      phylolm(p1 ~ m1, phy = restricted.tree),
      error = function(e) { warning(paste(e)); c(Estimate = NA, p.value = NA) })
    if (!is.na(x[1])) {
      summary(x)$coefficients["m1", c("Estimate", "p.value")]
    } else { x }
  })
  stopCluster(cl)
  r
}
parLapplyLB <- function(cl, x, fun, ...) {
  message("... ...initializing... ...")
  clusterCall(cl, LB.init, fun, ...)
  message("... ...executing... ...")
  r <- clusterApplyLB(cl, x, LB.worker)
  clusterEvalQ(cl, rm('.LB.fun', '.LB.args', pos=globalenv()))
  r
}
LB.init <- function(fun, ...) {
  assign('.LB.fun', fun, pos = globalenv())
  assign('.LB.args', list(...), pos = globalenv())
  NULL
}
LB.worker <- function(x) {
  do.call('.LB.fun', c(list(x), .LB.args))
}



matrix.wpgls <- function(tree, mtx, pheno, restrict, cores = 8, w, verbose = FALSE) {
  force(w)
  force(pheno)
  force(tree)
  force(mtx)
  force(verbose)
  if (verbose) message("variance-weighted phylogenetic linear model")
  restrict <- as.character(restrict)
  valid.cols <- Reduce(intersect, list(colnames(mtx), restrict, names(pheno)))
  restricted.tree <- keep.tips(tree, valid.cols)
  if (verbose) message("...preparing...")
  tree.cor <- corBrownian(phy = restricted.tree)
  if (cores > 1) {
    if (verbose) {
      cl <- makeCluster(cores, outfile = "")
    } else {
      cl <- makeCluster(cores)
    }
    clusterCall(cl, source, file="plm-functions.R")
    clusterExport(cl, c("restricted.tree", "pheno", "valid.cols", "tree.cor", "w"), 
      envir = environment())
  }
  message(sprintf("%d total rows", nrow(mtx)))
  mL <- lapply(seq_len(nrow(mtx)), function(i) list(n = i, v = mtx[i, ]))
  if (verbose) message("...running...")
  mL.fxn <- function(mL.i) {
    m <- mL.i$v
    if (verbose) { if ((mL.i$n %% 250) == 0) message(mL.i$n) }
    m1 <- m[valid.cols]
    p1 <- pheno[valid.cols]
    w1 <- w[valid.cols]
    names(m1) <- valid.cols
    names(p1) <- valid.cols
    x <- tryCatch(
      gls(p1 ~ m1,
        correlation = tree.cor,
        weights = ~w1,
        data = data.frame(p1 = p1, m1 = m1, w1 = w1),
        method = "ML"),
      error = function(e) { warning(paste(e)); c(Estimate = NA, p.value = NA) })
    if (!is.na(x[1])) {
      coefficients(summary(x))["m1", c("Value", "p-value")]
    } else { x }
  }
  if (cores > 1) {
    r <- parLapplyLB(cl, mL, mL.fxn)
    stopCluster(cl)
  } else {
    verbose <- FALSE
    r <- pblapply(mL, mL.fxn)
  }
  r2 <- simplify2array(r)
  colnames(r2) <- rownames(mtx)
  r2
}


matrix.censored.plm <- function(tree, mtx, pheno, restrict, cores = 8, subsample = 0.75, nsubsamples = 500) {
  message("censored phylogenetic linear model")
  restrict <- as.character(restrict)
  valid.cols <- Reduce(intersect, list(colnames(mtx), restrict, names(pheno)))
  restricted.tree <- keep.tips(tree, valid.cols)
  cl <- makeCluster(cores)
  clusterCall(cl, source, file="plm-functions.R")
  clusterCall(cl, source, file="test-induced-corr.R")
  clusterExport(cl, c("restricted.tree", "mtx", "pheno", "valid.cols"), 
    envir = environment())
  r <- parApply(cl, mtx, 1, function(m) {
    m1 <- m[valid.cols]
    p1 <- pheno[valid.cols]
    names(m1) <- valid.cols
    names(p1) <- valid.cols
    x <- tryCatch(
      cps.coefs(
        censored.phylolm.subsample(m1, p1, phy = restricted.tree, pct = subsample, B = nsubsamples))["sX", ],
      error = function(e) { warning(paste(e)); c(Estimate = NA, p.value = NA) })
    x
  })
  stopCluster(cl)
  r
}

matrix.interaction.plm <- function(tree, mtx, pheno, restrict, cores = 8) {
  message("phylogenetic linear model with interactions")
  restrict <- as.character(restrict)
  valid.cols <- Reduce(intersect, list(colnames(mtx), restrict, names(pheno)))
  restricted.tree <- keep.tips(tree, valid.cols)
  cl <- makeCluster(cores)
  clusterCall(cl, source, file="plm-functions.R")
  clusterExport(cl, c("restricted.tree", "mtx", "pheno", "valid.cols"), 
    envir = environment())
  r <- lapply(1:((nrow(mtx) - 1)), function(a) {
    q <- sapply((a + 1):nrow(mtx), function(b) {
      mA <- mtx[a, ]
      mB <- mtx[b, ]
      if (mB == mA) { return(c(Estimate = NA, p.value = NA)) } # if it's literally the same, can't do anything
      mA1 <- mA[valid.cols]
      mB1 <- mB[valid.cols]
      if (sum(mA1 * mB1) == 0) { return(c(Estimate = NA, p.value = NA)) } # also can't do anything here
      p1 <- pheno[valid.cols]
      names(mA1) <- valid.cols
      names(p1) <- valid.cols
      x <- tryCatch(
        phylolm(p1 ~ mA1 * mB1, phy = restricted.tree),
        error = function(e) { c(Estimate = NA, p.value = NA) })
      if (!is.na(x[1])) {
        summary(x)$coefficients[4, c("Estimate", "p.value")]
      } else { x }
      })
    simplify2array(q)
    q2 <- cbind(matrix(nr = 2, rep(NA, a * 2)), q)
    colnames(q2) <- rownames(mtx)
    q2
  })
  names(r) <- rownames(mtx[1:(nrow(mtx)-1)])
  simplify2array(r)
  stopCluster(cl)
  r
}

matrix.notree.lm <- function(tree, mtx, pheno, restrict, cores = 8) {
  message("sham phylogenetic linear model (regular linear model)")
  restrict <- as.character(restrict)
  valid.cols <- Reduce(intersect, list(colnames(mtx), restrict, names(pheno)))
  restricted.tree <- keep.tips(tree, valid.cols)
  print(restricted.tree$tip.label[1:5])
  cl <- makeCluster(cores)
  clusterCall(cl, source, file="plm-functions.R")
  clusterExport(cl, c("restricted.tree", "mtx", "pheno", "valid.cols"), 
    envir = environment())
  r <- parApply(cl, mtx, 1, function(m) {
    m1 <- m[valid.cols]
    p1 <- pheno[valid.cols]
    names(m1) <- valid.cols
    names(p1) <- valid.cols
    x <- tryCatch(
      lm(p1 ~ m1 + 1),
      error = function(e) { c(Estimate = NA, p.value = NA) })
    if (!is.na(x[1])) {
      z <- tryCatch({
        summary(x)$coefficients[2, c("Estimate", "Pr(>|t|)")]}, 
        error = function(e) {c(Estimate = NA, p.value = NA)})
      names(z) <- c("Estimate", "p.value")
      z
    } else { x }
  })
  stopCluster(cl)
  r
}


results.report <- function(results, sigs, signs) {
  pos.sigs <- lapply(c(strong="strong",med="med",weak="weak"), function(level) {
    make.pos.sig(sigs, signs, level)  
  })
  sapply(pos.sigs, function(x) sapply(x, function(y) length(y)))
  
}

get.top.N <- function(p, sigs, signs, results, level = "strong", exclude = NULL, N = 25, total.n.cutoff = 0, genomes.per.protein = NULL) {
  sig.up <- intersect(sigs[[p]][[level]], nw(signs[[p]] > 0))
  if (!is.null(genomes.per.protein)) {
    gn <- names(which(genomes.per.protein[[p]] >= total.n.cutoff))
    sig.up <- intersect(sig.up, gn)
  }
  sig.pv <- results[[p]][2, sig.up]
  setdiff(names(sort(sig.pv, dec = F)), exclude)[1:N]
}


# figure out which gene families are perfectly correlated so we can group them
# together

group.similar.genes <- function(gene.pres.abs, gene.lists) {
  mapply(gene.pres.abs, gene.lists, FUN=function(pa, gl) {
    bindist <- dist(pa[gl, ], 'binary')
    perfect.matches <- data.matrix((as.matrix(bindist) == 0) * 1)
    match.graph <- graph_from_adjacency_matrix(perfect.matches)
    cliques <- maximal.cliques(match.graph)
    if ((sapply(cliques, length) %>% sum) != length(gl)) {
      stop("error, perfect matches are somehow not symmetric")
    }
    cliques
  })
}

# given matrices for disease and healthy, get conditional probabilities

divide.by.boxplot <- function(x, y) { 
  i <- intersect(names(x), names(y))
  xi <- x[names(i)]
  yi <- y[names(i)]
  yn <- unique(yi)
  boxplot(lapply(yn, function(n) { xi[which(yi == n)] } ))
}

###

bootstrap.regularize.pET <- function(vec, env.ids, which.env = 1, prior = 0.05, b = 1, debug = FALSE, add.pc = FALSE,
  min.limit = -10, max.limit = 10, set.layout = TRUE, B = 200) {
  p.distro <- sapply(1:B, function(.) {
    resample.pos <- sample(vec[which(env.ids == which.env)], replace = TRUE)
    resample.neg <- sample(vec[which(env.ids != which.env)], replace = FALSE)
    resample.vec <- vec
    resample.vec[which(env.ids == which.env)] <- resample.pos
    resample.vec[which(env.ids != which.env)] <- resample.neg
    regularize.pET(resample.vec, env.ids, which.env, prior, b, debug, add.pc, min.limit, max.limit, set.layout)[2]
  })
}
# laplace regularization of P(T|E) and P(T|~E) via MAP estimation
# param is ~ p**k * (1-p)**(n-k) (the n choose k is constant wrt p)
# prior is ( (1/(2b)) * exp(-(|x-m|/b) ) where x, m are on logit scale

regularize.pET <- function(vec, env.ids, which.env = 1, prior = 0.05, b = 1, debug = FALSE, add.pc = FALSE, min.limit = -10, max.limit = 10, set.layout = TRUE) {
  # n and k in binomial are #(E) and #(E,T), respectively (used on P(T|E))
  n <- sum(env.ids == which.env) # #(E)
  k <- sum(vec[env.ids == which.env] > 0) # #(E & T)
  if (add.pc) { n <- n+1; k <- k+1 }
  if (!add.pc) {
    pT.E <- mean(c(vec[which(env.ids == which.env)] > 0))
    pT.nE <- mean(c(vec[which(env.ids != which.env)]) > 0)
  } else {
    pT.E <- mean(c(0, 1, vec[which(env.ids == which.env)] > 0))
    pT.nE <- mean(c(0, 1, vec[which(env.ids != which.env)]) > 0)
  }
  if (is.nan(pT.E)) { print("!"); print(vec); print(vec[which(env.ids == which.env)]) }
  pT <- pT.E * prior + pT.nE * (1 - prior)
  map <- function(p) {
    x <- (p * pT) / prior # P(T|C)
    dbinom(k, n, x) * ((1/2*b)) * exp(-(abs(logit(p)-logit(prior))/b))
  }
  map.logit <- function(logit.p) {
    p <- logistic(logit.p)
    x <- (p * pT) / prior # P(T|C)
    # return log probability also, better numerical stability
    dbinom(k, n, x, log = TRUE) + log(1 / (2*b)) - (abs(logit(p)-logit(prior))/b)
  }
  if (debug) {
    map.binomial.x <- function(x) {
      # show the actual value for x
      dbinom(k, n, x)
    }
    map.binomial.p <- function(logit.p) {
      # show the actual value for p
      p <- logistic(logit.p)
      x <- (p * pT) / prior # P(T|C)
      dbinom(k, n, x)
    }
    map.laplace <- function(logit.p) {
      p <- logistic(logit.p)
      ((1/2*b)) * exp(-(abs(logit(p)-logit(prior))/b))
    }
    y <- seq(0.001, 0.999, 0.001)
#    z <- seq(-10, 10, 0.01)
    max.p <- prior / pT
    z <- seq(-10, min(10, logit(max.p)), 0.01)
    bin <- sapply(y, map.binomial.x)
    bin2 <- sapply(z, map.binomial.p)
    lp <- sapply(z, map.laplace)
    overall <- sapply(z, map.logit)
    if (set.layout) layout(matrix(nr = 4, 1:4))
    plot(y, bin, type = 'l', xlab = "P(T|E)", ylab = paste0("Binomial (x): ", format(y[which.max(bin)], digits = 2)))
    plot(z, bin2, type = 'l', xlab = "logit(P(E|T))", ylab = paste0("Binomial (p): ", format(logistic(z[which.max(bin2)]), digits = 2)), xlim = c(min(z), max(z)))
    plot(z, lp, type = 'l', xlab = "logit(P(E|T))", ylab = paste0("Laplace: ", format(logistic(z[which.max(lp)]), digits = 2)), xlim =c(min(z), max(z)))
    plot(z, overall, type = 'l', xlab = "logit(P(E|T))", ylab = paste0("Posterior: ", format(logistic(z[which.max(overall)]), digits = 2)), xlim = c(min(z), max(z)))
    if (set.layout) layout(c(1))
  }
  max.p <- prior / pT # q is 1 here
  max.p <- min(max.p, logistic(max.limit))
  results <- optimize(map.logit, c(min.limit, logit(max.p)), maximum = TRUE)
  initial.x <- pT.E
  initial.p <- (pT.E * prior) / pT
  final.x <- logistic(results$maximum) * pT / prior
  final.p <- logistic(results$maximum)
  return(c(x = final.x, p = final.p, x.init = initial.x, p.init = initial.p, pT = pT))
}

bootstrap.pET.se <- function(vec, env.ids, which.env = 1, prior = 0.05, B = 250, add.pc = TRUE) {
  if (var(vec) == 0) {
    # save time when there's no variance in input vector
    pT.E <- mean(c(0, 1, vec[env.ids == which.env] > 0))
    pT.nE <- mean(c(0, 1, vec[env.ids != which.env] > 0))
    pT <- pT.E * prior + pT.nE * (1 - prior) # marginalize
    pET <- (pT.E * prior) / pT # Bayes' rule
    distro <- rep(pET, B)
    return(list(distro = distro, pET = pET, se = 0))
  }
  pET.distro <- sapply(1:B, function(.) {
    in.values <- sample(vec[env.ids == which.env], replace = TRUE)
    out.values <- sample(vec[env.ids != which.env], replace = TRUE)
    if (add.pc) {
      in.values <- c(in.values, 0, 1)
      out.values <- c(out.values, 0, 1)
    }
    pT.E <- mean(in.values > 0)
    pT.nE <- mean(out.values > 0)
    pT <- pT.E * prior + pT.nE * (1 - prior) # marginalize
    (pT.E * prior) / pT # Bayes' rule
  })
  return(list(distro = pET.distro, pET = mean(pET.distro), se = sd(pET.distro)))
}

iter.across <- function(list.of.lists, f, table = FALSE, reduce = TRUE) {
	rows <- Reduce(prod, sapply(list.of.lists, length))
	to.expand <- lapply(list.of.lists, function(x) {
		if (is.null(names(x))) seq_along(x) else names(x)
	})	
	input.grid <- do.call(expand.grid, to.expand)
	output <- apply(input.grid, 1, function(x) {
    i <- mapply(x, list.of.lists, SIMPLIFY = FALSE, FUN=function(y,z) {z[[y]]})
    names(i) <- names(formals(f))
    res <- do.call(f, i)
    #print(res)
    if (!table) { return(res) }
    if (!reduce) { return(res) }
    if (!is.null(nrow(res))) {
      cbind(sapply(x, function(y) rep(y, nrow(res))), res)
    } else {
    if (length(res) == 0) { c(x, "--")} else { c(x, res) }
    }
  })
  if (!table) {	
    output.grid <- data.frame(cbind(input.grid, t(output)))
  } else {
    if (reduce) {
      output.grid <- Reduce(rbind, Filter(f = function(x) { r <- nrow(x); if (is.null(r)) TRUE else if (r == 0) FALSE else TRUE}, output)) 
      rownames(output.grid) <- NULL
    }
    else {
      output.grid <- list(input = input.grid, output = output)
      names(output.grid$output) <- apply(input.grid, 1, function(x) paste(x, sep = '.', collapse = '.'))
    }
  }
  output.grid
}

reduce.ia <- function(ia.results, f) {
  input.grid <- ia.results$input
  output <- ia.results$output
  output2 <- lapply(1:length(output), function(i) {
    o2 <- f(output[[i]])
    if (!is.null(nrow(o2))) {
      cbind(sapply(input.grid[i, ], function(y) rep(y, nrow(o2))), o2)
    } else {
      cbind(input.grid[i, ], rbind(o2))
    }
  })
  all.ncols <- (sapply(output2, ncol))
  max.ncol <- max(all.ncols)
  output.filtered <- Filter(f = function(x) {
      r <- nrow(x)
      if (is.null(r)) FALSE else if (r == 0) FALSE else if (ncol(x) < max.ncol) FALSE else TRUE
  }, output2)
  output.grid <- Reduce(rbind, lapply(output.filtered, function(x) {colnames(x) <- 1:ncol(x); rownames(x) <- NULL; x} ))
  rownames(output.grid) <- NULL
  output.grid
}

fit.beta.list <-  function(mtx, ids) {
  lapply(unique(ids), function(i) {
    fitdistr(densfun = "beta", start = list(shape1 = 1, shape2 = 1), 
      apply(mtx[, which(ids == i)], 1, function(x) mean(c(x, 0, 1) > 0)))$estimate
  })
}

# e.g.: betas <- fit.beta.list(mgs.mtx.r, mgs.ids$ids)

# effect.size is in logit-space (sign randomly chosen with probability
# sign.pos.prob, so asymmetric effects could be modeled if desired)
# baseline.distro gives P(T) priors across matrix

bootstrap.logitP.se <- function() {

}

simulate.binom.mtx <- function(effect.size = 2, baseline.distro = c(shape1 = 0.66, shape2 = 2.62), samples = c(H = 38, D = 13), which.env = "D", taxa = 2000, tpr = 0.25, sign.pos.prob = 0.5) {
  tp.taxa <- round(tpr * taxa)
  neg.taxa <- taxa - (tp.taxa)
  n.classes <- length(samples)
  pT <- sapply(1:taxa, function(.) rbeta(1, baseline.distro[1], baseline.distro[2]))
  fx <- c(((2 * rbinom(n = tp.taxa, size = 1, sign.pos.prob)) - 1),
    rep(0, neg.taxa))
  pTbs <- sapply(1:taxa, function(i) { pTb <- logistic(logit(pT) + (fx * (effect.size))) })
  if (is.null(names(samples))) names(samples) <- 1:length(samples)
  sim.mtx <- t(sapply(1:taxa, function(i) {
    #c(rbinom(samples[1], size = 1, pT[i]), rbinom(samples[2], size = 1, pTbs[i]))
    Reduce(c, lapply.across.names(names(samples), function(smp) { 
        if (smp == which.env) { p <- pTbs[i] } else { p <- pT[i] }
        rbinom(samples[smp], size = 1, p) 
    }))
  }))
  ids <- Reduce(c, lapply(1:length(samples), function(i) rep(i, samples[i])))
  return(list(mtx = sim.mtx, pT = pT, pTbs = pTbs, fx = fx, ids = ids, input.params = list(effect.size = effect.size, baseline = baseline.distro, samples = samples, taxa = taxa, tpr = tpr, sign.pos.prob = sign.pos.prob)))
}

score.regularization <- function(mtx, ids, real.fx, which.env = 2, prior = 0.002, b = 0.1, add.pc = FALSE, tol = prior * 0.005, ...) {
  regularized <- apply(mtx, 1, function(x) regularize.pET(x, ids, which.env = which.env, prior = prior, b = b, debug = FALSE, add.pc = add.pc))
  posteriors <- regularized[2, ]
  signif <- 1 * (abs(posteriors - prior) > tol)
  predicted.signs <- signif * sign(posteriors - prior)
  fpr <- (signif[real.fx == 0] %>% mean)
  pwr.hi <- mean(predicted.signs[real.fx == 1] == 1)
  pwr.lo <- mean(predicted.signs[real.fx == -1] == -1)
  return(c(fpr = fpr, pwr.hi = pwr.hi, pwr.lo = pwr.lo))
}

optimize.b.wrapper.single <- function(real.mtx, real.ids, which.env = 2, which.shape = 1, prior = 0.002, effect.size = 2, ...) {
  shape.n <- which((real.ids %>% unique %>% sort) == which.shape)
  shapes <- fit.beta.list(real.mtx, real.ids)
  N <- count.each(real.ids)
  sim <- simulate.binom.mtx(effect.size = effect.size, baseline.distro = shapes[[shape.n]], samples = N, taxa = nrow(real.mtx), tpr = 0.25, sign.pos.prob = 0.5)
  optimize.b(sim, which(names(N) == which.env), prior, ...)
}

optimize.b.wrapper <- function(real.mtx, real.ids, which.real.env = 2, which.shape = 1, prior = 0.002, effect.size = 2, bounds = c(0.01, 5), add.pc = FALSE, tol = prior * 0.005, a = 0.05, verbose = FALSE, optim.fxn = s.opt.fxn.3, pos.prop = 0.5, ...) {
  shape.n <- which((real.ids %>% unique %>% sort) == which.shape)
  shapes <- suppressWarnings(fit.beta.list(real.mtx, real.ids))
  N <- count.each(real.ids)
  which.env <- which(names(N) == which.real.env)
  get.optim <- function(b) {
    sim <- simulate.binom.mtx(effect.size = effect.size, baseline.distro = shapes[[shape.n]], samples = N, taxa = nrow(real.mtx), tpr = 0.25, sign.pos.prob = pos.prop, which.env = which.real.env)
    s <- score.regularization(sim$mtx, sim$ids, sim$fx, which.env, prior, b, add.pc, tol, ...)
    if (verbose) print(c(s, optim.fxn(s, a)))
    optim.fxn(s, a)
  }
  optimize(get.optim, bounds, maximum = TRUE)
}

geommean <- function(x) exp(sum(log(x)) / length(x))
s.opt.fxn.3 <- function(s, a) {
  if (s["fpr"] <= a) (1 + geommean(s[c("pwr.hi", "pwr.lo")])) else ( (1 - s["fpr"]))
}
s.opt.fxn.2 <- function(s, a = 0.05) (mean(s[c("pwr.hi", "pwr.lo")])) * (a / max(s["fpr"], a))
s.opt.fxn <- function(s, a = 0.05) (geommean(s[c("pwr.hi", "pwr.lo")])) * (1 * (s["fpr"] <= a))

optimize.b <- function(sim, which.env, prior, add.pc = FALSE, tol = prior * 0.005, verbose = FALSE, ...) {
  get.optim <- function(b) {
    s <- score.regularization(sim$mtx, sim$ids, sim$fx, which.env, prior, b, add.pc, tol, ...)
    if (verbose) print(s)
    s.opt.fxn.2(s)
  }
  optimize(get.optim, c(0.01, 5), maximum = TRUE)  
}


# prior is P(E); vec is vector of taxon presence/absences; env.ids is a vector
# giving the class each sample belongs to (probably a factor but not
# necessarily); which.env is the class ID of E

make.class.ids <- function(mtx, sample.vecs, elim.zeros = TRUE) {
  if (any(is.null(names(sample.vecs)))) {
    names(sample.vecs) <- 1:length(sample.vecs)
  }
  u <- Reduce(union, sample.vecs)
  concat <- Reduce(c, sample.vecs)
  if (any(sort(u) != sort(concat))) { stop("class assignments must be unique") }
  mtx.reduced <- mtx[, u]
  prevs <- apply(mtx.reduced, 1, function(x) sum(x > 0))
  if (elim.zeros) {
    mtx.reduced <- mtx.reduced[which(prevs > 0), u]
  }
  ids <- vector(mode = "numeric", length = ncol(mtx.reduced))
  for (i in 1:length(sample.vecs)) {
    ids[which(colnames(mtx.reduced) %in% sample.vecs[[i]])] <- names(sample.vecs)[i]
  }
  names(ids) <- colnames(mtx.reduced)
  return(list(mtx = mtx.reduced, ids = ids))
}

# from flodel @ stackoverflow
list.depth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, list.depth)), 0L)
first.element.at.depth <- function(l, n) {
	if (n == 1) { l[[1]] } else { (first.element.at.depth(l[[1]], n - 1)) }
}

annotate.nested <- function(nested, summarize = NULL, n = NULL, n.names = NULL) {
	nestedness <- list.depth(nested)
	if (nestedness == 0) {
		if (is.null(summarize)) {
			values <- data.frame(value = nested)
			rownames(values) <- names(nested)	
		} else {
			values <- data.frame(value = summarize(nested))
		}
		n.mtx <- matrix(rep(n, length(values)), nc = length(n), byrow = TRUE)
		if (!is.null(n.names)) { colnames(n.mtx) <- n.names }
		final.df <- cbind(names = rownames(values), n.mtx, value = values)
		rownames(final.df) <- NULL
		return(final.df)
	} else {
    if (is.null(names(nested))) { names(nested) <- 1:length(nested) }
		Reduce(rbind, lapply.across.names(names(nested),
			function(x) annotate.nested(nested[[x]], summarize = summarize, n = c(n, x), n.names = n.names)))
	}
}

zipData <- function(iterover) {
  n <- names(iterover)
  names(n) <- n
  return(lapply(n, function(x) {
    list(name = n[x], data = iterover[[n[x]]])
  }))
}

zipLapply <- function(iterover, fxn) {
  lapply(zipData(iterover), fxn)
}

zipSapply <- function(iterover, fxn) {
  simplify2array(zipLapply(iterover, fxn))
}

get.paired.data <- function(phylum = "Bacteroidetes", gene = "FIG01376814", cTrait = species.prev.avg.logit, bTrait = phylum.ff.bin.restricted) {
  both <- intersect(colnames(bTrait[[phylum]]), names(cTrait))
  bTrait.1 <- bTrait[[phylum]][gene, both]
  cTrait.1 <- cTrait[both]
  return(list(bTrait = bTrait.1, cTrait = cTrait.1))
}

tree.to.dist <- function(tree) {
  ((1 - (vcv(tree, TRUE))) / 2) %>% as.dist
}


### calculate prevalence

prev.addw <- function(mtx, envir, meta) {
  if (!(envir %in% levels(meta$env))) {
    stop(paste0("environment ", envir, " not found in metadata"))
  }
  env.rows <- (meta$env == envir)
  dsets <- unique(meta[env.rows, "dataset"])
  if (length(dsets) > 1) {
    means.by.study <- lapply(dsets, function(d) {
      s <- intersect(colnames(mtx), 
                     meta$sample[(env.rows & (meta$dataset == d))] %>%
                       as.character)
      rm <- rowMeans(1 * (mtx[, s] > 0))
      list(rm = rm, s = s)
    })
    means.only <- sapply(means.by.study, function(x) x$rm)
    total.s <- Reduce(sum, lapply(means.by.study, function(x) x$s))
    avg.prev <- rowMeans(means.only)
    addw <- (1 + (avg.prev * total.s)) / (2 + total.s)
  } else {
    s <- intersect(colnames(mtx),
                   meta$sample[env.rows] %>% as.character)
    total.s <- length(s)
    rs <- rowSums(1 * (mtx[, s] > 0))
    addw <- (1 + rs) / (2 + total.s)
  }
  return(logit(addw))
}

calc.ess <- function(mtx, envir, meta, ptype, pdata = NULL, b.optim = NULL) {
  if (!(envir %in% levels(meta$env))) {
    stop(paste0("environment ", envir, " not found in metadata"))
  }
  env.rows <- (meta$env == envir)
  dsets <- unique(meta[env.rows, "dataset"])
  if (length(dsets) > 1) {
    stop("can't currently calculate specificity across more than one dataset")
  }
  meta.present <- meta[(meta$sample %>% as.character %in% colnames(mtx)), ]
  envirs <- unique(meta.present$env)
  if (ptype == "uninformative") {
    priors <- data.frame(
      env = envirs,
      prior = 1/length(envirs)
    )
  } else if (ptype == "alpha") {
    stop("not implemented yet, sorry")  
  } else if (ptype == "file") {
    priors <- pdata[pdata$env %in% envirs, ]
  } else {
    stop(paste0("don't know how to compute priors of type ", ptype))
  }
  ids <- sapply(colnames(mtx), function (sn) meta[(meta$sample == sn), 
                                                  "env"])
  names(ids) <- colnames(mtx)
  if (is.null(b.optim)) {
    b.optim <- optimize.b.wrapper(effect.size = 2,
                                  real.mtx = mtx,
                                  real.ids = ids,
                                  which.real.env = envir,
                                  which.shape = envir,
                                  prior = priors$prior[
                                    which(priors$env == envir)],
                                  a = 0.05,
                                  optim.fxn = s.opt.fxn.3,
                                  verbose = FALSE)$maximum
  } else {
    warning("skipping optimization of Laplace parameter (not recommended)")
  }
  regularized <- apply(mtx, 1, function(x) {
    regularize.pET(x, ids, which.env = envir, b = b.optim)
  })
  return(list(
    b.optim = b.optim,
    ess = logit(regularized[2, ]),
    regularized = regularized, 
    priors = priors))
}

gg.cont.tree <- function(phy,
                         ctrait,
                         cAnc = NULL,
                         model = "ARD",
                         cLimits = logit(c(0.025, 0.1)),
                         n = NULL,
                         reduced.phy = NULL,
                         colors = c(
                           low.col = "slateblue",
                           mid.col = "black",
                           high.col = "orange2"),
                         plot = T,
                         restrict = NULL,
                         cName = "prevalence",
                         reverse = F,
                         ladderize = T,
                         ...) {
  if (!is.null(restrict)) {
    ctrait <- ctrait[intersect(names(ctrait), restrict)]
  }
  if (is.null(n)) { n <- intersect(phy$tip.label, names(ctrait)) }
  if (is.null(reduced.phy)) { reduced.phy <- fix.tree(keep.tips(phy, n)) }
  message("getting continuous trait ancestry")
  if (is.null(cAnc)) cAnc <- fastAnc(reduced.phy, ctrait[n])
  # concatenate tip values and node values
  cDisplay <- truncated(c(ctrait[n], cAnc), cLimits)
  # make color scales
  if ("mid.col" %in% names(colors)) {
    message("found mid.col")
    cColors <- scale_color_gradient2(low = colors["low.col"],
                                     high = colors["high.col"],
                                     mid = colors["mid.col"],
                                     midpoint = mean(cLimits),
                                     guide = "colorbar",
                                     name = cName)
  } else {
    cColors <- scale_color_gradient(low = colors["low.col"],
                                    high = colors["high.col"],
                                    guide = "colorbar",
                                    name = cName)
  }
  # plot trees
  if (reverse) {
    ctree <- ggtree(reduced.phy,
                    ladderize = ladderize,
                    aes_string(color = cDisplay), ...) +
      cColors + scale_x_reverse() + theme(legend.position = "bottom")
  } else {
    ctree <- ggtree(reduced.phy,
                    ladderize = ladderize,
                    aes_string(color = cDisplay), ...) +
      cColors + theme(legend.position = "bottom")
  }
  if (plot) print(ctree)
  return(list(tree = ctree,
              cAnc = cAnc,
              rphy = reduced.phy,
              n = n,
              cols = colors,
              lims = cLimits,
              disp = cDisplay))
}

kable.recolor <- function(x, direction = 1, option = "D",
  na_color = "#BBBBBB", scale_from = NULL, colors = c("#000000","#FFFFFF"),
  limits = c(-Inf, Inf))
{
  if (x < limits[1]) { x[x < limits[1]] <- limits[1] }
  if (x > limits[2]) { x[x < limits[2]] <- limits[2] }
  if (is.null(scale_from)) {
    x <- round(rescale(x, c(1, 256)))
  }
  else {
    x <- round(rescale(x, to = c(1, 256), from = scale_from))
  }
  if (direction == -1) { x <- (257 - x) }
  cmap <- colorRampPalette(colors)(256)
  color_code <- cmap[x]
  color_code[is.na(color_code)] <- na_color
  return(color_code)
}
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
    {s <- substring(s, 2); if(strict) tolower(s) else s},
    sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
###

# Terrible XML Hack!
hack.tree.labels <- function(tree.obj,
  file,
  stroke.scale = 0.7,
  pheno = NULL,
  pheno.name = NULL,
  native.tooltip = FALSE,
  units = "",
  ...) {

  tip.labels <- with(tree.obj$data, label[isTip])
  xml <- xmlSVG(print(tree.obj), standalone = TRUE, ...)
  new.style.text <- " \n .faketip:hover ~ .realtip { \n stroke-width: 5; \n opacity: 1; \n  } \n .faketip:hover ~ .specieslabel { \n opacity: 1; \n } \n " 
  style.nodes <- xml_find_all(xml, "//*[local-name()='style']")
  xml_set_text(style.nodes[1], new.style.text)
  xml.text <- xml_find_all(xml, "//*[local-name()='text']")
  xml.text.contents <- sapply(xml.text, xml_text)
  xml.label.indices <- which(xml.text.contents %in% tip.labels)
  xml.label.heights <- sapply(xml.text, function(x) {
    xml_attrs(x)["y"]
  })
  xml.label.pair <- cbind(label = xml.text.contents,
    y = xml.label.heights)[xml.label.indices, ]
  ordered.labels <- xml.label.pair[order(xml.label.pair[, "y"] %>% as.numeric), "label"]
  xml.lines <- xml_find_all(xml, "//*[local-name()='line']")
  xml.line.props <- sapply(xml.lines, function(x) xml_attrs(x))
  xml.x2 <- xml.line.props["x2", ] %>% as.numeric
  xml.y2 <- xml.line.props["y2", ] %>% as.numeric
  # terminus <- max(xml.x2)
  uniq.x2 <- unique(xml.x2)
  count.x2 <- sapply(uniq.x2, function(x) sum(xml.x2 == x))
  terminus <- uniq.x2[which(count.x2 == length(tip.labels))]
  xml.y2.sorted <- sort(xml.y2[xml.x2 == terminus])
  # skootch over, remove tip, add title
  for (x in xml.text) {
    label <- xml_text(x)
    if (label %in% ordered.labels) {
      xml_set_attr(x, "x", as.character(terminus))
      xml_set_text(x, " ")
      if (native.tooltip) {
        xml_add_sibling(x, read_xml(paste0("<title>",
              label,
              "</title>")))
      }
    }
  }
  for (l in xml.lines) {
    l.attr <- xml_attrs(l)
    l.y2 <- as.numeric(l.attr["y2"])
    l.x2 <- as.numeric(l.attr["x2"])
    l.x1 <- as.numeric(l.attr["x1"])
    # This step is necessary because otherwise mouseover won't work
    style <- xml_attrs(l)["style"]
    s.parsed <- style.parse(style)
    for (n in 1:length(s.parsed)) {
      xml_set_attr(l, names(s.parsed)[n], s.parsed[n])
    }
    xml_set_attr(l, "style", "")
    if ("stroke-width" %in% names(s.parsed)) {
      xml_set_attr(l, "stroke-width", (as.numeric(s.parsed["stroke-width"]) * stroke.scale) %>% as.character)
    }
    if (!("stroke" %in% names(s.parsed))) {
      xml_set_attr(l, "stroke", "#000000")
    }
    if (l.x2 == terminus) {
      label <- ordered.labels[which(xml.y2.sorted == l.y2)]
      xml_set_attr(l, "id", label)
      xml_set_attr(l, "class", "realtip")
      new.group <- read_xml("<g class=\"tip\"> </g>")
      l2 <- xml_add_child(new.group, l)
      xml_add_child(new.group, l)
      xml_set_attr(l2, "opacity", "0")
      xml_set_attr(l2, "pointer-events", "all")
      xml_set_attr(l2, "stroke-width", 5)
      xml_set_attr(l2, "class", "faketip")
      xml_set_attr(l2, "x1",
        as.character(
          l.x1 - 100
          ))
      xml_set_attr(l2, "x2",
        as.character(
          terminus + 300
          ))
      extra.info <-  ""
      if (!is.null(pheno)) {
        if (is.null(pheno.name)) pheno.name <- "phenotype"
        if (label %in% names(pheno)) {
          phi <- format(pheno[label], digits = 3)
        } else {
          phi <- "NA" 
        }
        extra.info <- paste0("(", pheno.name, " = ", phi, units,  ")")
      }
      xml_add_child(new.group, read_xml(paste0(
            "<text x=\"",
            l.x2 + 5,
            "\" y = \"",
            l.y2 + 3, 
            "\" opacity=\"0\" pointer-events=\"all\" style=\"font-family: Arial; font-size: 10px;",
            " fill: ",
            xml_attr(l, "stroke"),
            ";",
            "\" class=\"specieslabel\"> ",
            label, 
            " ",
            extra.info,
            " ",
            "</text>")))
      if (native.tooltip) {
        xml_add_child(new.group, read_xml(paste0("<title>",
              label,
              "</title>")))
      }
      xml_replace(l, new.group)
    }
  }
  write_xml(x = xml, file) 
}

### minimal.example <- function(tree.obj, file, ...) {
minimal.example <- function(.) {
  library(xml2)
  library(ape)
  library(ggtree)
  library(svglite)
  random.tree <- rcoal(50)
  tree.xml <- xmlSVG(print(ggtree(random.tree) + geom_tiplab()), standalone = TRUE)
  tree.lines <- xml_find_all(tree.xml, "//*[local-name()='line']")
  for (l in tree.lines) {
    new.group <- read_xml("<g> </g>")
    xml_add_child(new.group, l)
#    xml_add_child(new.group, read_xml(paste0("<title>",
#          label,
#          "</title>")))
    #xml_add_child(tree.xml, new.group)
    #xml_remove(l)
    xml_replace(l, new.group)
  }
  return(tree.xml)
}

style.parse <- function(str) {
  semi.split <- strsplit(str, ";") %>% sapply(., trimws)
  if (is.null(dim(semi.split))) {
    c.split <- trimws(strsplit(semi.split, ":")[[1]])
    c.output <- c.split[2]
    names(c.output) <- c.split[1]
  } else {
    c.split <- apply(semi.split, 1, function(x) strsplit(x, ":")[[1]]) %>% trimws
    c.output <- c.split[2, ]
    names(c.output) <- c.split[1, ]
  }
  return(c.output)
}
