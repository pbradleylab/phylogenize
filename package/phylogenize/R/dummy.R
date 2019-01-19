# Functions for creating and writing dummy data

simulate.counts <- function(com, nr) {
    counts <- mapply(nr, data.frame(com), FUN = function(r, x) {
        if (all(x == 0)) { return(x) }
        m <- rmultinom(1, prob=x, size=r)
        names(m) <- rownames(com)
        m
    })
    colnames(counts) <- colnames(com)
    counts
}

simple.prevalence <- function(mtx) {
    apply(cbind(1 * (mtx > 0), 1, 0), 1, mean)
}

colnorm <- function(mtx) {
    apply(mtx, 2, function(x) {
        if (all(x == 0)) return(x)
        x/sum(x)
    })
}

make.sim <- function(n.taxa,
                     n.samples,
                     avg.reads, # can be a vector
                     prev.dist=c(-4,2),
                     effect.size.mtx=0,
                     alphas=NULL) {
    if (is.null(alphas)) {
        alphas <- round(runif(n.taxa, min = 1, max = 100000))
    }
    reads <- rnorm(n.samples, log(avg.reads), 0.5) %>% exp %>% round
    com <- gtools::rdirichlet(n.samples, alphas) %>% t
    colnames(com) <- paste0("sample", 1:n.samples)
    rownames(com) <- paste0("taxon", 1:n.taxa)
    ncells <- nrow(com) * ncol(com)
    random.logitprevs <- rnorm(n.taxa, prev.dist[1], prev.dist[2])
    logitprev.mtx <- rep.col(random.logitprevs, n.samples) + effect.size.mtx
    random.prevalences <- logistic(logitprev.mtx)
    rownames(random.prevalences) <- rownames(com)
    for (nt in 1:n.taxa) {
        prev.vec <- rbinom(n.samples, prob=random.prevalences[nt, ], size = 1)
        com[nt, (prev.vec == 0)] <- 0
    }
    # need to renormalize after zeroing out
    com <- colnorm(com)
    true.prevs <- rowMeans(com > 0)
    sim <- simulate.counts(com, reads) %>% colnorm
    return(list(sim=sim,
                com=com,
                prev=random.prevalences))
}


add.fx.to.matrix <- function(samples, taxa, which.taxa, which.samples, mtx=0, fx=1) {
    fx.mtx <- matrix(nc=length(samples),
                     nr=length(taxa),
                     0,
                     dimnames=list(taxa, samples))
    fx.mtx[which.taxa, which.samples] <- fx
    mtx + fx.mtx
}


build.fx.matrix <- function(samples,
                            taxa,
                            fx,
                            wtaxa,
                            wsamples) {
    fx.mtx <- 0
    if (!is.null(fx)) {
        if (length(fx) != length(wtaxa)) {
            stop("fx and wtaxa must be the same length")
        }
        if (length(fx) != length(wsamples)) {
            stop("fx and wsamples must be the same length")
        }
        for (e in 1:length(fx)) {
            fx.mtx <- add.fx.to.matrix(samples,
                                       taxa,
                                       wtaxa[[e]],
                                       wsamples[[e]], # from divided$env/dset
                                       fx.mtx,
                                       fx[[e]])
        }
    }
    fx.mtx
}

master.fx.matrix <- function(divided,
                             taxa,
                             env.fx,
                             env.taxa,
                             dset.fx,
                             dset.taxa) {
    env.mtx <- build.fx.matrix(divided$metadata$sample,
                               taxa,
                               env.fx,
                               env.taxa,
                               divided$env)
    dset.mtx <- build.fx.matrix(divided$metadata$sample,
                                taxa,
                                dset.fx,
                                dset.taxa,
                                divided$dset)
    env.mtx + dset.mtx
}

get.taxon.names <- function(n.taxa) { paste0("taxon", 1:n.taxa) }
get.sample.names <- function(n.taxa) { paste0("sample", 1:n.taxa) }

master.mtx.wrapper <- function(n.taxa,
                               divided,
                               env.n.affected,
                               dset.n.affected,
                               ...) {
    fx <- pick.env.dset.fx(n.taxa,
                           env.n.affected,
                           dset.n.affected,
                           divided,
                           ...)
    mtx <- master.fx.matrix(divided,
                            get.taxon.names(n.taxa),
                            fx$env.fx,
                            fx$env.taxa,
                            fx$dset.fx,
                            fx$dset.taxa)
    return(list(mtx=mtx, fx=fx))
}

pick.env.dset.fx <- function(n.taxa,
                             env.n.affected,
                             dset.n.affected,
                             divided,
                             env.mean = 0,
                             env.sd = 0,
                             dset.mean = 0,
                             dset.sd = 0) {
    taxa <- get.taxon.names(n.taxa)
    env <- divided$env
    dset <- divided$dset
    n.env <- length(env)
    n.dset <- length(dset)
    if (length(env.n.affected) == 1) {
        env.n.affected <- rep(env.n.affected, n.env)
        names(env.n.affected) <- names(env)
    }
    if (length(dset.n.affected) == 1) {
        dset.n.affected <- rep(dset.n.affected, n.dset)
        names(dset.n.affected) <- names(dset)
    }
    env.taxa <- lapply(env.n.affected, function(x) sample(taxa, x))
    dset.taxa <- lapply(dset.n.affected, function(x) sample(taxa, x))
    env.fx <- rnorm(n.env, env.mean, env.sd)
    dset.fx <- rnorm(n.dset, dset.mean, dset.sd)
    names(env.fx) <- names(env.n.affected)
    names(dset.fx) <- names(dset.n.affected)
    return(list(env.fx=env.fx, dset.fx=dset.fx, env.taxa=env.taxa, dset.taxa=dset.taxa))
}

divide.samples <- function(n.samples,
                           n.dsets,
                           n.envirs) {
    samples <- paste0("sample", 1:n.samples)
    dsets <- paste0("dset", 1:n.dsets)
    envirs <- paste0("env", 1:n.envirs)
    dset.map <- sample(dsets, n.samples, replace=TRUE)
    env.map <- sample(envirs, n.samples, replace=TRUE)
    dset.list <- split(samples, dset.map)
    env.list <- split(samples, env.map)
    metadata <- data.frame(
        sample=samples,
        dataset=dset.map,
        env=env.map)
    return(list(dset=dset.list,
                env=env.list,
                metadata=metadata))
}

rep.col <- function(vec, nc) {
    matrix(rep(vec, each=nc), ncol=nc, byrow=TRUE)
}

generate.fake.abd.meta <- function(n.samples=100,
                                   n.taxa=1000,
                                   n.envs=3,
                                   n.dsets=2,
                                   n.reads=1e6,
                                   env.fx.sizes=c(0,1,-1),
                                   dset.fx.sizes=c(0,0.1,-0.1),
                                   env.frac.affected=0.25,
                                   dset.frac.affected=0.1,
                                   env.sd=0,
                                   dset.sd=0,
                                   prev.dist=c(-4,2),
                                   make.16s=FALSE,
                                   tag.length=100,
                                   ...) {
    env.n.affected <- round(n.taxa * env.frac.affected)
    dset.n.affected <- round(n.taxa * dset.frac.affected)
    divided <- divide.samples(n.samples, n.dsets, n.envs)
    master.mtx <- master.mtx.wrapper(n.taxa,
                                     divided,
                                     env.n.affected,
                                     dset.n.affected,
                                     env.mean=env.fx.sizes,
                                     env.sd=env.sd,
                                     dset.mean=dset.fx.sizes,
                                     dset.sd=dset.sd)
    sim <- make.sim(n.taxa,
                    n.samples,
                    avg.reads=n.reads,
                    prev.dist=prev.dist,
                    master.mtx$mtx)
    abd.meta <- list(mtx=sim$sim,
                     metadata=divided$metadata,
                     fx=master.mtx,
                     com=sim$com,
                     prev=sim$prev)
    if (make.16s) {
        tagged <- make.simulated.denoised.data(abd.meta$mtx,
                                               tag.length=tag.length,
                                               ...)
        abd.meta$mtx <- tagged$mtx
        abd.meta$map <- tagged$map
        abd.meta$n <- tagged$n
    }
    sanity.check.abundance(abd.meta$mtx, ...)
    sanity.check.metadata(abd.meta$metadata, ...)
    abd.meta <- harmonize.abd.meta(abd.meta, ...)
    # binarize to save memory usage since we care about pres/abs
    abd.meta$mtx <- Matrix::Matrix(abd.meta$mtx > 0)
    abd.meta
}

write.test.tabular <- function(abd.meta,
                               abdfile="test-abundance.tab",
                               metafile="test-metadata.tab",
                               prep=TRUE,
                               ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    af <- file.path(opts('in_dir'), abdfile)
    mf <- file.path(opts('in_dir'), metafile)
    if (prep) abd.meta$mtx <- prep.mtx.for.write(abd.meta$mtx)
    write.table(as.matrix(abd.meta$mtx),
                af,
                sep='\t',
                row.names=FALSE,
                quote=FALSE)
    write.table(abd.meta$metadata,
                mf,
                sep='\t',
                row.names=FALSE,
                quote=FALSE)
}


prep.mtx.for.write <- function(mtx, initial.octo=FALSE) {
    mtx2 <- data.frame(as.matrix(mtx) * 1)
    if (initial.octo) {
        mtx2 <- cbind(data.frame(`#OTU ID`=rownames(mtx), check.names=FALSE),
                      mtx2)
    } else {
        mtx2 <- cbind(data.frame(OTU_ID=rownames(mtx), check.names=FALSE),
                      mtx2)
    }
    rownames(mtx2) <- NULL
    mtx2
}

write.test.biom <- function(abd.meta,
                            overwrite=FALSE,
                            ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    cn <- colnames(abd.meta$metadata)
    colnames(abd.meta$metadata)[which(cn=="sample")] <- "#SampleID"
    abd.meta$mtx <- prep.mtx.for.write(abd.meta$mtx, initial.octo=TRUE)
    td <- tempdir()
    af <- file.path(td, "test-abundance.tab")
    mf <- file.path(td, "test-metadata.tab")
    tmp.bf <- tempfile("test-", fileext=".biom")
    bf <- file.path(opts('in_dir'), opts('biom_file'))
    if (overwrite) file.remove(tmp.bf)
    if (overwrite) file.remove(bf)
    write.test.tabular(abd.meta,
                       in_dir=tempdir(),
                       abdfile="test-abundance.tab",
                       metafile="test-metadata.tab",
                       prep=FALSE)
    system2(file.path(opts('biom_dir'), "biom"),
            args = c("convert",
                     "-i",
                     af,
                     "-o",
                     tmp.bf,
                     "--table-type=\"OTU table\"",
                     "--to-hdf5"))
    system2(file.path(opts('biom_dir'), "biom"),
            args = c("add-metadata",
                     "-i",
                     tmp.bf,
                     "-o",
                     bf,
                     "--sample-metadata-fp",
                     mf))
}


#--- Dummy 16S ---#

random.species.from.file <- function(n.taxa,
                                     tag.length=100,
                                     ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    species.list <- seqinr::read.fasta(
                                file.path(opts('data_dir'),
                                          opts('burst_16sfile')),
                                seqtype='DNA',
                                as.string=TRUE)
    species.map <- sample(1:length(species.list), n.taxa, replace=FALSE)
    taxon.seqs <- sapply(species.map, function(x) {
        if (is.null(tag.length)) {
            species.list[x]
        } else {
            substr(species.list[x], 1, tag.length)
        }
    })
    #taxon.real.names <- sapply(species.map, function(x) {
    #    attr(species.list[x], 'name')
    #})
    taxon.real.names <- sapply(species.map, function(x) {
        strsplit(attr(species.list[[x]], 'Annot'), " ")[[1]][3]
    })
    return(list(seqs=taxon.seqs,
                names=taxon.real.names,
                map=species.map))
}

make.simulated.denoised.data <- function(mtx, tag.length=100, ...) {
    taxa <- rownames(mtx)
    rs <- random.species.from.file(n.taxa=length(taxa),
                                   tag.length=tag.length,
                                   ...)
    rownames(mtx) <- rs$seqs
    names(rs$seqs) <- taxa
    return(list(mtx=mtx,
                map=rs$seqs,
                n=rs$names))
}

