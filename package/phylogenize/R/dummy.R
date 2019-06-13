# Functions for creating and writing dummy data

#' Simulate microbe counts given multinomial probabilities and number of reads.
#'
#' @param com Multinomial probability distribution for all taxa.
#' @param nr Read depth.
#' @return Count matrix.
#' @keywords internal
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

#' Calculate prevalence (with additive smoothing) for taxa given an abundance matrix.
#'
#' @param mtx Abundance matrix (rows: taxa; columns: samples).
#' @return Prevalence with additive smoothing.
#' @keywords internal
simple.prevalence <- function(mtx) {
    apply(cbind(1 * (mtx > 0), 1, 0), 1, mean)
}

#' Normalize a matrix by columns (i.e., force them to sum to one).
#'
#' @param mtx Numeric matrix.
#' @return A numeric matrix where each column has been divided by its sum.
#' @keywords internal
colnorm <- function(mtx) {
    apply(mtx, 2, function(x) {
        if (all(x == 0)) return(x)
        x/sum(x)
    })
}

#' Construct a simulation of microbial sequencing data.
#'
#' @param n.taxa Integer; number of taxa to simulate.
#' @param n.samples Integer; number of samples to simulate.
#' @param avg.reads Number or numeric vector giving read depth(s) per sample.
#' @param prev.dist Two-element numeric vector giving beta distribution
#'     parameters for the prevalence distribution.
#' @param effect.size.mtx An optional numeric matrix of true effect sizes (in
#'     logit space).
#' @param alphas Optionally, a numeric vector of Dirichlet alpha parameters to
#'     use for Dirichlet-Multinomial simulation.
#' @return A list:
#'   \item{sim}{A simulated abundance matrix}
#'   \item{com}{Multinomial distribution used to simulate abundances}
#'   \item{prev}{Distribution of taxon prevalences.}
#' @keywords internal
make.sim <- function(n.taxa,
                     n.samples,
                     avg.reads, # can be a vector
                     prev.dist=c(-4, 2),
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
        prev.vec <- rbinom(n.samples, prob=random.prevalences[nt, , drop=FALSE], size = 1)
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

#' Add a true effect to a matrix of effect sizes (or generate a new one).
#'
#' @param samples Vector of strings representing the sample names.
#' @param taxa Vector of strings representing taxon names.
#' @param which.taxa Vector of strings giving subset of taxa to be affected.
#' @param which.samples Vector of strings giving which samples to be affected.
#' @param mtx Optionally, an existing matrix (if not provided, generate a new
#'     one).
#' @param fx Size of effect to be added.
#' @return A new effect size matrix.
#' @keywords internal
add.fx.to.matrix <- function(samples,
                             taxa,
                             which.taxa,
                             which.samples,
                             mtx=0,
                             fx=1) {
    fx.mtx <- matrix(nc=length(samples),
                     nr=length(taxa),
                     0,
                     dimnames=list(taxa, samples))
    fx.mtx[which.taxa, which.samples] <- fx
    mtx + fx.mtx
}


#' Assemble an effect-size matrix from a set of distinct effects.
#'
#' @param samples Vector of strings representing the sample names.
#' @param taxa Vector of strings representing taxon names.
#' @param wtaxa List of vectors of strings giving the subset of taxa to be
#'     affected in each distinct effect.
#' @param wsamples List of vectors of strings giving which samples to be
#'     affected, one for every distinct effect.
#' @param fx Numeric vector of sizes of effects to be added.
#' @return An effect size matrix.
#' @keywords internal
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
                                        # from divided$env/dset
                                       wsamples[[e]],
                                       fx.mtx,
                                       fx[[e]])
        }
    }
    fx.mtx
}

#' Construct a master effect size matrix given environmental and dataset
#' effects (wraps \code{master.fx.matrix}).
#'
#' @param divided Output of \code{divide.samples}.
#' @param taxa Vector of strings representing taxon names.
#' @param env.fx Named numeric vector giving effect sizes per environment.
#' @param env.taxa List of vectors of strings giving the subset of taxa to be
#'     affected in each distinct environment effect.
#' @param dset.fx Named numeric vector giving effect sizes per dataset.
#' @param dset.taxa List of vectors of strings giving the subset of taxa to be
#'     affected in each distinct dataset effect.
#' @return An effect size matrix incorporating both environmental and dataset
#'     effects.
#' @keywords internal
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

#' Generate fake names for taxa or samples.
#' @return A vector of fake names.
#' @name NameFaker
#' @keywords internal
NULL

#' @rdname NameFaker
#' @param n.taxa Number of taxa.
#' @keywords internal
get.taxon.names <- function(n.taxa) { paste0("taxon", 1:n.taxa) }

#' @rdname NameFaker
#' @param n.samples Number of samples.
#' @keywords internal
get.sample.names <- function(n.samples) { paste0("sample", 1:n.samples) }

#' Wrapper function for generating an effect size matrix automatically.
#'
#' @param n.taxa Number of taxa to simulate
#' @param divided Output of \code{divide.samples}
#' @param env.n.affected Number of taxa to be affected by environment effects.
#' @param dset.n.affected Number of taxa to be affected by dataset effects.
#' @return A list:
#'   \item{mtx}{An effect size matrix.}
#'   \item{fx}{The chosen, randomly generated effect sizes (see
#'     \code{pick.env.dset.fx}.}
#' @keywords internal
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

#' Function for picking taxa to simulate.
#'
#' @param n.taxa Number of taxa to simulate
#' @param divided Output of \code{divide.samples}
#' @param env.n.affected Number of taxa to be affected by environment effects.
#' @param dset.n.affected Number of taxa to be affected by dataset effects.
#' @return A list:
#'   \item{mtx}{An effect size matrix.}
#'   \item{fx}{The chosen, randomly generated effect sizes (see
#'     \code{pick.env.dset.fx}.}
#' @keywords internal
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

#' Divide a given number of samples into a given number of datasets and
#' environments.
#'
#' @param n.samples Number of samples to divide
#' @param n.dsets Number of datasets to divide samples into
#' @param n.envirs Number of environments to divide samples into
#' @return A list:
#'   \item{dset}{A list of vectors indicating which samples go with which
#'     datasets (see \code{?split}).}
#'   \item{env}{Like \code{dset} for environments.}
#'   \item{metadata}{A data frame of metadata associating samples with datasets
#'     and environments.}
#' @keywords internal
divide.samples <- function(n.samples,
                           n.dsets,
                           n.envirs) {
    samples <- get.sample.names(n.samples)
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

#' Convenience function to repeat a vector as columns.
#'
#' @param vec Vector to repeat.
#' @param nc Number of times to repeat vector.
#' @return A matrix consisting of \code{nc} repeats of the column-vector
#'     \code{vec}.
#' @keywords internal
rep.col <- function(vec, nc) {
    matrix(rep(vec, each=nc), ncol=nc, byrow=TRUE)
}

#' Wrapper script to generate fake data and metadata.
#'
#' @param n.samples Positive integer: how many samples should be simulated?
#' @param n.taxa Positive integer: how many taxa should be simulated?
#' @param n.envs Positive integer: how many environments should be simulated?
#' @param n.dsets Positive integer: how many datasets should be simulated?
#' @param n.reads Number or numeric vector giving read depth(s) per sample.
#' @param env.fx.sizes A numeric vector of effect sizes, one per environment.
#' @param dset.fx.sizes A numeric vector of effect sizes, one per dataset.
#' @param env.frac.affected Proportion of taxa that will be affected by an
#'     environment.
#' @param dset.frac.affected Proportion of taxa that will be affected by a
#'     dataset effect.
#' @param env.sd Number: when sampling environment effect sizes, use this
#'     standard deviation.
#' @param dset.sd Number: when sampling dataset effect sizes, use this standard
#'     deviation.
#' @param prev.dist Numeric vector of length 2 giving beta parameters for the
#'     distribution of taxon prevalences.
#' @param make.16s Simulate a 16S dataset instead of a WGS-MIDAS dataset.
#' @param tag.length If simulating 16S data, amplicon sequence variants will be
#'     this length.
#' @return A list (comparable to what is read in by *phylogenize*):
#'   \item{mtx}{Simulated abundance matrix}
#'   \item{metadata}{Simulated metadata data frame}
#' @keywords internal
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
                                   t.names=NULL,
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
    if (!is.null(t.names)) {
        rownames(abd.meta$mtx) <- sample(t.names,
                                         nrow(abd.meta$mtx),
                                         replace=FALSE)
    }
    abd.meta
}

#' Write and then read in tabular data.
#'
#' @param abd.meta List containing a matrix of abundance values, \code{mtx}, and
#'     a data frame of metadata, \code{metadata}.
#' @param abdfile String: write abundance matrix to this file.
#' @param metafile String: write metadata data frame to this file.
#' @param prep Boolean: convert matrix to a data frame and add a row name title.
#' @return Return value of \code{base::write.table}.
#' @keywords internal
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

#' Convert a matrix to a data frame of numeric values, then add a column for row
#' names.
#'
#' @param mtx Matrix to convert
#' @param initial.octo Boolean: if true, the name of the column with row names
#'     will be \code{#OTU ID}; if false, the name of the column will be
#'     \code{OTU_ID}. The former is necessary for the BIOM format.
#' @return A data frame derived from \code{mtx}.
#' @keywords internal
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

#' Write and then read in tabular data.
#'
#' Some global options that may be helpful include:
#' \describe{
#'   \item{biom_file}{String. Name of BIOM abundance-and-metadata file (in this
#'   case, to write to disk).}
#'   \item{in_dir}{String. Path to input directory (i.e., where to look for
#'   input files). In this case, this is the directory where the BIOM file will
#'   be *written*.}
#'   \item{biom_dir}{String. Path to BIOM command-line executables. These are
#'   necessary to write the BIOM file and add the metadata.}
#' }
#'
#' @param abd.meta List containing a matrix of abundance values, \code{mtx}, and
#'     a data frame of metadata, \code{metadata}.
#' @param overwrite Boolean: if a BIOM file exists at the location specified,
#'     overwrite it?
#' @keywords internal
write.test.biom <- function(abd.meta,
                            overwrite=FALSE,
                            ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    cn <- colnames(abd.meta$metadata)
    colnames(abd.meta$metadata)[which(cn==opts('sample_column'))] <- "#SampleID"
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

#' Simulate amplicon sequence variants (ASVs) from a database of 16S sequences.
#'
#' Some global options that may be helpful include:
#' \describe{
#'   \item{data_dir}{String. Path to directory containing the data files
#'   required to perform a \emph{phylogenize} analysis. Here, this is where the
#'   16S database is located.}
#'   \item{burst_16sfile}{String. Filename of the 16S FASTA database that maps
#'   back to MIDAS species.}
#' }
#'
#' @param n.taxa Positive integer: how many ASVs to simulate?
#' @param tag.length Positive integer: how long should the ASVs be?
#' @return A list:
#'   \item{seqs}{String vector: sequences of simulated ASVs.}
#'   \item{names}{String vector: names of the taxa to which the ASVs "should"
#'     map.}
#'   \item{species.map}{Numeric vector where the $i$'th element gives the number
#'     of the species that taxon $i$ was mapped to.}
#' @keywords internal
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

#' Wrapper function to simulate data "denoised" by an algorithm like DADA2 or
#' Deblur.
#'
#' Extra arguments get passed to \code{random.species.from.file}.
#'
#' @param mtx An input abundance matrix.
#' @param tag.length Positive integer: how long should amplicon sequence
#'     variants (ASVs) be?
#' @return A list:
#'   \item{mtx}{The input matrix, with renamed rows.}
#'   \item{map}{String vector of ASVs, named with the original row names of
#'     \code{mtx}.}
#'   \item{n}{Names of the MIDAS taxa to which the rows of \code{mtx} "should"
#'     map.}
#' @keywords internal
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

#' Make a dummy gene-to-species matrix by subsampling.
#'
#' For each matrix in the input list, find the top percentile most-variable
#' genes (using formula for binomial distribution), then subsample the matrix,
#' returning at most \code{n} genes. Fewer genes will be returned if the total
#' number of original genes or the number above the percentile cutoff is lower
#' than the desired \code{n}. Matrices with 1 or fewer genes after sampling will
#' be dropped.
#'
#' @param g2s.matrices A list of sparse binary matrices, one per phylum.
#' @param n Integer: how many genes per phylum to sample?
#' @param fp Double: what proportion of highest-variance genes to sample from?
#' @param minN Integer: minimum number of genes to return in a matrix
#' @keywords internal
dummy.g2s <- function(g2s.matrices, n, fp=0.1, minN=2) {
    new.g2s <- lapply(g2s.matrices, function(m) {
        ng <- ncol(m)
        if (ng < n) { s <- ng } else { s <- n }
        vars <- apply(m, 2, function(col) {
            length(col) * mean(col) * (1 - mean(col))
        })
        var.pctl <- order(vars, decreasing=TRUE) / length(vars)
        top.var <- which(var.pctl <= fp)
        ntv <- length(top.var)
        if (ntv < s) {
            s <- ntv
        }
        m[, sample(top.var, s), drop=FALSE]
    })
    Filter(function(x) {
        if (is.null(x)) return(FALSE)
        if (is.vector(x)) return(FALSE)
        if (ncol(x) >= minN) return(TRUE)
        return(FALSE)
    }, new.g2s)
}

#' Downsample trees.
#'
#' @param trees List of trees.
#' @param n Number of taxa to retain (maximum)
#' @keywords internal
dummy.trees <- function(trees, n=50) {
    lapply(trees, function(tr) {
        tips <- tr$tip.label
        s <- min(n, length(tips))
        keep.tips(tr, sample(tips, s))
    })
}

#' Wrapper to generate a low-memory, fast-running dummy database. Make sure
#' \code{Matrix} is loaded.
#'
#' @param nt Number of taxa to sample per phylum.
#' @param ng Number of genes to sample per phylum.
#' @param fp Fraction of most-variable genes to sample.
#' @param minN Only keep phyla with at least this many genes.
#' @keywords internal
generate.test.pzdb <- function(nt=75, ng=50, fp=0.1, minN=2, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    pz.db <- import.pz.db(...)
    dtr <- dummy.trees(pz.db$trees, nt)
    for (n in names(pz.db$gene.presence)) {
        d <- dtr[[n]]$tip.label
        if (methods::is(pz.db$gene.presence[[n]], "Matrix")) {
            dd <- intersect(d, colnames(pz.db$gene.presence[[n]]))
            if (length(dd) > 2) {
                pz.db$gene.presence[[n]] <- pz.db$gene.presence[[n]][, dd, drop=FALSE]
            }
        } else if (is.vector(pz.db$gene.presence[[n]])) {
            if (all(d %in% names(pz.db$gene.presence[[n]]))) {
                pz.db$gene.presence[[n]] <- pz.db$gene.presence[[n]][d]
            }
        } # otherwise don't change
    }
    gp <- dummy.g2s(pz.db$gene.presence, ng, fp, minN)
    saveRDS(gp, file.path(opts('data_dir'), "test-gene-presence-binary.rds"))
    saveRDS(dtr, file.path(opts('data_dir'), "test-trees.rds"))
    write.csv(pz.db$taxonomy, file.path(opts('data_dir'), "test-taxonomy.csv"))
}
