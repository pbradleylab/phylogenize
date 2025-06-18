#' Fit phylogenetic (or linear) models, or perform POMS.
#'
#' @param taxa Character vector giving the names of the taxa.
#' @param pheno Named numeric vector giving phenotype values per species. Can be
#'   NULL for POMS.
#' @param tree Either a single tree covering all species, or a list of per-taxon
#'   trees.
#' @param proteins Named list of gene presence/absence matrices, per taxon.
#' @param clusters Named list of character vectors of species IDs, per taxon.
#' @param method A function that returns a length-2 numeric vector of
#'   effect-size and p-value (see, e.g., \code{phylolm.fx.pv} or
#'   \code{lm.fx.pv}). Can be NULL for POMS.
#' @param restrict.figfams Optionally, a character vector giving a subset of
#'   genes to test.
#' @param drop.zero.var Boolean giving whether to drop genes that are always
#'   present or always absent in a particular taxon.
#' @param only.return.names Boolean giving whether to just return the names of
#'   genes to be tested (for debugging).
#' @param abd.meta List containing abundances and metadata (only required for
#'   POMS)
#' @param poms Boolean; whether to run POMS instead of phylolm.
#' @return Named list of p-value and effect-size matrices, one per taxon.
#' @export result.wrapper.plm
result.wrapper.plm <- function(taxa,
                               pheno,
                               tree,
                               proteins,
                               clusters,
                               method = phylolm.fx.pv,
                               restrict.figfams = NULL,
                               drop.zero.var = FALSE,
                               only.return.names = FALSE,
                               abd.meta = FALSE,
                               poms = FALSE,
                               ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    lapply.across.names(taxa, function(p) {
        message(p)
        if (class(tree) == "phylo") {
            tr <- tree
        } else if (class(tree) == "list") {
            tr <- tree[[p]]
        } else {
            stop("tree must be either an object of class phylo or a list")
        }
        valid <- Reduce(intersect,
                        list(colnames(proteins[[p]]),
                             clusters[[p]],
                             tr$tip.label))
        if (!is.null(pheno)) {
            valid <- intersect(valid, names(pheno))
        }
        
	if (is.null(restrict.figfams)) {
            restrict.figfams <- rownames(proteins[[p]])
        } else {
            restrict.figfams <- intersect(rownames(proteins[[p]]),
                                          restrict.figfams)
        }
        
	if (drop.zero.var) {
            fvar <- apply(proteins[[p]][, valid, drop=FALSE], 1, var)
            message(paste0(sprintf("%.01f",
                                   100 * mean(na.omit(fvar == 0))),
                           "% dropped [no variance]"))
            restrict.figfams <- intersect(restrict.figfams,
                                          rownames(proteins[[p]])[
                                              which(fvar > 0)])
        }
        
	if (!only.return.names) {
		if (poms) {
			matrix.POMS(tr,
				    proteins[[p]],
				    abd.meta,
				    restrict.taxa=valid,
				    restrict.ff=restrict.figfams,
				    ...)
		} else {
            matrix.plm(tr,
                       proteins[[p]],
                       pheno,
                       method = method,
                       restrict.taxa = valid,
                       restrict.ff = restrict.figfams,
                       ...)
		}
            
        } else {
            restrict.figfams
	}
})
}

#' Perform POMS modeling for a single clade.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{ncl}{Integer. Number of cores to use for parallel computation.
#'   Default: 1}
#'   \item{env_column}{String. Name of column in metadata file containing the
#'   environment annotations.}
#'   \item{sample_column}{String. Name of column in metadata file containing the
#'   sample names.}
#'   \item{which_envir}{String. Environment in which to calculate prevalence,
#'   specificity, or differential abundance. Must match annotations in metadata.}
#' }
#'
#' @param tree A tree relating taxa within a clade.
#' @param mtx Gene presence/absence matrix.
#' @param abd.meta A list giving an abundance matrix and metadata.
#' @param restrict.ff Optionally, a character vector giving a subset of genes to
#'   test.
#' @param restrict.taxa Optionally, a character vector giving a subset of taxa
#'   to test.
#' @param poms_min_tips Minimum number of tips (min_num_tips); default 2.
#' @param poms_min_func Minimum number of function instances
#'   (min_func_instances); default 2.
#' @param poms_pseudocount Pseudocount value for POMS; default 0.5.
#' @return Matrix of "effect-sizes" (row 1) and p-values (row 2) per gene
#'   (columns). Here, we "fake" effect sizes by taking the log2-ratio of
#'   num_FSNs_group1_enrich and num_FSNs_group2_enrich with a 0.5 pseudocount.
#' @export
matrix.POMS <- function(tree,
                        mtx,
                        abd.meta,
                        restrict.taxa=NULL,
                        restrict.ff=NULL,
                        poms_min_tips=2,
                        poms_min_func=2,
                        poms_pseudocount=0.5,
                        ...) {
    message("POMS")
    
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    E <- opts('env_column')
    S <- opts('sample_column')
    envir <- opts('which_envir')
    cores <- opts('ncl')
    
    # Note: dataset column is ignored
    poms_group1 <- abd.meta$meta[[S]][abd.meta$meta[[E]] == envir]
    poms_group2 <- abd.meta$meta[[S]][abd.meta$meta[[E]] != envir]
    
    if (is.null(restrict.taxa)) restrict.taxa <- colnames(mtx)
    if (is.null(restrict.ff)) restrict.ff <- rownames(mtx)
    
    phylotype_df <- data.frame(as.matrix(t(mtx[restrict.ff, ])))
    
    if (length(unique(as.numeric(abd.meta$mtx))) <= 2) {
        pz.error(paste0(
            "Abundance matrix has two or fewer unique values; ",
            "suggests matrix is binary (which will not work for POMS). ",
            "Please double check that you passed in the right data."
        ))
    }
    tree_nodes <- ape::makeNodeLabel(tree)
    # Don't check names because this will potentially change them, meaning 
    # they won't match the metadata anymore
    poms_output <- tryCatch({
			    POMS::POMS_pipeline(
        abun=data.frame(
            as.matrix(abd.meta$mtx),
            check.names=FALSE),
        func=phylotype_df,
        tree=tree_nodes,
        group1_samples=poms_group1,
        group2_samples=poms_group2,
        ncores=cores,
        min_num_tips=poms_min_tips,
        min_func_instances=poms_min_func,
        pseudocount=poms_pseudocount)},
			    error = function(e) {
				    pz.warning(e)
				    NA
			    })
    # handle error
    if ((length(poms_output)==0) && is.na(poms_output)) {
	    result_mtx <- rbind(Estimate = rep(NA, ncol(phylotype_df)),
			 p.value = NA,
			 StdErr = NA,
			 df = NA)
	    colnames(result_mtx) <- colnames(phylotype_df)
	    return(result_mtx)
    }
    poms_enrichments <- (log2(
        (poms_output$results[, "num_FSNs_group1_enrich"] + 0.5) /
            (poms_output$results[, "num_FSNs_group2_enrich"] + 0.5)))
    poms_pvals <- poms_output$results[, "multinomial_p"]
    result_mtx <- rbind(Estimate = poms_enrichments,
                 p.value = poms_pvals,
                 StdErr = NA,
                 df = NA)
    colnames(result_mtx) <- colnames(phylotype_df)
    result_mtx
}

#' Fit phylogenetic (or linear) models (single core version, single taxon).
#'
#' @param gene.matrix Gene presence/absence matrix.
#' @param tree Phylogeny relating taxa.
#' @param pheno Named numeric vector giving phenotype values per taxon.
#' @param taxon.name Name of taxon being considered.
#' @param method A function that returns a length-2 numeric vector of
#'     effect-size and p-value (see, e.g., \code{phylolm.fx.pv} or
#'     \code{lm.fx.pv}).
#' @param restrict.ff Optionally, a character vector giving a subset of
#'     genes to test.
#' @param remove.low.variance Boolean giving whether to drop genes that are
#'     always present or always absent in a particular taxon.
#' @param use.for.loop Boolean giving whether to use a for loop instead of a
#'     pbapply.
#' @return Named list of p-value and effect-size matrices, one per taxon.
#' @export
nonparallel.results.generator <- function(gene.matrix,
                                         tree,
                                         taxa,
                                         pheno,
                                         taxon.name="TestFamily",
                                         method=phylolm.fx.pv,
                                         restrict.ff=NULL,
                                         remove.low.variance=TRUE,
                                         use.for.loop=TRUE,
                                         ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    message(taxon.name)
    restrict.taxa <- Reduce(intersect, list(colnames(gene.matrix),
                                            tree$tip.label,
                                            names(pheno),
                                            taxa))
    if (is.null(restrict.ff)) restrict.ff <- rownames(gene.matrix)
    if (remove.low.variance) {
        restrict.lv <- names(which(
            apply(gene.matrix[restrict.ff, restrict.taxa, drop=FALSE],
                  1,
                  var) > 0))
        restrict.both <- intersect(restrict.ff, restrict.lv)
        message(sprintf("removing zero-variance genes (%d): %0.2f%% removed",
                        length(restrict.ff) - length(restrict.lv),
                        100 * (1 - (length(restrict.both) /
                                    length(restrict.ff)))
                        ))
        restrict.ff <- restrict.both
    }
    restrict.tree <- keep.tips(tree, restrict.taxa)
    if (use.for.loop) {
        ans <- matrix(nr = 4, nc = length(restrict.ff), NA)
        dimnames(ans) <- list(c("Estimate", "p.value", "StdErr", "df"),
                              restrict.ff)
        message(paste0(c('|',rep('-', 48), "|"), collapse=''))
        progress <- txtProgressBar(min = 1,
                                   max = length(restrict.ff),
                                   initial = 1,
                                   char = '=',
                                   style = 3,
                                   width = 50,
                                   file = stderr())
        pheno.restrict <- pheno[restrict.taxa]
        # replaces an apply to avoid copying
        for (fn in 1:length(restrict.ff)) {
            setTxtProgressBar(progress, fn)
            ans[, fn] <- method(gene.matrix[restrict.ff[fn], restrict.taxa],
                                pheno.restrict,
                                restrict.tree,
                                meas_err=opts('meas_err'))
        }
        close(progress)
        ans
    } else {
        pbapply::pbapply(gene.matrix[restrict.ff, restrict.taxa, drop=FALSE],
                         1,
                         function(m) {
                             method(m, pheno[restrict.taxa], restrict.tree,
                                    meas_err=opts('meas_err'))
                         })
    }
}

#' Wrapper around \emph{phylolm} that returns just the effect size and p-value.
#'
#' @param m Named numeric vector of gene presence/absences per taxon.
#' @param p Named numeric vector of phenotype values per taxon.
#' @param tr Phylogeny relating taxa (class \code{"phylo"}).
#' @param coefname Which coefficient from the phylolm to return?
#' @param restrict If not NULL, a character vector of taxa to consider.
#' @return Length-2 numeric vector with names \code{"Estimate"} and
#'     \code{"p.value"}. If there is an error in \code{phylolm}, the values of
#'     this vector will be \code{c(NA, NA)}.
#' @export
phylolm.fx.pv <- function(m, p, tr, coefname="mTRUE", restrict=NULL,
                          meas_err=FALSE) {
    # This seems redundant, but we can avoid touching the giant protein matrix
    # this way and therefore causing an expensive copy
    if (!is.null(restrict)) {
        p <- p[restrict]
        m <- m[restrict]
    }
    fx.pv <- tryCatch({
        plm <- phylolm::phylolm(p ~ m, phy=tr, measurement_error=meas_err)
        coef <- summary(plm)$coefficients
        c(coef[coefname, c("Estimate", "p.value", "StdErr")],
          df=(plm$n - plm$d))
    }, error = function(e) {
        pz.warning(paste(e))
        c(Estimate = NA, p.value = NA, StdErr = NA, df = NA)
    })
    fx.pv
}

#' Wrapper around \emph{lm} that returns just the effect size and p-value.
#'
#' @param m Named numeric vector of gene presence/absences per taxon.
#' @param p Named numeric vector of phenotype values per taxon.
#' @param tr Phylogeny relating taxa (class \code{"phylo"}).
#' @return Length-2 numeric vector with names \code{"Estimate"} and
#'     \code{"p.value"}. If there is an error in \code{phylolm}, the values of
#'     this vector will be \code{c(NA, NA)}.
#' @export
lm.fx.pv <- function(m, p, tr, coefname="mTRUE", restrict=NULL,
                     meas_err=FALSE) {
    if (!is.null(restrict)) {
        p <- p[restrict]
        m <- m[restrict]
    }
    fx.pv <- tryCatch({
        lm <- lm(p ~ m)
        coef <- summary(lm)$coefficients
        pair <- coef[coefname, c("Estimate", "Pr(>|t|)", "Std. Error")]
        names(pair) <- c("Estimate", "p.value", "StdErr")
        c(pair, df = df.residual(lm))
    }, error = function(e) {
        pz.warning(paste(e))
        c(Estimate = NA, p.value = NA, StdErr = NA, df = NA)
    })
    fx.pv
}

#' Wrapper to return random effect sizes and p-values; useful for testing.
#'
#' @param m Named numeric vector of gene presence/absences per taxon.
#' @param p Named numeric vector of phenotype values per taxon.
#' @param tr Phylogeny relating taxa (class \code{"phylo"}).
#' @param coefname Which coefficient from the phylolm to return (meaningless
#'     here)?
#' @param restrict If not NULL, a character vector of taxa to consider
#'     (meaningless here).
#' @return Length-2 numeric vector with names \code{"Estimate"} and
#'     \code{"p.value"}, with distributions N(0,1) and U(0,1), respectively.
#' @keywords internal
rnd.fx.pv <- function(m, p, tr, coefname="mTRUE", restrict=NULL,
                      meas_err=FALSE) {
    return(c(Estimate = rnorm(n=1, 0, 1),
             p.value = runif(n=1, 0, 1)))
}

#' Wrapper to return random effect sizes and p-values with a bias; useful for
#' testing.
#'
#' @param m Named numeric vector of gene presence/absences per taxon.
#' @param p Named numeric vector of phenotype values per taxon.
#' @param tr Phylogeny relating taxa (class \code{"phylo"}).
#' @param coefname Which coefficient from the phylolm to return (meaningless
#'     here)?
#' @param restrict If not NULL, a character vector of taxa to consider
#'     (meaningless here).
#' @param pos.pct Fraction of times to return a positive effect.
#' @param pos.fx Shift in mean for positive effects, in standard deviations.
#' @param pos.beta Length-2 numeric vector of beta parameters for true effect
#'     p-values.
#' @return Length-2 numeric vector with names \code{"Estimate"} and
#'     \code{"p.value"}, with distributions N(0,1) and U(0,1), respectively.
#' @keywords internal
rndpos.fx.pv <- function(m, p, tr, coefname="m", restrict=NULL,
                         pos.pct=0.1,
                         pos.fx=2,
                         pos.beta=c(shape1=1, shape2=20),
                         meas_err=FALSE,
                         se=0.1) {
    is.pos <- runif(n=1, min=0, max=1)
    if (is.pos <= pos.pct) {
        est <- rnorm(n=1, pos.fx, 1)
        pv <- rbeta(n=1, shape1=pos.beta[1], shape2=pos.beta[2])
    } else {
        est <- rnorm(n=1, 0, 1)
        pv <- runif(n=1, 0, 1)
    }
    return(c(Estimate = est, p.value = pv, StdErr = se))
}

#' Wrapper around \code{phylolm} that automatically discards taxa that are
#' absent from the tree tip labels, phenotype labels, or gene presence/absence
#' labels.
#'
#' @param m Named numeric vector of gene presence/absences per taxon.
#' @param p Named numeric vector of phenotype values per taxon.
#' @param phy Phylogeny relating taxa (class \code{"phylo"}).
#' @return Output of \code{phylolm}.
#' @keywords internal
phylolm.subset <- function(p, m, phy) {
    # For testing
    keep <- intersect(names(p), intersect(names(m), phy$tip.label))
    m <- m[keep]
    p <- p[keep]
    phys <- keep.tips(phy, keep)
    phylolm::phylolm(p ~ m, phy = phys)
}


#' Load the \code{phylogenize} package, either using \code{library} or
#' \code{devtools} as appropriate.
#'
#' @param cl A \code{parallel} cluster
#' @param devel Boolean; if TRUE, use \code{devtools}, otherwise use
#'     \code{library}
#' @param pkgdir Optional string giving path to package; used with
#'     \code{devtools} only
#' @keywords internal
cluster.load.pkg <- function(cl, devel, pkgdir="package/phylogenize") {
    if (devel) {
        parallel::clusterCall(cl,
                              "library",
                              "devtools",
                              character.only=TRUE)
        parallel::clusterCall(cl,
                              "load_all",
                              pkgdir)
    } else {
        parallel::clusterCall(cl,
                              "library",
                              "phylogenize",
                              character.only=TRUE)
    }
}

#' Perform phylogenetic (or linear) modeling for a single taxon.
#'
#' @param pheno Named numeric vector giving phenotype values per species.
#' @param tree A tree relating taxa within a taxon.
#' @param mtx Gene presence/absence matrix.
#' @param method A function that returns a length-2 numeric vector of
#'     effect-size and p-value (see, e.g., \code{phylolm.fx.pv} or
#'     \code{lm.fx.pv}).
#' @param restrict.ff Optionally, a character vector giving a subset of
#'     genes to test.
#' @param restrict.taxa Optionally, a character vector giving a subset of
#'     taxa to test.
#' @return Matrix of p-values (row 1) and effect-sizes (row 2) per gene
#'     (columns).
#' @export
matrix.plm <- function(tree,
                       mtx,
                       pheno,
                       method=phylolm.fx.pv,
                       restrict.taxa=NULL,
                       restrict.ff=NULL,
                       ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    cores <- opts('ncl')
    if (is.null(restrict.taxa)) restrict.taxa <- colnames(mtx)
    if (is.null(restrict.ff)) restrict.ff <- rownames(mtx)
    if (opts('separate_process') || (cores > 1)) {
        cl <- parallel::makeCluster(cores)
        cluster.load.pkg(cl, opts('devel'), opts('devel_pkgdir'))
        force(pheno)
        force(tree)
        force(restrict.taxa)
        force(restrict.ff)
        force(mtx)
    } else {
        cl <- NULL
    }
    r <- maybeParApply(mtx[restrict.ff, restrict.taxa, drop=FALSE],
                       1,
                       method,
                       cl,
                       p=pheno,
                       tr=tree,
                       restrict=restrict.taxa,
                       meas_err=opts('meas_err'))
    if (opts('separate_process') || (cores > 1)) {
        parallel::stopCluster(cl)
    }
    r
}

# laplace regularization of P(T|E) and P(T|~E) via MAP estimation
# param is ~ p**k * (1-p)**(n-k) (the n choose k is constant wrt p)
# prior is ( (1/(2b)) * exp(-(|x-m|/b) ) where x, m are on logit scale

#' Obtain a regularized estimate of specificity.
#'
#' @param vec A named numeric vector giving presence/absence across samples.
#' @param env.ids A named factor assigning an environment to each sample
#'   (names).
#' @param which.env A string: in which environment should specificity be
#'   calculated?
#' @param prior Prior probability of \code{which.env}
#' @param b Free parameter giving degree of regularization (see
#'   \code{optimize.b.wrapper}).
#' @param add.pc Boolean giving whether to add a pseudocount.
#' @param min.limit Will not optimize below this value.
#' @param max.limit Will not optimize above this value.
#' @return A list. \code{x}: regularized specificity estimate; \code{p}:
#'   regularized prevalence estimate; \code{x.init}: naive specificity estimate;
#'   \code{p.init}: naive prevalence estimate; \code{pT}: probability of
#'   encountering a particular taxon, marginalized across environments
#' @keywords internal
regularize.pET <- function(vec,
                          env.ids,
                          which.env = 1,
                          prior = 0.05,
                          b = 1,
                          add.pc = FALSE,
                          min.limit = -10,
                          max.limit = 10
                          ) {
    # n and k in binomial are #(E) and #(E,T), respectively (used on P(T|E))
    n <- sum(env.ids == which.env) # #(E)
    k <- sum(vec[env.ids == which.env] > 0) # #(E & T)
    if (add.pc) {
        n <- n + 1
        k <- k + 1
    }
    if (!add.pc) {
        pT.E <- mean(c(vec[which(env.ids == which.env)] > 0))
        pT.nE <- mean(c(vec[which(env.ids != which.env)]) > 0)
    } else {
        pT.E <- mean(c(0, 1, vec[which(env.ids == which.env)] > 0))
        pT.nE <- mean(c(0, 1, vec[which(env.ids != which.env)]) > 0)
    }
    pT <- pT.E * prior + pT.nE * (1 - prior)
    map <- function(p) {
        x <- (p * pT) / prior # P(T|C)
        dbinom(k, n, x) * ((1/2*b)) * exp(-(abs(logit(p)-logit(prior))/b))
    }
    map.logit <- function(logit.p) {
        p <- logistic(logit.p)
        x <- (p * pT) / prior # P(T|C)
        # return log probability also, better numerical stability
        dbinom(k, n, x, log = TRUE) +
            log(1 / (2*b)) -
            (abs(logit(p)-logit(prior))/b)
    }
    max.p <- prior / pT
    max.p <- min(max.p, logistic(max.limit))
    results <- optimize(map.logit, c(min.limit, logit(max.p)), maximum = TRUE)
    initial.x <- pT.E
    initial.p <- (pT.E * prior) / pT
    final.x <- logistic(results$maximum) * pT / prior
    final.p <- logistic(results$maximum)
    return(c(x = final.x,
             p = final.p,
             x.init = initial.x,
             p.init = initial.p,
             pT = pT))
}


#' Get a distribution of environment prevalences from a matrix.
#'
#' This is used to optimize the value of the free parameter $b$ in
#' \code{regularize.pET}.
#'
#' @param mtx A matrix of presence/absence values.
#' @param ids A named factor assigning samples (matrix columns) to environments.
#' @param fallback A two-element numeric vector, giving the beta parameters to
#'     use if fitting fails.
#' @return A best-fit of prevalences to a beta distribution.
#' @keywords internal
fit.beta.list <-  function(mtx, ids, fallback = c(NA, NA)) {
    lapply(unique(ids), function(i) {
        tryCatch(
            MASS::fitdistr(densfun = "beta",
                     start = list(shape1 = 1, shape2 = 1),
                     apply(mtx[, which(ids == i), drop=FALSE], 1, function(x) {
                         mean(c(x, 0, 1) > 0)
                     }))$estimate,
            error = function(e) fallback)
    })
}

#' Simulate a presence/absence matrix.
#'
#' @param effect.size A number giving the shift in logit-prevalence between
#'   environments for simulated true positives.
#' @param baseline.distro A two-element numeric vector of beta distribution
#'   parameters, giving the distribution of simulated prevalences.
#' @param which.env String; gives the environment in which an effect will be
#'   simulated.
#' @param samples Named integer vector giving the number of samples per
#'   environment.
#' @param taxa How many taxa to simulate?
#' @param tpr Fraction of taxa that should be simulated as true positives.
#' @param sign.pos.prob Fraction of effects that should be positive (vs.
#'   negative)?
#' @return A list. \code{mtx}: a simulated matrix; \code{pT}: true prevalences;
#'   \code{pTbs}: true prevalences after adding in effects; \code{fx}: effects;
#'   \code{ids}: mapping of samples to environments; \code{input.params}: input
#'   parameters used to call \code{simulate.binom.mtx}.
#' @keywords internal
simulate.binom.mtx <- function(effect.size = 2,
                              baseline.distro = c(shape1 = 0.66,
                                                  shape2 = 2.62),
                              samples = c(H = 38, D = 13),
                              which.env = "D",
                              taxa = 2000,
                              tpr = 0.25,
                              sign.pos.prob = 0.5) {
    tp.taxa <- round(tpr * taxa)
    neg.taxa <- taxa - (tp.taxa)
    n.classes <- length(samples)
    pT <- sapply(1:taxa, function(.) {
        rbeta(1, baseline.distro[1], baseline.distro[2])
    })
    fx <- c(((2 * rbinom(n = tp.taxa, size = 1, sign.pos.prob)) - 1),
            rep(0, neg.taxa))
    pTbs <- sapply(1:taxa, function(i) {
        pTb <- logistic(logit(pT) + (fx * (effect.size)))
    })
    if (is.null(names(samples))) names(samples) <- 1:length(samples)
    sim.mtx <- t(sapply(1:taxa, function(i) {
        Reduce(c, lapply.across.names(names(samples), function(smp) {
            if (smp == which.env) { p <- pTbs[i] } else { p <- pT[i] }
            rbinom(samples[smp], size = 1, p)
        }))
    }))
    ids <- Reduce(c, lapply(1:length(samples), function(i) rep(i, samples[i])))
    return(list(mtx = sim.mtx,
                pT = pT,
                pTbs = pTbs,
                fx = fx,
                ids = ids,
                input.params = list(effect.size = effect.size,
                                    baseline = baseline.distro,
                                    samples = samples,
                                    taxa = taxa,
                                    tpr = tpr,
                                    sign.pos.prob = sign.pos.prob)))
}

#' Score a simulated regularization by how well it recapitulates the ground truth.
#'
#' @param mtx A simulated matrix of presence/absences.
#' @param ids A factor mapping samples to environments.
#' @param real.fx A numeric vector giving "true" effect sizes.
#' @param which.env String or numeric: in which environment is there an effect?
#' @param prior Prior probability of encountering environment \code{which.env}.
#' @param b Free parameter governing strength of regularization. Typically, this
#'     function is called to evaluate different values of $b$.
#' @param add.pc Boolean: should \code{regularize.pET} add a pseudocount?
#' @param tol Numeric: values within \code{tol} of the prior will be considered
#'     to be shrunk back to the prior completely.
#' @return A vector. \code{fpr}: False positive rate; \code{pwr.hi}: power for
#'     positive effect sizes; \code{pwr.lo}: power for negative effect sizes.
#' @keywords internal
score.regularization <- function(mtx,
                                 ids,
                                 real.fx,
                                 which.env = 2,
                                 prior = 0.002,
                                 b = 0.1,
                                 add.pc = FALSE,
                                 tol = prior * 0.005,
                                 ...) {
    regularized <- apply(mtx, 1, function(x) {
        regularize.pET(x,
                       ids,
                       which.env = which.env,
                       prior = prior,
                       b = b,
                       add.pc = add.pc)
    })
    posteriors <- regularized[2, , drop=TRUE]
    signif <- 1 * (abs(posteriors - prior) > tol)
    predicted.signs <- signif * sign(posteriors - prior)
    fpr <- (signif[real.fx == 0] %>% mean)
    pwr.hi <- mean(predicted.signs[real.fx == 1] == 1)
    pwr.lo <- mean(predicted.signs[real.fx == -1] == -1)
    return(c(fpr = fpr, pwr.hi = pwr.hi, pwr.lo = pwr.lo))
}

#' Based on a real presence/absence matrix, optimize the value of the
#' regularization parameter $b$.
#'
#' @param real.mtx Presence/absence matrix.
#' @param real.ids A factor assigning environment labels to samples.
#' @param which.real.env A numeric or string giving the environment to calculate
#'     specificity within.
#' @param which.shape A numeric or string giving the environment in which the
#'     prevalence distribution should be fit.
#' @param prior Prior probability of encountering environment \code{which.real.env}.
#' @param effect.size Effect size to be simulated for "true" positives.
#' @param bounds Bounds for optimizing $b$.
#' @param add.pc Boolean giving whether to add a pseudocount when calculating
#'     prevalences and specificities.
#' @param tol Numeric: values within \code{tol} of the prior will be considered
#'     to be shrunk back to the prior completely.
#' @param a See ?\code{b.scorer}.
#' @param verbose Boolean: print debugging information?
#' @param optim.fxn Function summarizing the results of
#'     \code{score.regularization} into one metric to be optimized.
#' @param pos.prob Fraction of real effects that should be positive.
#' @return The output of \code{stats::optimize}.
#' @keywords internal
optimize.b.wrapper <- function(real.mtx,
                              real.ids,
                              which.real.env = 2,
                              which.shape = 1,
                              prior = 0.002,
                              effect.size = 2,
                              bounds = c(0.01, 5),
                              add.pc = FALSE,
                              tol = prior * 0.005,
                              a = 0.05,
                              verbose = FALSE,
                              optim.fxn = b.scorer,
                              pos.prop = 0.5,
                              ...
                              ) {
    shape.n <- which((real.ids %>% unique %>% sort) == which.shape)
    # Fall back to the overall beta if not enough information to fit a
    # distribution (or too wacky)
    overall.shape <- MASS::fitdistr(densfun = "beta",
                                    start = list(shape1 = 1, shape2 = 1),
                                    apply(real.mtx, 1, function(x) {
                                        mean(c(x, 0, 1) > 0)
                                    }))$estimate
    shapes <- suppressWarnings(fit.beta.list(real.mtx,
                                             real.ids,
                                             fallback = overall.shape))
    N <- count.each(real.ids)
    which.env <- which(names(N) == which.real.env)
    get.optim <- function(b) {
	sim <- simulate.binom.mtx(effect.size = effect.size,
                                  baseline.distro = shapes[[shape.n]],
                                  samples = N,
                                  taxa = nrow(real.mtx),
                                  tpr = 0.25,
                                  sign.pos.prob = pos.prop,
                                  which.env = which.real.env)
        s <- score.regularization(sim$mtx,
                                  sim$ids,
                                  sim$fx,
                                  which.env,
                                  prior,
                                  b,
                                  add.pc,
                                  tol,
                                  ...)
        if (verbose) print(c(s, optim.fxn(s, a)))
        optim.fxn(s, a)
    }
    optimize(get.optim, bounds, maximum = TRUE)
}

#' Summarize the statistics from \code{score.regularization} into a single metric.
#'
#' @param s Results of \code{score.regularization}
#' @param a Parameter controlling when to switch from optimizing the FPR to
#'     optimizing power. If proportion of false positives is above this value,
#'     return 1-FPR; otherwise, return 1 + the average (geometric mean) power on
#'     positive and negative effect sizes.
#' @return Summary statistic ranging between 0 and 2.
#' @keywords internal
b.scorer <- function(s, a) {
    if (s["fpr"] <= a) {
        (1 + geommean(s[c("pwr.hi", "pwr.lo")]))
    } else {
        ( (1 - s["fpr"]))
    }
}

### calculate prevalence

#' Main function to calculate taxon prevalences with additive smoothing.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{env_column}{String. Name of column in metadata file containing the
#'   environment annotations.}
#'   \item{dset_column}{String. Name of column in metadata file containing the
#'   dataset annotations.}
#'   \item{which_envir}{String. Environment in which to calculate prevalence or
#'   specificity. Must match annotations in metadata.}
#' }
#'
#' @param abd.meta A list giving an abundance matrix and metadata.
#' @return An additively-smoothed estimate of taxon prevalences.
#' @export
prev.addw <- function(abd.meta,
                      ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    envir <- opts('which_envir')
    E <- opts('env_column')
    D <- opts('dset_column')
    S <- opts('sample_column')
    if (!(envir %in% levels(abd.meta$metadata[[E]]))) {
        stop(paste0("environment ", envir, " not found in metadata"))
    }
    env.rows <- (abd.meta$metadata[[E]] == envir)
    dsets <- unique(abd.meta$metadata[env.rows, D, drop=TRUE])
    if (length(dsets) > 1) {
        means.by.study <- lapply(dsets, function(d) {
            s <- intersect(
                colnames(abd.meta$mtx),
                abd.meta$metadata[[S]][(env.rows &
                                        (abd.meta$metadata[[D]] == d))] %>%
                as.character)
            rm <- rowMeans(1 * (abd.meta$mtx[, s, drop=FALSE] > 0))
            list(rm = rm, s = s)
        })
        means.only <- sapply(means.by.study, function(x) x$rm)
        total.s <- Reduce(sum, lapply(means.by.study, function(x) length(x$s)))
        avg.prev <- rowMeans(means.only)
        addw <- (1 + (avg.prev * total.s)) / (2 + total.s)
    } else {
        s <- intersect(colnames(abd.meta$mtx),
                       abd.meta$metadata[[S]][env.rows] %>% as.character)
        total.s <- length(s)
        rs <- rowSums(1 * (abd.meta$mtx[, s, drop=FALSE] > 0))
        addw <- (1 + rs) / (2 + total.s)
    }
    return(logit(addw))
}


### calculate correlation

#' Main function to calculate taxon-to-phenotype correlations, using
#' clr-transformed abundances.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{env_column}{String. Name of column in metadata file containing (in
#'   this case) the correlation variable.}
#'   \item{dset_column}{String. Name of column in metadata file containing the
#'   dataset annotations.}
#'   \item{sample_column}{String. Name of column in metadata file containing the
#'   sample names.}
#' }
#'
#' @param abd.meta A list giving an abundance matrix and metadata.
#' @return An estimate of taxon abundance correlations with a phenotype.
#' @export
correl.clr <- function(abd.meta,
                      ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    R <- opts('env_column')
    S <- opts('sample_column')
    D <- opts('dset_column')
    if (!(R %in% colnames(abd.meta$metadata))) {
      stop(paste0("column ", R, " not found in metadata"))
    }
    R.rows <- !is.na(abd.meta$metadata[[R]])
    dsets <- unique(abd.meta$metadata[R.rows, D, drop=TRUE])
    if (length(dsets) > 1) {
        warning("datasets are ignored when calculating correlation")
    }
    clr_mtx <- clr(abd.meta$mtx)
    clr_cols <- colnames(clr_mtx)
    R_cols <- abd.meta$metadata[R.rows, S] %>% unlist %>% as.character
    R_values <- abd.meta$metadata[[R]][R.rows] %>% as.numeric
    names(R_values) <- R_cols
    keep_cols <- intersect(clr_cols, R_cols)
    clr_mtx <- clr_mtx[, keep_cols]
    trait <- R_values[keep_cols]
    if (sd(trait) == 0) { stop(
        "Trait has a constant value; can't take correlation") }
    if (all(trait %in% c(0, 1))) { warning(
        "Trait looks binary; are you sure you want to take correlation?") }
    corr_values <- apply(clr_mtx, 1, function(x) cor(x, trait))
    # take care of any perfect correlations before taking fisher transform
    corr_values[corr_values == 1] <- (
        (max(corr_values[corr_values < 1]) + 1) / 2
    )
    corr_values[corr_values == -1] <- (
        (min(corr_values[corr_values > 1]) + -1) / 2
    )
    return(atanh(corr_values))
}

#' Function to take the clr transform of a matrix.
#' @param pc Pseudocount to add to all values (default 0.5).
clr <- function(mtx, pc = 0.5) {
  apply(mtx + pc, 2, function(x) log(x) - mean(log(x)))
}

#' Main function to calculate environmental specificity scores.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{env_column}{String. Name of column in metadata file containing the
#'   environment annotations.}
#'   \item{dset_column}{String. Name of column in metadata file containing the
#'   dataset annotations.}
#'   \item{which_envir}{String. Environment in which to calculate prevalence or
#'   specificity. Must match annotations in metadata.}
#'   \item{prior_type}{String. What type of prior to use ("uninformative" or
#'   "file").}
#' }
#'
#' @param abd.meta A list giving an abundance matrix and metadata.
#' @param pdata Named numeric vector giving priors per environment.
#' @param b.optim If not NULL, use this value for the regularization parameter
#'     $b$, otherwise optimize it.
#' @return A list with the following components:
#' \describe{
#'   \item{b.optim}{Shrinkage parameter b obtained through optimization}
#'   \item{ess}{Logit-transformed, shrunken estimates of specificity.}
#'   \item{regularized}{Non-transformed regularized values.}
#'   \item{priors}{Values of priors.}
#'   \item{phenoP}{Prior for environment of interest.}
#' } 
#' @export
calc.ess <- function(abd.meta,
                     pdata = NULL,
                     b.optim = NULL,
                     ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    E <- opts('env_column')
    D <- opts('dset_column')
    S <- opts('sample_column')
    ptype <- opts('prior_type')
    envir <- opts('which_envir')
    if (!(envir %in% levels(abd.meta$metadata[[E]]))) {
        stop(paste0("environment ", envir, " not found in metadata"))
    }
    env.rows <- (abd.meta$metadata[[E]] == envir)
    dsets <- unique(abd.meta$metadata[env.rows, D, drop=TRUE])
    if (length(dsets) > 1) {
        warning("datasets are ignored when calculating specificity")
    }
    meta.present <- abd.meta$meta[(abd.meta$meta[[S]] %>%
                                   as.character %in%
                                   colnames(abd.meta$mtx)), ]
    envirs <- unique(meta.present[[E]])
    if (ptype == "uninformative") {
        priors <- data.frame(
            env = envirs,
            prior = 1/length(envirs)
        )
    } else if (ptype == "alpha") {
        stop("not implemented yet, sorry")
    } else if (ptype == "file") {
        priors <- pdata[pdata$env %in% envirs, , drop=FALSE]
    } else {
        stop(paste0("don't know how to compute priors of type ", ptype))
    }
    ids <- sapply(colnames(abd.meta$mtx),
                  function (sn) abd.meta$meta[(abd.meta$meta[[S]] == sn), E])
    if (length(unique(ids)) < 2) { 
        stop("error: only one environment found")
    }
    names(ids) <- colnames(abd.meta$mtx)
    ids <- simplify2array(ids)
    this.prior <- priors$prior[which(priors$env == envir)]
    tolerance <- this.prior * 0.01
    if (is.null(b.optim)) {
        b.optim <- optimize.b.wrapper(effect.size = 2,
                                      real.mtx = abd.meta$mtx,
                                      real.ids = ids %>% simplify2array,
                                      which.real.env = envir,
                                      which.shape = envir,
                                      prior = this.prior,
                                      tol = tolerance,
                                      a = 0.05,
                                      optim.fxn = b.scorer,
                                      verbose = FALSE,
                                      add.pc = TRUE)$maximum
    } else {
        warning("skipping optimization of Laplace parameter (not recommended)")
    }
    regularized <- apply(abd.meta$mtx, 1, function(x) {
        regularize.pET(x,
                       ids,
                       which.env = envir,
                       b = b.optim,
                       prior = this.prior,
                       add.pc = TRUE)
    })
    logist.pheno <- regularized[2, , drop=TRUE]
                                        # "hard-shrink" anything shrunk almost
                                        # to the prior to prevent these tiny
                                        # differences from affecting the result
                                        # in the absence of a strong change
    logist.pheno[which(logist.pheno %btwn%
                       c(this.prior - tolerance,
                         this.prior + tolerance))] <- this.prior
    phenoP <- logit(c(priors[priors$env == envir, "prior", drop=TRUE]))
    return(list(
        b.optim=b.optim,
        ess=logit(logist.pheno),
        regularized=regularized,
        priors=priors,
        phenoP=phenoP))
}

### Wrap ashr with default parameters

#' Wrapper function to perform adaptive shrinkage.
#'
#' @param m Estimates of parameter of interest.
#' @param s Standard errors of parameters of interest.
#' @param nw Null weight (default=10).
#' @param ashr_df Degrees of freedom (default=5)
#' @return A vector giving shrunken estimates of parameter.
ash_wrapper <- function(m, s, nw=10, ashr_df=5) {
    ashr::ash(m, s,
            mixcompdist="halfuniform",
            prior="nullbiased",
            nullweight=nw,
            df=ashr_df)
}
### calculate differential abundance

#' Main function to calculate differential abundances using either ANCOMBC2
#' or MaAsLin2, then smooth the results with adaptive shrinkage. Note that the
#' packages for ANCOMBC2 or MaAsLin2 must already be installed.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{env_column}{String. Name of column in metadata file containing the
#'   environment annotations.}
#'   \item{dset_column}{String. Name of column in metadata file containing the
#'   dataset annotations. Will be incorporated as a "nuisance" variable.}
#'   \item{diff_abund_method}{String. Either "ANCOMBC2" or "Maaslin2" (case 
#'   insensitive).}
#' }
#'
#' @param abd.meta A list giving an abundance matrix and metadata.
#' @return A vector giving shrunken estimates of differential abundance.
#' @export
ashr.diff.abund <- function(abd.meta,
                            ...) {
  opts <- clone_and_merge(PZ_OPTIONS, ...)
  categorical <- opts('categorical')
  envir <- opts('which_envir')
  E <- opts('env_column')
  D <- opts('dset_column')
  S <- opts('sample_column')
  M <- tolower(opts('diff_abund_method'))
  if (!(M %in% c('ancombc2', 'maaslin2'))) {
    pz.error(paste0("method ", M, " not recognized (see help)"))
  }
  if (categorical) {
    env_levels <- levels(abd.meta$metadata[[E]])
    if (!(envir %in% env_levels)) {
      pz.error(paste0("environment ", envir, " not found in metadata"))
    }
  } else {
    abd.meta$metadata[[E]] <- as.numeric(abd.meta$metadata[[E]])
    if (all(is.na(abd.meta$metadata[[E]]))) {
      pz.error(paste0(
          "Environment failed conversion to numeric; is this supposed to be ",
          "a categorical variable?"))
    }
  }
  # Both tools expect rownames to be sample IDs
  named_metadata <- data.frame(
    abd.meta$metadata[, setdiff(colnames(abd.meta$metadata), S)]
  )
  rownames(named_metadata) <- abd.meta$metadata[[S]]
  
  if (M=="ancombc2") { 
      
    named_metadata <- named_metadata %>%
	    dplyr::mutate(dplyr::across(c(E, D), as.factor))
    if (length(levels(named_metadata[[D]])) < 2) {
      ancom_formula = E
      } else {
      ancom_formula = paste0(E, "+", D)
    }
    if (length(levels(named_metadata[[E]])) < 2) {
      pz.error(
          "For abundance there must be at least two different environments")
    }
    named_metadata <- named_metadata %>%
	    filter(!is.na(E) | !is.na(D))
    
    ancom_tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
      assays=S4Vectors::SimpleList(counts=abd.meta$mtx),
      colData=S4Vectors::DataFrame(named_metadata))
    
    ancom_results <- ANCOMBC::ancombc2(ancom_tse,
                              assay_name="counts",
                              fix_formula=ancom_formula,
                              n_cl=opts('ncl'))
    ancom_results_tbl <- ancom_results$res %>% tibble::tibble()
    envir_stem <- paste0("env", envir)
    lfc_col <- paste0("lfc_", envir_stem)
    se_col <- paste0("se_", envir_stem)
    for (ctest in c(lfc_col, se_col)) {
        if (!(ctest %in% colnames(ancom_results_tbl))) {
            pz.error(paste0("The column ", ctest, " was not found in the ",
                            "ANCOMBC2 results. This should never happen and ",
                            "should be reported as a bug."))
        } 
    }
    ancom_ash <- ash_wrapper(dplyr::pull(ancom_results_tbl,
					 !!(lfc_col),
					 name=taxon),
			     dplyr::pull(ancom_results_tbl,
					 !!(se_col),
					 name=taxon))
    sample_ashr <- tibble::as_tibble(ancom_ash$result, rownames = "species")
    
  } else if (M == "maaslin2") {
      
      ref_env_level <- levels(named_metadata[[E]])[1]
      all_dset_level <- levels(named_metadata[[D]])
      ref_dset_level <- all_dset_level[1]
      if (length(all_dset_level) > 1) {
          m_fixed_effects <- c(E, D)
          m_reference <- paste0(ref_env_level, ",", ref_dset_level)
      } else {
          m_fixed_effects <- c(E)
          m_reference <- ref_env_level
      }
      maaslin_res <- Maaslin2::Maaslin2(
          input_data = as.data.frame(as.matrix(abd.meta$mtx)),
          input_metadata = named_metadata,
          # this is a vector of strings
          fixed_effects = m_fixed_effects,
          # this has to be a single string
          reference = m_reference,
          output = "temp/",
          plot_heatmap = FALSE,
          plot_scatter = FALSE
      )
      maaslin_results_tbl <- maaslin_res$results %>% tibble::tibble() %>%
          dplyr::filter(metadata == E, value == envir)
      if (any(is.na(maaslin_results_tbl$coef)) ||
          any(is.na(maaslin_results_tbl$stderr))) {
          pz.error("Error: Missing values found in coef or stderr columns.")
      }
      maaslin_ash <- ash_wrapper(dplyr::pull(maaslin_results_tbl,
                       coef,
                       name=feature),
                               dplyr::pull(maaslin_results_tbl,
                                           stderr,
                                           name=feature))
      sample_ashr <- tibble::as_tibble(maaslin_ash$result, rownames = "species")
  } else {
      pz.error(
          paste0(
              "This shouldn't be possible, but somehow an incorrect ",
              "differential abundance method was selected"
          )
      )
  }
  sample_pheno <- sample_ashr %>%
    dplyr::select(species, PosteriorMean) %>%
    tibble::deframe()
  
  # Fix automatic renaming of taxa that seem "numeric"
  spn <- names(sample_pheno)
  orign <- rownames(abd.meta$mtx)
  numeric_names <- tibble::tibble(old_name =
                                      as.character(orign[
                                          which(!is.na(as.numeric(orign)))
                                      ]),
                                  new_name = paste0("X", old_name))
  # if the names didn't change, just keep them:
  all_names <- dplyr::bind_rows(numeric_names,
                                tibble::tibble(old_name=orign,
                                               new_name=orign))
  sample_pheno_tbl <- dplyr::left_join(
      tibble::enframe(sample_pheno, name="new_name", value="pheno"),
      all_names)
  
  sample_pheno <- sample_pheno_tbl %>%
      dplyr::select(old_name, pheno) %>%
      tibble::deframe()

  return(sample_pheno)
}


#' Throw an error and optionally log it in errmsg.txt.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{error_to_file}{Boolean. Should pz.error, pz.warning, and pz.message
#'   output to an error message file?}
#'}
#'
#' @param errtext String: error message text.
#' @export
pz.error <- function(errtext, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    if (opts('error_to_file')) {
        tryCatch({
          cat(paste0(errtext, "\n"),
              file = file.path(opts('out_dir'), "errmsg.txt"),
              append = TRUE)
          }, error=function(e) NULL)
        }
    stop(errtext)
}

#' Report a message and optionally log it in errmsg.txt.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{error_to_file}{Boolean. Should pz.error, pz.warning, and pz.message
#'   output to an error message file?}
#'}
#'
#' @param errtext String: message text.
#' @export
pz.message <- function(msgtext, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    if (opts('error_to_file')) {
        tryCatch({
            cat(paste0(msgtext, '\n'),
                file = file.path(opts('out_dir'), "errmsg.txt"),
                append = TRUE)
        }, error=function(e) NULL)
    }
    message(msgtext)
}

#' Report a warning and optionally log it in errmsg.txt.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{error_to_file}{Boolean. Should pz.error, pz.warning, and pz.message
#'   output to an error message file?}
#'}
#'
#' @param errtext String: warning text.
#' @export
pz.warning <- function(msgtext, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    if (opts('error_to_file')) {
        tryCatch({
            cat(paste0(msgtext, '\n'),
                file = file.path(opts('out_dir'), "errmsg.txt"),
                append = TRUE)
        }, error=function(e) NULL)
    }
    warning(msgtext)
}

#' Convert results into a long (vs. wide) format.
#'
#' @param results Output of result.wrapper.plm.
#' @return A single data frame with entries from \code{results}.
#' @export
make.results.matrix <- function(results) {
    Reduce(dplyr::bind_rows, lapply(names(results), function(rn) {
        tibble::tibble(taxon = rn,
                       gene = results[[rn]] %>% colnames,
                       effect.size = results[[rn]][1,],
                       p.value = results[[rn]][2,],
                       std.err = results[[rn]][3,],
                       df = results[[rn]][4,])
    }))
}

#' Filter out genes that are almost always present or absent.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{minimum}{Integer. A particular gene must be observed, and also
#'   absent, at least this many times to be reported as a significant positive
#'   association with the phenotype.}
#' }
#' @param pz.db A database for use with *phylogenize* analyses.
#' @param phy.with.sigs A vector of strings giving which taxa had significant
#'     results.
#' @return A single data frame with entries from \code{results}.
#' @export
threshold.pos.sigs <- function(pz.db, phy.with.sigs, pos.sig, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    Min <- opts('minimum')
    mapply(pz.db$trees[phy.with.sigs],
           pos.sig[phy.with.sigs],
           pz.db$gene.presence[phy.with.sigs],
           FUN = function(tr, x, y) {
               i <- na.omit(intersect(tr$tip.label, colnames(y)))
               if (length(i) == 0) { return(character(0)) }
               if (length(x) == 0) { return(character(0)) }
               y2 <- y[x, i, drop=FALSE]
               r1 <- rowSums(y2)
               r2 <- rowSums(1 - y2)
               names(which((r2 >= Min) & (r1 >= Min)))
           },
           SIMPLIFY=FALSE)
}

#' Filter out genes that are almost always present or absent prior to running
#' regressions.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{minimum}{Integer. A particular gene must be observed, and also
#'   absent, at least this many times to be reported as a significant positive
#'   association with the phenotype.}
#' }
#' @param gene.presence A list of matrices of gene presence/absence, as included
#'   in a `pz.db` object.
#' @param trees A list of trees. Names should match names of `gene.presence`.
#'   Only taxa in these trees will be used for the analysis.
#' @return A revised `gene.presence` list. Note that some taxa may be dropped
#'   from the list (if they had zero genes left after filtering).
#' @export
above_minimum_genes <- function(gene.presence, trees, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    Min <- opts('minimum')
    taxa <- names(trees)
    to_remove <- rep(FALSE, length(taxa)) %>% setNames(taxa)
    for (tx in taxa) {
        tips <- trees[[tx]]$tip.label
        colns <- colnames(gene.presence[[tx]])
        i <- na.omit(intersect(tips, colns))
        if (length(i) > 0) {
            mtx <- gene.presence[[tx]][, i, drop=FALSE]
            g <- names(which((rowSums(mtx) >= Min) & (rowSums(!mtx) >= Min)))
            gene.presence[[tx]] <- mtx[g, , drop=FALSE]
	}
	if ((length(i) == 0) || (length(g) == 0)) {
            to_remove[tx] <- TRUE
	}
    }
    gene.presence[names(which(!to_remove))]
}

#' Add gene descriptions to significant results; return in a tibble.
#'
#' @param phy.with.sigs Character vector giving the taxa with significant hits.
#' @param pos.sig List of character vectors, one per species, of significant hits.
#' @param gene.to.fxn Data frame used to annotate genes to functions.
#' @return A single data frame of all significant results plus descriptions.
#' @export
add.sig.descs <- function(phy.with.sigs, pos.sig, gene.to.fxn) {
    pos.sig.tbl <- tibble::enframe(pos.sig, name="taxon", value="gene") %>%
        tidyr::unnest()
    
    column_names <- colnames(gene.to.fxn)
    na_columns <- which(is.na(column_names))
    if (length(na_columns) > 0) {
        for (i in seq_along(na_columns)) {
            column_names[na_columns[i]] <- paste0("NA_col_", i)
        }
    }
    colnames(gene.to.fxn) <- column_names
   
    # Ensure they are the same type and are not doubles
    gene.to.fxn <- gene.to.fxn %>%
	    dplyr::mutate(gene = as.character(gene)) 
    pos.sig.tbl <- pos.sig.tbl %>%
	    dplyr::mutate(gene = as.character(gene))
    pos.sig.descs <- dplyr::left_join(pos.sig.tbl,
                               gene.to.fxn,
                               by="gene") %>%
        dplyr::rename(description=`function`)
}
