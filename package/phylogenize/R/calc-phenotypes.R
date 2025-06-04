#--- Functions to calculate phenotypes ---#

#' Read in data and metadata and calculate a phenotype (differential abundance,
#' specificity, or prevalence).
#'
#' @param save_data Return processed data/metadata. Default is FALSE to save
#'   memory. Must set to TRUE if running POMS.
#' @param ... Parameters to override defaults.
#' @export
data_to_phenotypes <- function(save_data=FALSE, ...) {
    pz.options <- clone_and_merge(PZ_OPTIONS, ...)
    # Read in user-supplied data and metadata
    abd.meta <- read.abd.metadata(...)
    # Read in trees, gene presence/absence, taxonomy
    pz.db <- import.pz.db(...)
    # Figure out how many trees to retain
    pz.db <- adjust.db(pz.db, abd.meta, ...)
    if (pz.options('assume_below_LOD')) {
        abd.meta <- add.below.LOD(pz.db, abd.meta, ...)
        sanity.check.abundance(abd.meta$mtx, ...)
    }
    phenotype_results <- calculate_phenotypes(abd.meta, pz.db, ...)
    if (pz.options('which_phenotype') %in% c("specificity", "abundance")) {
        # only retain observed taxa
        pz.db$trees <- retain.observed.taxa(pz.db$trees,
                                            phenotype_results$phenotype,
                                            phenotype_results$phenoP,
                                            phenotype_results$mapped.observed)
        pz.db$species <- lapply(pz.db$trees, function(x) x$tip.label)
        pz.db$ntaxa <- length(pz.db$trees)
    }
    if (save_data) return(list(
        abd.meta=abd.meta,
        pz.db=pz.db,
        phenotype_results=phenotype_results
    ))
    # otherwise, don't need to save the original data beyond this point
    return(list(pz.db=pz.db, phenotype_results=phenotype_results))
}

#' Determine which phenotype to calculate and then calculate it.
#'
#' @param pz.db A database (typically obtained with \code{import.pz.db}).
#' @param abd.meta A list consisting of a species abundance matrix and the
#'     metadata (from read.abd.metadata or data_to_phenotypes).
#' @param ... Parameters to override defaults.
#' @export
calculate_phenotypes <- function(abd.meta, pz.db, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    mapped.observed <- names(which(Matrix::rowSums(abd.meta$mtx) > 0))
    if (tolower(opts('core_method')) == "poms") {
        pz.warning(paste0("POMS will ignore this phenotype, ",
                          "as it computes its own on balances"))
    }
    if (opts('which_phenotype') == "prevalence") {
        phenotype <- prev.addw(abd.meta, ...)
        phenoP <- NULL
    } else if (opts('which_phenotype') == "specificity") {
        if (opts('prior_type') == "file") {
            prior.data <- read.table(file.path(opts('input_dir'),
                                               opts('prior_file')))
        } else {
            prior.data <- NULL
        }
        ess <- calc.ess(abd.meta,
                        prior.data,
                        ...)
        phenotype <- ess$ess
        phenoP <- ess$phenoP
    } else if (pz.options("which_phenotype") == "provided") {
        p_tbl <- read_tsv(pz.options("phenotype_file"))
        if (ncol(p_tbl) == 2) { # assume we only have species IDs and values
            phenotype <- deframe(p_tbl)
        } else { # perform shrinkage on the provided values w/ their stderrs
            p_est <- as.numeric(p_tbl[["estimate"]])
            p_se <- as.numeric(p_tbl[["stderr"]])
            names(p_est) <- p_tbl[[1]]
            names(p_se) <- p_tbl[[1]]
            ashr_res <- ash_wrapper(p_est, p_se)
            phenotype <- ashr_res$result %>%
                as_tibble(rownames="species") %>%
                dplyr::select(species, PosteriorMean) %>%
                deframe
        }
        phenoP <- 0
    } else if (opts("which_phenotype") == "abundance") {
        phenotype <- ashr.diff.abund(abd.meta, ...)
        phenoP <- 0
    } else {
        pz.error(paste0("don't know how to calculate the phenotype ",
                        opts('which_phenotype')))
    }
    phenotype <- clean.pheno(phenotype, pz.db)
    if (pz.options("which_phenotype") != "prevalence") {
        # Except for prevalence, retain observed taxa
        pz.db$trees <- retain.observed.taxa(pz.db$trees,
                                            phenotype,
                                            phenoP,
                                            mapped.observed)
        pz.db$trees <- pz.db$trees[
            vapply(pz.db$trees, \(.) length(.$tip.label), 1L) >=
                pz.options("treemin")
        ]
        if (length(pz.db$trees) == 0) { pz.error("all trees dropped") }
        pz.db$species <- lapply(pz.db$trees, function(x) x$tip.label)
        pz.db$ntaxa <- length(pz.db$trees)
    }
    list(phenotype=phenotype, phenoP=phenoP, mapped.observed=mapped.observed)
}


#' Return a boolean telling whether a phenotype has nonzero variance in
#' different taxa.
#'
#' @param phenotype A quantitative phenotype.
#' @param taxa From `pz.db$taxa`.
#' @return A boolean vector with length equal to `length(taxa)`.
#' @export
pheno_nonzero_var <- function(phenotype,taxa) {
    vapply(taxa,
           function(tx) {
               p <- phenotype[intersect(names(phenotype),tx)]
               return((var(p) > 0))
           },
           TRUE)
}


#' Remove taxa from a phenotype that aren't in our trees and gene matrices
#' (usually only necessary for testing).
#'
#' @param phenotype A quantitative phenotype (named numeric vector).
#' @param pz.db A database.
#' @return A phenotype with only the measurements represented in the database.
#' @export
clean.pheno <- function(phenotype, pz.db) {
    tips <- Reduce(union, lapply(pz.db$trees, function(x) x$tip.label))
    cols <- Reduce(union, lapply(pz.db$gene.presence, colnames))
    valid.names <- intersect(tips, cols)
    phenotype[intersect(unlist(names(phenotype)), unlist(valid.names))]
}


#' Modify trees to retain only observed taxa (for use with specificity only).
#'
#' @param trees A list of tree objects.
#' @param phenotype A named vector giving the phenotype for each taxon ID.
#' @param phenoP The prior probability of the environment of interest.
#' @param mapped.observed A character vector giving which tips to retain.
#' @return An updated list of tree objects.
#' @export
retain.observed.taxa <- function(trees, phenotype, phenoP, mapped.observed) {
    trees <- lapply(trees, function(tr) {
        keep.tips(tr, intersect(tr$tip.label, mapped.observed))
    })
    n.not.prior <- sapply(trees, function(tr) {
        sum(phenotype[intersect(names(phenotype), tr$tip.label)] != phenoP)
    })
    if (all(n.not.prior == 0)) {
        pz.error(
            paste0("The phenotype for all taxa was shrunk to the prior. ",
                   "This probably means that there are no major taxonomic ",
                   "differences between groups. Without identifying some ",
                   "such taxonomic changes, phylogenize cannot continue."))
    }
    trees <- trees[which(n.not.prior > 0)]
    trees
}


#' Add in taxa that were not observed, assuming this means they were
#' zero-prevalence.
#'
#' @param pz.db A database (typically obtained with \code{import.pz.db}).
#' @param abd.meta A list consisting of a taxon abundance matrix and the
#'     metadata.
#' @return An updated version of \code{abd.meta}.
#' @export
add.below.LOD <- function(pz.db, abd.meta, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    species.observed <- rownames(abd.meta$mtx)
    all.possible.taxa <- Reduce(union,
                                lapply(pz.db$gene.presence, colnames))
    not.observed.taxa <- setdiff(all.possible.taxa, species.observed)
    if (length(not.observed.taxa) > 0) {
        taxa.zero <- matrix(rep(0, length(not.observed.taxa) *
                                    ncol(abd.meta$mtx)),
                            nr = length(not.observed.taxa),
                            byrow = TRUE,
                            dimnames = list(not.observed.taxa,
                                            colnames(abd.meta$mtx)))
        abd.meta$mtx <- rbind(abd.meta$mtx, taxa.zero)
    }
    # Make sure still binary
    if(opts('which_phenotype') != 'abundance'){
        abd.meta$mtx <- Matrix::Matrix(abd.meta$mtx > 0)
    }
    return(abd.meta)
}
