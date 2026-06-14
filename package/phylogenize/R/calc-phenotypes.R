#--- Functions to calculate phenotypes ---#

#' Read in data and metadata and calculate a phenotype (differential abundance,
#' specificity, or prevalence).
#'
#' @param save_data Return processed data/metadata. Default is FALSE to save
#'   memory. Must set to TRUE if running POMS.
#' @param ... Parameters to override defaults.
#' @export
data_to_phenotypes <- function(save_data=FALSE, ..., .opts=NULL) {
    opts <- pz.resolve.options(..., .opts=.opts)
    pz.message("  A) Reading data, metadata, and databases...")
    # Read in user-supplied data and metadata
    abd.meta <- read.abd.metadata(..., .opts=opts)
    pz.message(paste0(
        "  .....Input abundance matrix has ",
        nrow(abd.meta$mtx),
        " row(s) and ",
        ncol(abd.meta$mtx),
        " sample(s)"
    ))
    pz.message(paste0(
        "  .....Input metadata has ",
        nrow(abd.meta$metadata),
        " row(s)"
    ), level=2)
    pz.message("  B) Read in user-supplied data and metadata")
    # Read in trees, gene presence/absence, taxonomy
    pz.db <- import.pz.db(..., .opts=opts)
    pz.message(paste0(
        "  .....Database loaded with ",
        length(pz.db$trees),
        " tree(s) and ",
        length(pz.db$gene.presence),
        " gene presence matrix/matrices"
    ))
    pz.message("  C) Read in trees, gene presence/absence, taxonomy")
    # Figure out how many trees to retain
    pz.db <- adjust.db(pz.db, abd.meta, ..., .opts=opts)
    pz.message(paste0(
        "  .....Database adjusted to ",
        pz.db$ntaxa,
        " testable taxon/taxa"
    ))
    if (opts('assume_below_LOD')) {
        pz.message("  .....Adding unobserved taxa as below limit of detection")
        abd.meta <- add.below.LOD(pz.db, abd.meta, ..., .opts=opts)
        sanity.check.abundance(abd.meta$mtx, ...)
    }
    pz.message("  D) Calculating phenotypes...")
    phenotype_results <- calculate_phenotypes(abd.meta, pz.db, ..., .opts=opts)
    if (opts('quantile_normalize')) {
	    pz.message("  .....Quantile-normalizing phenotype")
	    phenotype_results <- quantile_normalize(phenotype_results)
    } 
    if (opts('which_phenotype') %in% c("specificity", "abundance")) {
        # only retain observed taxa
        pz.message("  .....Retaining observed taxa in trees")
        pz.db$trees <- retain.observed.taxa(pz.db$trees,
                                            phenotype_results$phenotype,
                                            phenotype_results$phenoP,
                                            phenotype_results$mapped.observed)
        pz.db$species <- lapply(pz.db$trees, function(x) x$tip.label)
        pz.db$ntaxa <- length(pz.db$trees)
        pz.message(paste0(
            "  .....Retained ",
            pz.db$ntaxa,
            " observed testable taxon/taxa"
        ))
    }
    if (tolower(opts('core_method')) == "poms") {
        pz.message("Saving abundance data to run POMS...", 2)
        save_data <- TRUE
    }
    if (save_data) return(list(
        abd.meta=abd.meta,
        pz.db=pz.db,
        phenotype_results=phenotype_results
    ))
    # otherwise, don't need to save the original data beyond this point
    return(list(pz.db=pz.db, phenotype_results=phenotype_results))
}

#' Quantile-normalize phenotype. Wraps the function `quant_norm()`.
#'
#' @param phenotype_results A list with named components "phenotype" and optionally "phenoP".
#' @export
quantile_normalize <- function(phenotype_results) {
  has_phenoP <- !is.null(phenotype_results$phenoP)
  if (has_phenoP) {
	  ph <- c(phenotype_results$phenoP, phenotype_results$phenotype)
  } else {
	  ph <- phenotype_results$phenotype
  }
  normed <- quant_norm(ph)
  if (has_phenoP) {
	  phenotype_results$phenotype <- normed[-1]
	  phenotype_results$phenoP <- normed[1]
  } else {
	  phenotype_results$phenotype <- normed
  }
  phenotype_results
}

#' Actually perform quantile normalization of a vector.
#'
#' @param x Vector to quantile-normalize to the normal distribution.
#' @export
quant_norm <- function(x) {
	l <- length(x)
	n <- qnorm((1:l)/(l+1))[rank(x)]
	names(n) <- names(x)
	n
}

#' Determine which phenotype to calculate and then calculate it.
#'
#' @param pz.db A database (typically obtained with \code{import.pz.db}).
#' @param abd.meta A list consisting of a species abundance matrix and the
#'     metadata (from read.abd.metadata or data_to_phenotypes).
#' @param ... Parameters to override defaults.
#' @export
calculate_phenotypes <- function(abd.meta, pz.db, ..., .opts=NULL) {
    opts <- pz.resolve.options(..., .opts=.opts)
    pz.message("  .....Collecting mapped observations greater than 0")
    mapped.observed <- names(which(Matrix::rowSums(abd.meta$mtx) > 0))
    pz.message(paste0(
        "  ..........Mapped observations retained: ",
        length(mapped.observed)
    ))
    pheno_sd <- NULL
    if (tolower(opts('core_method')) == "poms") {
        pz.warning(paste0("Generating an approximate phenotype just for ",
                          "plotting POMS output (logit-AUC)..."))
        phenotype <- logit_auc_pheno(abd.meta, ..., .opts=opts)
        phenoP <- 0
    } else {
        if (opts('which_phenotype') == "prevalence") {
            pz.message("  .....Running prevalence")
            phenotype <- prev.addw(abd.meta, ..., .opts=opts)
            phenoP <- NULL
        } else if (opts('which_phenotype') == "specificity") {
	    pz.message("  .....Running specificity")
            if (opts('prior_type') == "file") {
                prior.data <- read.table(file.path(opts('input_dir'),
                                                   opts('prior_file')))
            } else {
                prior.data <- NULL
            }
            ess <- calc.ess(abd.meta,
                            prior.data,
                            ...,
                            .opts=opts)
            phenotype <- ess$ess
            phenoP <- ess$phenoP
        } else if (opts("which_phenotype") == "provided") {
            pz.message("  .....Reading provided phenotype")
            phenotype_file <- opts("phenotype_file")
            if (!(file.exists(phenotype_file))) {
                pz.error(paste0("Phenotype file not found: ", phenotype_file),
                         .opts=opts)
            }
            p_tbl <- readr::read_tsv(phenotype_file, show_col_types = FALSE)
            pz.message(paste0(
                "  ..........Provided phenotype table has ",
                nrow(p_tbl),
                " row(s) and ",
                ncol(p_tbl),
                " column(s)"
            ))
            if (nrow(p_tbl) == 0 || ncol(p_tbl) < 2) {
                pz.error(paste0(
                    "Provided phenotype table must have at least one row and ",
                    "two columns"
                ), .opts=opts)
            }
            phenotype_ids <- trimws(as.character(p_tbl[[1]]))
            if (any(is.na(phenotype_ids)) || any(phenotype_ids == "")) {
                pz.error("Provided phenotype table contains blank or missing taxon IDs",
                         .opts=opts)
            }
            if (any(duplicated(phenotype_ids))) {
                duplicate_ids <- unique(phenotype_ids[duplicated(phenotype_ids)])
                pz.error(paste0(
                    "Provided phenotype table contains duplicate taxon IDs: ",
                    paste(duplicate_ids, collapse=", ")
                ), .opts=opts)
            }
            if (ncol(p_tbl) == 2) { # assume we only have species IDs and values
                phenotype_values <- suppressWarnings(as.numeric(p_tbl[[2]]))
                if (any(is.na(phenotype_values))) {
                    pz.error(
                        "Provided phenotype table contains nonnumeric phenotype values",
                        .opts=opts)
                }
                if (any(!is.finite(phenotype_values))) {
                    pz.error(
                        "Provided phenotype table contains non-finite phenotype values",
                        .opts=opts)
                }
                names(phenotype_values) <- phenotype_ids
                phenotype <- phenotype_values
            } else { # perform shrinkage on the provided values w/ their stderrs
                required_pheno_columns <- c("estimate", "stderr")
                missing_pheno_columns <- setdiff(required_pheno_columns,
                                                 names(p_tbl))
                if (length(missing_pheno_columns) > 0) {
                    pz.error(paste0(
                        "Provided phenotype table is missing required column(s): ",
                        paste(missing_pheno_columns, collapse=", ")
                    ), .opts=opts)
                }
                pz.message("  ..........Running shrinkage on provided phenotype")
                p_est <- suppressWarnings(as.numeric(p_tbl[["estimate"]]))
                p_se <- suppressWarnings(as.numeric(p_tbl[["stderr"]]))
                if (any(is.na(p_est))) {
                    pz.error(
                        "Provided phenotype table contains nonnumeric estimates",
                        .opts=opts)
                }
                if (any(is.na(p_se))) {
                    pz.error(
                        "Provided phenotype table contains nonnumeric standard errors",
                        .opts=opts)
                }
                if (any(!is.finite(p_est))) {
                    pz.error(
                        "Provided phenotype table contains non-finite estimates",
                        .opts=opts)
                }
                if (any(!is.finite(p_se)) || any(p_se <= 0)) {
                    pz.error(paste0(
                        "Provided phenotype table standard errors must be ",
                        "finite and greater than zero"
                    ), .opts=opts)
                }
                names(p_est) <- phenotype_ids
                names(p_se) <- phenotype_ids
                ashr_res <- ash_wrapper(p_est, p_se)
                phenotype <- ashr_res$result %>%
                    tibble::as_tibble(rownames="species") %>%
                    dplyr::select(species, PosteriorMean) %>%
                    tibble::deframe()
                pheno_sd <- ashr_res$result %>%
                    tibble::as_tibble(rownames = "species") %>%
                    dplyr::select(species, PosteriorSD) %>%
                    tibble::deframe()
            }
            phenoP <- 0
        } else if (opts("which_phenotype") == "abundance") {
            pz.message("  .....Running abundance")
            pheno_list <- ashr.diff.abund(abd.meta, ..., .opts=opts)
            phenotype <- pheno_list$pheno
            pheno_sd <- pheno_list$sd
            phenoP <- 0
        } else {
            pz.error(paste0("don't know how to calculate the phenotype ",
                            opts('which_phenotype')))
        }
    }
    pz.message("  .....cleaning phenotype")
    phenotype <- clean.pheno(phenotype, pz.db)
    if (length(phenotype) == 0) {
        pz.error("No phenotype values matched database taxa", .opts=opts)
    }
    pz.message(paste0(
        "  ..........Cleaned phenotype retains ",
        length(phenotype),
        " taxon/taxa"
    ))
    if (!is.null(pheno_sd)) pheno_sd <- pheno_sd[names(phenotype)]
    if (opts("which_phenotype") != "prevalence") {
        # Except for prevalence, retain observed taxa
        pz.message("  .....Filtering trees to observed taxa")
        pz.db$trees <- retain.observed.taxa(pz.db$trees,
                                            phenotype,
                                            phenoP,
                                            mapped.observed)
        pz.db$trees <- pz.db$trees[
            vapply(pz.db$trees, \(.) length(.$tip.label), 1L) >=
                opts("treemin")
        ]
        if (length(pz.db$trees) == 0) { pz.error("all trees dropped") }
        pz.db$species <- lapply(pz.db$trees, function(x) x$tip.label)
        pz.db$ntaxa <- length(pz.db$trees)
        pz.message(paste0(
            "  ..........",
            pz.db$ntaxa,
            " tree(s) remain after observed-taxa filtering"
        ))
    }
    list(phenotype=phenotype, phenoP=phenoP, mapped.observed=mapped.observed, pheno_sd=pheno_sd)
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
add.below.LOD <- function(pz.db, abd.meta, ..., .opts=NULL) {
    opts <- pz.resolve.options(..., .opts=.opts)
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
    # Make sure still binary, if appropriate
    if ((opts('which_phenotype') != 'abundance') &&
        tolower(opts('core_method')) != "poms") {
        pz.message("Binarizing input data...", 2)
        # binarize to save memory usage since we care about pres/abs
        abd.meta$mtx <- Matrix::Matrix(abd.meta$mtx > 0)
    }
    return(abd.meta)
}

#' Generate logit-AUC phenotype from Wilcoxon test on clr-transformed data.
#' Mostly used as an approximate phenotype for plotting POMS output.
#' 
#' @param abd.meta A list consisting of a taxon abundance matrix and the
#'   metadata.
#' @return An updated version of \code{abd.meta}.
#' @export
logit_auc_pheno <- function(abd.meta,
                            ...,
                            .opts=NULL) {
    opts <- pz.resolve.options(..., .opts=.opts)
    E <- opts('env_column')
    D <- opts('dset_column')
    S <- opts('sample_column')
    envir <- opts('which_envir')
    if (!(envir %in% levels(abd.meta$metadata[[E]]))) {
        stop(paste0("environment ", envir, " not found in metadata"))
    }
    env.rows <- (abd.meta$metadata[[E]] == envir)
    dsets <- unique(abd.meta$metadata[env.rows, D, drop=TRUE])
    if (length(dsets) > 1) {
        warning("datasets are ignored when calculating wilcox phenotype")
    }
    meta.present <- abd.meta$metadata[(abd.meta$metadata[[S]] %>%
                                       as.character %in%
                                       colnames(abd.meta$mtx)), ]
    envirs <- unique(meta.present[[E]])
    ids <- sapply(colnames(abd.meta$mtx),
                  function (sn) {
                      abd.meta$metadata[(abd.meta$metadata[[S]] == sn), E]
                  })
    if (length(unique(ids)) < 2) { 
        stop("error: only one environment found")
    }
    names(ids) <- colnames(abd.meta$mtx)
    ids <- simplify2array(ids)
    env.cols <- (ids == envir)
    nenv.cols <- (ids != envir)
    clr_mtx <- clr(abd.meta$mtx, pc = 1)
    logit_auc <- apply(clr_mtx, 1, \(x) {
        st <- wilcox.test(x[env.cols], x[nenv.cols])$statistic
        # convert to AUC, then take logit
        logit(st / (sum(env.cols) * sum(nenv.cols)))
    })
    return(logit_auc)
}
