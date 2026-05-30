#--- Functions to load and process data and databases ---#

#' Read in abundance and metadata file(s).
#'
#' Read in abundance and metadata, either as one BIOM-format file or as two
#' tab-delimited files.
#'
#' @details This function uses package-wide options (see \code{?pz.options}),
#'     which can be overridden using the \code{...} argument. Some particularly
#'     relevant options are:
#'
#' \describe{
#'   \item{env_column}{String. Name of column in metadata file containing the
#'     environment annotations.}
#'   \item{dset_column}{String. Name of column in metadata file containing the
#'     dataset annotations.}
#'   \item{input_format}{String. Whether to look for tabular or BIOM-formatted
#'     data ("tabular" or "biom").}
#'   \item{type}{String. Type of data to use, either "midas" (shotgun) or "16S"
#'     (amplicon).}
#' }
#'
#' @return A list with components \code{mtx} and \code{metadata}, corresponding
#'     to a sparse binary presence/absence matrix (see \code{Matrix} package)
#'     and a metadata data frame.
#' @export
read.abd.metadata <- function(...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    colns <- c(opts('env_column'), opts('dset_column'), opts('sample_column'))
    if (length(unique(colns)) < length(colns)) {
        pz.error(
            "Environment, dataset, and sample columns must all be different"
        )
    }
    if (opts('input_format') == "tabular") {
        pz.message("  .....Reading tabular abundance and metadata files")
        abd.meta <- read.abd.metadata.tabular(...)
    } else if (opts('input_format') == "biom") {
        pz.message("  .....Reading BIOM abundance and metadata file")
        abd.meta <- read.abd.metadata.biom(...)
    } else {
        pz.error(paste0("Invalid input format: ", opts('input_format')))
    }
    
    pz.message("  .....Checking abundance matrix")
    sanity.check.abundance(abd.meta$mtx, ...)
    pz.message("  .....Checking metadata")
    sanity.check.metadata(abd.meta$metadata, ...)
    if (opts('type_16S') == TRUE) {
        pz.message("  .....Processing 16S input")
        abd.meta <- process.16s(abd.meta, ...)
    }
    pz.message("  .....Harmonizing abundance matrix and metadata")
    abd.meta <- harmonize.abd.meta(abd.meta, ...)
    if ((opts('which_phenotype') != 'abundance') &&
        tolower(opts('core_method')) != "poms") {
        pz.message("Binarizing input data...", level=1)
        # binarize to save memory usage since we care about pres/abs
        abd.meta$mtx <- Matrix::Matrix(abd.meta$mtx > 0)
    }
    pz.message("  .....Harmonized data")
    pz.message(paste0(
        "  ..........Final abundance matrix: ",
        nrow(abd.meta$mtx),
        " row(s) x ",
        ncol(abd.meta$mtx),
        " sample(s)"
    ))
    pz.message(paste0(
        "  ..........Final metadata rows: ",
        nrow(abd.meta$metadata)
    ))
    gc()
    return(abd.meta)
}


#' Check and process metadata
#'
#' \code{check.process.metadata} is used to make sure that the metadata
#' satisfies the requirements specified by the global options and to make sure
#' that the metadata are of the correct type.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{env_column}{Name of metadata column containing environment
#'     annotations.}
#'   \item{envir}{Environment of interest.}
#'   \item{dset_column}{Name of metadata column containing dataset annotations.}
#'   \item{single_dset}{Boolean. If true, will assume that all samples come from
#'   a single dataset called \code{dset1} no matter what, if anything, is in
#'   \code{dset_column}.}
#' }
#'
#' @param metadata A data frame of metadata with environment, dataset, and
#'     sample columns corresponding to those in the global options (see
#'     \?pz.options).
#' @return A data frame of metadata, with environment and
#'     dataset columns converted to factors, *unless* calculating correlation
#'     in which case environment column will be cast as numeric. The environment
#'     factor will be automatically re-leveled so that \code{envir} is last, so
#'     that it is guaranteed to appear in the output of differential abundance 
#'     estimators.
#' @export
check.process.metadata <- function(metadata, ...) {
    opts <- clone_and_merge(phylogenize:::PZ_OPTIONS, ...)
    orig_md <- metadata
    E <- opts('env_column')
    S <- opts('sample_column')
    envir <- opts('which_envir')
    if (!(opts('env_column') %in% colnames(metadata))) {
        pz.error(paste0("environment column not found: ", opts('env_column')))
    }
    if (opts('single_dset')) {
        metadata[[opts('dset_column')]] <- rep("dset1", nrow(metadata))
    }
    if (!(opts('dset_column') %in% colnames(metadata))) {
        pz.error(paste0("dataset column not found: ", opts('dset_column')))
    }
    converted_rownames <- FALSE
    if (!(opts('sample_column') %in% colnames(metadata))) {
        if (!is.null(rownames(metadata))) {
            pz.warning(paste0("sample column not found: ",
                              opts('sample_column'),
                              "; assuming row names are sample IDs"))
            metadata[[opts('sample_column')]] <- rownames(metadata)
	    converted_rownames <- TRUE
        } else {
            pz.error(paste0("sample column not found: ", opts('sample_column')))
        }
    }
    if (opts('categorical')) {   # i.e., env not continuous
        env_factor <- factor(metadata[[E]])
        env_levels <- levels(env_factor)
        
        if (!(envir %in% env_levels)) {
            pz.error(paste0("environment ", envir, " not found in metadata"))
        }
        env_factor <- forcats::fct_relevel(env_factor, envir, after=Inf)
        metadata[[E]] <- env_factor
    } else {
        metadata[[E]] <- as.numeric(metadata[[E]])
        if (all(is.na(metadata[[E]]))) {
            pz.error(
                paste0("Environment failed conversion to numeric; is this ",
                       "supposed to be a categorical variable?"))
        }
    }
    metadata[[opts('dset_column')]] <- as.factor(metadata[[opts('dset_column')]])
    
    # One more sanity check before we return
    if (converted_rownames) orig_md <- tibble::as_tibble(orig_md, rownames=S)
    compare_md <- dplyr::full_join(orig_md[c(S, E)], metadata[c(S, E)], by=S)
    wrong_rows <- compare_md[
        compare_md[[paste0(E, ".x")]] != compare_md[[paste0(E, ".y")]],
    ]
    if (nrow(wrong_rows) > 0) {
        print(wrong_rows)
        pz.error(
            paste0("Something went wrong when loading metadata and ",
                   "environments no longer match; please report as a bug")
        )
    }
    
    return(metadata)
}

#' Read abundance matrix and metadata (BIOM)
#'
#' Read in taxon-by-sample matrix of abundances and metadata (sample
#' annotations) from a single BIOM-formatted file.
#'
#' @details This function uses package-wide options (see \code{?pz.options}),
#'     which can be overridden using the \code{...} argument. Some particularly
#'     relevant options are:
#'
#' \describe{
#'   \item{biom_file}{String. Name of BIOM abundance-and-metadata file. Default: "test.biom"}
#' }
#'
#' @return A list with components \code{mtx} (matrix of abundances) and
#'     \code{metadata} (data frame of metadata).
#' @keywords internal
read.abd.metadata.biom <- function(...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    bf <- opts('biom_file')
    pz.message(paste0("looking for file: ", normalizePath(bf)), level=2)
    if (!(file.exists(bf))) {
        pz.error(paste0("file not found: ", bf))
    } else { pz.message(paste0("located biom file: ", bf), level=2) }
    # biomf <- biomformat::read_biom(bf)
    pz.message("  ..........Reading BIOM file")
    biomf <- biomformat::read_biom(bf)
    abd.mtx <- biomformat::biom_data(biomf)
    pz.message(paste0(
        "  ..........BIOM abundance matrix: ",
        nrow(abd.mtx),
        " row(s) x ",
        ncol(abd.mtx),
        " sample(s)"
    ))
    if (!opts('separate_metadata')) {
        pz.message("  ..........Reading metadata from BIOM file")
        metadata <- biomformat::sample_metadata(biomf)
        metadata <- check.process.metadata(metadata, ...)
        # work around different naming convention
        if (is.null(rownames(metadata))) {
            pz.error(paste0("metadata had no sample names; should not be ",
                            "possible, check that your biom file is not ",
                            "corrupt"))
        }
    } else {
        mf <- opts('metadata_file')
        if (!(file.exists(mf))) {
            pz.error(paste0("file not found: ", mf))
        } else { pz.message(paste0("located metadata file: ", mf), level=2) }
        pz.message("  ..........Reading separate metadata file")
        metadata <- readr::read_tsv(mf, show_col_types = FALSE)
        metadata <- check.process.metadata(metadata, ...)
    }
    rm(biomf); gc()
    return(list(mtx=abd.mtx, metadata=metadata))
}

#' Read abundance matrix and metadata (tabular)
#'
#' Read in taxon-by-sample matrix of abundances and metadata (sample
#' annotations) from two tab-delimited files.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{dset_column}{Name of metadata column containing dataset annotations.}
#' }
#'
#' @return A list with components \code{mtx} (matrix of abundances) and
#'     \code{metadata} (data frame of metadata).
#' @keywords internal
read.abd.metadata.tabular <- function(...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    af <- opts('abundance_file')
    mf <- opts('metadata_file')
    if (!(file.exists(af))) {
        pz.error(paste0("file not found: ", af))
    } else { pz.message(paste0("located abundance file: ", af), level=2) }
    if (!(file.exists(mf))) {
        pz.error(paste0("file not found: ", mf))
    } else { pz.message(paste0("located metadata file: ", mf), level=2) }
    
    pz.message("  ..........Reading metadata table")
    metadata <- readr::read_tsv(mf, show_col_types = FALSE)
    metadata <- check.process.metadata(metadata, ...)
    pz.message(paste0(
        "  ..........Metadata rows: ",
        nrow(metadata)
    ))

    pz.message("  ..........Reading abundance table")
    abundance_header <- names(readr::read_tsv(
        af,
        n_max=0,
        show_col_types=FALSE
    ))
    if (length(abundance_header) < 2) {
        pz.error("Abundance table must contain a taxon column and at least one sample column")
    }
    metadata_samples <- trimws(as.character(metadata[[opts('sample_column')]]))
    abundance_keep_cols <- c(
        abundance_header[1],
        intersect(abundance_header[-1], metadata_samples)
    )
    if (length(abundance_keep_cols) < 2) {
        pz.error(paste0("No abundance table sample columns matched metadata ",
                        "sample IDs"))
    }
    pz.message(paste0(
        "  ..........Reading ",
        length(abundance_keep_cols) - 1,
        " abundance sample column(s) matching metadata"
    ))
    abundance_col_types <- do.call(
        readr::cols,
        c(
            stats::setNames(list(readr::col_character()), abundance_header[1]),
            list(.default=readr::col_double())
        )
    )
    abd.df <- readr::read_tsv(
        af,
        col_select=tidyselect::all_of(abundance_keep_cols),
        col_types=abundance_col_types,
        show_col_types=FALSE
    )
    abd.values <- abd.df[, -1, drop=FALSE]
    bad.cols <- names(abd.values)[
        vapply(abd.values, function(x) any(is.na(x)), logical(1))
    ]
    if (length(bad.cols) > 0) {
        pz.error(paste0(
            "Abundance table contains nonnumeric value(s) in column(s): ",
            paste(bad.cols, collapse=", ")
        ))
    }
    abd.mtx <- as.matrix(abd.values)
    rownames(abd.mtx) <- abd.df[[1]]
    pz.message(paste0(
        "  ..........Abundance matrix: ",
        nrow(abd.mtx),
        " row(s) x ",
        ncol(abd.mtx),
        " sample(s)"
    ))
    return(list(mtx=abd.mtx, metadata=metadata))
    
}

#--- Databases ---#

#' Import the data necessary for *phylogenize* analysis.
#'
#' \code{import.pz.db} decides based on global options which data files to
#' import.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{type_16S}{Boolean. If 16S data, TRUE, otherwise shotgun data is
#'     assumed. Default: FALSE}
#'   \item{db}{String. Which database to use ("gtdb" or "uhgp")).}
#'   \item{data_dir}{String. Path to directory containing the data files
#'   required to perform a \emph{phylogenize} analysis.}
#' }
#' @return A list of the data objects required to perform a *phylogenize*
#'     analysis, with components \code{gene.presence}, \code{trees},
#'     \code{taxonomy}, \code{g.mappings}, and \code{gene.to.fxn}.
#' @export
import.pz.db <- function(...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    db_csv <- file.path(opts('data_dir'), "databases.csv")
    pz.message(paste0("  .....Reading database index: ", db_csv), level=2)
    installed_dbs <- readr::read_delim(db_csv, show_col_types = FALSE)
    requested_db <- tolower(opts('db'))
    pz.message(paste0("  .....Requested database: ", requested_db))
    if (!(requested_db %in% installed_dbs[["database"]])) {
        pz.error(paste0("Database not installed in ", opts('data_dir'), ": ",
                        opts('db')))
    }
    found_db <- dplyr::filter(installed_dbs, database==requested_db)
    if (nrow(found_db) > 1) {
        pz.error(paste0("Duplicate data entries in db file: ", db_csv))
    }

    pz.message("  .....Read in phylogenize2 gene presence file")
    gene.presence <- readRDS(file.path(opts('data_dir'),
                                       found_db[["genes"]]))
    gene.presence <- gene.presence[names(gene.presence) != ""]
    pz.message(paste0(
        "  ..........Loaded ",
        length(gene.presence),
        " gene presence matrix/matrices"
    ))
    pz.message("  .....Read in phylogenize2 tree file")
    trees <- readRDS(file.path(opts('data_dir'),
                               found_db[["trees"]]))
    pz.message(paste0("  ..........Loaded ", length(trees), " tree(s)"))
    pz.message("  .....Read in phylogenize2 taxonomy file")
    taxonomy <- readr::read_csv(file.path(opts('data_dir'),
                                          found_db[["taxonomy"]]), show_col_types=FALSE)
    pz.message(paste0(
        "  ..........Loaded taxonomy with ",
        nrow(taxonomy),
        " row(s)"
    ))
    # Workaround for inconsistency with "_" vs. " " in pipeline...,
    names(trees) <- gsub(" ", "_", names(trees))
    names(gene.presence) <- gsub(" ", "_", names(gene.presence))
    # propagate cluster values up, as higher level taxonomic names may be missing
    taxonomy <- taxonomy %>% 
      rowwise() %>%
      mutate(species = ifelse(is.na(species), cluster, species)) %>%
      mutate(genus = ifelse(is.na(genus), species, genus)) %>%
      mutate(family = ifelse(is.na(family), genus, family)) %>%
      mutate(order = ifelse(is.na(order), family, order)) %>%
      mutate(class = ifelse(is.na(class), order, class)) %>%
      mutate(phylum = ifelse(is.na(phylum), class, phylum)) %>%
      mutate(domain = ifelse(is.na(domain), phylum, domain)) %>%
      mutate(across(domain:genus, \(x) gsub(" ", "_", x))) %>%
      ungroup()

    pz.message("  .....Read in phylogenize2 gene functions file")
    gene.to.fxn <- readr::read_csv(
        file.path(opts('data_dir'), found_db[["functions"]]),
        show_col_types = FALSE
    )
    pz.message(paste0(
        "  ..........Loaded ",
        nrow(gene.to.fxn),
        " gene function row(s)"
    ))
    # Check if the files exist instead of throwing a null error
    if (is.null(gene.presence) |
        is.null(trees) |
        is.null(taxonomy) |
        is.null(gene.to.fxn)) {
        pz.error(paste0("One or more of the database files used in phylogenize",
        " is not found. Please check these have been downloaded and are in the", 
        " expected location given by ", db_csv))
    }
    pz.message("  .....Remove extra info, make table a little more readable")
    # remove extra info, make table a little more readable
    gene.to.fxn <- gene.to.fxn %>%
        dplyr::rename(gene = node_head) %>%
        dplyr::select(gene, accession, `function`) %>%
        distinct()

    # Make the files at the user requested taxon level. Here we adjust the
    # classification level to look at i.e family, class, order, etc.
    if (opts('taxon_level') != "phylum") {
	pz.message("  .....Change taxonomy level for analysis")
        gene.presence <- change.presence.tax.level(gene.presence,
                                                   opts('taxon_level'),
                                                   taxonomy)
        trees <- change.tree.tax.level(trees, opts('taxon_level'), taxonomy)
        trees <- Filter(function(tr) length(tr$tip.label) > 1, trees)
        pz.message(paste0(
            "  ..........Taxon-level database now has ",
            length(trees),
            " tree(s)"
        ))
    }
    
    # filter based on the minimum number of observations
    pz.message("  .....Filter based on the minimum number of observations")
    gene.presence <- above_minimum_genes(gene.presence, trees)
    pz.message(paste0(
        "  ..........",
        length(gene.presence),
        " gene presence matrix/matrices remain"
    ))
   
    pz.message("  .....Drop any taxa that got culled in tree object") 
    # drop any taxa that got culled in above_minimum_genes
    trees <- trees[intersect(names(trees), names(gene.presence))]
    pz.message(paste0("  ..........", length(trees), " tree(s) remain"))
    
    return(list(gene.presence = gene.presence,
                trees = trees,
                taxonomy = taxonomy,
                gene.to.fxn = gene.to.fxn))
}

#' Clean up imported database.
#'
#' \code{adjust.db} removes any species with fewer than \code{opts('treemin')}
#' representatives, removes tips from trees that were not observed in the data,
#' and resolves polytomies. It also adds a couple of variables into the database
#' that give lists of species per tree and the total number of remaining taxa
#'
#' @param pz.db A database (typically obtained with \code{import.pz.db}).
#' @param abd.meta A list consisting of a species abundance matrix and the
#'     metadata.
#' @return An updated database.
#' @export
adjust.db <- function(pz.db, abd.meta, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    species.observed <- rownames(abd.meta$mtx)
    pz.message(paste0(
        "  .....Matching database to ",
        length(species.observed),
        " observed taxon/taxa"
    ))
    
    species.per.tree <- lapply(pz.db$trees, function(tr) {
        intersect(tr$tip.label, species.observed)
    })
    tL <- vapply(species.per.tree, length, 1L)
    if (all(tL < opts('treemin'))) {
        pz.error(paste0("Too few species found. Was the right database used? ",
                        "Are taxa named properly for this database?"))
    }
    passed.min <- nw(tL >= opts('treemin'))
    totalL <- vapply(pz.db$trees, function(tr) { length(tr$tip.label) }, 1L)
    pct.obs <- mapply(function(x, y) x / y, tL, totalL)
    passed.pct <- nw(pct.obs >= opts('pctmin'))
    pz.message("Determining which taxa to test...", level=1)
    for (tn in 1:length(pz.db$trees)) {
        pz.message(paste0(names(pz.db$trees)[tn],
                          " (pct): ", format(pct.obs[tn] * 100, digits=2),
                          "%; (number): ", tL[tn],
                          "; ", ifelse((pct.obs[tn] >= opts('pctmin') &&
                                            tL[tn] >= opts('treemin')),
                                       yes="kept",
                                       no="dropped")
        ), level=2)
    }
    saved.taxa <- intersect(passed.min, passed.pct)
    pz.message(paste0(
        "  .....",
        length(saved.taxa),
        " taxon/taxa passed tree-size and observation filters"
    ))
    pz.message(
        paste0(
            "  Taxa selected: ",
            paste(saved.taxa, collapse=", ")
        ),
        level=2
    )
    if (length(saved.taxa) == 0) {
        pz.error(paste0("All trees had less than ",
                        format(opts('pctmin') * 100, digits=2),
                        "% of species observed. Either a very small ASV ",
                        "table was provided, read depth was very shallow, ",
                        "the right database was not selected or very few ASVs ",
                        "mapped to entries in the database."))
    }
    pz.db$trees <- lapply(pz.db$trees[saved.taxa], function(tr) {
	tips <- intersect(tr$tip.label, species.observed)
	ape::keep.tip(tr, tips)
    })
    pz.message("  .....Filtering gene presence matrices to retained trees")
    pz.db$gene.presence <- above_minimum_genes(pz.db$gene.presence, pz.db$trees)
    pz.db$trees <- pz.db$trees[intersect(names(pz.db$trees), names(pz.db$gene.presence))]
    pz.db$species <- lapply(pz.db$trees, function(x) x$tip.label)
    pz.db$ntaxa <- length(pz.db$trees)
    pz.message(paste0(
        "  .....",
        pz.db$ntaxa,
        " taxon/taxa remain after database adjustment"
    ))
    pz.message("  .....Resolving tree edge cases")
    pz.db$trees <- lapply(pz.db$trees, fix.tree)
    # Re-run this to drop any genes that shouldn't be run
    pz.db
}



#' Changes the presence / absence binary for the desired taxonomic level given
#' by the user.
#'
#' \code{import.pz.db} filters down the binary file for database
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{binary}{Matrix. Object that represents a phylogenize-prepared
#'     internal database}
#'   \item{taxon}{String. Can be either be 'phylum', 'class', 'order', 'family',
#'     or 'genus'}
#'   \item{tax}{Dataframe. Object with all the taxonomy for the given database.
#'     Has column matching the taxon parameter}
#' }
#' @return list of binary matrices that are ready for use at the provided
#'   taxonomic level
#' @export
change.presence.tax.level <- function(binary, taxon, tax){
    # Make a mapping file that is at the taxonomic level selected from the tax
    # file.
    clean <- tax %>%
        dplyr::select(cluster, tidyselect::all_of(taxon), phylum) %>%
        dplyr::distinct()
    # Drop empty values from the taxonomic level selected if they are not
    clean <- clean[!(is.na(clean[[taxon]]) | clean[[taxon]] == ""), ]
    # Arrange them so that the runtime is slightly less in the lookup
    clean <- clean %>%
        dplyr::arrange(phylum, !!(rlang::sym(taxon))) %>%
        dplyr::mutate(cluster = as.character(cluster))
    # Split the dataframe by phylum so that we only need to lookup the taxon
    # associated with each one.
    clean <- clean %>%
        split(.$phylum)
    
    # For every species in the dataframe, if the species is selected, pop the
    # column and remove it from the sparse matrix.
    binary_matrices <- list()
    taxon <- rlang::sym(taxon)
    
    for (i in 1:length(names(binary))) {
        name <- names(binary)[i]
        b <- binary[[name]]
        t <- clean[[name]]
        if (is.null(t)) {
            pz.warning(paste("Warning: No data found for taxon classification",
                             name))
            next
        }
        
        # make a named list in one step
        named_tax_list <- t %>%
            dplyr::select(-phylum) %>%
            dplyr::group_by(!!taxon) %>%
            tidyr::nest() %>%
            tibble::deframe() %>%
            purrr::map(~ dplyr::pull(.x, cluster))
        
        columns <- colnames(binary[[name]])
        
        for (j in 1:length(named_tax_list)){
            shared <- intersect(named_tax_list[[j]], columns)
            if (length(shared) > 1) {
                binary_matrices[[names(named_tax_list)[j]]] <- b[,
                                                                 shared,
                                                                 drop = FALSE]
            }
        }
    }
    binary_matrices <- Filter(function(b) length(colnames(b)) > 1,
                              binary_matrices)
    return(binary_matrices)
}


#' Changes the internal database tree for the desired taxonomic level given by
#' the user.
#'
#' \code{import.pz.db} filters down the binary file for database
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{tree}{List. Object that represent a phylogenize-prepared internal
#'     database}
#'   \item{taxon}{String. Can be either be 'phylum', 'class', 'order', 'family',
#'     or 'genus'}
#'   \item{tax}{Dataframe. Gives full taxonomy for the given database.
#'     Must have column matching the taxon parameter}
#' }
#' @return list of tree objects that are ready for use at the user given
#'   tax_level
#' @export
change.tree.tax.level <- function(tree, taxon, tax) {
    # Make a mapping file that is at the taxonomic level selected from the tax
    # file.
    clean <- tax %>%
        dplyr::select(cluster, taxon, phylum) %>%
        dplyr::distinct()
    pz.message(head(clean), level=3)
    # Drop empty values from the taxonomic level selected if they are not 
    clean <- clean[!(is.na(clean[[taxon]]) | clean[[taxon]] == ""), ]
    # Arrange them so that the runtime is slightly less in the lookup
    clean <- clean %>%
        dplyr::group_by(
            dplyr::across(tidyselect::all_of(taxon))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(phylum, !!(rlang::sym(taxon))) %>%
        dplyr::mutate(cluster = as.character(cluster))
    # Split the dataframe by phylum so that we only need to lookup the taxon
    # associated with each one.
    clean <- clean %>% 
        split(.$phylum)
    
    tree_matrices <- list()
    names <- names(tree)
    for (i in seq_along(tree)) {
        name <- names[i]
        tr <- tree[[name]]
        t <- clean[[name]]
        
        if (is.null(clean[[name]])) {
            pz.message(names(clean), level=3)
            pz.warning(paste("Warning: No data found for taxon classification",
                             name), level=2)
            next 
        } else {
            pz.message(paste("Good news: Data found for taxon classification",
                             name), level=2)
        }	
        # Generate a list of unique split names based on taxon
        split_names <- t %>%
            dplyr::group_split(!!(rlang::sym(taxon))) %>%
            purrr::map(~ dplyr::pull(.x, !!(rlang::sym(taxon)))) %>%
            unlist() %>%
            unique()
        
        for (j in seq_along(split_names)) {
            tips <- tr$tip.label
            split_tips <- t %>%
                dplyr::filter(!!(rlang::sym(taxon)) == split_names[[j]])
            tips <- intersect(tips, split_tips[["cluster"]])
            subtree <- ape::keep.tip(tr, tips)
            if (!is.null(subtree) && length(subtree$tip.label) > 1) {
                tree_matrices[[split_names[j]]] <- subtree    
            }
        }
    }
    tree_matrices <- Filter(function(tr) length(tr$tip.label) > 1,
                            tree_matrices)
    return(tree_matrices)
}


#' Map denoised 16S sequence variants to MIDAS IDs using vsearch.
#'
#' \code{process.16s} is a wrapper for other functions that: output sequence
#' variants to a file; map them using vsearch against a reference database of
#' 16S sequences; then return a list of abundance and metadata values where the
#' rows of the abundance matrix are now MIDAS IDs.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{vsearch_infile}{String. File name of the sequences written to disk
#'     and then read into vsearch.}
#' }
#'
#' @param mtx A presence/absence or abundance matrix, with row names equal to
#'   amplicon sequence variant DNA sequences.
#' @return none
#' @export
process.16s <- function(abd.meta, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    if (!(all(is.dna(rownames(abd.meta$mtx))))) {
        pz.error(paste0("expected rows to be DNA sequences but found illegal ",
                        "characters"))
    }
    if (opts('which_16s_method')=="vsearch") {
        prepare.vsearch.input(abd.meta$mtx, ...)
        run.vsearch(...)
        results_16s <- get.vsearch.results(...)
    } else if (opts('which_16s_method')=="appspam") {
        prepare.vsearch.input(abd.meta$mtx, ...)
        run.appspam(...)
        results_16s <- get.appspam.results(...)
    } else if (opts('which_16s_method')=="jplace") {
        results_16s <- get.appspam.results(...)        
    } else {
        pz.error("which_16s_method must be vsearch, appspam, or jplace")
    }
    summed.uniq <- sum.nonunique.vsearch(results_16s, abd.meta$mtx, ...)
    csu <- colSums(summed.uniq)
    abd.meta$mtx <- summed.uniq[, which(csu > 0), drop=FALSE]
    # don't convert to relative abundance...
    # abd.meta$mtx <- apply(summed.uniq[, which(csu > 0), drop=FALSE],
    #                       2,
    #                       function(x) x / sum(x))
    abd.meta
    
}

#' Check metadata and abundance matrix against one another
#'
#' \code{harmonize.abd.meta} compares the abundance matrix to the metadata
#' matrix to make sure that enough samples are in common between the two to
#' perform an \emph{phylogenize} analysis, after dropping any singleton datasets
#' or environments (effects for these cannot be estimated).
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{env_column}{Name of metadata column containing environment
#'   annotations.}
#'   \item{dset_column}{Name of metadata column containing dataset annotations.}
#' }
#'
#' @param abd.meta A list with components \code{mtx} and \code{metadata},
#'     corresponding to a sparse binary presence/absence matrix (see Matrix
#'     package) and a metadata data frame.
#' @return A list of the same form as \code{abd.meta}.
#' @export
harmonize.abd.meta <- function(abd.meta, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    align.metadata.to.matrix <- function(abd.meta) {
        metadata_order <- match(colnames(abd.meta$mtx),
                                abd.meta$metadata[[opts('sample_column')]])
        abd.meta$metadata <- abd.meta$metadata[metadata_order, , drop=FALSE]
        abd.meta
    }
    abd.meta$metadata[[opts('sample_column')]] <- trimws(
        abd.meta$metadata[[opts('sample_column')]]
    )
    colnames(abd.meta$mtx) <- trimws(colnames(abd.meta$mtx))
    samples.present <- intersect(
        abd.meta$metadata[[opts('sample_column')]],
        colnames(abd.meta$mtx))
    pz.message(paste0(
        "  ..........Samples present in both metadata and matrix: ",
        length(samples.present)
    ))
    if (length(samples.present) == 0) {
        pz.error(paste0("No samples found in both metadata and ",
                        "abundance matrix; check for illegal characters ",
                        "in sample ID column"))
    }
    abd.meta$mtx <- abd.meta$mtx[, samples.present, drop=FALSE]
    abd.meta$metadata <- abd.meta$metadata[
        abd.meta$metadata[[opts('sample_column')]] %in%
            samples.present, ]
    abd.meta <- align.metadata.to.matrix(abd.meta)
    
    if (opts('which_phenotype') %in% c("specificity",
                                       "prevalence",
                                       "abundance")) {
        env.number <- table(abd.meta$metadata[[opts('env_column')]])
        singleton.envs <- names(which(env.number == 1))
        if (length(singleton.envs) > 0) {
            pz.warning(paste0(
                "Warning: each environment requires at least two samples."))
            pz.warning(
                "Dropped the following environment(s) from the analysis: ")
            for (s in singleton.envs) {
                pz.warning(s)
            }
        }
        nonsingleton.envs <- names(which(env.number > 1))
        pz.message(paste0(length(nonsingleton.envs),
                          " non-singleton environment(s) found"), level=1)
        if ((length(nonsingleton.envs) < 2) &&
            (opts('which_phenotype') == 'specificity')) {
            pz.error(paste0(
                "In order to calculate specificity, there must be at least",
                " two environments with two samples each."))
        }
        if ((length(nonsingleton.envs) < 1) &&
            (opts('which_phenotype') == 'prevalence')) {
            pz.error(paste0(
                "In order to calculate prevalence, there must be at least",
                " one environment with two samples."))
        }
    } else {
        # Must be correlation or provided
        n_present <- sum(!is.na(as.numeric(
            abd.meta$metadata[[opts('env_column')]])))
        if (n_present < 3) pz.error(paste0(
            "There must be at least three non-missing samples in an ",
            "environment"))
    }
    
    pz.message("  .....collecting all unique dsets")
    dset.number <- table(abd.meta$metadata[[opts('dset_column')]])
    nonsingleton.dsets <- names(which(dset.number > 1))
    pz.message(paste0(length(nonsingleton.dsets),
                      " non-singleton dataset(s) found"), level=2)
    f_dsets <- (abd.meta$metadata[[opts('dset_column')]] %in%
                    nonsingleton.dsets)
    if (!(opts('which_phenotype') %in% c("correlation", "provided"))) {
        f_envs <- (abd.meta$metadata[[opts('env_column')]] %in%
                       nonsingleton.envs)
        wrows <- which(f_envs & f_dsets)
    } else {
        wrows <- which(f_dsets)
    }
    if (length(wrows) < 2) {
        pz.error(paste0("Too few rows found in metadata matrix after ",
                        "dropping singletons and matching with abundance ",
                        "matrix (need at least 2)"))
    }
    pz.message("  .....Finding metadata and matrix intersections")
    abd.meta$metadata <- abd.meta$metadata[wrows, , drop=FALSE]
    wcols <- intersect(colnames(abd.meta$mtx),
                       abd.meta$metadata[[opts('sample_column')]])
    pz.message(paste0(
        "  ..........Samples retained after singleton filtering: ",
        length(wcols)
    ))
    if (length(wcols) < 2) {
        pz.error(paste0("Too few columns found in abundance matrix ",
                        "after dropping singletons and matching with metadata ",
                        "(need at least 2)"))
    }
    abd.meta$mtx <- abd.meta$mtx[, wcols, drop=FALSE]
    abd.meta <- align.metadata.to.matrix(abd.meta)
    pz.message("  .....Removing all rows and columns that are completely 0")
    abd.meta$mtx <- remove.allzero.abundances(abd.meta$mtx)
    abd.meta <- align.metadata.to.matrix(abd.meta)
    return(abd.meta)
}


#' Sanity-check abundance data
#'
#' \code{sanity.check.abundance} is used to make sure that the abundance matrix
#' satisfies the requirements specified by the \emph{phylogenize} application.
#'
#' @param abd.mtx A matrix or Matrix of abundance or presence values (double or
#'     logical).
#' @return Always returns TRUE, but will throw errors if the
#'     abundance data is the wrong type or class.
#' @export
sanity.check.abundance <- function(abd.mtx, ...) {
    if ((!methods::is(abd.mtx, "matrix")) & (!methods::is(abd.mtx, "Matrix"))) {
        pz.error(paste0(
            "Abundance matrix must be a matrix of logicals/doubles/integers ",
            "and instead was a ",
            class(abd.mtx)))
    }
    if (methods::is(abd.mtx, "Matrix")) {
        if (!(typeof(abd.mtx@x[1]) %in% c("logical", "double", "integer"))) {
            pz.error(paste0(
                "Abundance matrix must be a matrix of logicals/doubles/",
                "integers and instead contained ",
                typeof(abd.mtx@x[1])))
        }
    } else if (methods::is(abd.mtx, "matrix")) {
        if (!(typeof(abd.mtx) %in% c("logical", "double", "integer"))) {
            pz.error(paste0(
                "Abundance matrix must be a matrix of logicals/doubles/",
                "integers and instead contained ",
                typeof(abd.mtx)))
        }
    }
    if (ncol(abd.mtx) < 2) {
        pz.error("Abundance matrix had fewer than 2 columns")
    }
    if (nrow(abd.mtx) < 2) {
        pz.error("Abundance matrix had fewer than 2 rows")
    }
    if (is.null(dimnames(abd.mtx))) {
        pz.error(
            "Abundance matrix is lacking row (taxon) and column (sample) names")
    }
    if (is.null(rownames(abd.mtx))) {
        pz.error("Abundance matrix is lacking row (taxon) names")
    }
    if (is.null(colnames(abd.mtx))) {
        pz.error("Abundance matrix is lacking column (sample) names")
    }
    return(TRUE)
}


#' Check that dataset, environment, and sample columns all present
#'
#' \code{sanity.check.abundance} is used to make sure that the metadata data
#' frame satisfies the requirements specified by the \emph{phylogenize}
#' application.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{env_column}{Name of metadata column containing environment
#'     annotations.}
#'   \item{dset_column}{Name of metadata column containing dataset
#'     annotations.}
#' }
#'
#' @param metadata A data frame giving sample annotations.
#' @return Always returns TRUE, but will throw errors if the metadata does not
#'   match specifications.
#' @export
sanity.check.metadata <- function(metadata, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    if (!(opts('env_column') %in% colnames(metadata))) {
        pz.error(
            paste0("When looking for environment, no column found labeled ",
                   opts('env_column'))
        )
    }
    if (!(opts('dset_column') %in% colnames(metadata))) {
        pz.error(
            paste0("When looking for dataset, no column found labeled ",
                   opts('dset_column'))
        )
    }
    if (!(opts('sample_column') %in% colnames(metadata))) {
        pz.error(paste0(
            "When looking for dataset, no column found labeled ",
            opts('sample_column'))
        )
    }
    if (nrow(metadata) < 2) {
        pz.error("Fewer than two rows found in metadata")
    }
    return(TRUE)
}


#' Remove rows and columns of a matrix that are all zero.
#'
#' \code{remove.allzero.abundances} removes all rows and columns of a matrix
#' where every observation is zero, starting with columns and then proceeding to
#' rows.
#'
#' @param abd.mtx A matrix of abundance values (double or logical).
#' @return A matrix of abundance values (double), with all-zero columns and rows
#'     removed.
#' @export
remove.allzero.abundances <- function(abd.mtx, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    cs <- Matrix::colSums(abd.mtx)
    nz.cols <- which(cs > 0)
    z.col.logical <- (cs == 0)
    
    pz.message("  ..........calculating column with zero")
    if (sum(z.col.logical) > 0) {
        pz.warning(paste0("Dropping ", sum(z.col.logical),
                          " column(s), since no mapped taxa had observations"))
        if (!is.null(names(z.col.logical))) {
            pz.warning("Columns dropped: ")
            for (n in names(which(z.col.logical))) {
                pz.warning(n)
            }
        }
    }
    if (length(nz.cols) < 2) {
        pz.error("Too few columns with at least one non-zero entry")
    }

    pz.message("  ..........Final drop of columns")
    abd.mtx <- abd.mtx[, nz.cols, drop = FALSE]

    pz.message("  ..........calculating rows with zero")
    rs <- Matrix::rowSums(abd.mtx)
    nz.rows <- which(rs > 0)

    #nz.rows <- which(Matrix::rowSums(abd.mtx[, nz.cols, drop=FALSE]) > 0)
    if (length(nz.rows) < 2) {
        pz.error("Too few rows with at least one non-zero entry")
    } else {
	pz.message("  ..........Final drop of rows")
        n_drop <- nrow(abd.mtx) - length(nz.rows)
        pz.message(paste0("  ..........dropping ", n_drop, " all-zero row(s) out of ", nrow(abd.mtx)))
        abd.mtx <- abd.mtx[nz.rows, , drop = FALSE]
    }
    
    return(abd.mtx)
}

#' Copy and/or unzip downloaded data into the Phylogenize2 database directory.
#'
#' @param path Provide a path to a local file containing a .zip file (which will
#'   be unzipped) or a .csv file (which will be copied).
#' @param force Boolean; overwrite existing files? (Default: FALSE)
#' @export
install_data <- function(path, force=FALSE) {
    fn <- basename(path)
    extd_path <- system.file("extdata/", package="phylogenize")
    if (stringr::str_ends(basename(path), "\\.csv")) {
        if (!file.exists(file.path(extd_path, fn)) || force) {
            print(paste0(
                "Copying ",
                path,
                " to ",
                system.file("", package="phylogenize")))
            file.copy(path, extd_path, overwrite = force)
        }
    } else {
        print(paste0(
            "Unzipping ",
            path,
            " to ",
            system.file("", package="phylogenize")))
        unzip(path, overwrite = force, exdir = extd_path)
    }
}
