#' @import Matrix
NULL

#--- I/O and initial processing ---#

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
        pz.error("Environment, dataset, and sample columns must all be different")
    }
    if (opts('input_format') == "tabular") {
        abd.meta <- read.abd.metadata.tabular(...)
    } else if (opts('input_format') == "biom") {
        abd.meta <- read.abd.metadata.biom(...)
    } else {
        pz.error(paste0("Invalid input format: ", opts('input_format')))
    }
    sanity.check.abundance(abd.meta$mtx, ...)
    sanity.check.metadata(abd.meta$metadata, ...)
    if (opts('type_16S') == TRUE) {
        abd.meta <- process.16s(abd.meta, ...)
    }
    abd.meta <- harmonize.abd.meta(abd.meta, ...)
    # binarize to save memory usage since we care about pres/abs
    abd.meta$mtx <- Matrix::Matrix(abd.meta$mtx > 0)
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
#'   \item{env_column}{Name of metadata column containing environment annotations.}
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
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    E <- opts('env_column')
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
    if (!(opts('sample_column') %in% colnames(metadata))) {
        if (!is.null(rownames(metadata))) {
            pz.warning(paste0("sample column not found: ",
                              opts('sample_column'),
                              "; assuming row names are sample IDs"))
            metadata[[opts('sample_column')]] <- rownames(metadata)
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
      levels(env_factor) <- c(setdiff(env_levels, envir), envir)
      metadata[[E]] <- env_factor
    } else {
      metadata[[E]] <- as.numeric(metadata[[E]])
      if (all(is.na(abd.meta$metadata[[E]]))) {
        pz.error("Environment failed conversion to numeric; is this supposed to be a categorical variable?")
      }
    }
    metadata[[opts('dset_column')]] <- as.factor(metadata[[opts('dset_column')]])
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
    pz.message(paste0("looking for file: ", normalizePath(bf)))
    if (!(file.exists(bf))) {
        pz.error(paste0("file not found: ", bf))
    } else { pz.message(paste0("located biom file: ", bf)) }
    # biomf <- biomformat::read_biom(bf)
    biomf <- biomformat::read_biom(bf)
    abd.mtx <- biomformat::biom_data(biomf)
    if (!opts('separate_metadata')) {
      metadata <- biomformat::sample_metadata(biomf)
      metadata <- check.process.metadata(metadata, ...)
      # work around different naming convention
      if (is.null(rownames(metadata))) {
          pz.error(paste0("metadata had no sample names; should not be possible,",
                          " check that your biom file is not corrupt"))
      }
    } else {
        mf <- opts('metadata_file')
        if (!(file.exists(mf))) {
            pz.error(paste0("file not found: ", mf))
        } else { pz.message(paste0("located metadata file: ", mf)) }
        metadata <- data.frame(data.table::fread(mf))
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
    } else { pz.message(paste0("located abundance file: ", af)) }
    if (!(file.exists(mf))) {
        pz.error(paste0("file not found: ", mf))
    } else { pz.message(paste0("located metadata file: ", mf)) }
    abd.mtx <- fastread(af, cn=FALSE)
    metadata <- data.frame(data.table::fread(mf))
    metadata <- check.process.metadata(metadata, ...)
    return(list(mtx=abd.mtx, metadata=metadata))
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
                "Abundance matrix must be a matrix of logicals/doubles/integers ",
                "and instead contained ",
                typeof(abd.mtx@x[1])))
        }
    } else if (methods::is(abd.mtx, "matrix")) {
        if (!(typeof(abd.mtx) %in% c("logical", "double", "integer"))) {
            pz.error(paste0(
                "Abundance matrix must be a matrix of logicals/doubles/integers ",
                "and instead contained ",
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
        pz.error("Abundance matrix is lacking row (taxon) and column (sample) names")
    }
    if (is.null(rownames(abd.mtx))) {
        pz.error("Abundance matrix is lacking row (taxon) names")
    }
    if (is.null(colnames(abd.mtx))) {
        pz.error("Abundance matrix is lacking column (sample) names")
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
    nz.rows <- which(Matrix::colSums(abd.mtx[, nz.cols, drop=FALSE]) > 0)
    if (length(nz.cols) < 2) {
        pz.error("Too few rows with at least one non-zero entry")
    }
    # pass
    return(abd.mtx)
}

#' Check that dataset, environment, and sample columns all present
#'
#' \code{sanity.check.abundance} is used to make sure that the metadata data frame satisfies the requirements specified by the \emph{phylogenize} application.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{env_column}{Name of metadata column containing environment annotations.}
#'   \item{dset_column}{Name of metadata column containing dataset annotations.}
#' }
#'
#' @param metadata A data frame giving sample annotations.
#' @return Always returns TRUE, but will throw errors if the metadata does not match specifications.
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
    samples.present <- intersect(abd.meta$metadata[[opts('sample_column')]],
                                 colnames(abd.meta$mtx))
    if (length(samples.present) == 0) {
        pz.error(paste0("No samples found in both metadata and ",
                        "abundance matrix; check for illegal characters ",
                        "in sample ID column"))
    }
    abd.meta$mtx <- abd.meta$mtx[, samples.present, drop=FALSE]
    abd.meta$metadata <- abd.meta$metadata[
                            abd.meta$metadata[[opts('sample_column')]] %in%
                            samples.present, ]

    if (opts('which_phenotype') %in% c("specificity", "prevalence", "abundance")) {
        all.envs <- unique(abd.meta$metadata[[opts('env_column')]])
        env.number <- sapply(all.envs, function(e) {
            sum(abd.meta$metadata[[opts('env_column')]] == e)
        })
        names(env.number) <- all.envs
        singleton.envs <- names(which(env.number == 1))
	if (length(singleton.envs) > 0) {
            pz.warning(paste0("Warning: each environment requires at least two samples."))
            pz.warning("Dropped the following environment(s) from the analysis: ")
            for (s in singleton.envs) {
                pz.warning(s)
            }
        }
        nonsingleton.envs <- names(which(env.number > 1))
        pz.message(paste0(length(nonsingleton.envs),
                        " non-singleton environment(s) found"))
        if ((length(nonsingleton.envs) < 2) &&
            (opts('which_phenotype') == 'specificity')) {
            pz.error(paste0("In order to calculate specificity, there must be at least",
                            " two environments with two samples each."))
        }
        if ((length(nonsingleton.envs) < 1) &&
            (opts('which_phenotype') == 'prevalence')) {
            pz.error(paste0("In order to calculate prevalence, there must be at least",
                            " one environment with two samples."))
        }
    } else {
      # Must be correlation
      n_present <- sum(!is.na(as.numeric(abd.meta$metadata[[opts('env_column')]])))
      if (n_present < 3) pz.error("In order to calculate correlation, there must be at least 3 non-missing values.")
    }

    all.dsets <- unique(abd.meta$metadata[[opts('dset_column')]])
    dset.number <- sapply(all.dsets, function(d) {
        sum(abd.meta$metadata[[opts('dset_column')]] == d)
    })
    names(dset.number) <- all.dsets
    nonsingleton.dsets <- names(which(dset.number > 1))
    pz.message(paste0(length(nonsingleton.dsets),
                      " non-singleton dataset(s) found"))
    if (opts('which_phenotype') != "correlation") {
      wrows <- which(
      (abd.meta$metadata[[opts('env_column')]] %in% nonsingleton.envs) &
      (abd.meta$metadata[[opts('dset_column')]] %in% nonsingleton.dsets))
    } else {
      wrows <- which(abd.meta$metadata[[opts('dset_column')]] %in% nonsingleton.dsets)
    }
    if (length(wrows) < 2) {
        pz.error(paste0("Too few rows found in metadata matrix after ",
                        "dropping singletons and matching with abundance ",
                        "matrix (need at least 2)"))
    }
    abd.meta$metadata <- abd.meta$metadata[wrows, , drop=FALSE]
    wcols <- intersect(colnames(abd.meta$mtx),
                       abd.meta$metadata[[opts('sample_column')]])
    if (length(wcols) < 2) {
        pz.error(paste0("Too few columns found in abundance matrix ",
                        "after dropping singletons and matching with metadata ",
                        "(need at least 2)"))
    }
    abd.meta$mtx <- abd.meta$mtx[, wcols, drop=FALSE]
    abd.meta$mtx <- remove.allzero.abundances(abd.meta$mtx)
    abd.meta
}

#--- Process 16S data---#

#' Prepare input file for alignment
#'
#' \code{prepare.vsearch.input} outputs a FASTA file of the sequences in the input
#' 16S data for analysis using vsearch or vsearch.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   input files).}
#'   \item{vsearch_infile}{String. File name of the sequences written to disk and
#'   then read into vsearch/vsearch.}
#' }
#'
#' @param mtx A presence/absence or abundance matrix, with row names equal to
#'     amplicon sequence variant DNA sequences.
#' @export
prepare.vsearch.input <- function(mtx, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    binary = basename(opts('vsearch_bin'))
    check.dna <- is.dna(rownames(mtx))
    if (!(all(check.dna))) {
        pz.error(paste0(
            "Rownames in 16S data must correspond to the actual sequences of denoised ",
            "amplicon sequence variants, but some rownames contained invalid characters ",
            "(e.g. ", rownames(mtx)[which(!check.dna)[1]], ")."))
    }
    asvnames = paste0("Row", (1:nrow(mtx)))
    asvs = rownames(mtx)
    if (binary %in% c("vsearch")) {
      pz.message("Assuming aligner CAN do reverse complement by itself...")
    } else {
      pz.message("Assuming aligner CANNOT do reverse complement by itself...")
      asvnames = c(asvnames, asvnames)
      revcomp = function(x) seqinr::c2s(rev(seqinr::comp(seqinr::s2c(x))))
      asvs = c(asvs, vapply(asvs, revcomp, ""))
    }
    seqinr::write.fasta(as.list(asvs),
                        asvnames,
                        file.out=opts('vsearch_infile'),
                        nbchar=99999,
                        as.string=TRUE)
    return(TRUE)
}

#' Run vsearch analysis on a FASTA file of sequences.
#'
#' \code{run.vsearch} runs, by default, the optimal sequence aligner vsearch
#' (doi.org/10.5281/zenodo.806850) on a FASTA file, typically one generated by
#' \code{prepare.vsearch.input}. By changing the \code{vsearch_bin}
#' option, the aligner vsearch can also be used. If the provided binary name
#' is not"vsearch", \emph{phylogenize} will assume
#' an old version of vsearch is being called that doesn't have the ability to 
#' search both strands.
#'
#' Some particularly relevant global options are:
#' \describe{
#' input files).}
#'   \item{vsearch_infile}{String. File name of the sequences to be read into
#' vsearch.}
#'   \item{vsearch_outfile}{String. File name where vsearch writes output which is
#' then read back into \emph{phylogenize}.}
#'   \item{vsearch_dir}{String. Path where the binary of vsearch is found.}
#'   \item{vsearch_bin}{String. File name of the binary of vsearch.}
#'   \item{vsearch_cutoff}{Float. Between 0.95 and 1.00; percent identity minimum
#'   for alignment results.}
#'   \item{vsearch_16sfile}{String. Path to the 16S FASTA database that maps back
#' to MIDAS species.}
#'   \item{data_dir}{String. Path to directory containing the data files
#' required to perform a \emph{phylogenize} analysis. }
#' }
#'
#' @return Returns TRUE unless an error is thrown.
#' @export run.vsearch
run.vsearch <- function(...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    binary = basename(opts('vsearch_bin'))
    pid = opts('vsearch_cutoff')
    if (binary %in% c("vsearch")) {
      vsearch_args = c("--db",
                     file.path(opts('data_dir'),
                               opts('vsearch_16sfile')),
                     "--usearch_global",
                     opts('vsearch_infile'),
                     "--strand both",
		     "--id",
                     pid,
                     "--blast6out",
                     opts('vsearch_outfile'))
    } else if (binary == "vsearch") {
      vsearch_args = c("--usearch_global",
        opts('vsearch_infile'),
        "--db",
        file.path(opts('data_dir'), opts('vsearch_16sfile')),
        "--strand both",
        "--blast6out",
        opts('vsearch_outfile'),
        "--id",
        pid)
    } else {
      pz.warning(paste0("Aligner not recognized, calling as old version of ",
         "vsearch that does not support reverse complements"))
      vsearch_args = c("--db",
                     file.path(opts('data_dir'),
                               opts('vsearch_16sfile')),
                     "--usearch_global",
                     opts('vsearch_infile'),
                     "--id",
                     pid,
                     "--blast6out",
                     opts('vsearch_outfile'))

    }
    pz.message(paste0("Calling aligner ", binary, " with arguments: ",
                      paste(vsearch_args, sep=" ", collapse=" ")))
    r <- system2(file.path(opts('vsearch_dir'), opts('vsearch_bin')),
                 args = vsearch_args)
    if (r != 0) {
        pz.error(paste0("Aligner failed with error code ", r))
    }
    return(TRUE)
}

#' Read in results from vsearch.
#'
#' \code{get.vsearch.results} reads and parses the output of vsearch to get the
#' best-hit MIDAS species identifier for any 16S hit. Note that the reference
#' 16S FASTA database file must describe entries in the format: ">gene;;
#' species_or_genus_ID;;MIDAS_ID". Only MIDAS_ID is used so the contents of
#' "gene" and "species_or_genus_ID" can be arbitrary.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   input files).}
#'   \item{vsearch_outfile}{String. File name where vsearch writes output which is
#'   then read back into \emph{phylogenize}.}
#' }
#'
#' @return List containing a vector of hits, a vector of MIDAS ID targets, and a
#'     data frame of the assignments as they came out of vsearch.
#' @export
get.vsearch.results <- function(...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    # map to MIDAS IDs using vsearch
    assignments <- data.frame(
        data.table::fread(opts('vsearch_outfile')))
    row.hits <- as.numeric(gsub("Row", "", assignments[, 1]))
    row.targets <- sapply(assignments[, 2],
                          function(x) strsplit(x, ";;")[[1]][3])
    return(list(hits=row.hits, targets=row.targets, assn=assignments))
}

#' Sum non-unique rows after vsearch mapping.
#'
#' \code{sum.nonunique.vsearch} takes vsearch results and an abundance or
#' presence/absence matrix, drops any rows that mapped to multiple MIDAS IDs
#' (i.e. that couldn't confidently be assigned to a MIDAS species),
#' then sums any rows that mapped to the same MIDAS ID.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{out_dir}{String. Path to output directory. Default: "output"}
#'   input files).}
#'   \item{data_dir}{String. Path to directory containing the data files
#'   required to perform a \emph{phylogenize} analysis. Default: on package
#'   load, this default is set to the result of \code{system.file("extdata",
#'   package="phylogenize")}.}
#' }
#'
#' @param vsearch A list obtained by running \code{get.vsearch.results}.
#' @param mtx A presence/absence or abundance matrix, with row names equal to
#'     amplicon sequence variant DNA sequences.
#' @return A new matrix with MIDAS IDs as rows.
#' @export sum.nonunique.vsearch
sum.nonunique.vsearch <- function(vsearch, mtx, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    uniq.hits <- which(count.each(vsearch$hits) < 2)
    rh <- vsearch$hits[uniq.hits]
    rt <- vsearch$targets[uniq.hits]
    subset.abd <- mtx[rh, , drop=FALSE]
    urt <- unique(rt)
    summed.uniq <- sapply(urt, function(r) {
        w <- which(rt == r)
        apply(subset.abd[w, , drop=FALSE], 2, sum)
    }) %>% t
    rownames(summed.uniq) <- urt
    summed.uniq
}

#' Map denoised 16S sequence variants to MIDAS IDs using vsearch.
#'
#' \code{process.16s} is a wrapper for other functions that: output sequence
#' variants to a file; map them using vsearch against a reference database of 16S
#' sequences; then return a list of abundance and metadata values where the rows
#' of the abundance matrix are now MIDAS IDs.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   input files).}
#'   \item{vsearch_infile}{String. File name of the sequences written to disk and
#'   then read into vsearch.}
#' }
#'
#' @param mtx A presence/absence or abundance matrix, with row names equal to
#'     amplicon sequence variant DNA sequences.
#' @return none
#' @export
process.16s <- function(abd.meta, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    if (!(all(is.dna(rownames(abd.meta$mtx))))) {
        pz.error("expected rows to be DNA sequences but found illegal characters")
    }
    prepare.vsearch.input(abd.meta$mtx, ...)
    run.vsearch(...)
    vsearch <- get.vsearch.results(...)
    summed.uniq <- sum.nonunique.vsearch(vsearch, abd.meta$mtx, ...)
    csu <- colSums(summed.uniq)
    abd.meta$mtx <- apply(summed.uniq[, which(csu > 0), drop=FALSE],
                          2,
                          function(x) x / sum(x))
    abd.meta
}

create_matrix <- function(sub_df) {
  # Pivot the dataframe to wide format using pivot_wider
  pivoted <- sub_df %>%
    tidyr::pivot_wider(names_from = species, values_from = presence, values_fill = list(presence = 0))

  mat <- as.matrix(pivoted %>% select(-all_of("max_count")) %>% select(-all_of("function")) %>% select(-all_of("protein")))
  rownames(mat) <- pivoted$protein
  return(mat)
}

#--- Import data necessary for analyses ---#

#' Import the data necessary for *phylogenize* analysis.
#'
#' \code{import.pz.db} decides based on global options which data files to import.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{type_16S}{Boolean. If 16S data, TRUE, otherwise shotgun data is assumed. Default: FALSE}
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
    
    if (opts('db') == "gtdb") {
            # Read in gene presence and the functions file 
            gene.presence <- readRDS(file.path(opts('data_dir'), "gtdb-gene-presence-binary.rds"))
            # Read in the phylogenetic trees
            trees <- readRDS(file.path(opts('data_dir'), "gtdb_214-trees.rds"))
            # Add in the species column from the cluster column *This can be removed later on so that the external db is format fully
            taxonomy <- read_csv(file.path(opts('data_dir'),"gtdb_214-taxonomy.csv"))
	    # Read in the KO annotations file
            gene.to.fxn <- read_csv(file.path(opts('data_dir'), "gtdb.functions"))
    } else if (opts('db') == "uhgp") {
	    gene.presence <- readRDS(file.path(opts('data_dir'), "uhgp-gene-presence-binary.rds"))
	    trees <- readRDS(file.path(opts('data_dir'), "uhgp-trees.rds"))
	    taxonomy <- read_csv(file.path(opts('data_dir'),"uhgp-taxonomy.csv"))
	    gene.to.fxn <- read_csv(file.path(opts('data_dir'), "uhgp.functions"))
    } else {
            pz.error(paste0("Unknown data type ", opts('db')))
    }
    colnames(gene.to.fxn) <- c("gene", "function")

    # finished
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
    species.per.tree <- lapply(pz.db$trees, function(tr) {
	intersect(tr$tip.label, species.observed)
    })
    tL <- vapply(species.per.tree, length, 1L)
    if (all(tL < opts('treemin'))) {
        pz.error("Too few species found. Was the right database used?")
    }
    passed.min <- nw(tL >= opts('treemin'))
    totalL <- vapply(pz.db$trees, function(tr) { length(tr$tip.label) }, 1L)
    pct.obs <- mapply(function(x, y) x / y, tL, totalL)
    passed.pct <- nw(pct.obs >= opts('pctmin'))
    pz.message("Determining which taxa to test...")
    for (tn in 1:length(pz.db$trees)) {
        pz.message(paste0(names(pz.db$trees)[tn],
                          " (pct): ", format(pct.obs[tn] * 100, digits=2),
                          "%; (number): ", tL[tn],
                          "; ", ifelse((pct.obs[tn] >= opts('pctmin') &&
                                        tL[tn] >= opts('treemin')),
                                       yes="kept",
                                       no="dropped")
                          ))
    }
    saved.taxa <- intersect(passed.min, passed.pct)
    if (length(saved.taxa) == 0) {
        pz.error(paste0("All trees had less than ",
                        format(opts('pctmin') * 100, digits=2),
                        "% of species observed. Either a very small ASV ",
                        "table was provided, read depth was very shallow, ",
                        "the right database was not selected or very few ASVs ",
                        "mapped to entries in the database."))
    }
    pz.db$trees <- pz.db$trees[saved.taxa]
    pz.db$species <- lapply(pz.db$trees, function(x) x$tip.label)
    pz.db$ntaxa <- length(pz.db$trees)
    pz.db$trees <- lapply(pz.db$trees, fix.tree)
    pz.db
}

#--- Tools for processing data ---#

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
    abd.meta$mtx <- Matrix::Matrix(abd.meta$mtx > 0)
    abd.meta
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
        sum(phenotype[tr$tip.label] != phenoP)
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
    phenotype[intersect(names(phenotype), valid.names)]
}


#--- Plotting ---#

#' Get phenotype plotting scales.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence"
#'   or "specificity" or "correlation").}
#'   \item{prev_color_low}{String. When graphing prevalence on a tree, this
#'   color is the lowest value.} \item{prev_color_high}{String. When graphing
#'   prevalence on a tree, this color is the highest value.}
#'   \item{spec_color_high}{String. When graphing specificity on a tree, this
#'   color is the lowest value (most anti-specific).}
#'   \item{spec_color_med}{String. When graphing specificity on a tree, this
#'   color denotes the prior (no association).} \item{spec_color_high}{String.
#'   When graphing specificity on a tree, this color is the highest value (most
#'   specific).}
#' }
#' @param phenotype A named vector giving the phenotype for each taxon ID.
#' @param trees A list of tree objects.
#' @param phenoP An optional value giving the prior probability for the
#'     environment of interest.
#' @return A list of overall limits (\code{limits}), taxon-specific limits
#'     (\code{phy.limits}), a color scale (\code{colors}), and the zero point
#'     (\code{zero}).
#' @export
get.pheno.plotting.scales <- function(phenotype, trees, phenoP=0, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    if (opts('which_phenotype') == 'prevalence') {
        get.pheno.plotting.scales.prevalence(phenotype,
                                             trees,
                                             phenoP,
                                             ...)
    } else if (opts('which_phenotype') %in% c('specificity', 'abundance')) {
        get.pheno.plotting.scales.specificity(phenotype,
                                              trees,
                                              phenoP,
                                              ...)
    }
}


#' Get phenotype plotting scales (prevalence-specific).
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{prev_color_low}{String. When graphing prevalence on a tree, this
#'   color is the lowest value.}
#'   \item{prev_color_high}{String. When graphing prevalence on a tree, this
#'   color is the highest value.}
#' }
#' @param phenotype A named vector giving the phenotype for each taxon ID.
#' @param trees A list of tree objects.
#' @param phenoP An optional value giving the prior probability for the
#'     environment of interest.
#' @return A list of overall limits (\code{limits}), taxon-specific limits
#'     (\code{phy.limits}), a color scale (\code{colors}), and the zero point
#'     (\code{zero}).
#' @keywords internal
get.pheno.plotting.scales.prevalence <- function(phenotype,
                                                 trees,
                                                 phenoP=0,
                                                 ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    phenoLimits <- quantile(unique(phenotype), c(0.2, 0.8))
    phenoLimitsTaxon <- lapply(trees, function(tr) {
        phi <- phenotype[intersect(names(phenotype), tr$tip.label)] %>%
            na.omit
        quantile(unique(phi), c(0.2, 0.8))
    })
    phenoColors <- c(low.col=opts('prev_color_low'),
                     high.col=opts('prev_color_high'))
    return(list(limits=phenoLimits,
                phy.limits=phenoLimitsTaxon,
                colors=phenoColors,
                zero=0))
}


#' Get phenotype plotting scales (specificity or abundance).
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{spec_color_high}{String. When graphing specificity on a tree, this
#'   color is the lowest value (most anti-specific).}
#'   \item{spec_color_med}{String. When graphing specificity on a tree, this
#'   color denotes the prior (no association).} \item{spec_color_high}{String.
#'   When graphing specificity on a tree, this color is the highest value (most
#'   specific).}
#' }
#' @param phenotype A named vector giving the phenotype for each taxon ID.
#' @param trees A list of tree objects.
#' @param phenoP An optional value giving the prior probability for the
#'     environment of interest.
#' @return A list of overall limits (\code{limits}), taxon-specific limits
#'     (\code{phy.limits}), a color scale (\code{colors}), and the zero point
#'     (\code{zero}).
#' @keywords internal
get.pheno.plotting.scales.specificity <- function(phenotype,
                                                  trees,
                                                  phenoP=0,
                                                  ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    phenoLimits <- c(phenoP - 1 * sd(phenotype),
                     phenoP + 1 * sd(phenotype))
    phenoLimitsTaxon <- lapply(trees, function(tr) {
        phi <- phenotype[intersect(names(phenotype), tr$tip.label)] %>%
            na.omit %>%
            unique
        c(phenoP - (1 * sd(phi)), phenoP + (1 * sd(phi)))
    })
    phenoColors <- c(low.col=opts('spec_color_low'),
                     mid.col=opts('spec_color_mid'),
                     high.col=opts('spec_color_high'))
    return(list(limits=phenoLimits,
                phy.limits=phenoLimitsTaxon,
                colors=phenoColors,
                zero=phenoP))
}

#' Plot a phenotype along a list of trees.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence"
#'   or "specificity").}
#' }
#'
#' @param phenotype A named vector with the phenotype values for each taxon.
#' @param trees A list of trees.
#' @param scale A list returned from \code{get.pheno.plotting.scales}.
#' @return A list of ggtree objects in which the phenotype has been plotted
#'     across each tree in \code{trees}.
#' @export plot.phenotype.trees
plot.phenotype.trees <- function(phenotype,
                                 trees,
                                 scale,
                                 ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    if (any(!(names(trees) %in% names(scale$phy.limits)))) {
        pz.error("taxon-specific limits must be calculated for every tree")
    }
    plotted.pheno.trees <- lapply(names(trees), function(tn) {
        tryCatch({gg.cont.tree(trees[[tn]],
                               phenotype,
                               cLimits=scale$phy.limits[[tn]],
                               colors=scale$colors,
                               cName=paste0(tn,
                                            ": ",
                                            opts('which_phenotype')),
                               plot=FALSE)},
                 error=function(e) {
                     pz.warning(e)
                     NA
                 })
    })
    names(plotted.pheno.trees) <- names(trees)
    plotted.pheno.trees <- plotted.pheno.trees[vapply(plotted.pheno.trees, is.list, TRUE)]
    if (length(plotted.pheno.trees) == 0) {
        pz.warning("No trees were plotted: too few taxa for any taxon?")
    }
    plotted.pheno.trees
}

#' Plot distributions of a phenotype across taxa.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence"
#'   or "specificity").}
#' }
#' @param phenotype A named vector with the phenotype values for each species.
#' @param pz.db A database containing a \code{taxonomy} and \code{trees}.
#' @return A ggplot object with the phenotype distribution plotted per taxon.
#' @export plot.pheno.distributions
plot.pheno.distributions <- function(phenotype,
                                     pz.db,
                                     ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    kept.species <- Reduce(c, lapply(pz.db$trees, function(x) x$tip.label))
    pheno.taxon <- pz.db$taxonomy[match(names(phenotype),
                                         pz.db$taxonomy$cluster),
                                   "phylum",
                                   drop=TRUE]
    pheno.characteristics <- data.frame(pheno=phenotype,
                                        taxon=pheno.taxon,
                                        cluster=names(phenotype))
    sub.pheno <- subset(pheno.characteristics, cluster %in% kept.species)
    distros <- ggplot2::ggplot(sub.pheno,
                               ggplot2::aes(pheno,
                                            color = taxon,
                                            fill = taxon)) +
        ggplot2::geom_density() +
        ggplot2::xlab(opts('which_phenotype')) +
        ggplot2::ggtitle(paste0("Distributions of phenotype (",
                                opts('which_phenotype'), ")"))
    if (length(unique(sub.pheno$taxon)) > 1) { # don't assume >1 taxon
        distros <- distros + ggplot2::facet_grid(taxon ~ .)
    }
    return(distros)
}


#' Edit a list of plotted trees to add fancy highlight labels.
#'
#' This function adds fancy SVG highlight labels to ggtree objects and then
#' plots them. If there's an error, it will fall back to a regular plot.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence"
#'   or "specificity").}
#' }
#'
#' @param plotted.pheno.trees A named list of ggtree plots (per taxon).
#' @param phenotype Phenotype to plot/label.
#' @param label Label to give to the phenotype.
#' @param stroke.scale How thick to make the highlight.
#' @param units A string appended to each label, used to give units of phenotype.
#' @export plot.labeled.phenotype.trees
plot.labeled.phenotype.trees <- function(plotted.pheno.trees,
                                         phenotype,
                                         label='prevalence',
                                         stroke.scale=0.3,
                                         units='%') {
    if (is.null(plotted.pheno.trees)) {
        pz.message("warning: no trees found")
        return(NULL)
    }
    if (length(plotted.pheno.trees) == 0) {
        pz.message("warning: no trees found")
        return(NULL)
    }

    for (pn in 1:length(plotted.pheno.trees)) {
	    p <- plotted.pheno.trees[[pn]]
	    rp <- p$rphy
	    tr <- p$tree
	    rp2 <- rp
	    rp2$tip.label <- paste0(rp2$tip.label, " (phenotype = ...", units, ")")
	    xlim <- plot(rp2, plot=FALSE, no.margin=TRUE)$x.lim
	    new.tr <- tr +
		    ggtree::geom_tiplab() +
		    ggplot2::xlim(xlim[1], xlim[2])
	    fn <- knitr::fig_path('svg', number = pn)

	    #new.tr + geom_tiplab() + xlim(NA, 5)
            # Make the labels as species not the midas id


	    non.interactive.plot(new.tr, fn)
    }
    #tree_list <- lapply(plotted.pheno.trees, function(tree) {
    #    pz.error(tree)        
    #})
    
    

    #interactive_plot_list <- lapply(plotted.pheno.trees, plot.interactive.tree, phenotype)
    #for (i in seq_along(interactive_plot_list)) {
    #    plot.tree.with.phenotype(tree_list[[i]], phenotype_data, paste0("tree_plot_", i, ".svg"))
    #}
}


#' Make a hybrid tree-heatmap plot showing the taxon distribution of significant
#' hits.
#'
#' This function wraps single.cluster.plot, running it in a separate process.
#' This is because there can be problems with memory leaks. For relevant global
#' options, see the doumentation for that function.
#'
#' @param gene.presence Gene presence/absence matrix.
#' @param sig.genes Character vector of the significant genes.
#' @param tree A tree object.
#' @param plotted.tree A ggtree plot of \code{tree}.
#' @param taxon Name of the taxon represented by \code{tree}
#' @param verbose Whether to report debugging information (boolean).
#' @export do.clust.plot
do.clust.plot <- function(gene.presence,
                          sig.genes,
                          tree,
                          plotted.tree,
                          taxon,
                          verbose=FALSE,
                          ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    # Run these on a separate process to avoid memory leak
    cl <- parallel::makeCluster(1)
    if (verbose) message("exporting data...")
    parallel::clusterExport(cl,
                  c("gene.presence",
                    "sig.genes",
                    "tree",
                    "plotted.tree",
                    "taxon",
                    "verbose",
                    "PZ_OPTIONS"),
                  envir=environment())
    if (verbose) message("importing source...")
    # parallel::clusterCall(cl, library, phylogenize)
    cluster.load.pkg(cl, opts("devel"), opts("devel_pkgdir"))
    if (verbose) message("performing call...")
    tmpL <- parallel::clusterCall(cl,
                        single.cluster.plot,
                        gene.presence,
                        sig.genes,
                        tree,
                        plotted.tree,
                        taxon,
                        verbose,
                        ...)
    print(tmpL[[1]])
    rm(tmpL)
    parallel::stopCluster(cl)
    gc()
    return(NULL) # Avoid wasting memory since we never touch these
}

#' Make a hybrid tree-heatmap plot showing the taxon distribution of significant
#' hits.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence"
#'   or "specificity").}
#'   \item{gene_color_absent}{String. When graphing gene presence/absence, this color indicates absence.}
#'   \item{gene_color_present}{String. When graphing gene presence/absence, this color indicates presence.}
#' }
#'
#' @param gene.presence Gene presence/absence matrix.
#' @param sig.genes Character vector of the significant genes.
#' @param tree A tree object.
#' @param plotted.tree A ggtree plot of \code{tree}.
#' @param taxon Name of the taxon represented by \code{tree}
#' @param verbose Whether to report debugging information (boolean).
#' @return A faceted ggplot object.
#' @export
single.cluster.plot <- function(gene.presence,
                                sig.genes,
                                tree,
                                plotted.tree,
                                taxon,
                                verbose=FALSE,
                                ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    sig.bin <- gene.presence[intersect(rownames(gene.presence),
                                       sig.genes), , drop=FALSE]
    if (is.null(dim(sig.bin))) {
        sig.bin <- as.matrix(sig.bin) %>% t
        rownames(sig.bin) <- sig.genes
    }
    names(dimnames(sig.bin)) <- c("gene", "id")
    col.df <- data.frame(id=names(plotted.tree$disp),
                         disp=plotted.tree$disp)
    p <- ggtree::ggtree(tree, ladderize = TRUE) %<+%
        col.df +
        ggtree::geom_tippoint(ggplot2::aes(color=disp, shape='.'))
    if (opts('which_phenotype') == "prevalence") {
        p <- p + ggplot2::scale_color_gradient(
                              low=plotted.tree$cols["low.col"],
                              high=plotted.tree$cols["high.col"])
    } else if (opts('which_phenotype') %in% c("specificity", "correlation")) {
        p <- p + ggplot2::scale_color_gradient2(
                              low=plotted.tree$cols["low.col"],
                              mid=plotted.tree$cols["mid.col"],
                              high=plotted.tree$cols["high.col"],
                              midpoint=mean(plotted.tree$lims))
    }
    if (length(sig.genes) == 0) { return(p) }
    if (length(sig.genes) > 1) {
        clust <- hclust(dist(sig.bin, method = "binary"))
        sig.ord <- sparseMelt(t(sig.bin)[, clust$order, drop=FALSE])
    } else {
        sig.ord <- sparseMelt(t(sig.bin))
        sig.ord <- sig.ord[order(sig.ord[, 3]), , drop=FALSE]
    }
    tmp <- ggtree::facet_plot(p,
                              panel=paste0('heatmap: ', taxon),
                              data=sig.ord,
                              geom=ggplot2::geom_tile,
                              mapping=ggplot2::aes(x = as.numeric(as.factor(gene)),
                                                   fill = as.numeric(as.factor(value)))) +
        ggplot2::scale_fill_gradient(low = opts('gene_color_absent'),
                                     high = opts('gene_color_present'),
                                     na.value = opts('gene_color_absent'))
    tmp
}

#--- Report generation ---#

#' Run *phylogenize* start to finish.
#'
#' @param output_file Path giving what to name the resulting HTML file.
#' @param report_input Optionally override which notebook to knit (useful for
#'     testing).
#' @param do_cache Turn on or off Rmarkdown's caching.
#' @param reset_after Reset global options to package defaults after running?
#' @param ... Parameters to override defaults.
#' @export
render.report <- function(output_file='report_output.html',
                          report_input='phylogenize-report.Rmd',
                          do_cache=TRUE,
                          reset_after=TRUE,
                          ...) {
    prev.options <- pz.options()
    do.call(pz.options, list(...))
    pz.options(working_dir=normalizePath(getwd()))
    pz.options(out_dir=normalizePath(pz.options("out_dir")))
    if (!dir.exists(pz.options('out_dir'))) {dir.create(pz.options('out_dir'))}
    p <- pz.options()
    for (n in names(p)) {
        message(paste0(n, ": ", p[[n]]))
    }
    file.copy(system.file("rmd",
                          report_input,
                          package="phylogenize"),
              pz.options("out_dir"),
              overwrite=TRUE)
    r <- rmarkdown::render(file.path(pz.options("out_dir"),
                                     report_input),
                           output_file=basename(output_file),
                           output_dir=pz.options("out_dir"),
                           output_options=list(
                               cache.path=pz.options("out_dir"),
                               cache=do_cache
                           ),
                           intermediates_dir=pz.options("out_dir"),
                           knit_root_dir=pz.options("out_dir"))
    if (reset_after) {
        do.call(pz.options, prev.options)
        set_data_internal()
    }
    r
}

#' Make a pretty enrichment table.
#'
#' @param enr.table Input enrichment table.
#' @return Mutated enrichment table with better-labeled columns and significance
#'     coloring.
#' @export
output.enr.table <- function(enr.table) {
    enr.table %>%
        dplyr::mutate(
            Gene_significance=capwords(Gene_significance),
            Taxon=capwords(Taxon)
        ) %>%
        dplyr::mutate(
              q_value=kableExtra::cell_spec(
                prettyNum(q_value, digits=2),
                "html",
                background=kable.recolor(
                    -log10(q_value),
                    colors=c("#FFFFFF", "#FF8888"),
                    limits=c(0, 10)))) %>%
        dplyr::mutate(
              Enrichment_odds_ratio=kableExtra::cell_spec(
                prettyNum(Enrichment_odds_ratio, digits=3),
                "html",
                background=kable.recolor(
                    Enrichment_odds_ratio,
                    colors=c("#FFFFFF", "#FFFF44"),
                    limits=c(1, 10)))) %>%
        knitr::kable("html", escape=FALSE, align="l") %>%
        kableExtra::kable_styling(c("striped", "condensed"))
}


#' Check if a particular string is likely to be DNA.
#'
#' @param seq String to check for illegal characters.
#' @return TRUE if it contains no illegal characters, FALSE otherwise.
#' @keywords internal
is.dna <- function(seq) {
    !(grepl("[^actguwsmkrybdhvn]", tolower(seq)))
}

#' Download data from figshare (or provide it locally) and un-gzip it into the
#' package directory so that it can be imported.
#'
#' @param data_path Optional: provide a path to a local file containing a
#'     compressed .tar archive of data. Must extract to the subdirectory
#'     \code{extdata/}.
#' @param figshare_url Optional: override the URL from which to obtain the data.
#' @export
install.data.figshare <- function(data_path=NULL,
                                  figshare_url="https://ndownloader.figshare.com/files/43692576?private_link=987aeecdfebd2da02302") {
    if (is.null(data_path)) {
        data_path = tempfile()
        curl::curl_download(figshare_url, data_path)
    }
    print(system.file("", package="phylogenize"))
    untar(data_path, exdir = system.file("", package="phylogenize"))
}

#' Return a boolean telling whether a phenotype has nonzero variance in different taxa.
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
