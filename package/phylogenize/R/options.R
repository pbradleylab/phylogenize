# Set the default parameters
default_params <- list(
    abundance_file = "",
    assume_below_LOD = TRUE,
    biom_file = "",
    db = "gtdb",
    data_dir = "",
    devel = FALSE,
    devel_pkgdir = 'package/phylogenize',
    dset_column = "dataset",
    env_column = "env",
    taxon_level = "phylum",
    error_to_file = TRUE,
    gene_color_absent = 'black',
    gene_color_present = 'slateblue2',
    input_format = "tabular",
    metadata_file = "",
    meas_err = TRUE,
    minimum = 3,
    min_fx = 0,
    ncl = 1,
    output_file="results.html",
    out_dir="./",
    pctmin = 0.025,
    phenotype_file = "",
    prior_type = "uninformative",
    prior_file = "",
    prev_color_high = 'orange2',
    prev_color_low = 'black',
    relative_out_dir = NULL,
    sample_column = "sample",
    separate_process = TRUE,
    single_dset = FALSE,
    skip_graphs = FALSE,
    spec_color_high = 'tomato',
    spec_color_low = 'slateblue',
    spec_color_mid = 'gray50',
    treemin = 5,
    type_16S = FALSE,
    tax_level = "family",
    use_rmd_params = FALSE,
    vsearch_16sfile = "16s_gtdb.frn",
    vsearch_cutoff = 0.985,
    vsearch_dir = "",
    vsearch_infile = "input_seqs.txt",
    vsearch_outfile = "output_assignments.txt",
    which_envir = "stool",
    which_phenotype = "prevalence",
    phenotype_file = "phenotype.tsv",
    categorical = TRUE,
    diff_abund_method = "maaslin2",
    working_dir = '.',
    core_method = "phylogenize",
    rds_output_file = "core_output.rds"
)


# Options given by user
PZ_OPTIONS <- options_manager(.list=default_params)

#' Set and get options for phylogenize.
#'
#' Function to set and get global options for the \emph{phylogenize} package.
#'
#' These options are global because they affect how most of the functions in
#' \emph{phylogenize} work. Descriptions of these options follow.
#'
#' @section File input/output and paths:
#' \describe{
#'   \item{abundance_file}{String. Name of abundance tabular file. Default: "test-abundance.tab"}
#'   \item{biom_file}{String. Name of BIOM abundance-and-metadata file. Default: "test.biom"}
#'   \item{data_dir}{String. Path to directory containing the data files required to perform a \emph{phylogenize} analysis. Default: empty string, but on package load, this default is set to the result of \code{system.file("extdata", package="phylogenize")}.}
#'   \item{error_to_file}{Boolean. Should pz.error, pz.warning, and pz.message output to an error message file? Default: FALSE}
#'   \item{input_format}{String. Whether to look for tabular or BIOM-formatted data ("tabular" or "biom"). Default: "tabular"}
#'   \item{metadata_file}{String. Name of metadata tabular file. Default: "test-metadata.tab"}
#'   \item{rds_output_file}{String. Name of output RDS file containing the full results of applying `phylogenize_core()`. Set to empty string to disable. Default: "core_output.rds"}
#'   \item{output_file}{String. Name of output file: "results.html"}
#'   \item{phenotype_file}{String. Name of input file for optional pre-calculated phenotype. Default: ""}
#'   \item{prior_file}{String. File name of optional pre-computed prior. Default: ""}
#'   \item{separate_metadata}{Boolean. For BIOM data, is there a separate tabular abundance table? Default: FALSE}
#'   \item{vsearch_16sfile}{String. Path to the 16S FASTA database that maps back to MIDAS species. Default: "16s_gtdb.frn"}
#'   \item{vsearch_dir}{String. Path where the binary of the aligner is found. Default: "/usr/local/bin/"}
#'   \item{vsearch_infile}{String. File name of the sequences written to disk and then read into the aligner. Default: "input_seqs.txt"}
#'   \item{vsearch_outfile}{String. File name where the aligner writes output which is then read back into \emph{phylogenize}. Default: "output_assignments.txt"}
#' }
#'
#' @section Computing phenotypes and results:
#' \describe{
#'   \item{assume_below_LOD}{Boolean. If TRUE, MIDAS species that are not present are assumed to have a prevalence of zero; if FALSE, they are dropped from the analysis. Default: TRUE}
#'   \item{db}{String. Which database to use. Can be "gtdb" or "uhgp." Default: "gtdb"}
#'   \item{dset_column}{String. Name of column in metadata file containing the dataset annotations. Default: "dataset"}
#'   \item{env_column}{String. Can either be set to 'phylum', 'class', 'order', 'family', or 'genus'. Default: "phylum"}
#'   \item{taxon_level}{String. Can either be set to 'phylum', 'class', 'order', 'family', or 'genus'. Default: "phylum"}
#'   \item{linearize}{Boolean. If TRUE, use a regular linear model instead of a phylogenetic linear model. Mostly useful for testing report generation, since the linear model is much faster but returns many more false positives. Default: FALSE}
#'   \item{meas_err}{Boolean. Separately estimate measurement error from phenotype variation in the phylogenetic linear model. Default: TRUE}
#'   \item{min_fx}{Positive double. Effects that are significantly equivalent to this effect size will be excluded from significant positive hits. If zero, the equivalence test will be skipped. Default: 0}
#'   \item{minimum}{Integer. A particular gene must be observed, and also absent, at least this many times to be reported as a significant positive association with the phenotype. Default: 3}
#'   \item{ncl}{Integer. Number of cores to use for parallel computation. Default: 1}
#'   \item{pctmin}{Float. A taxon must have at least this percent of observed representatives in order to be processed. Default: 0.01}
#'   \item{prior_type}{String. What type of prior to use ("uninformative" or "file"). Default: "uninformative"}
#'   \item{single_dset}{Boolean. If true, will assume that all samples come from a single dataset called \code{dset1} no matter what, if anything, is in \code{dset_column}. Default: FALSE}
#'   \item{treemin}{Integer. A taxon must have at least this many representatives in order to be processed. Default: 5}
#'   \item{type_16S}{Boolean. If 16S data, TRUE, otherwise shotgun data is assumed. Default: FALSE}
#'   \item{vsearch_cutoff}{Float. Value between 0.95 and 1.00 giving the percent ID cutoff to use when assigning denoised sequence variants to MIDAS species using vsearch. Default: 0.985}
#'   \item{which_envir}{String. Environment for which prevalence, specificity, or differential abundance scores will be the phenotype of interest. Must match annotations in metadata. Default: "Stool"}
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence", "specificity", "abundance", "provided"). Default: "prevalence"}
#'   \item{phenotype_file}{String. If phenotype is provided, what is the path to the file? Default: "phenotype.tsv"}
#'   \item{categorical}{Boolean. For abundance estimates, is the environment in env_column a categorical variable (TRUE) or continuous (FALSE)? Default: TRUE}
#'   \item{diff_abund_method}{String. Which method to use to calculate differential abundance. Either "ANCOMBC2" or "Maaslin2" (case insensitive). Default: "Maaslin2"}
#'   \item{core_method}{String. Which method to use to associate genes with phenotypes. Either "phylogenize" or "POMS" (case insensitive). Default: "phylogenize"}
#' }
#'
#' @section Graphing:
#' \describe{
#'   \item{gene_color_absent}{String. When graphing gene presence/absence, this color indicates absence. Default: "black"}
#'   \item{gene_color_present}{String. When graphing gene presence/absence, this color indicates presence. Default: "black"}
#'   \item{prev_color_high}{String. When graphing prevalence on a tree, this color is the highest value. Default: "orange2"}
#'   \item{prev_color_low}{String. When graphing prevalence on a tree, this color is the lowest value. Default: "black"}
#'   \item{skip_graphs}{Boolean. If TRUE, skip making graphs in the report, which can be time- and memory-consuming. Default: FALSE}
#'   \item{spec_color_high}{String. When graphing specificity on a tree, this color is the highest value (most specific). Default: "tomato"}
#'   \item{spec_color_med}{String. When graphing specificity on a tree, this color denotes the prior (no association). Default: "gray50"}
#'   \item{spec_color_low}{String. When graphing specificity on a tree, this color is the lowest value (most anti-specific). Default: "slateblue"}
#' }
#'
#' @section Memory management:
#' \describe{
#'   \item{separate_process}{Boolean. When displaying clustered top gene associations alongside a tree colored by phenotype, this flag indicates whether to use a separate subprocess. This allows memory used by clustering to be released back to the operating system immediately. Default: TRUE}
#' }
#'
#' @param ... Names of options (to retrieve) or \code{[key]=[value]} pairs (to set).
#'
#' @export
pz.options <- function(...) {
    settings::stop_if_reserved(...)
    PZ_OPTIONS(...)
}

#' Set data directory to internal
#'
#' @param fail Boolean. If TRUE, set_data_internal will not attempt to download
#'     and install data from Figshare if it is missing.
#' @param startup Boolean. Is this function being called by .onLoad?
#' @export
set_data_internal <- function(fail=FALSE, startup=FALSE) {
    if (startup) {
        M <- packageStartupMessage
    } else {
        M <- message
    }
    dd <- system.file("extdata", package="phylogenize")
    instdir <- system.file("", package="phylogenize")
    success <- FALSE
    if (grepl("00LOCK-", instdir)) {
        success <- TRUE
        M("Skipping check for data during staged install")
    } else {
        if (dd == "") {
            if (!fail) {
                M("Installing data from Figshare...")
                success <- tryCatch(install.data.figshare(), error = function(e) {
                    warning(paste("Installing from Figshare failed: ", e))
                    return(FALSE)
                })
            }
            if (!success) {
                warning(paste("Data not found; *phylogenize* will not run",
                              "properly. Please try",
                              "phylogenize::install.data.figshare()",
                              "later or install the data manually into",
                              file.path(system.file("", package="phylogenize"),
                                        "extdata"),
                              "."))
            }
        } else {
            success <- TRUE
        }
    }
    if (success) pz.options(data_dir=dd)
}

#' Test whether data is installed and warn user if not.
#'
#' @param startup Boolean. Is this function being called by .onLoad?
#' @export
check_data_found <- function(fail=FALSE, startup=FALSE) {
    if (startup) {
        M <- packageStartupMessage
    } else {
        M <- message
    }
    
    dd <- system.file("extdata", package="phylogenize")
    instdir <- system.file("", package="phylogenize")
    success <- FALSE
    if (dd == "") {
        M(paste("Data not found; *phylogenize* will not run",
                "properly. Please try",
                "phylogenize::install.data.figshare()",
                "later or install the data manually into",
                file.path(system.file("", package="phylogenize"),
                          "extdata"),
                "."))
    } else {
        success <- TRUE
    }
    if (success && pz.options('data_dir')=="") pz.options(data_dir=dd)
}

.onLoad <- function(libname, pkgname) {
    check_data_found(startup=TRUE)
}

