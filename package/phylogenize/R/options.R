# Set the default parameters
default_params <- list(
    dset_column = "dataset",
    env_column = "env",
    diff_abund_method = "ANCOMBC2",
    core_method = "permutrate-rlm",
    verbosity = 1,
    error_to_file = TRUE,
    error_file = "errmsg.txt",
    abundance_file = "",
    assume_below_LOD = TRUE,
    biom_file = "",
    db = "human-gut",
    working_dir = '.',
    data_dir = "",
    devel = FALSE,
    devel_pkgdir = 'package/phylogenize',
    gene_color_absent = 'white',
    gene_color_present = 'coral',
    input_format = "tabular",
    metadata_file = "",
    meas_err = TRUE,
    minimum = 3,
    gene_min_frac = 0.5,
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
    separate_metadata = FALSE,
    separate_process = TRUE,
    single_dset = FALSE,
    skip_graphs = FALSE,
    spec_color_high = 'tomato',
    spec_color_low = 'slateblue',
    spec_color_mid = 'gray80',
    treemin = 10,
    type_16S = FALSE,
    taxon_level = "family",
    use_rmd_params = FALSE,
    which_envir = "stool",
    which_phenotype = "prevalence",
    phenotype_file = "phenotype.tsv",
    categorical = TRUE,
    rds_output_file = "core_output.rds",
    fdr_method = "BH",
    quantile_normalize = FALSE,
    only_specific_taxa = NULL,
    # 16S stuff (not officially supported yet)
    aln_path_16s="",
    tree_path_16s="",
    jplace_file="",
    min_frac_16s=0.8,		       
    appspam_path="/usr/local/bin/appspam",
    which_16s_method="appspam",
    vsearch_16sfile = "16s_gtdb.frn",
    vsearch_cutoff = 0.985,
    vsearch_dir = "",
    named_asv_file = "input_seqs.txt",
    vsearch_outfile = "output_assignments.txt"
)


# Options given by user
PZ_OPTIONS <- settings::options_manager(.list=default_params)

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
#'   \item{error_file}{String. What should the filename of the error message file be? Default: "errmsg.txt"}
#'   \item{verbosity}{Number. What level of pz.message (out of 3) should be reported? Default: 1}
#'   \item{input_format}{String. Whether to look for tabular or BIOM-formatted data ("tabular" or "biom"). Default: "tabular"}
#'   \item{metadata_file}{String. Name of metadata tabular file. Default: "test-metadata.tab"}
#'   \item{rds_output_file}{String. Name of output RDS file containing the full results of applying `phylogenize_core()`. Set to empty string to disable. Default: "core_output.rds"}
#'   \item{output_file}{String. Name of output file: "results.html"}
#'   \item{phenotype_file}{String. Name of input file for optional pre-calculated phenotype. Default: ""}
#'   \item{prior_file}{String. File name of optional pre-computed prior. Default: ""}
#'   \item{separate_metadata}{Boolean. For BIOM data, is there a separate tabular abundance table? Default: FALSE}
#'   \item{which_16s_method}{String. Can be "vsearch" (best-hit alignment), "appspam" (perform phylogenetic placement), or "jplace" (bring-your-own .jplace file). Default: "appspam"}
#'   \item{vsearch_16sfile}{String. Path to the 16S FASTA database that maps back to MIDAS species. Default: "16s_gtdb.frn"}
#'   \item{vsearch_dir}{String. Path where the binary of the aligner is found. Default: "/usr/local/bin/"}
#'   \item{named_asv_file}{String. Path where sequences will be written to disk and then read into the aligner/AppSpam. Default: "input_seqs.txt"}
#'   \item{vsearch_outfile}{String. File name where the aligner writes output which is then read back into \emph{phylogenize}. Default: "output_assignments.txt"}
#'   \item{jplace_file}{String. Path to write .jplace file (if which_16s_method is "appspam") or to read user-provided .jplace file (if which_16s_method is "jplace".) \emph{phylogenize}. Default: ""}
#'   \item{aln_path_16s}{String. Path to the multiple alignment of 16S sequences used for phylogenetic placement. Default: ""}
#'   \item{tree_path_16s}{String. Path to the tree of 16S sequences used for phylogenetic placement. \emph{phylogenize}. Default: ""}
#'   \item{min_frac_16s}{Numeric. Should be between 0.5 and 1. Only keep ASVs where at least this fraction of assignments are to the same species. Allows some tolerance for mislabeled or phylogenetically-inconsistent 16S sequences in the database. Default: 0.8}
#' }
#'
#' @section Computing phenotypes and results:
#' \describe{
#'   \item{assume_below_LOD}{Boolean. If TRUE, MIDAS species that are not present are assumed to have a prevalence of zero; if FALSE, they are dropped from the analysis. Default: TRUE}
#'   \item{quantile_normalize}{Boolean. If TRUE, all phenotypes will be quantile-normalized to the normal distribution. Default: FALSE}
#'   \item{db}{String. Which database to use. Can be "gtdb" or "uhgp." Default: "gtdb"}
#'   \item{dset_column}{String. Name of column in metadata file containing the dataset annotations. Default: "dataset"}
#'   \item{env_column}{String. Name of column in metadata file containing the environment annotations. Default: "env"}
#'   \item{taxon_level}{String. Can either be set to 'phylum', 'class', 'order', 'family', or 'genus'. Default: "phylum"}
#'   \item{linearize}{Boolean. If TRUE, use a regular linear model instead of a phylogenetic linear model. Mostly useful for testing report generation, since the linear model is much faster but returns many more false positives. Default: FALSE}
#'   \item{meas_err}{Boolean. Separately estimate measurement error from phenotype variation in the phylogenetic linear model. Default: TRUE}
#'   \item{min_fx}{Positive double. Effects that are significantly equivalent to this effect size will be excluded from significant positive hits. If zero, the equivalence test will be skipped. Default: 0}
#'   \item{minimum}{Integer. A particular gene must be observed, and also absent, at least this many times to be reported as a significant positive association with the phenotype. Default: 3}
#'   \item{ncl}{Integer. Number of cores to use for parallel computation. Default: 1}
#'   \item{pctmin}{Float. A taxon must have at least this percent of observed representatives in order to be processed. Default: 0.025}
#'   \item{prior_type}{String. What type of prior to use ("uninformative" or "file"). Default: "uninformative"}
#'   \item{single_dset}{Boolean. If true, will assume that all samples come from a single dataset called \code{dset1} no matter what, if anything, is in \code{dset_column}. Default: FALSE}
#'   \item{treemin}{Integer. A taxon must have at least this many representatives in order to be processed. Default: 10}
#'   \item{type_16S}{Boolean. If 16S data, TRUE, otherwise shotgun data is assumed. Default: FALSE}
#'   \item{vsearch_cutoff}{Float. Value between 0.95 and 1.00 giving the percent ID cutoff to use when assigning denoised sequence variants to MIDAS species using vsearch. Default: 0.985}
#'   \item{which_envir}{String. Environment for which prevalence, specificity, or differential abundance scores will be the phenotype of interest. Must match annotations in metadata. Default: "Stool"}
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence", "specificity", "abundance", "provided"). Default: "prevalence"}
#'   \item{phenotype_file}{String. If phenotype is provided, what is the path to the file? Default: "phenotype.tsv"}
#'   \item{categorical}{Boolean. For abundance estimates, is the environment in env_column a categorical variable (TRUE) or continuous (FALSE)? Default: TRUE}
#'   \item{diff_abund_method}{String. Which method to use to calculate differential abundance. Either "ANCOMBC2" or "Maaslin2" (case insensitive). Default: "ANCOMBC2"}
#'   \item{core_method}{String. Which method to use to associate genes with phenotypes. Either "permutrate-rlm", "permutrate-lm", "permulate-lm", "permulate-rlm", "phylolm", "lm", or "POMS" (case insensitive). Default: "permutrate-rlm"}
#'   \item{fdr_method}{String. Which method to correct FDR for significant results? Either "BH", "BY", or "qvalue". Default: "qvalue"}
#' }
#'
#' @section Graphing:
#' \describe{
#'   \item{gene_color_absent}{String. When graphing gene presence/absence, this color indicates absence. Default: "white"}
#'   \item{gene_color_present}{String. When graphing gene presence/absence, this color indicates presence. Default: "coral"}
#'   \item{prev_color_high}{String. When graphing prevalence on a tree, this color is the highest value. Default: "orange2"}
#'   \item{prev_color_low}{String. When graphing prevalence on a tree, this color is the lowest value. Default: "black"}
#'   \item{skip_graphs}{Boolean. If TRUE, skip making graphs in the report, which can be time- and memory-consuming. Default: FALSE}
#'   \item{spec_color_high}{String. When graphing specificity or abundance on a tree, this color is the highest value (most specific). Default: "tomato"}
#'   \item{spec_color_med}{String. When graphing specificity or abundance on a tree, this color denotes the prior (no association). Default: "gray50"}
#'   \item{spec_color_low}{String. When graphing specificity or abundance on a tree, this color is the lowest value (most anti-specific). Default: "slateblue"}
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
    dd <- system.file("extdata/databases.csv", package="phylogenize")
    phd <- system.file("", package="phylogenize")
    success <- FALSE
    #if (grepl("00LOCK-", instdir)) {
    #    success <- TRUE
    #    M("Skipping check for data during staged install")
    #} else {
        if (dd == "") {
            M(sprintf("Note: databases.csv was not found under directory '%s'. You will need to manually set the directory later with the option 'data_dir=<PATH>'.", phd))
        } else {
            success <- TRUE
            db <- readr::read_delim(dd)
            M(sprintf("Databases listed:\n\t - %s", paste0(db[["database"]], collapse="\n\t - ")))
        }
    #}
    if (success && pz.options('data_dir') == "") {
        pz.options(data_dir = dd)
    }
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
    dd <- system.file("extdata/databases.csv", package = "phylogenize")
    phd <- system.file("", package = "phylogenize")
    success <- FALSE
    #if (grepl("00LOCK-", instdir)) {
    #    success <- TRUE
    #    M("Skipping check for data during staged install")
    #} else {
        if (dd == "") {
            M(sprintf(
                "Note: databases.csv was not found under directory '%s'. You will need to manually set the directory later with the option 'data_dir=<PATH>'.",
                phd
            ))
        } else {
            success <- TRUE
            db <- readr::read_delim(dd)
            M(sprintf(
                "Databases listed:\n\t - %s",
                paste0(db[["database"]], collapse = "\n\t - ")
            ))
        }
    #}
    if (success && pz.options('data_dir') == "") {
        pz.options(data_dir = dd)
    }
}

.onLoad <- function(libname, pkgname) {
    check_data_found(startup=TRUE)
}

