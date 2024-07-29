# Set the default parameters
default_params <- list(
    abundance_file = "",
    assume_below_LOD = TRUE,
    biom_file = "",
    db = "gtdb",
    data_dir = "./data",
    devel = FALSE,
    devel_pkgdir = 'package/phylogenize',
    dset_column = "dataset",
    env_column = "env",
    error_to_file = TRUE,
    gene_color_absent = 'black',
    gene_color_present = 'slateblue2',
    input_format = "tabular",
    linearize = FALSE,
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
    working_dir = '.'
)


# Options given by user
PZ_OPTIONS <- options_manager(
    abundance_file = default_params["abundance_file"][[1]],
    assume_below_LOD = default_params["assume_below_LOD"][[1]],
    biom_file = default_params["biom_file"][[1]],
    db = default_params["db"][[1]],
    data_dir = default_params["data_dir"][[1]],
    devel = default_params["devel"][[1]],
    devel_pkgdir = default_params["devel_pkgdir"][[1]],
    dset_column = default_params["dset_column"][[1]],
    env_column = default_params["env_column"][[1]],
    error_to_file = default_params["error_to_file"][[1]],
    gene_color_absent = default_params["gene_color_absent"][[1]],
    gene_color_present = default_params["gene_color_present"][[1]],
    input_format = default_params["input_format"][[1]],
    linearize = default_params["linearize"][[1]],
    meas_err = default_params["meas_err"][[1]],
    metadata_file = default_params["metadata_file"][[1]],
    min_fx = default_params["min_fx"][[1]],
    minimum = default_params["minimum"][[1]],
    ncl = default_params["ncl"][[1]],
    output_file = default_params["output_file"][[1]],
    out_dir = default_params["out_dir"][[1]],
    pctmin = default_params["pctmin"][[1]],
    phenotype_file = default_params["phenotype_file"][[1]],
    prev_color_high = default_params["prev_color_high"][[1]],
    prev_color_low = default_params["prev_color_low"][[1]],
    prior_file = default_params["prior_file"][[1]],
    prior_type = default_params["prior_type"][[1]],
    relative_out_dir = default_params["relative_out_dir"][[1]],
    sample_column = default_params["sample_column"][[1]],
    separate_process = default_params["separate_process"][[1]],
    single_dset = default_params["single_dset"][[1]],
    skip_graphs = default_params["skip_graphs"][[1]],
    spec_color_high = default_params["spec_color_high"][[1]],
    spec_color_low = default_params["spec_color_low"][[1]],
    spec_color_mid = default_params["spec_color_mid"][[1]],
    treemin = default_params["treemin"][[1]],
    type_16S = default_params["type_16S"][[1]],
    tax_level = default_params["tax_level"][[1]],
    use_rmd_params = default_params["use_rmd_params"][[1]],
    vsearch_16sfile = default_params["vsearch_16sfile"][[1]],
    vsearch_cutoff = default_params["vsearch_cutoff"][[1]],
    vsearch_dir = default_params["vsearch_dir"][[1]],
    vsearch_infile = default_params["vsearch_infile"][[1]],
    vsearch_outfile = default_params["vsearch_outfile"][[1]],
    which_envir = default_params["which_envir"][[1]],
    which_phenotype = default_params["which_phenotype"][[1]],
    working_dir = default_params["working_dir"][[1]]
)

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
#'   \item{data_dir}{String. Path to directory containing the data files required to perform a \emph{phylogenize} analysis. Default: "./data", but on package load, this default is set to the result of \code{system.file("extdata", package="phylogenize")}.}
#'   \item{error_to_file}{Boolean. Should pz.error, pz.warning, and pz.message output to an error message file? Default: FALSE}
#'   \item{input_format}{String. Whether to look for tabular or BIOM-formatted data ("tabular" or "biom"). Default: "tabular"}
#'   \item{metadata_file}{String. Name of metadata tabular file. Default: "test-metadata.tab"}
#'   \item{output_file}{String. Name of output file: "results.html"}
#'   \item{phenotype_file}{String. Name of input file for optional pre-calculated phenotype. Default: ""}
#'   \item{prior_file}{String. File name of optional pre-computed prior. Default: ""}
#'   \item{separate_metadata}{Boolean. For BIOM data, is there a separate tabular abundance table? Default: FALSE}
#'   \item{vsearch_16sfile}{String. Path to the 16S FASTA database that maps back to MIDAS species. Default: "16s_gtdb.frn"}
#'   \item{vsearch_bin}{String. File name of the binary of the aligner. Default: "vsearch"}
#'   \item{vsearch_dir}{String. Path where the binary of the aligner is found. Default: "/usr/local/bin/"}
#'   \item{vsearch_infile}{String. File name of the sequences written to disk and then read into the aligner. Default: "input_seqs.txt"}
#'   \item{vsearch_outfile}{String. File name where the aligner writes output which is then read back into \emph{phylogenize}. Default: "output_assignments.txt"}
#' }
#'
#' @section Computing phenotypes:
#' \describe{
#'   \item{assume_below_LOD}{Boolean. If TRUE, MIDAS species that are not present are assumed to have a prevalence of zero; if FALSE, they are dropped from the analysis. Default: TRUE}
#'   \item{db}{String. Type of data to use, gtdb or uhgp. Default: "gtdb"}
#'   \item{dset_column}{String. Name of column in metadata file containing the dataset annotations. Default: "dataset"}
#'   \item{env_column}{String. Name of column in metadata file containing the environment annotations. Default: "env"}
#'   \item{linearize}{Boolean. If TRUE, use a regular linear model instead of a phylogenetic linear model. Mostly useful for testing report generation, since the linear model is much faster but returns many more false positives. Default: FALSE}
#'   \item{meas_err}{Boolean. Separately estimate measurement error from phenotype variation in the phylogenetic linear model. Default: TRUE}
#'   \item{min_fx}{Positive double. Effects that are significantly equivalent to this effect size will be excluded from significant positive hits. If zero, the equivalence test will be skipped. Default: 0}
#'   \item{minimum}{Integer. A particular gene must be observed, and also absent, at least this many times to be reported as a significant positive association with the phenotype. Default: 3}
#'   \item{ncl}{Integer. Number of cores to use for parallel computation. Default: 1}
#'   \item{pctmin}{Float. A phylum must have at least this percent of observed representatives in order to be processed. Default: 0.01}
#'   \item{prior_type}{String. What type of prior to use ("uninformative" or "file"). Default: "uninformative"}
#'   \item{single_dset}{Boolean. If true, will assume that all samples come from a single dataset called \code{dset1} no matter what, if anything, is in \code{dset_column}. Default: FALSE}
#'   \item{treemin}{Integer. A phylum must have at least this many representatives in order to be processed. Default: 5}
#'   \item{type_16S}{Boolean. If 16S data other wise shotgun data is assumed. Default: False}
#'   \item{tax_level}{String. A classification of taxonomy which is either "phylum", "class", "order", "family", or "genus" Default: "family"}
#'   \item{vsearch_cutoff}{Float. Value between 0.95 and 1.00 giving the percent ID cutoff to use when assigning denoised sequence variants to MIDAS species using vsearch. Default: 0.985}
#'   \item{which_envir}{String. Environment in which to calculate prevalence or specificity. Must match annotations in metadata. Default: "Stool"}
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence" or "specificity"). Default: "prevalence"}
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

.onLoad <- function(libname, pkgname) {
    set_data_internal(startup=TRUE)
}

