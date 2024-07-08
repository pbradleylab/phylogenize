# Does a quick sanity check to make sure that the repository is set up correctly
options_manager_check <- function(type, db_version, in_dir) {
    # Check if 'type' is valid
    valid_types <- c("midas", "16S")
    if (!type[[1]] %in% valid_types) {
        stop("Invalid value for 'type'. It should be either 'midas' for shotgun or '16S' for amplicon.")
    }
    # Check if 'db_version' is valid
    valid_db_versions <- c("midas_v1.0", "midas_v1.2", "gtdb_v214", "midas2_uhgg_fam")
    if (!db_version[[1]] %in% valid_db_versions) {
        stop("Invalid value for 'db_version'. It should be either 'midas_v1.0','midas_v1.2','gtdb_v214','midas2_uhgg_fam'.")
    }
    # Check if the output directory exists
    if (!dir.exists(in_dir[[1]])) {
        stop("The directory used by 'in_dir' does not exist. Please check your permissions and that the path is correct")
    }
    return(TRUE)
}

# Set the default parameters
default_params <- list(
    ncl=1,type="midas",
    out_dir="./output",
    in_dir=".",
    data_dir="./data",
    abundance_file="hmp-shotgun-bodysite.tab",
    metadata_file="hmp-shotgun-bodysite-metadata.tab",
    biom_file="test.biom",
    input_format="tabular",
    separate_metadata=FALSE,
    env_column="env",
    dset_column="dataset",
    sample_column="sample",
    phenotype_file="",
    db_version="midas_v1.2",
    which_phenotype="prevalence",
    which_envir="Stool",
    prior_type="uninformative",
    prior_file="",
    minimum=3,
    treemin=5,
    pctmin=0.025,
    assume_below_LOD=TRUE,
    skip_graphs=FALSE,
    vsearch_dir="/usr/local/bin",
    linearize=FALSE,
    check_mem_usage=FALSE,
    prev_color_low='black',
    prev_color_high='orange2',
    spec_color_low='slateblue',
    spec_color_mid='gray50',
    spec_color_high='tomato',
    gene_color_absent='black',
    gene_color_present='slateblue2',
    separate_process=TRUE,
    biom_dir='/usr/local/bin/',
    error_to_file=TRUE,
    vsearch_16sfile="16s_renamed.frn",
    vsearch_infile="input_seqs.txt",
    vsearch_outfile="output_assignments.txt",
    vsearch_cutoff=0.985,
    vsearch_bin='vsearch',
    use_rmd_params=FALSE,
    devel=FALSE,
    devel_pkgdir='package/phylogenize',
    relative_out_dir=NULL,
    single_dset=FALSE,
    working_dir='.',
    meas_err=TRUE,min_fx=0
)

# # Check that the inputs are valid. Note, this section coudl be expanded
# if (options_manager_check(default_params["type"], default_params["db_version"],
#                             default_params["in_dir"], default_params["abundance_file"])
# } else {
#     error("There was a problem with one or more of your inputs. Please see previous error messages")
# }

# Options given by user
PZ_OPTIONS <- options_manager(
  ncl = default_params["ncl"][[1]],
  type = default_params["type"][[1]],
  out_dir = default_params["out_dir"][[1]],
  in_dir = default_params["in_dir"][[1]],
  data_dir = default_params["data_dir"][[1]],
  abundance_file = default_params["abundance_file"][[1]],
  metadata_file = default_params["metadata_file"][[1]],
  biom_file = default_params["biom_file"][[1]],
  input_format = default_params["input_format"][[1]],
  separate_metadata = default_params["separate_metadata"][[1]],
  env_column = default_params["env_column"][[1]],
  dset_column = default_params["dset_column"][[1]],
  sample_column = default_params["sample_column"][[1]],
  phenotype_file = default_params["phenotype_file"][[1]],
  db_version = default_params["db_version"][[1]],
  which_phenotype = default_params["which_phenotype"][[1]],
  which_envir = default_params["which_envir"][[1]],
  prior_type = default_params["prior_type"][[1]],
  prior_file = default_params["prior_file"][[1]],
  minimum = default_params["minimum"][[1]],
  treemin = default_params["treemin"][[1]],
  pctmin = default_params["pctmin"][[1]],
  assume_below_LOD = default_params["assume_below_LOD"][[1]],
  skip_graphs = default_params["skip_graphs"][[1]],
  vsearch_dir = default_params["vsearch_dir"][[1]],
  linearize = default_params["linearize"][[1]],
  check_mem_usage = default_params["check_mem_usage"][[1]],
  prev_color_low = default_params["prev_color_low"][[1]],
  prev_color_high = default_params["prev_color_high"][[1]],
  spec_color_low = default_params["spec_color_low"][[1]],
  spec_color_mid = default_params["spec_color_mid"][[1]],
  spec_color_high = default_params["spec_color_high"][[1]],
  gene_color_absent = default_params["gene_color_absent"][[1]],
  gene_color_present = default_params["gene_color_present"][[1]],
  separate_process = default_params["separate_process"][[1]],
  biom_dir = default_params["biom_dir"][[1]],
  error_to_file = default_params["error_to_file"][[1]],
  vsearch_16sfile = default_params["vsearch_16sfile"][[1]],
  vsearch_infile = default_params["vsearch_infile"][[1]],
  vsearch_outfile = default_params["vsearch_outfile"][[1]],
  vsearch_cutoff = default_params["vsearch_cutoff"][[1]],
  vsearch_bin = default_params["vsearch_bin"][[1]],
  use_rmd_params = default_params["use_rmd_params"][[1]],
  devel = default_params["devel"][[1]],
  devel_pkgdir = default_params["devel_pkgdir"][[1]],
  relative_out_dir = default_params["relative_out_dir"][[1]],
  single_dset = default_params["single_dset"][[1]],
  working_dir = default_params["working_dir"][[1]],
  meas_err = default_params["meas_err"][[1]],
  min_fx = default_params["min_fx"][[1]]
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
#'   \item{out_dir}{String. Path to output directory. Default: "output"}
#'   \item{in_dir}{String. Path to input directory (i.e., where to look for input files). Default: "."}
#'   \item{data_dir}{String. Path to directory containing the data files required to perform a \emph{phylogenize} analysis. Default: "./data", but on package load, this default is set to the result of \code{system.file("extdata", package="phylogenize")}.}
#'   \item{working_dir}{String. Path to directory where relative paths should originate from. Default: \code{"."}}
#'   \item{abundance_file}{String. Name of abundance tabular file. Default: "test-abundance.tab"}
#'   \item{metadata_file}{String. Name of metadata tabular file. Default: "test-metadata.tab"}
#'   \item{biom_file}{String. Name of BIOM abundance-and-metadata file. Default: "test.biom"}
#'   \item{separate_metadata}{Boolean. For BIOM data, is there a separate tabular abundance table? Default: FALSE}
#'   \item{input_format}{String. Whether to look for tabular or BIOM-formatted data ("tabular" or "biom"). Default: "tabular"}
#'   \item{phenotype_file}{String. Name of input file for optional pre-calculated phenotype. Default: ""}
#'   \item{prior_file}{String. File name of optional pre-computed prior. Default: ""}
#'   \item{error_to_file}{Boolean. Should pz.error, pz.warning, and pz.message output to an error message file? Default: FALSE}
#'   \item{biom_dir}{String. Path to BIOM executables. Only used during testing. Default: "/usr/local/bin/"}
#'   \item{vsearch_dir}{String. Path where the binary of the aligner is found. Default: "/usr/local/bin/"}
#'   \item{vsearch_bin}{String. File name of the binary of the aligner. Default: "vsearch"}
#'   \item{vsearch_16sfile}{String. Path to the 16S FASTA database that maps back to MIDAS species. Default: "16s_renamed.frn"}
#'   \item{vsearch_infile}{String. File name of the sequences written to disk and then read into the aligner. Default: "input_seqs.txt"}
#'   \item{vsearch_outfile}{String. File name where the aligner writes output which is then read back into \emph{phylogenize}. Default: "output_assignments.txt"}
#' }
#'
#' @section Computing phenotypes:
#' \describe{
#'   \item{ncl}{Integer. Number of cores to use for parallel computation. Default: 1}
#'   \item{type}{String. Type of data to use, either "midas" (shotgun) or "16S" (amplicon). Default: "midas"}
#'   \item{env_column}{String. Name of column in metadata file containing the environment annotations. Default: "env"}
#'   \item{dset_column}{String. Name of column in metadata file containing the dataset annotations. Default: "dataset"}
#'   \item{sample_column}{String. Name of column in metadata file containing the sample IDs. Default: "sample_id"}
#'   \item{single_dset}{Boolean. If true, will assume that all samples come from a single dataset called \code{dset1} no matter what, if anything, is in \code{dset_column}. Default: FALSE}
#'   \item{db_version}{String. Which version of the MIDAS database to use ("midas_v1.2", "midas_v1.0", "gtdb_v214", "midas2_uhgg_fam"). Default: "midas_v1.2"}
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence" or "specificity"). Default: "prevalence"}
#'   \item{which_envir}{String. Environment in which to calculate prevalence or specificity. Must match annotations in metadata. Default: "Stool"}
#'   \item{prior_type}{String. What type of prior to use ("uninformative" or "file"). Default: "uninformative"}
#'   \item{minimum}{Integer. A particular gene must be observed, and also absent, at least this many times to be reported as a significant positive association with the phenotype. Default: 3}
#'   \item{assume_below_LOD}{Boolean. If TRUE, MIDAS species that are not present are assumed to have a prevalence of zero; if FALSE, they are dropped from the analysis. Default: TRUE}
#'   \item{linearize}{Boolean. If TRUE, use a regular linear model instead of a phylogenetic linear model. Mostly useful for testing report generation, since the linear model is much faster but returns many more false positives. Default: FALSE}
#'   \item{vsearch_cutoff}{Float. Value between 0.95 and 1.00 giving the percent ID cutoff to use when assigning denoised sequence variants to MIDAS species using vsearch. Default: 0.985}
#'   \item{meas_err}{Boolean. Separately estimate measurement error from phenotype variation in the phylogenetic linear model. Default: TRUE}
#'   \item{min_fx}{Positive double. Effects that are significantly equivalent to this effect size will be excluded from significant positive hits. If zero, the equivalence test will be skipped. Default: 0}
#'   \item{treemin}{Integer. A phylum must have at least this many representatives in order to be processed. Default: 5}
#'   \item{pctmin}{Float. A phylum must have at least this percent of observed representatives in order to be processed. Default: 0.01}
#' }
#'
#' @section Graphing:
#' \describe{
#'   \item{skip_graphs}{Boolean. If TRUE, skip making graphs in the report, which can be time- and memory-consuming. Default: FALSE}
#'   \item{prev_color_low}{String. When graphing prevalence on a tree, this color is the lowest value. Default: "black"}
#'   \item{prev_color_high}{String. When graphing prevalence on a tree, this color is the highest value. Default: "orange2"}
#'   \item{spec_color_high}{String. When graphing specificity on a tree, this color is the lowest value (most anti-specific). Default: "slateblue"}
#'   \item{spec_color_med}{String. When graphing specificity on a tree, this color denotes the prior (no association). Default: "gray50"}
#'   \item{spec_color_high}{String. When graphing specificity on a tree, this color is the highest value (most specific). Default: "tomato"}
#'   \item{gene_color_absent}{String. When graphing gene presence/absence, this color indicates absence. Default: "black"}
#'   \item{gene_color_present}{String. When graphing gene presence/absence, this color indicates presence. Default: "black"}
#' }
#'
#' @section Memory management:
#' \describe{
#'   \item{check_mem_usage}{Boolean. If TRUE, report memory usage when generating the report. Default: FALSE}
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
