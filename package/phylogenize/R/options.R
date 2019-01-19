# Options

PZ_OPTIONS <- options_manager(
  ncl=1,
  type="midas",
  out_dir="output",
  in_dir=".",
  data_dir="./data",
  abundance_file="test-abundance.tab",
  metadata_file="test-metadata.tab",
  biom_file="test.biom",
  input_format="tabular",
  env_column="env",
  dset_column="dataset",
  phenotype_file="",
  db_version="midas_v1.2",
  which_phenotype="prevalence",
  which_envir="Stool",
  prior_type="uninformative",
  prior_file="",
  minimum=3,
  treemin=5,
  assume_below_LOD=TRUE,
  skip_graphs=FALSE,
  burst_dir="/usr/local/bin",
  linearize=FALSE,
  pryr=TRUE,
  prev_color_low='black',
  prev_color_high='orange2',
  spec_color_low='slateblue',
  spec_color_mid='gray50',
  spec_color_high='tomato',
  gene_color_absent='black',
  gene_color_present='slateblue2',
  separate_process=TRUE,
  biom_dir='/usr/local/bin/',
  error_to_file=FALSE,
  burst_16sfile="16s_renamed.frn",
  burst_infile="input_seqs.txt",
  burst_outfile="output_assignments.txt",
  burst_cutoff=0.985
)

#' Set and get options for phylogenize.
#'
#' Function to set and get global options for the \emph{phylogenize} package.
#'
#' These options are global because they affect how most of the functions in \emph{phylogenize} work. Descriptions of these options follow.
#' 
#' @section File input/output and paths:
#' \itemize{
#'   \item \strong{out_dir} String. Path to output directory. Default: "output"
#'   \item \strong{in_dir} String. Path to input directory (i.e., where to look for input files). Default: "."
#'   \item \strong{data_dir} String. Path to directory containing the data files required to perform a \emph{phylogenize} analysis. Default: on package load, this default is set to the result of \code{system.file("extdata", package="phylogenize")}.
#'   \item \strong{abundance_file} String. Name of abundance tabular file. Default: "test-abundance.tab"
#'   \item \strong{metadata_file} String. Name of metadata tabular file. Default: "test-metadata.tab"
#'   \item \strong{biom_file} String. Name of BIOM abundance-and-metadata file. Default: "test.biom"
#'   \item \strong{input_format} String. Whether to look for tabular or BIOM-formatted data ("tabular" or "biom"). Default: "tabular"
#'   \item \strong{phenotype_file} String. Name of input file for optional pre-calculated phenotype. Default: ""
#'   \item \strong{prior_file} String. File name of optional pre-computed prior. Default: ""
#'   \item \strong{error_to_file} Boolean. Should pz.error, pz.warning, and pz.message output to an error message file? Default: FALSE
#'   \item \strong{biom_dir} String. Path to BIOM executables. Only used during testing. Default: "/usr/local/bin/"
#'   \item \strong{burst_dir} String. Path to the binary of BURST. Default: "/usr/local/bin/"
#'   \item \strong{burst_16sfile} String. Path to the 16S FASTA database that maps back to MIDAS species. Default: "16s_renamed.frn"
#'   \item \strong{burst_infile} String. File name of the sequences written to disk and then read into BURST. Default: "input_seqs.txt"
#'   \item \strong{burst_outfile} String. File name where BURST writes output which is then read back into \emph{phylogenize}. Default: "output_assignments.txt"
#' }
#' @section Computing phenotypes:
#' \itemize{
#'   \item \strong{ncl} Integer. Number of cores to use for parallel computation. Default: 1
#'   \item \strong{type} String. Type of data to use, either "midas" (shotgun) or "16S" (amplicon). Default: "midas"
#'   \item \strong{env_column} String. Name of column in metadata file containing the environment annotations. Default: "env"
#'   \item \strong{dset_column} String. Name of column in metadata file containing the dataset annotations. Default: "dataset"
#'   \item \strong{db_version} String. Which version of the MIDAS database to use ("midas_v1.2" or "midas_v1.0"). Default: "midas_v1.2"
#'   \item \strong{which_phenotype} String. Which phenotype to calculate ("prevalence" or "specificity"). Default: "prevalence"
#'   \item \strong{which_envir} String. Environment in which to calculate prevalence or specificity. Must match annotations in metadata. Default: "Stool"
#'   \item \strong{prior_type} String. What type of prior to use ("uninformative" or "file"). Default: "uninformative"
#'   \item \strong{minimum} Integer. A particular gene must be observed, and also absent, at least this many times to be reported as a significant positive association with the phenotype. Default: 3
#'   \item \strong{assume_below_LOD} Boolean. If TRUE, MIDAS species that are not present are assumed to have a prevalence of zero; if FALSE, they are dropped from the analysis. Default: TRUE
#'   \item \strong{linearize} Boolean. If TRUE, use a regular linear model instead of a phylogenetic linear model. Mostly useful for testing report generation, since the linear model is much faster but returns many more false positives. Default: FALSE
#'   \item \strong{burst_cutoff} Float. Value between 0.95 and 1.00 giving the percent ID cutoff to use when assigning denoised sequence variants to MIDAS species using BURST. Default: 0.985
#' }
#' @section Graphing:
#' \itemize{
#'   \item \strong{treemin} Integer. A phylum must have at least this many representatives in order to be graphed in the report. Default: 5
#'   \item \strong{skip_graphs} Boolean. If TRUE, skip making graphs in the report, which can be time- and memory-consuming. Default: FALSE
#'   \item \strong{prev_color_low} String. When graphing prevalence on a tree, this color is the lowest value. Default: "black"
#'   \item \strong{prev_color_high} String. When graphing prevalence on a tree, this color is the highest value. Default: "orange2"
#'   \item \strong{spec_color_high} String. When graphing specificity on a tree, this color is the lowest value (most anti-specific). Default: "slateblue"
#'   \item \strong{spec_color_med} String. When graphing specificity on a tree, this color denotes the prior (no association). Default: "gray50"
#'   \item \strong{spec_color_high} String. When graphing specificity on a tree, this color is the highest value (most specific). Default: "tomato"
#'   \item \strong{gene_color_absent} String. When graphing gene presence/absence, this color indicates absence. Default: "black"
#'   \item \strong{gene_color_present} String. When graphing gene presence/absence, this color indicates presence. Default: "black"
#' }
#' @section Memory management:
#' \itemize{
#'   \item \strong{pryr} Boolean. If TRUE, report memory usage when generating the report. Default: FALSE
#'   \item \strong{separate_process} Boolean. When displaying clustered top gene associations alongside a tree colored by phenotype, this flag indicates whether to use a separate subprocess. This allows memory used by clustering to be released back to the operating system immediately. Default: TRUE
#' }
#'
#' @param ... Names of options (to retrieve) or \code{[key]=[value]} pairs (to set).
#'
#' @export
pz.options <- function(...) {
    stop_if_reserved(...)
    PZ_OPTIONS(...)
}

pz.options(data_dir=system.file("extdata", package="phylogenize"))
