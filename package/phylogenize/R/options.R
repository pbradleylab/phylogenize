# Options

PZ_OPTIONS <- options_manager(
  ncl=1,
  type="midas",
  source_dir=".",
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
  burst_dir="./bin",
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
  burst_outfile="output_assignments.txt"
)

#' Set and get options for phylogenize
#' 
#' @param ... Names of options (to retrieve) or \code{[key]=[value]} pairs (to set).
#'
#' @export
pz.options <- function(...) {
    stop_if_reserved(...)
    PZ_OPTIONS(...)
}
