devtools::load_all("/fs/project/bradley.720/users/kananen.13/tools/phylogenize/package/phylogenize/")

out_dir <- "/fs/project/bradley.720/users/kananen.13/tools/phylogenize/output/cirrhosis/uhgp/ancombc2/family"
 
core <- phylogenize_core(
    abundance_file = "/fs/project/bradley.720/users/kananen.13/tools/phylogenize/comms/phylogenize2/cirrhosis_uhgg2.tsv",
    metadata_file  = "/fs/project/bradley.720/users/kananen.13/tools/phylogenize/comms/phylogenize2/cirrhosis_uhgg2_metadata.tsv",
    input_format = "tabular",
    db = "human-gut",
    taxon_level = "family",
    which_envir = "case",
    sample_column = "sample",
    env_column = "env",
    dset_column = "dataset",
    only_specific_taxa="Lachnospiraceae",
    min_fx = 0,
    ncl = 4,
    out_dir = out_dir
  )

saveRDS(core, file.path(out_dir, "cirrhosis-ancombc2-family.rds"))

# To load this output back into memory after writing to disk:
# cirrhosis_family_abundance <- readRDS("output/cirrhosis_uhgp_abd_family/cirrhosis-fam-abd.rds")

render_core_report(
  core,
  output_file="cirrhosis-fam-abun.html",
  out_dir=out_dir)
