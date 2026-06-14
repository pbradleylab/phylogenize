devtools::load_all("/fs/project/bradley.720/users/kananen.13/tools/phylogenize/package/phylogenize")

out_dir <- "/fs/project/bradley.720/users/kananen.13/tools/phylogenize/output/cirrhosis/uhgp/poms/family"

core <- phylogenize_core(
    abundance_file = "/fs/project/bradley.720/users/kananen.13/tools/phylogenize/comms/phylogenize2/cirrhosis_uhgg2.tsv",                                                                          
    metadata_file  = "/fs/project/bradley.720/users/kananen.13/tools/phylogenize/comms/phylogenize2/cirrhosis_uhgg2_metadata.tsv",                                                                 
    input_format = "tabular",
    db = "human-gut",
    taxon_level = "family",
    core_method = "poms",
    which_envir = "case",
    sample_column = "sample",
    env_column = "env",
    only_specific_taxa="Lachnospiraceae",
    dset_column = "dataset",
    min_fx = 0,
    ncl = 4,
    out_dir = out_dir
  )

saveRDS(core, file.path(out_dir, "cirrhosis-poms-family.rds"))
