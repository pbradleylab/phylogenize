# Script to run correlations with phylogenize

library(tidyverse)
library(readxl)
library(phylogenize)

min_pos <- function(x) min(x[x>0])
floor_data <- function(x) { x[x==0] <- min_pos(x)/2; x }
floor_mtx <- function(mtx) { apply(mtx, 1, floor_data) %>% t }
mean_ctr_mtx <- function(mtx) { apply(mtx, 1, mean_ctr) %>% t }
mean_ctr <- function(x) { x - mean(x) }

lipidomics <- read_xlsx("mmc3.xlsx", na="n.d.") %>%
    mutate(metab_id = map2_chr(`Metabolite names`,
                               1:nrow(.),
                               ~ paste0(.x, '__', .y)))

lipid_pg <- data.matrix(
    lipidomics %>%
    filter(Unit=="Picogram") %>%
    select(-Platforms, -`Metabolite names`, -`Unit`, -metab_id)
)

lipid_pa <- data.matrix(
    lipidomics %>%
    filter(Unit=="Peak area") %>%
    select(-Platforms, -`Metabolite names`, -`Unit`, -metab_id)
)

rownames(lipid_pg) <- lipidomics %>% filter(Unit=="Picogram") %>% `$`("metab_id")
rownames(lipid_pa) <- lipidomics %>% filter(Unit=="Peak area") %>% `$`("metab_id")

lipid_pg[is.na(lipid_pg)] <- 0
lipid_pg_floored <- apply(lipid_pg, 1, floor_data) %>% t
lipid_pa_floored <- apply(lipid_pa, 1, floor_data) %>% t

lipid_pg_logctr <- lipid_pg_floored %>% log2 %>% mean_ctr_mtx
lipid_pa_logctr <- lipid_pa_floored %>% log2 %>% mean_ctr_mtx
lipid_logctr <- rbind(lipid_pg_logctr, lipid_pa_logctr)

lipid_svd <- svd(lipid_logctr)

manifest <- read_tsv("DRP006201/manifest.tsv") %>%
  select(`sample-id`:AvgSpotLen) %>%
  mutate(`sample` = map_chr(`sample-id`, ~ {
    gsub("2w-", "", .x) %>%
      gsub("(.*)(\\d)$", "\\1_\\2", .) %>%
      gsub("H", "High", .) %>%
      gsub("L", "Low", .) %>%
      gsub("Vanc", "Van", .) %>%
      gsub("C", "Control", .) %>%
      gsub("X", "Abx", .) %>%
      gsub("-", "_", .)
  }))

lipid_logctr_t <- as_tibble(t(lipid_logctr),
                            rownames="sample")

lipid_pc_mtx <- lipid_svd$v
colnames(lipid_pc_mtx) <- paste0("PC_", 1:ncol(lipid_pc_mtx))
rownames(lipid_pc_mtx) <- colnames(lipid_logctr)

lipid_pcs <- lipid_pc_mtx %>% as_tibble(rownames="sample")

metadata <- left_join(manifest, lipid_logctr_t) %>% left_join(., lipid_pcs)
write_tsv(metadata, "lipidomics_metadata.tsv")

metadata_justpcs <- left_join(manifest, lipid_pcs)
write_tsv(metadata_justpcs, "lipidomics_metadata_pcs.tsv")


# HMP TEST
md_hmp <- read_tsv("../hmp/hmp-16s-phylogenize-metadata-full.tab") %>%
  mutate(corr_trait = map_dbl(env, ~ {
    if (.x %in% c("Buccal mucosa", "Palatine Tonsils", "Subgingival plaque", "Supragingival plaque", "Throat",
                  "Tongue dorsum", "Attached/Keratinized gingiva", "Saliva", "Hard palate")) { return(0.5) }
    if (.x == "Stool") return(1)
    return(0)
  } ))
write_tsv(md_hmp, "../hmp/test-hmp-md-silly.tab")
bs <- read_tsv("../hmp/hmp-shotgun-bodysite.tab")
colnames(bs) <- gsub("\\.\\.", "__", colnames(bs))
write_tsv(bs, "../hmp/shotgun-renamed.tab")


    phylogenize::render.report(
                     output_file=file.path(hmp_dir, "shotgun-corr-results-phylo.html"),
                     in_dir=hmp_dir,
                     out_dir=file.path(hmp_dir, "shotgun-corr-output"),
                     linearize=FALSE,
                     type="midas",
                     sample_column="sampleid",
                     db_version="midas_v1.0",
                     which_phenotype="correlation",
                     which_envir="Stool",
                     env_column="corr_trait",
                     abundance_file="shotgun-renamed.tab",
                     metadata_file="test-hmp-md-silly.tab",
                     burst_bin="vsearch",
                     data_dir=system.file(package="phylogenize", "extdata"),
                     input_format="tabular",
                     burst_dir=BURST_DIR,
                     ncl=NCL,
                     meas_err=TRUE,
                     pryr=FALSE)
