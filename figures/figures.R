library(rlang)
library(tidyverse)
library(cowplot)
library(furrr)
library(phylogenize)
library(hexbin)
library(curl)
library(nlme)

source("figure-functions.R")

plan(multiprocess, workers=6)
hmp_dir <- normalizePath(file.path("..", "hmp"))
emp_dir <- normalizePath(file.path("..", "emp"))
NCL <- 8
BURST_DIR <- path.expand("~/bin/")

# uncomment to actually run associations (time consuming)
perform_associations <- strsplit(
#    "hmp16s hmp16s-linear hmpshotgun emp emp-linear",
    "",
    " ")[[1]]

if (!file.exists(file.path(hmp_dir, "hmp-16s-dada2-full.tab"))) {
    if (file.exists(file.path(hmp_dir, "hmp-16s-dada2-full.tab.xz"))) {
        message("unzipping HMP 16S data")
        system2("xz", paste0("-d ",
                       file.path(hmp_dir, "hmp-16s-dada2-full.tab.xz")),
                     "-k")
    } else {
        stop("HMP 16S files not found")
    }
}

if (!file.exists(file.path(hmp_dir, "hmp-shotgun-bodysite.tab"))) {
    if (file.exists(file.path(hmp_dir, "hmp-shotgun-bodysite.tab.xz"))) {
        message("unzipping HMP shotgun data")
        system2("xz", paste0("-d ",
                       file.path(hmp_dir, "hmp-shotgun-bodysite.tab.xz")),
                     "-k")
    } else {
        stop("HMP shotgun files not found")
    }
}

if (!file.exists(file.path(emp_dir, "emp_deblur_orig_metadata.biom"))) {
    message("downloading EMP data from Figshare")
    curl_download("https://ndownloader.figshare.com/files/17348873",
                  file.path(emp_dir, "emp_deblur_orig_metadata.biom"))
}

if ("hmp16s" %in% perform_associations) {
    phylogenize::render.report(
                     output_file=file.path(hmp_dir, "16S-results.html"),
                     in_dir=hmp_dir,
                     out_dir=file.path(hmp_dir, "16S-results"),
                     type="16S",
                     db_version="midas_v1.2",
                     which_phenotype="prevalence",
                     which_envir="Stool",
                     abundance_file="hmp-16s-dada2-full.tab",
                     metadata_file="hmp-16s-phylogenize-metadata-full.tab",
                     data_dir=system.file(package="phylogenize", "extdata"),
                     input_format="tabular",
                     burst_dir=BURST_DIR,
                     ncl=NCL,
                     meas_err=TRUE,
                     pryr=FALSE)
}


if ("hmp16s-linear" %in% perform_associations) {
    phylogenize::render.report(
                     output_file=file.path(hmp_dir, "16S-linear-results.html"),
                     in_dir=hmp_dir,
                     out_dir=file.path(hmp_dir, "16S-linear-output"),
                     type="16S",
                     db_version="midas_v1.2",
                     which_phenotype="prevalence",
                     which_envir="Stool",
                     abundance_file="hmp-16s-dada2-full.tab",
                     metadata_file="hmp-16s-phylogenize-metadata-full.tab",
                     data_dir=system.file(package="phylogenize", "extdata"),
                     input_format="tabular",
                     burst_dir=BURST_DIR,
                     ncl=NCL,
                     linearize=TRUE,
                     pryr=FALSE)
}

if ("hmpshotgun" %in% perform_associations) {
    phylogenize::render.report(
                     output_file=file.path(hmp_dir, "shotgun-results.html"),
                     in_dir=hmp_dir,
                     out_dir=file.path(hmp_dir, "shotgun-output"),
                     type="midas",
                     db_version="midas_v1.0",
                     which_phenotype="prevalence",
                     which_envir="Stool",
                     abundance_file="hmp-shotgun-bodysite.tab",
                     metadata_file="hmp-shotgun-bodysite-metadata.tab",
                     data_dir=system.file(package="phylogenize", "extdata"),
                     input_format="tabular",
                     burst_dir=BURST_DIR,
                     ncl=NCL,
                     meas_err=TRUE,
                     pryr=FALSE)
}

if ("emp" %in% perform_associations) {
   phylogenize::render.report(
                     output_file=file.path(emp_dir, "emp-plant-rhizosphere.html"),
                     out_dir = file.path(emp_dir, "plant-rhizosphere-phylo"),
                     in_dir = emp_dir,
                     type = "16S",
                     db_version = "midas_v1.2",
                     which_phenotype = "specificity",
                     single_dset = TRUE,
                     which_envir = "Plant rhizosphere",
                     env_column = "empo_3",
                     abundance_file = "",
                     metadata_file = "",
                     biom_file = "emp_deblur_orig_metadata.biom",
                     input_format = "biom",
                     burst_dir=BURST_DIR,
                     meas_err=TRUE,
                     ncl=NCL,
                     use_rmd_params = FALSE)
}

if ("emp-linear" %in% perform_associations) {
    phylogenize::render.report(
                     output_file=file.path(emp_dir, "emp-plant-rhizosphere-linear.html"),
                     out_dir = file.path(emp_dir, "plant-rhizosphere-linear"),
                     in_dir = emp_dir,
                     type = "16S",
                     db_version = "midas_v1.2",
                     which_phenotype = "specificity",
                     single_dset = TRUE,
                     which_envir = "Plant rhizosphere",
                     env_column = "empo_3",
                     abundance_file = "",
                     metadata_file = "",
                     biom_file = "emp_deblur_orig_metadata.biom",
                     input_format = "biom",
                     burst_dir=BURST_DIR,
                     meas_err=TRUE,
                     ncl=NCL,
                     linearize = TRUE,
                     use_rmd_params = FALSE)
}

## Compare phenotypes (prevalence)
## Compare effect sizes and significance by phylum (for top 4 phyla)

pz.db <- import.pz.db(db_version="midas_v1.2")
pz.db.0 <- import.pz.db(db_version="midas_v1.0")


tax <- mutate(pz.db$taxonomy,
              id=map_chr(cluster, ~ last_elem(strsplit(.x, "_")[[1]]))) %>%
    as_tibble %>%
    select(taxon_id, kingdom, phylum,class,order,family,genus,species,cluster,id)

shotgun_pheno <- read_tsv(file.path(hmp_dir, "shotgun-output", "phenotype.tab"),
                          col_names=c("id", "pheno_sh"),
                          col_types=c("cd"), skip=1)
sixteen_pheno <- read_tsv(file.path(hmp_dir, "16S-results", "phenotype.tab"),
                          col_names=c("id", "pheno_16"), skip=1)
sixteen_rn <- mutate(sixteen_pheno,
                     id=map_chr(strsplit(id, "_"),
                                function(x) x[length(x)]))
pheno_cmp <- full_join(sixteen_rn, shotgun_pheno, by="id")
pheno_mg <- left_join(pheno_cmp, select(tax, phylum, id), by="id") %>% distinct

per_phylum_cor <- pheno_mg %>%
    group_by(phylum) %>%
    summarize(cor=cor(pheno_16, pheno_sh, use="pairwise.complete.obs")) %>%
    filter(!is.na(cor))

detected_phy <- per_phylum_cor$phylum

ggplot(pheno_mg %>%
       filter(phylum %in% detected_phy),
       aes(x=pheno_sh, y=pheno_16)) +
    geom_point() +
    facet_wrap(~phylum) +
    geom_text(data=mutate(per_phylum_cor,
                          label=map_chr(cor, ~ sprintf("r = %0.2f", .))),
              x=-5, y=1, hjust=0, size=5,
              aes(label=label))

pdf(height=4.33,width=4.33,file="shotgun_vs_16s_abundances.pdf")
ggplot(pheno_mg %>%
       filter(phylum %in% detected_phy),
       aes(x=pheno_sh, y=pheno_16)) +
    geom_hex(bins=25) +
    scale_fill_continuous(low="#666666", high="black") +
    annotate("text",
             x=-5, y=0, hjust=0, size=5,
             label=sprintf("r = %0.2f",
                           with(pheno_cmp, cor(pheno_16,
                                               pheno_sh,
                                               use="pairwise.complete.obs")))) +
    theme(legend.justification=c(0,1), legend.position=c(0,1))
dev.off()

shotgun_results <- read_csv(file.path(hmp_dir, "shotgun-output", "all-results.csv")) %>%
    select(-X1)
sixteen_results <- read_csv(file.path(hmp_dir, "16S-results", "all-results.csv")) %>%
    mutate(phylum=phylogenize:::capwords(phylum)) %>%
    select(-X1)
linear_hmp16s_results <- read_csv(file.path(hmp_dir, "16S-linear-output",
                                        "all-results.csv")) %>%
    mutate(phylum=phylogenize:::capwords(phylum)) %>%
    select(-X1)

sixteen_equiv_qvs <- process_data(sixteen_results)
shotgun_equiv_qvs <- process_data(shotgun_results)
linear_hmp16s_equiv_qvs <- process_data(linear_hmp16s_results)

genes_cmp <- inner_join(shotgun_equiv_qvs,
                        sixteen_equiv_qvs,
                        by=c("gene","phylum"),
                        suffix=c(".sh",".16"))

hmp_linear_cmp <- inner_join(sixteen_equiv_qvs,
                             linear_hmp16s_equiv_qvs,
                             by=c("gene","phylum"),
                             suffix=c(".phylo",".linear"))

genes_signif <- mutate(genes_cmp,
                       signif.sh=-log10(q.value.sh) * sign(effect.size.sh),
                       signif.16=-log10(q.value.16) * sign(effect.size.16),
                       mean_fx=(effect.size.sh + effect.size.16) / 2) %>%
    mutate(fx_range=cut(abs(mean_fx), breaks=c(0, 0.5, 1, Inf))) %>%
    mutate(best_qv=pmap_dbl(list(q.value.sh, q.value.16),
                            ~ min(.x, .y)))


restricted_cor <- genes_signif %>%
    filter(q.value.sh <= 0.05 | q.value.16 <= 0.05) %>%
    group_by(phylum) %>%
    summarize(cor=cor(effect.size.sh, effect.size.16, use="pairwise.complete.obs"))

pdf(width=4.33,height=4.33,file="hmp-compare-genefx.pdf")
ggplot(genes_signif,
      aes(x=effect.size.sh, y=effect.size.16)) +
    geom_hex(bins=40) + scale_fill_gradient(low="darkgray",high="black") +
    facet_wrap(~phylum) +
    geom_text(data=mutate(restricted_cor,
                          label=sprintf("cor_restricted = %0.2f", cor)),
              x=-5, y=4, aes(label=label), hjust=0) +
    geom_hline(yintercept=0, lty=2, col="gray") +
    geom_vline(xintercept=0, lty=2, col="gray") +
    theme_cowplot(font_size=12) + 
    theme(legend.position=c(1,0), legend.justification=c(1,0)) 
dev.off()



for (p in pz.db.0$phyla) {
    pz.db.0$trees[[p]] <- ap_tree(pz.db.0$trees[[p]])
}
for (p in pz.db$phyla) {
    pz.db$trees[[p]] <- ap_tree(pz.db$trees[[p]])
}

rhizo_phylo_enr <- read_csv(file.path(emp_dir,
                                      "plant-rhizosphere-phylo/enr-table.csv")) %>%
    select(-X1)
rhizo_linear_enr <- read_csv(file.path(emp_dir,
                                       "plant-rhizosphere-linear/enr-table.csv")) %>%
    select(-X1)

rhizo_phylo_results <- read_csv(file.path(emp_dir,
                                          "plant-rhizosphere-phylo/all-results.csv"),
                                col_types="cccdddd") %>% select(-X1) %>%
    process_data
rhizo_linear_results <- read_csv(file.path(emp_dir,
                                           "plant-rhizosphere-linear/all-results.csv"),
                                 col_types="cccdddd") %>% select(-X1) %>%
    process_data

rhizo_cmp_results <- full_join(rhizo_phylo_results, rhizo_linear_results, 
  by = c("phylum", "gene"),
  suffix = c(".phylo", ".linear"))

rhizo_best_phyla <- rhizo_phylo_results %>%
    filter(!is.na(effect.size)) %>%
    group_by(phylum) %>%
    filter(q.value <= 0.05) %>%
    summarize(n=n()) %>%
    filter(n > 250) %>%
    select(phylum) %>%
    unlist

if (!(file.exists("rhizo_phylo_enr.rds") && file.exists("rhizo_linear_enr.rds"))) {
  rhizo_phylo_enr <- alt.multi.enrich(rhizo_phylo_results %>%
                                      filter(phylum %in% rhizo_best_phyla),
                                      pz.db$g.mappings,
                                      future=TRUE)
  rhizo_linear_enr <- alt.multi.enrich(rhizo_linear_results %>%
                                      filter(phylum %in% rhizo_best_phyla),
                                      pz.db$g.mappings,
                                      future=TRUE)
  saveRDS(file="rhizo_phylo_enr.rds", rhizo_phylo_enr)
  saveRDS(file="rhizo_linear_enr.rds", rhizo_linear_enr)
} else {
  rhizo_phylo_enr <- readRDS("rhizo_phylo_enr.rds")
  rhizo_linear_enr <- readRDS("rhizo_linear_enr.rds")
}

rhizo_cmp_enr <- left_join(rhizo_phylo_enr, rhizo_linear_enr,
                           by=c("phylum", "cutoff", "termset", "term"),
                           suffix=c(".phylo",".linear")) %>%
    mutate(conditional.phylo=ifelse(enr.estimate.phylo > 1,
                                    enr.qval.phylo,
                                    1),
           conditional.linear=ifelse(enr.estimate.linear > 1,
                                     enr.qval.linear,
                                     1)) %>%
    mutate(enr.overlap.phylo=map_chr(enr.overlap.phylo, ~ paste0(.x, collapse=","))) %>%
    mutate(enr.overlap.linear=map_chr(enr.overlap.linear, ~ paste0(.x, collapse=","))) %>%
    mutate(nitrogen = factor(grepl("[Nn]itrogen", term)))

pdf(height=4.33, width=8.66, file="rhizo-compare-new.pdf")
ggplot(rhizo_cmp_enr %>% filter(cutoff=="strong"),
       aes(x=conditional.phylo, y=conditional.linear, size=nitrogen)) +
    geom_abline(slope=1, intercept=0, col="gray") +
    geom_hline(yintercept=0, col="gray") +
    geom_vline(xintercept=0, col="gray") +
    geom_hline(yintercept=0.25, col="gray", lty=2) +
    geom_vline(xintercept=0.25, col="gray", lty=2) +
    geom_hline(yintercept=1, col="gray") +
    geom_vline(xintercept=1, col="gray") +
    geom_point(aes(size=nitrogen), alpha=0.3) +
    facet_wrap(~phylum) +
    scale_size_manual(values=c(2,8)) +
    theme_cowplot(font_size=12) +
    theme(legend.position = "none") 
dev.off()

write_tsv(rhizo_cmp_enr %>% filter(cutoff=="strong") %>%
          filter(conditional.phylo <= 0.25, conditional.linear > 0.25) %>%
          select(phylum, termset, term,
                 conditional.phylo, conditional.linear,
                 enr.estimate.phylo, enr.estimate.linear),
          "phylo_only.tsv")
write_tsv(rhizo_cmp_enr %>% filter(cutoff=="strong") %>%
          filter(conditional.phylo > 0.25, conditional.linear <= 0.25) %>%
          select(phylum, termset, term,
                 conditional.phylo, conditional.linear,
                 enr.estimate.phylo, enr.estimate.linear),
          "linear_only.tsv")
write_tsv(rhizo_cmp_enr %>% filter(cutoff=="strong") %>%
          filter(conditional.phylo <= 0.25, conditional.linear <= 0.25) %>%
          select(phylum, termset, term,
                 conditional.phylo, conditional.linear,
                 enr.estimate.phylo, enr.estimate.linear),
          "phylo_linear_both.tsv")


write_tsv(rhizo_cmp_enr %>% filter(cutoff=="strong") %>%
          filter(conditional.phylo <= 0.25 | conditional.linear <= 0.25) %>%
          select(phylum, termset, term,
            conditional.phylo, conditional.linear,
            enr.estimate.phylo, enr.estimate.linear),
  "phylo_linear_all_enr.tsv")

rhizo_cmp_sig_alpha <- rhizo_cmp_results %>%
    filter(q.value.phylo <= 0.05 | q.value.linear <= 0.05) %>%
    mutate(alpha = future_map2_dbl(phylum,
                                   gene,
                                   ~ get_alpha(.x, .y, lab=8),
                                   .progress=TRUE))
rhizo_cmp2_sig_alpha <- rhizo_cmp_sig_alpha %>%
    filter(phylum %in% setdiff(rhizo_best_phyla, 'fusobacteria')) %>%
    mutate(where_signif = factor(ifelse(q.value.phylo <= 0.05,
                                 ifelse(q.value.linear <= 0.05,
                                        '2_both',
                                        '1_phylo'), '3_linear'),
                                 ordered=TRUE))

ggplot(rhizo_cmp2_sig_alpha, aes(x=where_signif, y=alpha)) +
    facet_wrap(~phylum) +
geom_boxplot()

rhizo_sig_alpha_lvp <- rhizo_cmp2_sig_alpha %>%
    mutate(where_signif = ifelse(where_signif=='3_linear', TRUE, FALSE)) %>%
    mutate(log_alpha = log(alpha))


## Make contrasts, relative to the mean

phyl_factor <- rhizo_sig_alpha_lvp$phylum %>% as.factor %>% relevel("firmicutes")
rhizo_sig_alpha_lvp$phylum <- phyl_factor
nphyla <- length(levels(phyl_factor))
phyl_L <- vapply(2:nphyla, function(x) {
    col <- rep(-1/nphyla, nphyla)
    col[x] <- col[x] + 1
    col
}, rep(1, nphyla)) %>% t
rownames(phyl_L) <- paste0(levels(phyl_factor)[2:nphyla],
                                   ".vs.mean")
colnames(phyl_L) <- levels(phyl_factor)
phyl_L <- rbind(grande = 1/nphyla, phyl_L)
phyl_C <- solve(phyl_L)

test_lm <- lm(log_alpha ~ where_signif * phylum,
              data = rhizo_sig_alpha_lvp,
              contrasts=list(phylum = phyl_C[, -1],
                             where_signif = cbind(c(-1, 1))))
print(summary(test_lm))

pdf("lm-fit-new.pdf", useDingbats=FALSE)
rhizo_sig_alpha_lvp %>%
    mutate(fit=predict(test_lm)) %>%
    nest(-phylum, -where_signif, -fit) %>%
    select(-data) %>%
    bind_rows(rhizo_sig_alpha_lvp, .) %>%
    ggplot(., aes(where_signif, log_alpha)) +
      facet_wrap(~phylum) +
      geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) +
    geom_point(aes(group=1, y=fit), lwd=3) +
    stat_summary(fun.y="mean", aes(group=1, y=fit), lwd=1.5, geom="line")
dev.off()


## Make reduced dataset

assn <- read_tsv(file.path(hmp_dir, "output_assignments.txt"),
                 col_names=c("row", "match", "pctid", "len", "mismatch",
                             "gap", "startpos", "endpos", "starttgt", "endtgt",
                             "eval", "bitscore")) %>%
    mutate(midas_id = map_chr(match, ~ strsplit(.x, " ")[[1]][3]))

assn <- left_join(assn,
                  pz.db$taxonomy %>% select(cluster, phylum) %>% distinct,
                  by=c("midas_id"="cluster"))

bac_assn <- assn %>%
    nest(-row) %>%
    filter(map_lgl(data, ~ all(.$phylum == "Bacteroidetes"))) %>%
    unnest

bac_rows <- map_int(bac_assn$row, ~ gsub("Row", "", .x) %>% as.integer)

hmp_16s <- read_tsv(file.path(hmp_dir, "hmp-16s-dada2-full.tab"))
hmp_16s_metadata <- read_tsv(file.path(hmp_dir, "hmp-16s-phylogenize-metadata-full.tab"))

retain_cols <- hmp_16s_metadata %>%
    filter(env %in% c("Stool", "Saliva", "Hard palate", "Tongue dorsum",
                      "Supragingival plaque", "Subgingival plaque",
                      "Buccal mucosa")) %>%
    select(sample) %>%
    deframe %>%
    intersect(., colnames(hmp_16s))

small_hmp <- hmp_16s[bac_rows, c("rowname", retain_cols)]
                                        # remove zero-counts
small_hmp <- small_hmp[-which(rowSums(small_hmp[, -1]) == 0), ]
small_hmp <- small_hmp[, -(1+which(colSums(small_hmp[, -1]) == 0))]
write_tsv(small_hmp, file.path(hmp_dir, "hmp-16s-dada2-bacteroidetes.tab"))
write_tsv(hmp_16s_metadata %>%
          filter(sample %in% colnames(small_hmp)),
          file.path(hmp_dir, "hmp-16s-metadata-bacteroidetes.tab"))

## This is how you would analyze the smaller dataset
if (FALSE) {
  phylogenize::render.report(
                  output_file=file.path(hmp_dir, "16S-bacteroidetes-results.html"),
                  in_dir=hmp_dir,
                  out_dir=file.path(hmp_dir, "16S-bacteroidetes-output"),
                  type="16S",
                  db_version="midas_v1.2",
                  which_phenotype="specificity",
                  which_envir="Stool",
                  abundance_file="hmp-16s-dada2-bacteroidetes.tab",
                  metadata_file="hmp-16s-metadata-bacteroidetes.tab",
                  data_dir=system.file(package="phylogenize", "extdata"),
                  input_format="tabular",
                  burst_dir="/home/pbradz/bin/",
                  ncl=1,
                  linearize=FALSE,
                  pryr=FALSE)
}
