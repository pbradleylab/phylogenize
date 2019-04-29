devtools::document("package/phylogenize")
devtools::load_all("package/phylogenize")
library(rlang)
library(tidyverse)
library(cowplot)
library(furrr)
plan(multiprocess, workers=6)
setwd("~/projects/phylogenize")
hmp_dir <- normalizePath(file.path(".", "hmp"))
emp_dir <- normalizePath(file.path(".", "emp"))

perform_associations <- FALSE
if (perform_associations) {
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
                     burst_dir="/home/pbradz/bin/",
                     ncl=10,
                     meas_err=TRUE,
                     devel=TRUE,
                     devel_pkgdir=file.path(getwd(), "package/phylogenize"),
                     pryr=FALSE)

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
                     burst_dir="/home/pbradz/bin/",
                     ncl=1,
                     devel=TRUE,
                     devel_pkgdir=file.path(getwd(), "package/phylogenize"),
                     linearize=TRUE,
                     pryr=FALSE)

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
                     burst_dir="/home/pbradz/bin/",
                     ncl=10,
                     meas_err=TRUE,
                     devel=TRUE,
                     devel_pkgdir=file.path(getwd(), "package/phylogenize"),
                     pryr=FALSE)

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
                     burst_dir = "/home/pbradz/bin/",
                     meas_err=TRUE,
                     ncl = 10,
                     data_dir = "/home/pbradz/projects/phylogenize/data/",
                     devel = TRUE,
                     devel_pkgdir = "/home/pbradz/projects/phylogenize/package/phylogenize/",
                     use_rmd_params = FALSE)

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
                     burst_dir = "/home/pbradz/bin/",
                     meas_err=TRUE,
                     ncl = 10,
                     data_dir = "/home/pbradz/projects/phylogenize/data/",
                     linearize = TRUE,
                     devel = TRUE,
                     devel_pkgdir = "/home/pbradz/projects/phylogenize/package/phylogenize/",
                     use_rmd_params = FALSE)

    phylogenize::render.report(
                     output_file=file.path(emp_dir, "emp-plant-specificity.html"),
                     out_dir = file.path(emp_dir, "plant-specificity-output"),
                     in_dir = "/home/pbradz/projects/phylogenize",
                     type = "16S",
                     db_version = "midas_v1.2",
                     which_phenotype = "specificity",
                     which_envir = "Plant",
                     abundance_file = "",
                     metadata_file = "",
                     biom_file = "emp_deblur_freeliving.biom",
                     input_format = "biom",
                     burst_dir = "/home/pbradz/bin/",
                     meas_err=TRUE,
                     ncl = 10,
                     data_dir = "/home/pbradz/projects/phylogenize/data/",
                     devel = TRUE,
                     devel_pkgdir = "/home/pbradz/projects/phylogenize/package/phylogenize/",
                     use_rmd_params = FALSE)

    phylogenize::render.report(
                     output_file=file.path(emp_dir, "emp-plant-specificity-linear.html"),
                     out_dir = file.path(emp_dir, "plant-specificity-linear-output"),
                     in_dir = "/home/pbradz/projects/phylogenize",
                     type = "16S",
                     db_version = "midas_v1.2",
                     which_phenotype = "specificity",
                     which_envir = "Plant",
                     abundance_file = "",
                     metadata_file = "",
                     biom_file = "emp_deblur_freeliving.biom",
                     burst_dir = "/home/pbradz/bin/",
                     input_format = "biom",
                     meas_err=TRUE,
                     ncl = 10,
                     data_dir = "/home/pbradz/projects/phylogenize/data/",
                     devel = TRUE,
                     devel_pkgdir = "/home/pbradz/projects/phylogenize/package/phylogenize/",
                     linearize=TRUE,
                     use_rmd_params = FALSE)
}

## Compare phenotypes (prevalence)
## Compare effect sizes and significance by phylum (for top 4 phyla)

pz.db <- import.pz.db(db_version="midas_v1.2")
pz.db.0 <- import.pz.db(db_version="midas_v1.0")

last_elem <- function(x) x[length(x)]
tax <- mutate(pz.db$taxonomy,
              id=map_chr(cluster, ~ last_elem(strsplit(.x, "_")[[1]]))) %>%
    as_tibble %>%
    select(taxon_id, kingdom, phylum,class,order,family,genus,species,cluster,id)

shotgun_pheno <- read_tsv(file.path(hmp_dir, "shotgun-output", "phenotype.tab"),
                          col_names=c("id", "pheno_sh"),
                          col_types=c("cd"), skip=1)
sixteen_pheno <- read_tsv(file.path(hmp_dir, "16S-results", "phenotype.tab"),
                          col_names=c("id", "pheno_16"), skip=1)
plant_pheno <- read_tsv(file.path(emp_dir, "plant-specificity-output", "phenotype.tab"),
                          col_names=c("id", "pheno_16"),
                          col_types=c("cd"), skip=1)
sixteen_rn <- mutate(sixteen_pheno,
                     id=map_chr(strsplit(id, "_"),
                                function(x) x[length(x)]))
pheno_cmp <- full_join(sixteen_rn, shotgun_pheno, by="id")
pheno_mg <- left_join(pheno_cmp, select(tax, phylum, id), by="id") %>% distinct

per_phylum_cor <- pheno_mg %>% group_by(phylum) %>% summarize(cor=cor(pheno_16, pheno_sh, use="pairwise.complete.obs")) %>% filter(!is.na(cor))
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

pdf(height=4.33,width=4.33,file="shotgun_vs_16s.pdf")
ggplot(pheno_mg %>%
       filter(phylum %in% detected_phy),
       aes(x=pheno_sh, y=pheno_16)) +
                                        #geom_point() +
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
    mutate(phylum=capwords(phylum)) %>%
    select(-X1)
linear_hmp16s_results <- read_csv(file.path(hmp_dir, "16S-linear-output",
                                        "all-results.csv")) %>%
    mutate(phylum=capwords(phylum)) %>%
    select(-X1)

make_qvs <- function(tbl,

                     .nestby="phylum",
                     .pvcol="p.value",
                     .qvcol="q.value") {
    .qv <- sym(.qvcol)
    .pv <- sym(.pvcol)
    .nb <- sym(.nestby)
    tbl %>%
        group_by(!!(.nb)) %>%
        nest %>%
        mutate(data=map(data, function(x) {
            mutate(x,
                   !!(.qv) := qvals(!!(.pv), error_to_file=FALSE))
        })) %>%
        unnest
}

make_equivs <- function(tbl, mfx=0.5) {
    mutate(tbl,
           equiv.pv = pmap_dbl(list(effect.size, std.err, df),
                               equiv_test,
                               min_fx=mfx))
}

process_data <- function(tbl) {
    tbl %>% make_equivs %>% make_qvs %>%
        make_qvs(., .pvcol="equiv.pv", .qvcol="equiv.qv")
}

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



total_cor <- genes_signif %>%
    group_by(phylum) %>%
    summarize(cor=cor(shotgun_fx, sixteen_fx))

restricted_cor <- genes_signif %>%
    filter(q.value.sh <= 0.05 | q.value.16 <= 0.05) %>%
    filter(equiv.qv.sh > 0.25 | equiv.qv.16 > 0.25) %>%
    group_by(phylum) %>%
    summarize(cor=cor(effect.size.sh, effect.size.16, use="pairwise.complete.obs"))

pos_sig_jaccard <- genes_signif %>%
    group_by(phylum) %>%
    filter(effect.size.sh > 0 | effect.size.16 > 0) %>%
    filter(q.value.sh <= 0.05 | q.value.16 <= 0.05) %>%
    filter(!is.na(q.value.sh), !is.na(q.value.16)) %>%
    summarize(jaccard=mean(q.value.sh <= 0.05 & q.value.16 <= 0.05 &
                           effect.size.sh > 0 & effect.size.16 > 0))

#g_test <- genes_signif %>%
#    group_by(phylum) %>%
#    nest %>%
#    mutate(data=map(data,
#                    ~ fisher.test(matrix(nr=2,
#                                         c(nrow(.),
#                                           sum(.$shotgun_qv <= 0.25),
#                                           sum(.$sixteen_qv <= 0.25),
#                                           sum(.$shotgun_qv <= 0.25 &
#                                               .$sixteen_qv <= 0.25)))))) %>%
#    mutate(pv=map_dbl(data, ~.$p.value),
#           est=map_dbl(data, ~.$estimate))
#
pdf(width=4.33,height=4.33,file="hmp-compare-genefx.pdf")
ggplot(genes_signif %>% filter(best_qv <= 0.05), #%>%
#      filter(equiv.qv.sh > 0.25 | equiv.qv.16 > 0.25),
      aes(x=effect.size.sh, y=effect.size.16)) +
#    geom_point(aes(color=pmax(std.err.sh, std.err.16))) + scale_color_gradient(low="black",high="white",limits=c(0, 0.25)) +
#    geom_point(aes(color=best_qv)) +
                                        #    scale_color_gradient(low="black",high="white", limits=c(0, 0.05)) +
    geom_hex(bins=40) + scale_fill_gradient(low="darkgray",high="black") +
    facet_wrap(~phylum) +
    geom_text(data=mutate(restricted_cor,
                          label=sprintf("cor_restricted = %0.2f", cor)),
              x=-5, y=4, aes(label=label), hjust=0) +
    geom_hline(yintercept=0, lty=2, col="gray") +
    geom_vline(xintercept=0, lty=2, col="gray") +
    theme(legend.position=c(1,0), legend.justification=c(1,0))
dev.off()


pz.db.0 <- import.pz.db(db_version="midas_v1.0")

sixteen_thresh_results <- group_by(sixteen_results, phylum) %>%
    nest %>%
    mutate(data=map(data, function(x) {
        y=rbind(x[, 2] %>% unlist,
                x[, 3] %>% unlist)
        colnames(y) <- x[, 1] %>% unlist
        rownames(y) <- c("fx", "pv")
        y
    })) %>%
    deframe
sixteen_signif <- make.sigs(sixteen_thresh_results, min.fx=0.5)
sixteen_signs <- make.signs(sixteen_thresh_results)
shotgun_thresh_results <- group_by(shotgun_results, phylum) %>%
    nest %>%
    mutate(data=map(data, function(x) {
        y=rbind(x[, 2] %>% unlist,
                x[, 3] %>% unlist)
        colnames(y) <- x[, 1] %>% unlist
        rownames(y) <- c("fx", "pv")
        y
    })) %>%
    deframe
shotgun_signif <- make.sigs(shotgun_thresh_results, min.fx=0.5)
shotgun_signs <- make.signs(shotgun_thresh_results)

sixteen_enr <- multi.enrich(sixteen_signif, sixteen_signs, pz.db$g.mappings)
shotgun_enr <- multi.enrich(shotgun_signif, shotgun_signs, pz.db$g.mappings)
merged_enr <- inner_join(rename(shotgun_enr, enr.estimate.shotgun=enr.estimate),
                         rename(sixteen_enr, enr.estimate.sixteen=enr.estimate),
                         by=c("phylum","cutoff","termset","term"))
ggplot(merged_enr, aes(x=log10(enr.estimate.shotgun), y=log10(enr.estimate.sixteen))) +
    geom_point() +
    geom_hline(yintercept=log10(2)) +
    geom_vline(xintercept=log10(2))


emp_dir=normalizePath(file.path(".", "emp"))
emp_phylo_results <- read_csv(file.path(emp_dir,
                                         "plant-specificity-output",
                                         "all-results.csv"),
                               col_types='cccdddd') %>%
    select(-X1) %>%
    process_data
emp_linear_results <- read_csv(file.path(emp_dir,
                                         "plant-specificity-linear-output",
                                        "all-results.csv"),
                              col_types='cccdddd') %>%
    select(-X1) %>%
    process_data

plant_cmp <- inner_join(emp_linear_results, emp_phylo_results, by=c("gene","phylum"), suffix=c(".phylo",".linear")) %>% filter(!is.na(effect.size.phylo), !is.na(effect.size.linear))

plant_signif <- mutate(plant_cmp,
                       phylo_sigq=-log10(q.value.phylo) * sign(effect.size.phylo),
                       linear_sigq=-log10(q.value.linear) * sign(effect.size.linear),
                       phylo_sigp=-log10(p.value.phylo) * sign(effect.size.phylo),
                       linear_sigp=-log10(p.value.linear) * sign(effect.size.linear),
                       mean_fx=(effect.size.linear + effect.size.phylo) / 2) %>%
    mutate(fx_range=cut(abs(mean_fx), c(0, 0.5, 1, Inf))) %>%
    mutate(best_qv=pmap_dbl(list(q.value.phylo, q.value.linear),
                        ~ min(.x, .y)))

plant_total_cor <- plant_signif %>%
    group_by(phylum) %>%
    summarize(cor=cor(effect.size.phylo, effect.size.linear))

plant_restricted_cor <- plant_signif %>%
    group_by(phylum) %>%
    filter(q.value.phylo <= 0.05 | q.value.linear <= 0.05) %>%
    summarize(cor=cor(effect.size.phylo, effect.size.linear))

plant_pos_sig_jaccard <- plant_signif %>%
    group_by(phylum) %>%
    filter(q.value.phylo <= 0.05 | q.value.linear <= 0.05) %>%
    summarize(jaccard=mean(q.value.phylo <= 0.05 & q.value.linear <= 0.05))

hmp_pos_sig_jaccard <- hmp_phylo_mod%>%
    group_by(phylum) %>%
    filter(q.value.phylo <= 0.05 | q.value.linear <= 0.05) %>%
    summarize(jaccard=mean(q.value.phylo <= 0.05 & q.value.linear <= 0.05))

#ggplot(plant_signif,
#       aes(x=phylo_sigp, y=linear_sigp)) +
#    geom_point(alpha=0.05) +
#    scale_color_continuous(low="blue", high="orange", limits=c(0, 2)) +
#    facet_grid(facets=fx_range~phylum) +
#    geom_vline(xintercept=-log10(0.25), lty=2, col="gray") +
#    geom_vline(xintercept=log10(0.25), lty=2, col="gray") +
#    geom_hline(yintercept=-log10(0.25), lty=2, col="gray") +
#    geom_hline(yintercept=log10(0.25), lty=2, col="gray") +
#    geom_text(data=mutate(plant_pos_sig_jaccard,
#                          label=sprintf("jaccard=%0.2f", jaccard)),
#              x=-10, y=20, aes(label=label), hjust=0)

plant_signif_bestphy <- plant_signif %>%
    filter(best_qv <= 0.05) %>%
    group_by(phylum) %>%
    summarize(n=length(gene)) %>%
    filter(n > 250) %>%
    select(phylum) %>%
    unlist

ggplot(plant_signif %>% filter(best_qv <= 0.05, phylum %in% plant_signif_bestphy),
       aes(x=effect.size.linear, y=effect.size.phylo)) +
                                        #    geom_point(aes(color=best_qv)) +
                                        #    scale_color_gradient(low="black",high="white", limits=c(0, 0.05)) +
    geom_hex() + scale_fill_gradient(low="darkgray",high="black") +
    facet_wrap(~phylum) +
    geom_text(data=mutate(plant_restricted_cor %>%
                          filter(phylum %in% plant_signif_bestphy),
                          label=sprintf("cor_restricted = %0.2f", cor)),
              x=-2, y=2, aes(label=label), hjust=0) +
    geom_hline(yintercept=0, lty=2, col="gray") +
    geom_vline(xintercept=0, lty=2, col="gray")

get_alpha <- function(phylum, gene, db=pz.db, lab=4) {
                                        #    attr(p, "class") <- "phylo"
    p <- db$trees[[phylum]]
    tips <- intersect(colnames(db$gene.presence[[phylum]]),
                      p$tip.label)
    df <- data.frame(g=db$gene.presence[[phylum]][gene, tips])
    if (var(df$g) == 0) { warning(paste0("zero-variance: ", gene)); return(NA) }
    tryCatch({phyloglm(g ~ 1,
                       data=df,
                       phy=keep.tip(p, tips),
                       log.alpha.bound=lab)$alpha },
             error=function(e) { warning(e); NA })
}

get_alpha <- function(phylum, gene, db=pz.db, lab=4) {
                                        #    attr(p, "class") <- "phylo"
    p <- db$trees[[phylum]]
    tips <- intersect(colnames(db$gene.presence[[phylum]]),
                      p$tip.label)
    df <- data.frame(g=db$gene.presence[[phylum]][gene, tips])
    if (var(df$g) == 0) { warning(paste0("zero-variance: ", gene)); return(NA) }
    tryCatch({phyloglm(g ~ 1,
                       data=df,
                       phy=keep.tip(p, tips),
                       log.alpha.bound=lab)$alpha },
             error=function(e) { warning(e); NA })
}

                                        # time-consuming...


ap_tree <- function(x) {
    attr(x, "class") <- "phylo"
    x
}

for (p in pz.db.0$phyla) {
    pz.db.0$trees[[p]] <- ap_tree(pz.db.0$trees[[p]])
}
for (p in pz.db$phyla) {
    pz.db$trees[[p]] <- ap_tree(pz.db$trees[[p]])
}

hmp_phylo_signal <- hmp_linear_cmp %>%
    filter((q.value.phylo <= 0.05) | (q.value.linear <= 0.05)) %>%
    filter(phylum %in% names(pz.db.0$trees)) %>%
    mutate(signal=future_pmap_dbl(list(phylum, gene), get_alpha, db=pz.db.0,
                                  .progress=TRUE))
write_tsv(select(hmp_phylo_signal, phylum, gene, signal), path="signal_hmp.tsv")
hmp_phylo_signal <- mutate(left_join(select(hmp_phylo_signal, phylum, gene, signal),
                                     hmp_linear_cmp,
                                     by=c("phylum", "gene")),
                           phylo=(q.value.phylo <= 0.05),
                           linear=(q.value.linear <= 0.05))
hmp_phylo_mod <- filter(hmp_phylo_signal, phylo | linear) %>%
    mutate(type=ifelse(phylo, ifelse(linear, yes="2_both", no="3_phylo"), no="1_linear"))

signal_cmp <- plant_signif %>%
    filter(best_qv <= 0.05) %>%
    filter(phylum %in% names(pz.db$trees)) %>%
    mutate(signal=future_pmap_dbl(list(phylum, gene), get_alpha),
           .progress=TRUE)

write_tsv(select(signal_cmp, phylum, gene, signal), path="signal_cmp_r.tsv")
# signal_cmp_r <- read_tsv("signal_cmp_r.tsv")
# signal_cmp <- left_join(signal_cmp_r, plant_signif)
signal_cmp_mod <- mutate(signal_cmp,
                         phylo=(q.value.phylo <= 0.05),
                         linear=(q.value.linear <= 0.05)) %>%
    filter(phylo | linear) %>%
    mutate(type=ifelse(phylo, ifelse(linear, yes="2_both", no="3_phylo"), no="1_linear"))

#emp_best_phyla <- signal_cmp %>%
#    group_by(phylum) %>%
#    summarize(n=length(gene)) %>%
#    filter(n > 250) %>%
#    select(phylum) %>%
#    unlist

nsig_each <- signal_cmp_mod %>%
    filter(q.value.phylo <= 0.05 | q.value.linear <= 0.05) %>%
    group_by(phylum) %>%
    summarize(l_only=sum(linear & !phylo),
              p_only=sum(phylo & !linear),
              nb=sum(phylo & linear))

emp_best_phyla <- nsig_each %>%
    filter(p_only + nb >= 250, l_only + nb >= 250) %>%
    select(phylum) %>%
    unlist

pdf(height=4.33, width=4.33, file="emp-compare-signal.pdf")
ggplot(signal_cmp_mod %>%
       filter(phylum %in% emp_best_phyla) %>%
       gather(key="type", value="value", phylo, linear) %>%
       filter(value) %>%
       rename(ives_garland_alpha=signal),
       aes(y=ives_garland_alpha, x=type)) +
    geom_violin() +
    facet_wrap(~phylum)
dev.off()

hmp_nsig_each <- hmp_phylo_mod %>%
    filter(q.value.phylo <= 0.05 | q.value.linear <= 0.05) %>%
    group_by(phylum) %>%
    summarize(l_only=sum(linear & !phylo),
              p_only=sum(phylo & !linear),
              nb=sum(phylo & linear))

hmp_best_phyla <- hmp_nsig_each %>%
    filter(l_only + nb >= 250, p_only + nb >= 250) %>%
    select(phylum) %>%
    unlist

pdf(height=4.33, width=4.33, file="hmp-compare-signal.pdf")
ggplot(hmp_phylo_mod %>%
       filter(phylum %in% hmp_best_phyla) %>%
       gather(key="type", value="value", phylo, linear) %>%
       filter(value) %>%
       rename(ives_garland_alpha=signal),
       aes(y=ives_garland_alpha, x=type)) +
    geom_violin() +
    facet_wrap(~phylum)

dev.off()

pheno_phy <- left_join(sixteen_pheno, pz.db$taxonomy %>% as_tibble %>% select(cluster, phylum) %>% distinct, by=c("id"="cluster"))
pheno_phy2 <- left_join(shotgun_pheno, pz.db$taxonomy %>% as_tibble %>% select(cluster, phylum) %>% distinct, by=c("id"="cluster"))


a <- rcoal(500)
b <- rbinTrait(n=1, beta=1, alpha=10, phy=a)
d <- rTrait(n=1, phy=a, model="BM") + (b * rnorm(50, mean=1, sd=0.5))
d_noise <- vapply(d, function(k) rnorm(n=100, mean=k, sd=sqrt(abs(k))), rep(1.0, 100)) %>% t
colnames(d_noise) <- paste0("rep", 1:100)
d_noise_long <- gather(cbind(tips=rownames(d_noise), as_tibble(d_noise)), rep, value, rep1:rep100, factor_key=TRUE) %>% as_tibble
expected_sd <- apply(d_noise, 1, sd)
d_full <- left_join(d_noise_long, enframe(b), by=c("tips"="name"))
tcor <- corBrownian(phy=a)
w <- apply(d_noise, 1, function(x) 1/var(x))

min_fx_pv <- function(fx, se, df, min_fx=0.25) {
    if (abs(fx) < abs(min_fx)) return(1)
    t_stat = abs((fx-min_fx)/se)
    2 * pt(-t_stat, df)
}


plant_pheno_byp <- left_join(plant_pheno, select(tax, phylum, cluster),
                             by=c("id"="cluster")) %>%
    mutate(phylum=make.names(tolower(phylum))) %>%
    group_by(phylum) %>%
    nest %>%
    left_join(., pz.db$trees %>% enframe %>% rename(phylum=name))

plant_pheno_signal <- plant_pheno_byp %>%
    mutate(signal=pmap(list(data, value),
                       ~ {
                           pheno=deframe(.x)
                           tips=intersect(names(pheno), .y$tip.label)
                           phylolm(pheno[tips] ~ 1, phy=keep.tips(.y, tips),
                                   measurement_error=TRUE)
                       }))

plant_pheno_r2 <- plant_pheno_signal %>%
    mutate(r2 = map_dbl(signal,
                        ~ {
                            rss = sum(.x$residuals ** 2)
                            tss = sum(.x$y ** 2)
                            r2 = 1 - (rss/tss)
                            }))


hmp_linear_enr <- read_csv("./hmp/16S-linear-output/enr-table.csv") %>% select(-X1)
hmp_phylo_enr <- read_csv("./hmp/16S-results/enr-table.csv") %>% select(-X1)
hmp_enr_cmp <- left_join(hmp_linear_enr, hmp_phylo_enr, by=c("Phylum",
                                                             "Gene_significance",
                                                             "Subsystem_level",
                                                             "Subsystem"),
                         suffix=c(".linear", ".phylo"))

plan(multiprocess, workers=10)

emp_phylo_enr <- alt.multi.enrich(emp_phylo_results, pz.db$g.mappings, future=TRUE)
emp_linear_enr <- alt.multi.enrich(emp_linear_results, pz.db$g.mappings, future=TRUE)

hmp_phylo_enr <- alt.multi.enrich(sixteen_equiv_qvs %>%
                                    filter(phylum %in% hmp_best_phyla),
                                  pz.db$g.mappings,
                                  future=TRUE)
hmp_linear_enr <- alt.multi.enrich(linear_hmp16s_equiv_qvs %>%
                                     filter(phylum %in% hmp_best_phyla),
                                   pz.db$g.mappings,
                                   future=TRUE)

emp_compare_enr <- left_join(emp_linear_enr, emp_phylo_enr, by=c("phylum","cutoff","termset","term"), suffix=c(".linear", ".phylo"))


rhizo_phylo_enr <- read_csv("./emp/plant-rhizosphere-phylo/enr-table.csv") %>% select(-X1)
rhizo_linear_enr <- read_csv("./emp/plant-rhizosphere-linear/enr-table.csv") %>% select(-X1)

rhizo_phylo_results <- read_csv("./emp/plant-rhizosphere-phylo/all-results.csv",
                                col_types="cccdddd") %>% select(-X1) %>%
    process_data
rhizo_linear_results <- read_csv("./emp/plant-rhizosphere-linear/all-results.csv",
                                 col_types="cccdddd") %>% select(-X1) %>%
    process_data
rhizo_best_phyla <- rhizo_phylo_results %>%
    filter(!is.na(effect.size)) %>%
    group_by(phylum) %>%
    filter(q.value <= 0.05) %>%
    summarize(n=n()) %>%
    filter(n > 250) %>%
    select(phylum) %>%
    unlist
#rhizo_phylo_enr <- alt.multi.enrich(rhizo_phylo_results %>%
#                                    filter(phylum %in% rhizo_best_phyla),
#                                    pz.db$g.mappings,
#                                    future=TRUE)
#rhizo_linear_enr <- alt.multi.enrich(rhizo_linear_results %>%
#                                    filter(phylum %in% rhizo_best_phyla),
#                                    pz.db$g.mappings,
#                                    future=TRUE)
#saveRDS(file="rhizo_phylo_enr.rds", rhizo_phylo_enr)
#saveRDS(file="rhizo_linear_enr.rds", rhizo_linear_enr)
rhizo_phylo_enr <- readRDS("rhizo_phylo_enr.rds")
rhizo_linear_enr <- readRDS("rhizo_linear_enr.rds")

rhizo_cmp_enr <- left_join(rhizo_phylo_enr, rhizo_linear_enr,
                           by=c("phylum", "cutoff", "termset", "term"),
                           suffix=c(".phylo",".linear")) %>%
    mutate(conditional.phylo=ifelse(enr.estimate.phylo > 1,
                                    enr.qval.phylo,
                                    1),
           conditional.linear=ifelse(enr.estimate.linear > 1,
                                     enr.qval.linear,
                                     1))
pdf(height=4.33, width=8.66, file="rhizo-compare.pdf")
ggplot(rhizo_cmp_enr %>% filter(cutoff=="strong"),
       aes(x=conditional.phylo, y=conditional.linear)) +
    geom_abline(slope=1, intercept=0, col="gray") +
    geom_hline(yintercept=0, col="gray") +
    geom_vline(xintercept=0, col="gray") +
    geom_hline(yintercept=0.25, col="gray", lty=2) +
    geom_vline(xintercept=0.25, col="gray", lty=2) +
    geom_hline(yintercept=1, col="gray") +
    geom_vline(xintercept=1, col="gray") +
    geom_point(size=3, alpha=0.3) +
    facet_wrap(~phylum)
dev.off()

write_tsv(rhizo_cmp_enr %>% filter(cutoff=="strong") %>%
          filter(conditional.phylo <= 0.25, conditional.linear > 0.25) %>%
         select(phylum, termset, term, conditional.phylo, conditional.linear, enr.estimate.phylo, enr.estimate.linear),
          "phylo_only.tsv")
write_tsv(rhizo_cmp_enr %>% filter(cutoff=="strong") %>%
          filter(conditional.phylo > 0.25, conditional.linear <= 0.25) %>%
          select(phylum, termset, term, conditional.phylo, conditional.linear, enr.estimate.phylo, enr.estimate.linear),
          "linear_only.tsv")
write_tsv(rhizo_cmp_enr %>% filter(cutoff=="strong") %>%
          filter(conditional.phylo <= 0.25, conditional.linear <= 0.25) %>%
          select(phylum, termset, term, conditional.phylo, conditional.linear, enr.estimate.phylo, enr.estimate.linear),
          "phylo_linear_both.tsv")
