devtools::document("package/phylogenize")
devtools::load_all("package/phylogenize")
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
                     input_format = "biom",
                     burst_dir = "/home/pbradz/bin/",
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
library(tidyverse)
library(cowplot)
last_elem <- function(x) x[length(x)]
tax <- mutate(pz.db$taxonomy,
              id=map_chr(cluster, ~ last_elem(strsplit(.x, "_")[[1]]))) %>%
    as_tibble %>%
    select(taxon_id, kingdom, phylum,class,order,family,genus,species,cluster,id)
shotgun_pheno <- read_tsv(file.path(hmp_dir, "shotgun-output", "phenotype.tab"),
                          col_names=c("id", "pheno_sh"),
                          col_types=c("cd"), skip=1)
sixteen_pheno <- read_tsv(file.path(hmp_dir, "16S-output", "phenotype.tab"),
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
    rename(shotgun_fx=effect.size, shotgun_pv=p.value) %>%
    select(-X1)
sixteen_results <- read_csv(file.path(hmp_dir, "16S-output", "all-results.csv")) %>%
    mutate(phylum=capwords(phylum)) %>%
    rename(sixteen_fx=effect.size, sixteen_pv=p.value) %>%
    select(-X1)
linear_hmp16s_results <- read_csv(file.path(hmp_dir, "16S-linear-output",
                                        "all-results.csv")) %>%
    mutate(phylum=capwords(phylum)) %>%
    rename(linear_fx=effect.size, linear_pv=p.value) %>%
    select(-X1)

genes_cmp <- inner_join(shotgun_results, sixteen_results, by=c("gene","phylum"))
hmp_linear_cmp <- inner_join(sixteen_results %>%
                             rename(phylo_fx=sixteen_fx,
                                    phylo_pv=sixteen_pv),
                             linear_hmp16s_results,
                             by=c("gene","phylum"))

hmp_cmp_qvs <- group_by(hmp_linear_cmp, phylum) %>%
    nest %>%
    mutate(data=map(data, function(x) {
        mutate(x,
               phylo_qv=qvals(phylo_pv, error_to_file=FALSE),
               linear_qv=qvals(linear_pv, error_to_file=FALSE))
    })) %>%
    unnest

genes_qvs <- group_by(genes_cmp, phylum) %>%
    nest %>%
    mutate(data=map(data, function(x) {
        mutate(x,
               shotgun_qv=qvals(shotgun_pv),
               sixteen_qv=qvals(sixteen_pv))
    })) %>%
    unnest
genes_1sig <- filter(genes_qvs, (sixteen_qv <= 0.25) | (shotgun_qv <= 0.25))
genes_cmp <- inner_join(shotgun_results, sixteen_results, by=c("gene","phylum"))

logit_pseud <- function(x, fudge=1e-10) logit((x + fudge) / (1 + (2 * fudge)))

genes_signif <- mutate(genes_qvs,
                       shotgun_sig=-log10(shotgun_qv) * sign(shotgun_fx),
                       sixteen_sig=-log10(sixteen_qv) * sign(sixteen_fx),
                       mean_fx=(sixteen_fx + shotgun_fx) / 2) %>%
    mutate(fx_range=cut(abs(mean_fx), breaks=c(0, 0.5, 1, Inf))) %>%
    mutate(best_qv=pmap_dbl(list(shotgun_qv, sixteen_qv),
                            ~ min(.x, .y)))

total_cor <- genes_signif %>%
    group_by(phylum) %>%
    summarize(cor=cor(shotgun_fx, sixteen_fx))

restricted_cor <- genes_signif %>%
    filter(sixteen_qv <= 0.05 | shotgun_qv <= 0.05) %>%
    group_by(phylum) %>%
    summarize(cor=cor(shotgun_fx, sixteen_fx, use="pairwise.complete.obs"))

pos_sig_jaccard <- genes_signif %>%
    group_by(phylum) %>%
    filter(shotgun_qv <= 0.05 | sixteen_qv <= 0.05) %>%
    summarize(jaccard=mean(shotgun_qv <= 0.05 & sixteen_qv <= 0.05))

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
ggplot(genes_signif %>% filter(best_qv <= 0.05),
       aes(x=shotgun_fx, y=sixteen_fx)) +
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
emp_linear_results <- read_csv(file.path(emp_dir,
                                         "plant-linear",
                                         "all-results.csv")) %>%
        rename(linear_fx=effect.size, linear_pv=p.value) %>%
        select(-X1)
emp_phylo_results <- read_csv(file.path(emp_dir,
                                         "plant-phylo",
                                         "all-results.csv")) %>%
    rename(phylo_fx=effect.size, phylo_pv=p.value) %>%
    select(-X1)

plant_cmp <- inner_join(emp_linear_results, emp_phylo_results, by=c("gene","phylum"))
plant_qvs <- group_by(plant_cmp, phylum) %>%
    nest %>%
    mutate(data=map(data, function(x) {
        mutate(x,
               linear_qv=qvals(linear_pv, error_to_file=FALSE),
               phylo_qv=qvals(phylo_pv, error_to_file=FALSE))
    })) %>%
    unnest

plant_signif <- mutate(plant_qvs,
                       phylo_sigq=-log10(phylo_qv) * sign(phylo_fx),
                       linear_sigq=-log10(linear_qv) * sign(linear_fx),
                       phylo_sigp=-log10(phylo_pv) * sign(phylo_fx),
                       linear_sigp=-log10(linear_pv) * sign(linear_fx),
                       mean_fx=(linear_fx + phylo_fx) / 2) %>%
    mutate(fx_range=cut(abs(mean_fx), c(0, 0.5, 1, Inf))) %>%
    mutate(best_qv=pmap_dbl(list(phylo_qv, linear_qv),
                        ~ min(.x, .y)))

plant_total_cor <- plant_signif %>%
    group_by(phylum) %>%
    summarize(cor=cor(phylo_fx, linear_fx))

plant_restricted_cor <- plant_signif %>%
    group_by(phylum) %>%
    summarize(cor=cor(phylo_fx, linear_fx))

plant_pos_sig_jaccard <- plant_signif %>%
    group_by(phylum, fx_range) %>%
    filter(phylo_qv <= 0.25 | linear_qv <= 0.25) %>%
    summarize(jaccard=mean(phylo_qv <= 0.25 & linear_qv <= 0.25))

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
       aes(x=linear_fx, y=phylo_fx)) +
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

get_alpha <- function(phylum, gene) {
    tips <- intersect(colnames(pz.db$gene.presence[[phylum]]),
                      pz.db$trees[[phylum]]$tip.label)
    df <- data.frame(g=pz.db$gene.presence[[phylum]][gene, tips])
    phyloglm(g ~ 1, data=df, phy=keep.tip(pz.db$trees[[phylum]], tips))$alpha
}

get_alpha <- function(phylum, gene, db=pz.db) {
    tips <- intersect(colnames(db$gene.presence[[phylum]]),
                      db$trees[[phylum]]$tip.label)
    df <- data.frame(g=db$gene.presence[[phylum]][gene, tips])
    phyloglm(g ~ 1, data=df, phy=keep.tip(db$trees[[phylum]], tips))$alpha
}


                                        # time-consuming...
library(furrr)
plan(multiprocess)
hmp_phylo_signal <- hmp_cmp_qvs %>%
    filter((phylo_qv <= 0.05) | (linear_qv <= 0.05)) %>%
    mutate(phylum=tolower(phylum)) %>%
    mutate(signal=future_pmap_dbl(list(phylum, gene), get_alpha, db=pz.db.0,
                                  .progress=TRUE))
write_tsv(select(hmp_phylo_signal, phylum, gene, signal), path="signal_hmp.tsv")
hmp_phylo_signal <- mutate(hmp_phylo_signal,
                           phylo=phylo_qv <= 0.05,
                           linear=linear_qv <= 0.05)

# signal_cmp <- plant_signif %>%
#     filter(best_qv <= 0.05) %>%
#     mutate(signal=pmap_dbl(list(phylum, gene), get_alpha))
# write_tsv(select(signal_cmp, phylum, gene, signal), path="signal_cmp_r.tsv")
signal_cmp_r <- read_tsv("signal_cmp_r.tsv")
signal_cmp <- left_join(signal_cmp_r, plant_signif)
signal_cmp_mod <- mutate(signal_cmp,
                         phylo=phylo_qv <= 0.05,
                         linear=linear_qv <= 0.05)
pdf(height=4.33, width=4.33, file="emp-compare-signal.pdf")
ggplot(signal_cmp_mod %>%
       filter(phylum %in% plant_signif_bestphy) %>%
       rename(ives_garland_alpha=signal),
       aes(y=ives_garland_alpha, x=phylo)) +
    geom_violin() +
    facet_wrap(~phylum)
dev.off()

hmp_best_phyla <- hmp_phylo_signal %>%
    group_by(phylum) %>%
    summarize(n=length(gene)) %>%
    filter(n > 250) %>%
    select(phylum) %>%
    unlist

pdf(height=4.33, width=4.33, file="hmp-compare-signal.pdf")
ggplot(hmp_phylo_signal %>%
       filter(phylum %in% hmp_best_phyla) %>%
       rename(ives_garland_alpha=signal),
       aes(y=ives_garland_alpha, x=phylo)) +
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

