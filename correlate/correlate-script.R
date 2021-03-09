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
