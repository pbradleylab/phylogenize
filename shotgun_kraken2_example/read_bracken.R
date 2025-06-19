library(pbapply)
library(Matrix)
library(tidyverse)
select <- dplyr::select

wide_to_mtx <- function(tbl) {
  mtx <- data.matrix(tbl[,-1])
  rownames(mtx) <- tbl[[1]]
  mtx
}

### Part 1: combine Bracken outputs

data_dir = "./bracken/"
bracken_files <- list.files(path = file.path(data_dir, .), pattern = ".*bracken")
if (length(bracken_files) == 0) { stop("bracken files not found") }

bracken_data <- map(bracken_files, ~ { read_tsv(file.path(data_dir, .), show_col_types=FALSE) })
bracken_samples <- gsub("\\.bracken", "", bracken_files)
bracken_data_newcol <- map2(bracken_data, bracken_samples, ~ { mutate(.x, sample = .y)})
bracken_data_tidy <- bind_rows(bracken_data_newcol)
bracken_data_wide <- pivot_wider(bracken_data_tidy, names_from = sample, values_from = fraction_total_reads, id_cols = name, values_fill = 0)
bracken_data_wide_counts <- pivot_wider(bracken_data_tidy, names_from = sample, values_from = new_est_reads, id_cols = name, values_fill = 0)

### Part 2: rename to match Phylogenize2 database

# Import species taxonomy from MIDAS2/UHGG database and make them match what's in Kraken database
genome_metadata <- read_tsv("metadata.tsv.gz") %>%
  separate(
    col = Lineage,
    into = c(
      "domain",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    ),
    sep = ";"
  )
genome_metadata_fixed <- genome_metadata %>%
  mutate(species_name = map2_chr(species, MGnify_accession, ~ {
    if (.x == "s__") {
      paste0("s__", .y)
    } else {
      .x
    }
  }))
# There is one more species that needs to be renamed before everything matches:
#   setdiff(ge, nome_metadata_fixed$species_name, bracken_data_tidy$name %>% unique)
#   [1] "s__Agathobacter rectalis"
#   setdiff(bracken_data_tidy$name %>% unique, genome_metadata_fixed$species_name)
#   [1] "s__Agathobacter rectale"
genome_metadata_fixed2 <-
  genome_metadata_fixed %>% mutate(species_name = str_replace(
    species_name,
    "s__Agathobacter rectalis",
    "s__Agathobacter rectale"
  ))

### Part 3: collapse reads for the same individual
# Load in and match sample/run metadata

sample_metadata <- read_tsv("cirrhosis-igg-samples.tsv")
sra_metadata <- read_csv("SraRunTable.txt")
sra_metadata_with_env <- left_join(sra_metadata, sample_metadata, by=c("sample_acc"="sample_id"))

indiv_md <- sra_metadata_with_env %>%
    separate("Sample_Name", "_", into=c("individual", "seq_run")) %>%
    distinct

bracken_count_setup <- bracken_data_tidy %>% 
  left_join(indiv_md, by = c("sample" = "Run")) %>%
  ungroup %>%
  select(name, sample, individual, new_est_reads) %>%
  group_by(name, individual) %>%
  summarize(total_reads = sum(new_est_reads)) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = individual,
    values_from = total_reads,
    id_cols = name,
    values_fill = 0
  )

### Done!
write_tsv(bracken_count_setup, "bracken_output.tsv")

