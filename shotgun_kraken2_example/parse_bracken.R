#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(optparse)
})

### utilities

wide_to_mtx <- function(tbl) {
  mtx <- data.matrix(tbl[,-1])
  rownames(mtx) <- tbl[[1]]
  mtx
}

min_nonzero <- function(x) min(x[x>0])

### parse arguments

opt_list <- list(
  make_option(c("-i", "--input"), type="character", help="path to input directory with bracken output files"),
  make_option(c("-t", "--taxonomy"), type="character", help="path to Phylogenize2 taxonomy .csv table (for making species names match species IDs)"),
  make_option(c("-m", "--metadata"), type="character", help="optional path to a metadata file that describes which runs should be merged (must have columns 'run' and 'sample')"),
  make_option(c("-o", "--output"), type="character", help="path to output .tsv file"),
  make_option(c("-u", "--underscores"), type="logical", default=TRUE, help="remove underscores from taxonomy file when matching")
)
prs <- OptionParser(option_list = opt_list)
p <- parse_args(prs)

### script

# Read in and process raw data
message("Reading raw data...")

#me    taxonomy_id     taxonomy_lvl    kraken_assigned_reads   added_reads     new_est_reads   fraction_total_reads
#phocaeicola dorei       3484    S       1077476 187677  1265153 0.03612

bracken_files <- list.files(path = p$input, pattern = ".*bracken$")
if (length(bracken_files) <= 0) { stop("Error: bracken files not found") }

bracken_data <- map(bracken_files, ~ {
  read_tsv(file.path(p$input, .), col_types="cccdddd")
}, .progress=TRUE)

bracken_run_names <- gsub("\\.bracken", "", bracken_files)
bracken_data_newcol <- map2(bracken_data, bracken_run_names, ~ {
  mutate(.x, run = .y)
}, .progress=TRUE)
bracken_data_tidy <- bind_rows(bracken_data_newcol) %>% ungroup()

if (p$underscores) {
  bracken_data_tidy <- bracken_data_tidy %>%
    mutate(name = gsub("_", " ", name))
}

# Import species taxonomy from taxonomy.csv file
# Note that in the Kraken2 databases, species are named by the representative genome (the cluster ID in Phylogenize2) if they have no name of their own
message("Reading taxonomy data...")

taxonomy <- read_csv(p$taxonomy, col_types="c") %>%
  mutate(species_name = map2_chr(species, cluster, ~ {
    if (is.na(.x) || (.x == "")) {
      .y
    } else {
      if (p$underscores) gsub("_", " ", .x) else .x
    }  
  }))

# Show which names, if any, do not match

# message("Names found only in taxonomy file but not Bracken data: ")
# message(paste(setdiff(taxonomy$species_name, bracken_data_tidy$name %>% unique), sep=", ", collapse=", "))
warning("Names found only in Bracken data but not taxonomy file: ")
warning(paste(setdiff(bracken_data_tidy$name %>% unique, taxonomy$species_name), sep=", ", collapse=", "))

# Load in metadata, if provided

if (("metadata" %in% names(p)) && (!(is.na(p$metadata) || (p$metadata == "")))) {
  message("Reading metadata...")
  metadata <- read_tsv(p$metadata, col_types="c")
  bracken_data_tidy <- bracken_data_tidy %>% 
    inner_join(metadata, by = "run") %>%
    select(name, run, sample, new_est_reads) %>%
    group_by(name, sample) %>%
    summarize(total_reads = sum(new_est_reads), .groups="keep") %>% 
    ungroup() %>%
    rename(individual = sample)
} else {
  bracken_data_tidy <- bracken_data_tidy %>% rename(individual = run, total_reads = new_est_reads)
}

# Now pivot wider
message("Correcting species IDs and pivoting data to table...")

bracken_count_table <- bracken_data_tidy %>% 
  pivot_wider(
    names_from = individual,
    values_from = total_reads,
    id_cols = name,
    values_fill = 0
  ) %>%
  inner_join(., select(taxonomy, name=species_name, species_id=cluster), by="name") %>%
  select(-name) %>%
  relocate(species_id)

write_tsv(bracken_count_table, p$output)

