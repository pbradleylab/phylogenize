devtools::load_all("/fs/project/bradley.720/users/kananen.13/tools/phylogenize/package/phylogenize/")

out_dir <- file.path("output", "bug_wrong_overlap_size")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pz.options(out_dir = out_dir)

cirrhosis_family_abundance <- phylogenize_core(
  db = "globdb",
  taxon_level = "family",
  type_16S = FALSE,
  data_dir = "/fs/project/bradley.720/users/kananen.13/tools/phylogenize/package/phylogenize/inst/extdata/",
  which_phenotype = "abundance",
  diff_abund_method = "ancombc2",
  which_envir = "control",
  abundance_file = "/fs/project/bradley.720/users/kananen.13/analysis/260427/results/sylph/matrices/sylph_est_read_counts_genomes_by_individuals.tsv",
  metadata_file = "/fs/project/bradley.720/users/kananen.13/analysis/260427/case_control_with_dataset.csv",
  input_format = "tabular",
  sample_column = "sampleid",
  env_column = "env",
  ncl = 4
)

saveRDS(cirrhosis_family_abundance, file.path(out_dir, "core.rds"))

overlap_file <- file.path(out_dir, "enr-overlaps.csv")
if (!file.exists(overlap_file)) {
  stop("Expected enrichment overlap file was not created: ", overlap_file)
}

overlaps <- read.csv(overlap_file, check.names = FALSE)
if (!("effectsize" %in% colnames(overlaps))) {
  stop("enr-overlaps.csv is missing the effectsize column")
}

if (nrow(overlaps) == 0) {
  stop("enr-overlaps.csv was created but contains no rows")
}

matched_overlaps <- overlaps[!is.na(overlaps$effectsize), , drop = FALSE]
if (nrow(matched_overlaps) == 0) {
  stop("All enrichment overlap effectsize values are NA")
}

results_matrix <- cirrhosis_family_abundance$list_signif$results.matrix
check_rows <- utils::head(matched_overlaps, 25)

for (i in seq_len(nrow(check_rows))) {
  row <- check_rows[i, , drop = FALSE]
  expected <- results_matrix |>
    dplyr::filter(taxon == row$taxon, gene == row$gene) |>
    dplyr::pull(effect.size)

  if (length(expected) == 0) {
    stop("No matching result found for taxon=", row$taxon, ", gene=", row$gene)
  }

  if (!isTRUE(all.equal(row$effectsize, expected[1]))) {
    stop(
      "Effect size mismatch for taxon=", row$taxon,
      ", gene=", row$gene,
      ": overlap=", row$effectsize,
      ", expected=", expected[1]
    )
  }
}

message(
  "PASS: checked ", nrow(check_rows),
  " non-NA enrichment overlap effect sizes against results.matrix"
)
