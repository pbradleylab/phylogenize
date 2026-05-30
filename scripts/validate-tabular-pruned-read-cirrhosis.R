#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grepl(file_arg, args)])
repo_root <- if (length(script_path) > 0) {
    normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
    normalizePath(".", mustWork = TRUE)
}

pkg_dir <- file.path(repo_root, "package", "phylogenize")
abundance_file <- file.path(
    repo_root,
    "comms",
    "phylogenize2",
    "cirrhosis_uhgg2.tsv"
)
metadata_file <- file.path(
    repo_root,
    "comms",
    "phylogenize2",
    "cirrhosis_uhgg2_metadata.tsv"
)

if (!dir.exists(pkg_dir)) {
    stop("Could not find package directory at: ", pkg_dir)
}
if (!file.exists(abundance_file)) {
    stop("Could not find abundance file at: ", abundance_file)
}
if (!file.exists(metadata_file)) {
    stop("Could not find metadata file at: ", metadata_file)
}

suppressMessages(
    suppressPackageStartupMessages(devtools::load_all(pkg_dir, quiet = TRUE))
)

metadata_values_identical <- function(x, y) {
    x <- as.data.frame(x)
    y <- as.data.frame(y)
    rownames(x) <- NULL
    rownames(y) <- NULL
    identical(x, y)
}

legacy_unpruned_read_abd_metadata_tabular <- function(...) {
    opts <- settings::clone_and_merge(phylogenize:::PZ_OPTIONS, ...)
    af <- opts("abundance_file")
    mf <- opts("metadata_file")

    abundance_header <- names(readr::read_tsv(
        af,
        n_max = 0,
        show_col_types = FALSE
    ))
    if (length(abundance_header) < 2) {
        stop("Abundance table must contain a taxon column and at least one sample column")
    }
    abundance_col_types <- do.call(
        readr::cols,
        c(
            stats::setNames(list(readr::col_character()), abundance_header[1]),
            list(.default = readr::col_double())
        )
    )
    abd.df <- readr::read_tsv(
        af,
        col_types = abundance_col_types,
        show_col_types = FALSE
    )
    abd.values <- abd.df[, -1, drop = FALSE]
    bad.cols <- names(abd.values)[
        vapply(abd.values, function(x) any(is.na(x)), logical(1))
    ]
    if (length(bad.cols) > 0) {
        stop(
            "Legacy reader found nonnumeric value(s) in column(s): ",
            paste(bad.cols, collapse = ", ")
        )
    }
    abd.mtx <- as.matrix(abd.values)
    rownames(abd.mtx) <- abd.df[[1]]

    metadata <- readr::read_tsv(mf, show_col_types = FALSE)
    metadata <- check.process.metadata(metadata, ...)
    list(mtx = abd.mtx, metadata = metadata)
}

make_augmented_abundance <- function(input_file, n_extra = 500) {
    cat("Preparing augmented abundance file with", n_extra,
        "unused sample column(s)...\n")
    abd <- readr::read_tsv(
        input_file,
        col_types = readr::cols(.default = readr::col_character()),
        show_col_types = FALSE
    )
    sample_cols <- names(abd)[-1]
    if (length(sample_cols) == 0) {
        stop("Abundance file has no sample columns")
    }
    source_cols <- sample_cols[((seq_len(n_extra) - 1L) %% length(sample_cols)) + 1L]
    for (i in seq_len(n_extra)) {
        abd[[paste0("unused_extra_sample_", i)]] <- abd[[source_cols[i]]]
    }
    out_file <- file.path(tempdir(), "cirrhosis_uhgg2_with_unused_samples.tsv")
    readr::write_tsv(abd, out_file)
    out_file
}

extra_arg <- commandArgs(trailingOnly = TRUE)
n_extra <- if (length(extra_arg) > 0) as.integer(extra_arg[1]) else 500L
if (is.na(n_extra) || n_extra < 1) {
    stop("Optional argument n_extra must be a positive integer")
}

augmented_abundance_file <- make_augmented_abundance(abundance_file, n_extra)

old_opts <- pz.options()
on.exit(do.call(pz.options, old_opts), add = TRUE)
pz.options(
    abundance_file = augmented_abundance_file,
    metadata_file = metadata_file,
    input_format = "tabular",
    db = "human-gut",
    taxon_level = "family",
    core_method = "poms",
    which_envir = "case",
    sample_column = "sample",
    env_column = "env",
    dset_column = "dataset",
    min_fx = 0,
    ncl = 4,
    error_to_file = FALSE,
    verbosity = 0
)

cat("\nReal-data tabular read pruning validation using poms_cirrhosis.R inputs\n")
cat("  Base abundance:", abundance_file, "\n")
cat("  Augmented abundance:", augmented_abundance_file, "\n")
cat("  Metadata:", metadata_file, "\n")
cat("  Extra unused sample columns:", n_extra, "\n\n")

legacy_time <- system.time({
    legacy <- legacy_unpruned_read_abd_metadata_tabular()
})[["elapsed"]]

current_time <- system.time({
    current <- phylogenize:::read.abd.metadata.tabular()
})[["elapsed"]]

cat("Raw read dimensions:\n")
cat("  Legacy unpruned matrix:",
    paste(dim(legacy$mtx), collapse = " x "), "\n")
cat("  Current pruned matrix:",
    paste(dim(current$mtx), collapse = " x "), "\n\n")

cat("Read runtime:\n")
cat("  Legacy unpruned read elapsed:", legacy_time, "sec\n")
cat("  Current metadata-pruned read elapsed:", current_time, "sec\n")
cat("  Speedup:", round(legacy_time / current_time, 2), "x\n\n")

legacy_harmonized_time <- system.time({
    legacy_harmonized <- harmonize.abd.meta(legacy)
})[["elapsed"]]
current_harmonized_time <- system.time({
    current_harmonized <- harmonize.abd.meta(current)
})[["elapsed"]]

same_metadata <- metadata_values_identical(
    current_harmonized$metadata,
    legacy_harmonized$metadata
)
same_colnames <- identical(
    colnames(current_harmonized$mtx),
    colnames(legacy_harmonized$mtx)
)
same_rownames <- identical(
    rownames(current_harmonized$mtx),
    rownames(legacy_harmonized$mtx)
)
same_matrix <- identical(
    current_harmonized$mtx,
    legacy_harmonized$mtx
)

cat("Post-harmonization runtime:\n")
cat("  Legacy harmonize elapsed:", legacy_harmonized_time, "sec\n")
cat("  Current harmonize elapsed:", current_harmonized_time, "sec\n\n")

cat("Post-harmonization output equivalence:\n")
cat("  Metadata values identical:", same_metadata, "\n")
cat("  Matrix colnames identical:", same_colnames, "\n")
cat("  Matrix rownames identical:", same_rownames, "\n")
cat("  Matrix values identical:", same_matrix, "\n")
cat("  Final matrix dimensions:",
    paste(dim(current_harmonized$mtx), collapse = " x "), "\n")
cat("  Final metadata rows:", nrow(current_harmonized$metadata), "\n")

if (!all(same_metadata, same_colnames, same_rownames, same_matrix)) {
    stop("Current pruned reader changes post-harmonization output")
}

cat("\nPASS: metadata-pruned reader preserves final output on augmented real cirrhosis data.\n")
