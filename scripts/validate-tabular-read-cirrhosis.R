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

legacy_read_abd_metadata_tabular <- function(...) {
    opts <- settings::clone_and_merge(phylogenize:::PZ_OPTIONS, ...)
    af <- opts("abundance_file")
    mf <- opts("metadata_file")

    abd.df <- readr::read_tsv(
        af,
        col_types = readr::cols(.default = readr::col_character()),
        show_col_types = FALSE
    )
    abd.values <- as.data.frame(
        lapply(abd.df[, -1, drop = FALSE], function(x) {
            suppressWarnings(as.numeric(x))
        }),
        check.names = FALSE
    )
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

old_opts <- pz.options()
on.exit(do.call(pz.options, old_opts), add = TRUE)
pz.options(
    abundance_file = abundance_file,
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

cat("Real-data tabular read validation using poms_cirrhosis.R inputs\n")
cat("  Abundance:", abundance_file, "\n")
cat("  Metadata:", metadata_file, "\n\n")

legacy_time <- system.time({
    legacy <- legacy_read_abd_metadata_tabular()
})[["elapsed"]]

current_time <- system.time({
    current <- phylogenize:::read.abd.metadata.tabular()
})[["elapsed"]]

same_metadata_strict <- identical(current$metadata, legacy$metadata)
same_metadata_values <- metadata_values_identical(current$metadata, legacy$metadata)
same_colnames <- identical(colnames(current$mtx), colnames(legacy$mtx))
same_rownames <- identical(rownames(current$mtx), rownames(legacy$mtx))
same_matrix <- identical(current$mtx, legacy$mtx)

cat("Read runtime:\n")
cat("  Legacy character-read + as.numeric elapsed:", legacy_time, "sec\n")
cat("  Current direct numeric-read elapsed:", current_time, "sec\n")
cat("  Speedup:", round(legacy_time / current_time, 2), "x\n\n")

cat("Raw read output equivalence:\n")
cat("  Metadata strictly identical:", same_metadata_strict, "\n")
cat("  Metadata values identical:", same_metadata_values, "\n")
cat("  Matrix colnames identical:", same_colnames, "\n")
cat("  Matrix rownames identical:", same_rownames, "\n")
cat("  Matrix values identical:", same_matrix, "\n")
cat("  Matrix dimensions:", paste(dim(current$mtx), collapse = " x "), "\n")
cat("  Metadata rows:", nrow(current$metadata), "\n\n")

if (!same_metadata_values) {
    cat("Metadata all.equal output:\n")
    print(all.equal(current$metadata, legacy$metadata))
    cat("\n")
}

if (!all(same_metadata_values, same_colnames, same_rownames, same_matrix)) {
    stop("Current tabular reader output differs from legacy output")
}

legacy_harmonized <- harmonize.abd.meta(legacy)
current_harmonized <- harmonize.abd.meta(current)

same_harmonized_metadata_strict <- identical(
    current_harmonized$metadata,
    legacy_harmonized$metadata
)
same_harmonized_metadata_values <- metadata_values_identical(
    current_harmonized$metadata,
    legacy_harmonized$metadata
)
same_harmonized_colnames <- identical(
    colnames(current_harmonized$mtx),
    colnames(legacy_harmonized$mtx)
)
same_harmonized_rownames <- identical(
    rownames(current_harmonized$mtx),
    rownames(legacy_harmonized$mtx)
)
same_harmonized_matrix <- identical(
    current_harmonized$mtx,
    legacy_harmonized$mtx
)

cat("Post-harmonization equivalence:\n")
cat("  Metadata strictly identical:", same_harmonized_metadata_strict, "\n")
cat("  Metadata values identical:", same_harmonized_metadata_values, "\n")
cat("  Matrix colnames identical:", same_harmonized_colnames, "\n")
cat("  Matrix rownames identical:", same_harmonized_rownames, "\n")
cat("  Matrix values identical:", same_harmonized_matrix, "\n")
cat("  Harmonized matrix dimensions:",
    paste(dim(current_harmonized$mtx), collapse = " x "), "\n")
cat("  Harmonized metadata rows:", nrow(current_harmonized$metadata), "\n")

if (!all(
    same_harmonized_metadata_values,
    same_harmonized_colnames,
    same_harmonized_rownames,
    same_harmonized_matrix
)) {
    stop("Current post-harmonization output differs from legacy output")
}

cat("\nPASS: current direct numeric reader matches legacy output on real cirrhosis data.\n")
