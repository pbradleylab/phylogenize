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

legacy_harmonize_abd_meta <- function(abd.meta, ...) {
    opts <- settings::clone_and_merge(phylogenize:::PZ_OPTIONS, ...)
    align.metadata.to.matrix <- function(abd.meta) {
        metadata_order <- match(colnames(abd.meta$mtx),
                                abd.meta$metadata[[opts("sample_column")]])
        abd.meta$metadata <- abd.meta$metadata[metadata_order, , drop = FALSE]
        abd.meta
    }
    abd.meta$metadata[[opts("sample_column")]] <- trimws(
        abd.meta$metadata[[opts("sample_column")]]
    )
    colnames(abd.meta$mtx) <- trimws(colnames(abd.meta$mtx))
    samples.present <- intersect(
        abd.meta$metadata[[opts("sample_column")]],
        colnames(abd.meta$mtx)
    )
    if (length(samples.present) == 0) {
        stop("Legacy harmonizer found no shared samples")
    }
    abd.meta$mtx <- abd.meta$mtx[, samples.present, drop = FALSE]
    abd.meta$metadata <- abd.meta$metadata[
        abd.meta$metadata[[opts("sample_column")]] %in% samples.present,
        ,
        drop = FALSE
    ]
    abd.meta <- align.metadata.to.matrix(abd.meta)

    if (opts("which_phenotype") %in% c("specificity", "prevalence", "abundance")) {
        all.envs <- unique(abd.meta$metadata[[opts("env_column")]])
        env.number <- sapply(all.envs, function(e) {
            sum(abd.meta$metadata[[opts("env_column")]] == e)
        })
        names(env.number) <- all.envs
        nonsingleton.envs <- names(which(env.number > 1))
    } else {
        nonsingleton.envs <- NULL
    }

    all.dsets <- unique(abd.meta$metadata[[opts("dset_column")]])
    dset.number <- sapply(all.dsets, function(d) {
        sum(abd.meta$metadata[[opts("dset_column")]] == d)
    })
    names(dset.number) <- all.dsets
    nonsingleton.dsets <- names(which(dset.number > 1))
    f_dsets <- abd.meta$metadata[[opts("dset_column")]] %in% nonsingleton.dsets
    if (!(opts("which_phenotype") %in% c("correlation", "provided"))) {
        f_envs <- abd.meta$metadata[[opts("env_column")]] %in% nonsingleton.envs
        wrows <- which(f_envs & f_dsets)
    } else {
        wrows <- which(f_dsets)
    }
    if (length(wrows) < 2) {
        stop("Legacy harmonizer retained fewer than two metadata rows")
    }

    abd.meta$metadata <- abd.meta$metadata[wrows, , drop = FALSE]
    wcols <- intersect(colnames(abd.meta$mtx),
                       abd.meta$metadata[[opts("sample_column")]])
    if (length(wcols) < 2) {
        stop("Legacy harmonizer retained fewer than two matrix columns")
    }
    abd.meta$mtx <- abd.meta$mtx[, wcols, drop = FALSE]
    abd.meta <- align.metadata.to.matrix(abd.meta)
    abd.meta$mtx <- remove.allzero.abundances(abd.meta$mtx)
    align.metadata.to.matrix(abd.meta)
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

cat("Real-data harmonization validation using poms_cirrhosis.R inputs\n")
cat("  Abundance:", abundance_file, "\n")
cat("  Metadata:", metadata_file, "\n\n")

read_time <- system.time({
    abd.meta <- phylogenize:::read.abd.metadata.tabular()
})[["elapsed"]]

cat("Read input:\n")
cat("  Elapsed:", read_time, "sec\n")
cat("  Matrix dimensions before harmonization:",
    paste(dim(abd.meta$mtx), collapse = " x "), "\n")
cat("  Metadata rows before harmonization:", nrow(abd.meta$metadata), "\n\n")

legacy_time <- system.time({
    legacy <- legacy_harmonize_abd_meta(abd.meta)
})[["elapsed"]]

new_time <- system.time({
    current <- harmonize.abd.meta(abd.meta)
})[["elapsed"]]

same_metadata <- identical(current$metadata, legacy$metadata)
same_colnames <- identical(colnames(current$mtx), colnames(legacy$mtx))
same_rownames <- identical(rownames(current$mtx), rownames(legacy$mtx))
same_matrix <- identical(as.matrix(current$mtx), as.matrix(legacy$mtx))

cat("Runtime:\n")
cat("  Legacy harmonize elapsed:", legacy_time, "sec\n")
cat("  Current harmonize elapsed:", new_time, "sec\n")
cat("  Speedup:", round(legacy_time / new_time, 2), "x\n\n")

cat("Output equivalence:\n")
cat("  Metadata identical:", same_metadata, "\n")
cat("  Matrix colnames identical:", same_colnames, "\n")
cat("  Matrix rownames identical:", same_rownames, "\n")
cat("  Matrix values identical:", same_matrix, "\n")
cat("  Final matrix dimensions:", paste(dim(current$mtx), collapse = " x "), "\n")
cat("  Final metadata rows:", nrow(current$metadata), "\n")

if (!all(same_metadata, same_colnames, same_rownames, same_matrix)) {
    stop("Current harmonization output differs from legacy output")
}

cat("\nPASS: current harmonization matches legacy output on real cirrhosis data.\n")
