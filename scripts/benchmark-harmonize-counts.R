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
if (!dir.exists(pkg_dir)) {
    stop("Could not find package directory at: ", pkg_dir)
}

suppressMessages(
    suppressPackageStartupMessages(devtools::load_all(pkg_dir, quiet = TRUE))
)

legacy_counts <- function(metadata, column) {
    values <- unique(metadata[[column]])
    counts <- sapply(values, function(value) {
        sum(metadata[[column]] == value)
    })
    names(counts) <- values
    counts
}

table_counts <- function(metadata, column) {
    table(metadata[[column]])
}

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

    abd.meta$metadata <- abd.meta$metadata[wrows, , drop = FALSE]
    wcols <- intersect(colnames(abd.meta$mtx),
                       abd.meta$metadata[[opts("sample_column")]])
    abd.meta$mtx <- abd.meta$mtx[, wcols, drop = FALSE]
    abd.meta <- align.metadata.to.matrix(abd.meta)
    abd.meta$mtx <- remove.allzero.abundances(abd.meta$mtx)
    align.metadata.to.matrix(abd.meta)
}

make_abd_meta <- function(n_samples = 250000,
                          n_taxa = 20,
                          n_envs = 200,
                          n_dsets = 500) {
    set.seed(20260601)
    samples <- paste0("sample", seq_len(n_samples))
    env <- paste0("env", ((seq_len(n_samples) - 1L) %% n_envs) + 1L)
    dataset <- paste0("dataset", ((seq_len(n_samples) - 1L) %% n_dsets) + 1L)

    # Include singleton categories so the benchmark exercises filtering behavior.
    env[1] <- "singleton_env"
    dataset[2] <- "singleton_dataset"

    metadata <- data.frame(
        sample = samples,
        env = env,
        dataset = dataset,
        stringsAsFactors = FALSE
    )

    mtx <- matrix(1, nrow = n_taxa, ncol = n_samples)
    rownames(mtx) <- paste0("taxon", seq_len(n_taxa))
    colnames(mtx) <- samples

    list(mtx = Matrix::Matrix(mtx, sparse = TRUE), metadata = metadata)
}

time_it <- function(expr) {
    gc()
    elapsed <- system.time(force(expr))[["elapsed"]]
    gc()
    elapsed
}

same_named_counts <- function(old, new) {
    old <- old[order(names(old))]
    new <- new[order(names(new))]
    identical(names(old), names(new)) &&
        identical(unname(as.integer(old)), unname(as.integer(new)))
}

abd.meta <- make_abd_meta()

cat("Synthetic input:\n")
cat("  Samples:", ncol(abd.meta$mtx), "\n")
cat("  Taxa:", nrow(abd.meta$mtx), "\n")
cat("  Environments:", length(unique(abd.meta$metadata$env)), "\n")
cat("  Datasets:", length(unique(abd.meta$metadata$dataset)), "\n\n")

env_old <- legacy_counts(abd.meta$metadata, "env")
env_new <- table_counts(abd.meta$metadata, "env")
dset_old <- legacy_counts(abd.meta$metadata, "dataset")
dset_new <- table_counts(abd.meta$metadata, "dataset")

stopifnot(same_named_counts(env_old, env_new))
stopifnot(same_named_counts(dset_old, dset_new))

cat("Count correctness:\n")
cat("  Environment counts identical: TRUE\n")
cat("  Dataset counts identical: TRUE\n\n")

n_iter <- 5
legacy_time <- time_it({
    for (i in seq_len(n_iter)) {
        legacy_counts(abd.meta$metadata, "env")
        legacy_counts(abd.meta$metadata, "dataset")
    }
})
table_time <- time_it({
    for (i in seq_len(n_iter)) {
        table_counts(abd.meta$metadata, "env")
        table_counts(abd.meta$metadata, "dataset")
    }
})

cat("Counting benchmark over", n_iter, "iteration(s):\n")
cat("  Legacy unique+sapply+sum elapsed:", legacy_time, "sec\n")
cat("  New table() elapsed:", table_time, "sec\n")
cat("  Speedup:", round(legacy_time / table_time, 2), "x\n\n")

old_opts <- pz.options()
on.exit(do.call(pz.options, old_opts), add = TRUE)
pz.options(
    sample_column = "sample",
    env_column = "env",
    dset_column = "dataset",
    which_phenotype = "prevalence",
    error_to_file = FALSE,
    verbosity = 0
)

gc()
legacy_harmonize_time <- system.time({
    legacy_harmonized <- legacy_harmonize_abd_meta(abd.meta)
})[["elapsed"]]
gc()
new_harmonize_time <- system.time({
    harmonized <- harmonize.abd.meta(abd.meta)
})[["elapsed"]]
gc()

stopifnot(identical(colnames(harmonized$mtx), colnames(legacy_harmonized$mtx)))
stopifnot(identical(rownames(harmonized$mtx), rownames(legacy_harmonized$mtx)))
stopifnot(identical(harmonized$metadata, legacy_harmonized$metadata))
stopifnot(identical(as.matrix(harmonized$mtx), as.matrix(legacy_harmonized$mtx)))

cat("Ending output check:\n")
cat("  Legacy harmonize elapsed:", legacy_harmonize_time, "sec\n")
cat("  New harmonize elapsed:", new_harmonize_time, "sec\n")
cat("  Matrix dimensions:", paste(dim(harmonized$mtx), collapse = " x "), "\n")
cat("  Metadata rows:", nrow(harmonized$metadata), "\n")
cat("  Harmonized output is identical to legacy counting output: TRUE\n")
