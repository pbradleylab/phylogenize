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

legacy_change_tree_tax_level <- function(tree, taxon, tax) {
    clean <- tax %>%
        dplyr::select(cluster, tidyselect::all_of(taxon), phylum) %>%
        dplyr::distinct()
    clean <- clean[!(is.na(clean[[taxon]]) | clean[[taxon]] == ""), ]
    clean <- clean %>%
        dplyr::group_by(
            dplyr::across(tidyselect::all_of(taxon))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(phylum, !!(rlang::sym(taxon))) %>%
        dplyr::mutate(cluster = as.character(cluster))
    clean <- clean %>%
        split(.$phylum)

    tree_matrices <- list()
    tree_names <- names(tree)
    for (i in seq_along(tree)) {
        name <- tree_names[i]
        tr <- tree[[name]]
        t <- clean[[name]]

        if (is.null(t)) {
            next
        }

        split_names <- t %>%
            dplyr::group_split(!!(rlang::sym(taxon))) %>%
            purrr::map(~ dplyr::pull(.x, !!(rlang::sym(taxon)))) %>%
            unlist() %>%
            unique()

        for (j in seq_along(split_names)) {
            tips <- tr$tip.label
            split_tips <- t %>%
                dplyr::filter(!!(rlang::sym(taxon)) == split_names[[j]])
            tips <- intersect(tips, split_tips[["cluster"]])
            subtree <- ape::keep.tip(tr, tips)
            if (!is.null(subtree) && length(subtree$tip.label) > 1) {
                tree_matrices[[split_names[j]]] <- subtree
            }
        }
    }
    Filter(function(tr) length(tr$tip.label) > 1, tree_matrices)
}

make_tree_taxonomy_case <- function(n_phyla = 8,
                                    families_per_phylum = 150,
                                    species_per_family = 8) {
    set.seed(20260602)
    tax_parts <- vector("list", n_phyla)
    tree <- vector("list", n_phyla)
    names(tree) <- paste0("phylum_", seq_len(n_phyla))

    for (i in seq_len(n_phyla)) {
        phylum <- names(tree)[i]
        families <- paste0(phylum, "_family_", seq_len(families_per_phylum))
        clusters <- paste0(
            phylum,
            "_s",
            seq_len(families_per_phylum * species_per_family)
        )
        tax_parts[[i]] <- data.frame(
            cluster = clusters,
            phylum = phylum,
            family = rep(families, each = species_per_family),
            stringsAsFactors = FALSE
        )
        tree[[i]] <- ape::rtree(length(clusters), tip.label = clusters)
    }

    list(tree = tree, tax = do.call(rbind, tax_parts))
}

canonical_tree_output <- function(x) {
    lapply(x[sort(names(x))], function(tr) sort(tr$tip.label))
}

trailing <- commandArgs(trailingOnly = TRUE)
n_phyla <- if (length(trailing) >= 1) as.integer(trailing[1]) else 8L
families_per_phylum <- if (length(trailing) >= 2) as.integer(trailing[2]) else 150L
species_per_family <- if (length(trailing) >= 3) as.integer(trailing[3]) else 8L
if (any(is.na(c(n_phyla, families_per_phylum, species_per_family))) ||
    any(c(n_phyla, families_per_phylum, species_per_family) < 1)) {
    stop("Arguments must be positive integers: n_phyla families_per_phylum species_per_family")
}

old_opts <- pz.options()
on.exit(do.call(pz.options, old_opts), add = TRUE)
pz.options(error_to_file = FALSE, verbosity = 0)

case <- make_tree_taxonomy_case(
    n_phyla = n_phyla,
    families_per_phylum = families_per_phylum,
    species_per_family = species_per_family
)

cat("Synthetic change.tree.tax.level validation\n")
cat("  Phyla:", n_phyla, "\n")
cat("  Families per phylum:", families_per_phylum, "\n")
cat("  Species per family:", species_per_family, "\n")
cat("  Input trees:", length(case$tree), "\n")
cat("  Input taxonomy rows:", nrow(case$tax), "\n\n")

legacy_time <- system.time({
    legacy <- legacy_change_tree_tax_level(case$tree, "family", case$tax)
})[["elapsed"]]

current_time <- system.time({
    current <- change.tree.tax.level(case$tree, "family", case$tax)
})[["elapsed"]]

legacy_canonical <- canonical_tree_output(legacy)
current_canonical <- canonical_tree_output(current)
same_names <- identical(names(current_canonical), names(legacy_canonical))
same_tips <- identical(current_canonical, legacy_canonical)

cat("Runtime:\n")
cat("  Legacy repeated-filter elapsed:", legacy_time, "sec\n")
cat("  Current precomputed-split elapsed:", current_time, "sec\n")
cat("  Speedup:", round(legacy_time / current_time, 2), "x\n\n")

cat("Output equivalence:\n")
cat("  Subtree names identical:", same_names, "\n")
cat("  Subtree tip labels identical:", same_tips, "\n")
cat("  Output subtree count:", length(current), "\n")

if (!all(same_names, same_tips)) {
    stop("Current change.tree.tax.level output differs from legacy output")
}

cat("\nPASS: current change.tree.tax.level matches legacy output.\n")
