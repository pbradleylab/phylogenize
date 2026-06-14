test_that("null association p-values are not inflated under independent data", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    set.seed(20240614)

    n_taxa <- 80
    n_genes <- 120
    taxa <- paste0("taxon", seq_len(n_taxa))
    genes <- paste0("gene", seq_len(n_genes))
    phenotype <- stats::rnorm(n_taxa)
    names(phenotype) <- taxa

    gene_rows <- replicate(
        n_genes,
        sample(rep(c(FALSE, TRUE), each=n_taxa / 2)),
        simplify=FALSE
    )
    gene_presence <- do.call(rbind, gene_rows)
    rownames(gene_presence) <- genes
    colnames(gene_presence) <- taxa
    tree <- ape::rtree(n_taxa, tip.label=taxa)
    opts <- pz.resolve.options(
        core_method="lm",
        separate_process=FALSE,
        meas_err=FALSE,
        error_to_file=FALSE
    )

    result <- matrix.plm(
        tree,
        gene_presence,
        phenotype,
        method=lm.fx.pv,
        restrict.taxa=taxa,
        restrict.ff=genes,
        .opts=opts
    )
    p_values <- as.numeric(result["p.value", ])
    q_values <- p.adjust(p_values, method="BH")

    expect_true(all(is.finite(p_values)))
    expect_true(all(p_values >= 0 & p_values <= 1))
    expect_lte(mean(p_values < 0.05), 0.10)
    expect_gt(stats::ks.test(p_values, "punif")$p.value, 0.01)
    expect_equal(sum(q_values < 0.05), 0)
})

test_that("association statistics are invariant to tree and matrix ordering", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    set.seed(20240615)

    n_taxa <- 40
    n_genes <- 12
    taxa <- paste0("taxon", seq_len(n_taxa))
    genes <- paste0("gene", seq_len(n_genes))
    phenotype <- stats::rnorm(n_taxa)
    names(phenotype) <- taxa
    gene_presence <- matrix(
        stats::rbinom(n_taxa * n_genes, 1, 0.45),
        nrow=n_genes,
        dimnames=list(genes, taxa)
    )
    keep <- rowSums(gene_presence) > 3 &
        rowSums(gene_presence) < (n_taxa - 3)
    gene_presence <- gene_presence[keep, , drop=FALSE]
    tree <- ape::rtree(n_taxa, tip.label=taxa)
    opts <- pz.resolve.options(
        core_method="lm",
        separate_process=FALSE,
        meas_err=FALSE,
        error_to_file=FALSE
    )
    run_assoc <- function(tr, genes_by_taxa, pheno) {
        matrix.plm(
            tr,
            genes_by_taxa,
            pheno,
            method=lm.fx.pv,
            restrict.taxa=intersect(colnames(genes_by_taxa), tr$tip.label),
            restrict.ff=rownames(genes_by_taxa),
            .opts=opts
        )
    }

    baseline <- run_assoc(tree, gene_presence, phenotype)
    permuted_gene_presence <- gene_presence[
        sample(rownames(gene_presence)),
        sample(colnames(gene_presence)),
        drop=FALSE
    ]
    permuted_phenotype <- phenotype[sample(names(phenotype))]
    permuted_tree <- ape::keep.tip(tree, sample(tree$tip.label))
    permuted <- run_assoc(
        permuted_tree,
        permuted_gene_presence,
        permuted_phenotype
    )

    expect_equal(rownames(permuted), rownames(baseline))
    expect_setequal(colnames(permuted), colnames(baseline))
    expect_equal(
        permuted[, colnames(baseline), drop=FALSE],
        baseline,
        tolerance=1e-12,
        ignore_attr=TRUE
    )
})
