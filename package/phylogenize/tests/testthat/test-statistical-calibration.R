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
