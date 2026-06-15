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

test_that("known synthetic signals rank above null genes", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    set.seed(20240616)

    n_taxa <- 80
    n_null <- 30
    taxa <- paste0("taxon", seq_len(n_taxa))
    latent <- seq(-2, 2, length.out=n_taxa)
    phenotype <- latent + stats::rnorm(n_taxa, sd=0.25)
    names(phenotype) <- taxa
    high_phenotype <- phenotype >= stats::median(phenotype)
    null_rows <- replicate(
        n_null,
        sample(rep(c(0, 1), each=n_taxa / 2)),
        simplify=FALSE
    )
    gene_presence <- rbind(
        causal_positive=as.integer(high_phenotype),
        causal_negative=as.integer(!high_phenotype),
        do.call(rbind, null_rows)
    )
    rownames(gene_presence)[-(1:2)] <- paste0("null", seq_len(n_null))
    colnames(gene_presence) <- taxa
    tree <- ape::rtree(n_taxa, tip.label=taxa)
    opts <- pz.resolve.options(
        core_method="lm",
        separate_process=FALSE,
        meas_err=FALSE,
        fdr_method="BH",
        error_to_file=FALSE
    )

    result <- matrix.plm(
        tree,
        gene_presence,
        phenotype,
        method=lm.fx.pv,
        restrict.taxa=taxa,
        restrict.ff=rownames(gene_presence),
        .opts=opts
    )
    p_values <- result["p.value", ]
    q_values <- phylogenize:::qvals(p_values, .opts=opts)
    effects <- result["Estimate", ]
    causal_genes <- c("causal_positive", "causal_negative")
    null_genes <- paste0("null", seq_len(n_null))

    expect_setequal(names(sort(p_values))[1:2], causal_genes)
    expect_gt(effects[["causal_positive"]], 0)
    expect_lt(effects[["causal_negative"]], 0)
    expect_true(all(q_values[causal_genes] < 0.05))
    expect_equal(sum(q_values[null_genes] < 0.05), 0)
    expect_lt(max(q_values[causal_genes]), min(q_values[null_genes]))
})

test_that("deterministic association methods agree on core signal invariants", {
    testthat::skip_if_not_installed("phylolm")
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    set.seed(20240617)

    n_taxa <- 60
    n_null <- 8
    taxa <- paste0("taxon", seq_len(n_taxa))
    latent <- seq(-2, 2, length.out=n_taxa)
    phenotype <- latent + stats::rnorm(n_taxa, sd=0.2)
    names(phenotype) <- taxa
    high_phenotype <- phenotype >= stats::median(phenotype)
    null_rows <- replicate(
        n_null,
        sample(rep(c(0, 1), each=n_taxa / 2)),
        simplify=FALSE
    )
    gene_presence <- rbind(
        causal_positive=as.integer(high_phenotype),
        do.call(rbind, null_rows)
    )
    rownames(gene_presence)[-1] <- paste0("null", seq_len(n_null))
    colnames(gene_presence) <- taxa
    tree <- ape::rtree(n_taxa, tip.label=taxa)
    opts <- pz.resolve.options(
        separate_process=FALSE,
        meas_err=FALSE,
        error_to_file=FALSE
    )
    methods <- list(
        lm=lm.fx.pv,
        phylolm=phylolm.fx.pv
    )

    results <- lapply(methods, function(method) {
        matrix.plm(
            tree,
            gene_presence,
            phenotype,
            method=method,
            restrict.taxa=taxa,
            restrict.ff=rownames(gene_presence),
            .opts=opts
        )
    })

    expect_equal(lapply(results, rownames), stats::setNames(
        rep(list(c("Estimate", "p.value", "StdErr", "df")),
            length(results)),
        names(results)
    ))
    expect_equal(lapply(results, colnames), stats::setNames(
        rep(list(rownames(gene_presence)), length(results)),
        names(results)
    ))
    for (method_name in names(results)) {
        result <- results[[method_name]]
        expect_true(all(is.finite(result[, "causal_positive"])))
        expect_gt(result["Estimate", "causal_positive"], 0)
        expect_lt(result["p.value", "causal_positive"], 0.05)
        expect_equal(names(sort(result["p.value", ]))[1], "causal_positive")
    }
})
