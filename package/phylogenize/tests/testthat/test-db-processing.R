test_that("change.presence.tax.level splits sparse matrices by requested taxon", {
    tax <- tibble::tibble(
        cluster=paste0("s", 1:4),
        phylum="p1",
        family=c("f1", "f1", "f2", "f2")
    )
    binary <- list(
        p1=Matrix::Matrix(
            matrix(
                c(1, 0, 1, 0,
                  0, 1, 0, 1),
                nrow=2,
                byrow=TRUE,
                dimnames=list(c("gene1", "gene2"), paste0("s", 1:4))
            ),
            sparse=TRUE
        )
    )

    expect_warning(
        split_binary <- change.presence.tax.level(binary, "family", tax),
        "Using an external vector"
    )

    expect_setequal(names(split_binary), c("f1", "f2"))
    expect_equal(colnames(split_binary$f1), c("s1", "s2"))
    expect_equal(colnames(split_binary$f2), c("s3", "s4"))
    expect_equal(rownames(split_binary$f1), c("gene1", "gene2"))
})

test_that("change.tree.tax.level returns matching subtrees", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(verbosity=0, error_to_file=FALSE)

    tax <- tibble::tibble(
        cluster=paste0("s", 1:4),
        phylum="p1",
        family=c("f1", "f1", "f2", "f2")
    )
    tree <- list(p1=ape::read.tree(text="((s1:1,s2:1):1,(s3:1,s4:1):1);"))

    expect_warning(
        split_trees <- change.tree.tax.level(tree, "family", tax),
        "Using an external vector"
    )

    expect_setequal(names(split_trees), c("f1", "f2"))
    expect_equal(split_trees$f1$tip.label, c("s1", "s2"))
    expect_equal(split_trees$f2$tip.label, c("s3", "s4"))
})

test_that("above_minimum_genes keeps genes with enough present and absent tips", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(minimum=1, gene_min_frac=0.5, verbosity=0, error_to_file=FALSE)

    gene.presence <- list(
        TaxonA=Matrix::Matrix(
            matrix(
                c(1, 0, 1,
                  1, 1, 1,
                  0, 0, 0),
                nrow=3,
                byrow=TRUE,
                dimnames=list(c("g_keep", "g_all", "g_none"), c("s1", "s2", "s3"))
            ),
            sparse=TRUE
        )
    )
    trees <- list(TaxonA=ape::read.tree(text="(s1:1,s2:1,s3:1);"))

    filtered <- above_minimum_genes(gene.presence, trees)

    expect_equal(names(filtered), "TaxonA")
    expect_equal(rownames(filtered$TaxonA), "g_keep")
    expect_equal(colnames(filtered$TaxonA), c("s1", "s2", "s3"))
})
