test_that("enrichment overlap exports include effect sizes", {
    out_dir <- file.path(tempdir(), "phylogenize-enrichment-export-test")
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(out_dir=out_dir, error_to_file=FALSE)

    testthat::local_mocked_bindings(
        multi.kegg.enrich=function(...) {
            tibble::tibble(
                taxon="TaxonA",
                cutoff="strong",
                ID="ko00001",
                Description="Mock pathway",
                qvalue=0.01,
                enr.estimate=2,
                geneID="K00001"
            )
        },
        .package="phylogenize"
    )

    pz.db <- list(
        gene.to.fxn=tibble::tibble(
            gene="gene_1",
            accession="K00001",
            `function`="Mock gene"
        )
    )
    results.matrix <- tibble::tibble(
        taxon="TaxonA",
        gene="gene_1",
        effect.size=1.5,
        p.value=0.001
    )
    signif <- list(TaxonA=list(strong="gene_1", med="gene_1", weak="gene_1"))
    signs <- list(TaxonA=c(gene_1=1))

    get_enrichment_tbls(
        signif,
        signs,
        pz.db,
        results.matrix,
        export=TRUE,
        kegg_pw_data=list(),
        kegg_mod_data=list(),
        use_kegg_cache=FALSE
    )

    overlap_path <- file.path(out_dir, "enr-overlaps.csv")
    sorted_path <- file.path(out_dir, "enr-overlaps-sorted.csv")
    expect_true(file.exists(overlap_path))
    expect_true(file.exists(sorted_path))

    overlap <- read.csv(overlap_path, check.names=FALSE)
    expect_true("effectsize" %in% names(overlap))
    expect_false(all(is.na(overlap$effectsize)))
    expect_equal(overlap$effectsize, 1.5)
})
