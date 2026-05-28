test_that("make.results.matrix converts result matrices to long rows", {
    results <- list(
        TaxonA=matrix(
            c(1.5, -2.0,
              0.01, 0.20,
              0.10, 0.30,
              8.00, 8.00),
            nrow=4,
            byrow=TRUE,
            dimnames=list(c("Estimate", "p.value", "StdErr", "df"), c("g1", "g2"))
        )
    )

    result_tbl <- make.results.matrix(results)

    expect_equal(result_tbl$taxon, c("TaxonA", "TaxonA"))
    expect_equal(result_tbl$gene, c("g1", "g2"))
    expect_equal(as.numeric(result_tbl$effect.size), c(1.5, -2.0))
    expect_equal(as.numeric(result_tbl$p.value), c(0.01, 0.20))
})

test_that("qvals uses configured p.adjust methods", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    p <- c(0.01, 0.02, 0.50)

    pz.options(fdr_method="BH")
    expect_equal(phylogenize:::qvals(p), p.adjust(p, "BH"))

    pz.options(fdr_method="BY")
    expect_equal(phylogenize:::qvals(p), p.adjust(p, "BY"))
})
