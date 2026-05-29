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

    pz.options(fdr_method="BH")
    expect_equal(phylogenize:::qvals(p, fdr_method="BY"), p.adjust(p, "BY"))
})

test_that("calc.alpha.power compares named rejected tests to named truth sets", {
    pvs <- c(g_null=0.01, g_alt=0.02, g_miss=0.50, g_alt2=0.20)

    out <- calc.alpha.power(
        pvs,
        null=c("g_null", "g_miss"),
        alt=c("g_alt", "g_alt2"),
        alpha=0.05
    )

    expect_equal(unname(out["r"]), 0.5)
    expect_equal(unname(out["p"]), 0.5)
    expect_equal(unname(out["a"]), 0.5)

    filtered <- calc.alpha.power(
        pvs,
        null=c("g_null", "g_miss"),
        alt=c("g_alt", "g_alt2"),
        alpha=0.05,
        filter=c("g_alt", "g_miss")
    )

    expect_equal(unname(filtered["r"]), 0.25)
    expect_equal(unname(filtered["p"]), 1)
    expect_equal(unname(filtered["a"]), 0)
})

test_that("negative-only significant taxa are retained after thresholding", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(fdr_method="BH",
               min_fx=0,
               minimum=3,
               error_to_file=FALSE)

    results <- list(
        TaxonNeg=matrix(
            c(-1.5,
              0.001,
              0.10,
              8.00),
            nrow=4,
            dimnames=list(c("Estimate", "p.value", "StdErr", "df"), "g_neg")
        )
    )
    pz.db <- list(
        trees=list(TaxonNeg=ape::read.tree(text="(s1:1,s2:1,s3:1,s4:1,s5:1,s6:1);")),
        gene.presence=list(
            TaxonNeg=Matrix::Matrix(
                matrix(
                    c(1, 1, 1, 0, 0, 0),
                    nrow=1,
                    dimnames=list("g_neg", paste0("s", 1:6))
                ),
                sparse=TRUE
            )
        ),
        gene.to.fxn=tibble::tibble(
            gene="g_neg",
            accession="K00001",
            `function`="Mock negative gene"
        )
    )

    out <- get_signif_associated_genes(pz.db, results)

    expect_equal(out$phy.with.pos.sigs, character(0))
    expect_equal(out$phy.with.neg.sigs, "TaxonNeg")
    expect_equal(out$neg.sig.thresh$TaxonNeg, "g_neg")
    expect_equal(out$phy.with.sigs, "TaxonNeg")
})
