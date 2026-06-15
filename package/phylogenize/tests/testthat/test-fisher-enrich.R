test_that("legacy enrichment helpers unnest without tidyr warnings", {
    sigs <- list(
        TaxonA=list(strong=c("g1", "g2"))
    )
    signs <- list(
        TaxonA=c(g1=1, g2=1, g3=-1, g4=1)
    )
    mappings <- list(
        pathway=tibble::tibble(
            gene=c("g1", "g3", "g4"),
            term=c("term-a", "term-b", "term-a")
        )
    )

    expect_warning(
        enriched <- multi.enrich(sigs, signs, mappings),
        NA
    )

    expect_s3_class(enriched, "data.frame")
    expect_equal(enriched$taxon, c("TaxonA", "TaxonA"))
    expect_equal(enriched$cutoff, c("strong", "strong"))
    expect_named(enriched, c(
        "taxon", "cutoff", "termset", "term", "data", "enr",
        "enr.pval", "enr.estimate", "enr.overlap", "enr.qval"
    ))
})

test_that("enrichment mapping preparation is shared and stable", {
    mappings <- list(
        pathway=tibble::tibble(
            gene=c("g1", "g3", "g4"),
            term=c("term-a", "term-b", "term-a")
        )
    )

    prepared <- phylogenize:::prepare.enrichment.mappings(mappings)

    expect_named(prepared, c("tbl.mappings", "map.bg"))
    expect_named(prepared$tbl.mappings, c("termset", "terms"))
    expect_equal(prepared$tbl.mappings$termset, "pathway")
    expect_equal(prepared$tbl.mappings$terms[[1]]$term,
                 c("term-a", "term-b"))
    expect_equal(prepared$map.bg, c("g1", "g4", "g3"))
})

test_that("tbl.result.qvs unnests result rows without tidyr warnings", {
    results <- tibble::tibble(
        taxon=c("TaxonA", "TaxonA", "TaxonB"),
        gene=c("g1", "g2", "g3"),
        p.value=c(0.01, 0.20, 0.03)
    )

    expect_warning(
        out <- tbl.result.qvs(results, method=function(x) p.adjust(x, "BH")),
        NA
    )

    expect_equal(nrow(out), 3)
    expect_equal(out$q.value[1:2], p.adjust(results$p.value[1:2], "BH"))
    expect_equal(out$q.value[3], results$p.value[3])
})
