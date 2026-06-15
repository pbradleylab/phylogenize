test_that("KEGG enrichment formatting handles null and non-null results", {
    if (!methods::isClass("MockKeggEnrichResult")) {
        methods::setClass("MockKeggEnrichResult",
                          slots=c(result="data.frame"))
    }
    pwy <- methods::new(
        "MockKeggEnrichResult",
        result=data.frame(
            ID="ko00001",
            Description="Mock pathway",
            GeneRatio="1/2",
            BgRatio="3/4",
            stringsAsFactors=FALSE
        )
    )
    mod <- methods::new(
        "MockKeggEnrichResult",
        result=data.frame(
            ID="M00001",
            Description="Mock module",
            GeneRatio="2/3",
            BgRatio="4/5",
            stringsAsFactors=FALSE
        )
    )

    both <- phylogenize:::format.kegg.enrichment.result(
        list(pwy, mod),
        cutoff="strong",
        taxon="TaxonA"
    )
    one <- phylogenize:::format.kegg.enrichment.result(
        list(pwy, NULL),
        cutoff="weak",
        taxon="TaxonB"
    )
    none <- phylogenize:::format.kegg.enrichment.result(
        list(NULL, NULL),
        cutoff="strong",
        taxon="TaxonA"
    )

    expect_s3_class(both, "tbl_df")
    expect_equal(both$ID, c("ko00001", "M00001"))
    expect_equal(both$cutoff, c("strong", "strong"))
    expect_equal(both$taxon, c("TaxonA", "TaxonA"))
    expect_equal(one$ID, "ko00001")
    expect_equal(one$cutoff, "weak")
    expect_equal(one$taxon, "TaxonB")
    expect_null(none)
})
