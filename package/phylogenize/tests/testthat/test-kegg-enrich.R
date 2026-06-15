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

test_that("KEGG background is restricted to tested genes", {
    pid_to_ko <- tibble::tibble(
        gene=c("g1", "g2", "g3", "untested"),
        accession=c("K00001", "K00002", "K00003", "K99999")
    )
    signs <- c(g1=1, g2=-1, g3=NA)

    observed <- phylogenize:::tested.kegg.background(signs, pid_to_ko)

    expect_equal(observed, c("K00001", "K00002"))
})

test_that("KEGG enrichment passes taxon-specific tested backgrounds", {
    if (!methods::isClass("MockKeggEnrichResult")) {
        methods::setClass("MockKeggEnrichResult",
                          slots=c(result="data.frame"))
    }
    seen <- new.env(parent=emptyenv())
    seen$backgrounds <- list()

    testthat::local_mocked_bindings(
        enrich_downloaded_KEGG=function(KOs, background, qCut, downloaded, ...) {
            seen$backgrounds[[length(seen$backgrounds) + 1L]] <- background
            methods::new(
                "MockKeggEnrichResult",
                result=data.frame(
                    ID="ko00001",
                    Description="Mock pathway",
                    GeneRatio="1/1",
                    BgRatio=paste0("1/", length(background)),
                    pvalue=0.01,
                    p.adjust=0.01,
                    qvalue=0.01,
                    geneID=paste(KOs, collapse="/"),
                    stringsAsFactors=FALSE
                )
            )
        },
        .package="phylogenize"
    )

    sigs <- list(
        TaxonA=list(strong="g1"),
        TaxonB=list(strong="g3")
    )
    signs <- list(
        TaxonA=c(g1=1, g2=-1),
        TaxonB=c(g3=1, g4=-1)
    )
    pid_to_ko <- tibble::tibble(
        gene=c("g1", "g2", "g3", "g4", "global_only"),
        accession=c("K00001", "K00002", "K00003", "K00004", "K99999")
    )

    out <- multi.kegg.enrich(
        sigs,
        signs,
        pid_to_ko,
        kegg_pw=list(),
        kegg_mod=list()
    )

    expect_s3_class(out, "tbl_df")
    expect_equal(
        seen$backgrounds,
        list(
            c("K00001", "K00002"),
            c("K00001", "K00002"),
            c("K00003", "K00004"),
            c("K00003", "K00004")
        )
    )
})
