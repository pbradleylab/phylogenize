test_that("phylogenize_core respects explicit do_POMS", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    seen <- new.env(parent=emptyenv())

    testthat::local_mocked_bindings(
        data_to_phenotypes=function(save_data=FALSE, ...) {
            dots <- list(...)
            seen$phenotype_core_method <- dots[["core_method"]]
            list(
                phenotype_results=list(phenotype=c(sp1=1)),
                pz.db=list(ntaxa=1)
            )
        },
        get_all_associated_genes=function(list_pheno, p.method, ...) {
            dots <- list(...)
            seen$association_core_method <- dots[["core_method"]]
            NULL
        },
        .package="phylogenize"
    )

    phylogenize_core(
        do_enr=FALSE,
        do_POMS=TRUE,
        core_method="lm",
        error_to_file=FALSE
    )

    expect_equal(seen$phenotype_core_method, "poms")
    expect_equal(seen$association_core_method, "poms")
})

test_that("phylogenize_core skips enrichment when associations return NULL", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    seen <- new.env(parent=emptyenv())
    seen$enrichment_called <- FALSE

    testthat::local_mocked_bindings(
        data_to_phenotypes=function(save_data=FALSE, ...) {
            list(
                phenotype_results=list(phenotype=c(sp1=1)),
                pz.db=list(ntaxa=1)
            )
        },
        get_all_associated_genes=function(list_pheno, p.method, ...) {
            NULL
        },
        get_enrichment_tbls=function(...) {
            seen$enrichment_called <- TRUE
            tibble::tibble()
        },
        .package="phylogenize"
    )

    expect_warning(
        core <- phylogenize_core(
            do_enr=TRUE,
            error_to_file=FALSE
        ),
        "Association testing returned no results"
    )

    expect_null(core$list_signif)
    expect_null(core$enr_tbls)
    expect_false(seen$enrichment_called)
})

test_that("phylogenize_core returns a stable result object shape", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    set.seed(20240618)
    tips <- paste0("sp", seq_len(30))
    phenotype <- seq(-1.5, 1.5, length.out=length(tips))
    names(phenotype) <- tips
    high_phenotype <- phenotype >= stats::median(phenotype)
    gene_presence <- matrix(
        c(as.integer(high_phenotype), as.integer(!high_phenotype)),
        nrow=2,
        byrow=TRUE,
        dimnames=list(c("g_pos", "g_neg"), tips)
    )
    pz_db <- list(
        ntaxa=1,
        species=list(TestTaxon=tips),
        trees=list(TestTaxon=ape::rtree(length(tips), tip.label=tips)),
        gene.presence=list(TestTaxon=gene_presence),
        gene.to.fxn=tibble::tibble(
            gene=rownames(gene_presence),
            accession=c("K00001", "K00002"),
            `function`=c("positive fixture gene", "negative fixture gene")
        )
    )
    out_dir <- file.path(tempdir(), "phylogenize-core-shape-test")

    testthat::local_mocked_bindings(
        data_to_phenotypes=function(save_data=FALSE, ...) {
            list(
                phenotype_results=list(phenotype=phenotype),
                pz.db=pz_db
            )
        },
        get_enrichment_tbls=function(...) {
            tibble::tibble(
                taxon="TestTaxon",
                term="mock_term",
                enr.qval=0.01,
                enr.pval=0.001
            )
        },
        .package="phylogenize"
    )

    core <- phylogenize_core(
        do_enr=TRUE,
        p.method=lm.fx.pv,
        core_method="lm",
        min_fx=0,
        minimum=3,
        out_dir=out_dir,
        error_to_file=FALSE,
        verbosity=0
    )

    expect_named(core, c("list_pheno", "list_signif", "enr_tbls", "options"))
    expect_named(core$list_pheno, c("phenotype_results", "pz.db"))
    expect_equal(names(core$list_pheno$phenotype_results$phenotype), tips)
    expect_equal(core$list_pheno$pz.db$ntaxa, 1)
    expect_named(core$list_signif, c(
        "results", "signif", "signs", "pos.sig", "results.matrix",
        "phy.with.sigs", "pos.sig.descs", "pos.sig.thresh",
        "pos.sig.thresh.descs", "neg.sig", "neg.sig.descs",
        "neg.sig.thresh", "neg.sig.thresh.descs", "phy.with.pos.sigs",
        "phy.with.neg.sigs"
    ))
    expect_named(core$list_signif$results, "TestTaxon")
    expect_equal(rownames(core$list_signif$results$TestTaxon),
                 c("Estimate", "p.value", "StdErr", "df"))
    expect_equal(colnames(core$list_signif$results$TestTaxon),
                 c("g_pos", "g_neg"))
    expect_equal(nrow(core$list_signif$results.matrix), 2)
    expect_equal(
        names(core$list_signif$results.matrix),
        c("taxon", "gene", "effect.size", "p.value", "std.err", "df",
          "q.value")
    )
    expect_setequal(core$list_signif$results.matrix$gene, c("g_pos", "g_neg"))
    expect_equal(core$list_signif$phy.with.pos.sigs, "TestTaxon")
    expect_equal(core$list_signif$phy.with.neg.sigs, "TestTaxon")
    expect_s3_class(core$enr_tbls, "tbl_df")
    expect_equal(nrow(core$enr_tbls), 1)
    expect_equal(core$options("out_dir"),
                 normalizePath(out_dir, mustWork=FALSE))
    expect_equal(core$options("core_method"), "lm")
})

test_that("phylogenize_core preserves deterministic golden association results", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    set.seed(20240619)
    tips <- paste0("sp", seq_len(30))
    phenotype <- seq(-1.5, 1.5, length.out=length(tips))
    names(phenotype) <- tips
    high_phenotype <- phenotype >= stats::median(phenotype)
    gene_presence <- matrix(
        c(as.integer(high_phenotype), as.integer(!high_phenotype)),
        nrow=2,
        byrow=TRUE,
        dimnames=list(c("g_pos", "g_neg"), tips)
    )
    pz_db <- list(
        ntaxa=1,
        species=list(TestTaxon=tips),
        trees=list(TestTaxon=ape::rtree(length(tips), tip.label=tips)),
        gene.presence=list(TestTaxon=gene_presence),
        gene.to.fxn=tibble::tibble(
            gene=rownames(gene_presence),
            accession=c("K00001", "K00002"),
            `function`=c("positive fixture gene", "negative fixture gene")
        )
    )

    testthat::local_mocked_bindings(
        data_to_phenotypes=function(save_data=FALSE, ...) {
            list(
                phenotype_results=list(phenotype=phenotype),
                pz.db=pz_db
            )
        },
        .package="phylogenize"
    )

    core <- phylogenize_core(
        do_enr=FALSE,
        p.method=lm.fx.pv,
        core_method="lm",
        min_fx=0,
        minimum=3,
        out_dir=tempdir(),
        error_to_file=FALSE,
        verbosity=0
    )
    results_matrix <- core$list_signif$results.matrix

    expect_equal(results_matrix$taxon, c("TestTaxon", "TestTaxon"))
    expect_equal(results_matrix$gene, c("g_pos", "g_neg"))
    expect_equal(
        unname(results_matrix$effect.size),
        c(1.55172413793103, -1.55172413793103),
        tolerance=1e-12
    )
    expect_equal(
        unname(results_matrix$p.value),
        rep(6.065437457762578e-10, 2),
        tolerance=1e-15
    )
    expect_equal(
        unname(results_matrix$std.err),
        rep(0.168930327088495, 2),
        tolerance=1e-12
    )
    expect_equal(unname(results_matrix$df), c(28, 28))
    expect_equal(
        unname(results_matrix$q.value),
        rep(6.065437457762578e-10, 2),
        tolerance=1e-15
    )
    expect_equal(core$list_signif$signif$TestTaxon$strong,
                 c("g_pos", "g_neg"))
    expect_equal(core$list_signif$pos.sig$TestTaxon, "g_pos")
    expect_equal(core$list_signif$neg.sig$TestTaxon, "g_neg")
    expect_equal(core$list_signif$pos.sig.thresh$TestTaxon, "g_pos")
    expect_equal(core$list_signif$neg.sig.thresh$TestTaxon, "g_neg")
    expect_equal(core$list_signif$phy.with.pos.sigs, "TestTaxon")
    expect_equal(core$list_signif$phy.with.neg.sigs, "TestTaxon")
})

test_that("serial association honors core_method instead of p.method", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    set.seed(20240620)
    tips <- paste0("sp", seq_len(30))
    phenotype <- seq(-1.5, 1.5, length.out=length(tips))
    names(phenotype) <- tips
    high_phenotype <- phenotype >= stats::median(phenotype)
    gene_presence <- matrix(
        c(as.integer(high_phenotype), as.integer(!high_phenotype)),
        nrow=2,
        byrow=TRUE,
        dimnames=list(c("g_pos", "g_neg"), tips)
    )
    list_pheno <- list(
        phenotype_results=list(phenotype=phenotype),
        pz.db=list(
            ntaxa=1,
            species=list(TestTaxon=tips),
            trees=list(TestTaxon=ape::rtree(length(tips), tip.label=tips)),
            gene.presence=list(TestTaxon=gene_presence),
            gene.to.fxn=tibble::tibble(
                gene=rownames(gene_presence),
                accession=c("K00001", "K00002"),
                `function`=c("positive fixture gene", "negative fixture gene")
            )
        )
    )
    failing_method <- function(...) {
        stop("p.method should not be used when core_method selects lm")
    }

    observed <- get_all_associated_genes(
        list_pheno,
        p.method=failing_method,
        core_method="lm",
        ncl=1,
        min_fx=0,
        minimum=3,
        error_to_file=FALSE,
        verbosity=0
    )

    expect_named(observed$results, "TestTaxon")
    expect_equal(colnames(observed$results$TestTaxon), c("g_pos", "g_neg"))
    expect_equal(observed$pos.sig$TestTaxon, "g_pos")
    expect_equal(observed$neg.sig$TestTaxon, "g_neg")
})

test_that("phylogenize_core creates configured output directory", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    out_dir <- file.path(tempdir(), "phylogenize-core-output-dir-test")
    if (dir.exists(out_dir)) unlink(out_dir, recursive=TRUE)

    testthat::local_mocked_bindings(
        data_to_phenotypes=function(save_data=FALSE, ...) {
            list(
                phenotype_results=list(phenotype=c(sp1=1)),
                pz.db=list(ntaxa=1)
            )
        },
        get_all_associated_genes=function(list_pheno, p.method, ...) {
            NULL
        },
        .package="phylogenize"
    )

    phylogenize_core(
        do_enr=FALSE,
        out_dir=out_dir,
        error_to_file=FALSE
    )

    expect_true(dir.exists(out_dir))
})

test_that("phylogenize wrapper passes resolved options without mutating globals", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    out_dir <- file.path(tempdir(), "phylogenize-wrapper-options-test")
    seen <- new.env(parent=emptyenv())

    testthat::local_mocked_bindings(
        phylogenize_core=function(..., .opts=NULL) {
            seen$core_opts <- .opts
            list(
                list_pheno=list(),
                list_signif=NULL,
                enr_tbls=NULL,
                options=.opts
            )
        },
        render_core_report=function(core, ..., .opts=NULL) {
            seen$render_opts <- .opts
            file.path(.opts("out_dir"), basename(.opts("output_file")))
        },
        .package="phylogenize"
    )

    pz.options(out_dir="global-out", output_file="global.html")
    phylogenize(
        out_dir=out_dir,
        output_file="local.html",
        rds_output_file="",
        error_to_file=FALSE,
        reset_after=TRUE
    )

    expect_equal(seen$core_opts("out_dir"), normalizePath(out_dir, mustWork=FALSE))
    expect_equal(seen$core_opts("output_file"), "local.html")
    expect_equal(seen$render_opts("out_dir"), normalizePath(out_dir, mustWork=FALSE))
    expect_equal(pz.options("out_dir"), "global-out")
    expect_equal(pz.options("output_file"), "global.html")
})

test_that("augment_with_enrichments uses options stored in core", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    seen <- new.env(parent=emptyenv())
    core_opts <- pz.resolve.options(
        out_dir="core-out",
        fdr_method="BY",
        error_to_file=FALSE
    )
    core <- list(
        list_signif=list(
            signif=list(TaxonA=character()),
            signs=list(TaxonA=character()),
            results.matrix=tibble::tibble()
        ),
        list_pheno=list(pz.db=list()),
        options=core_opts
    )

    testthat::local_mocked_bindings(
        get_enrichment_tbls=function(..., .opts=NULL) {
            seen$opts <- .opts
            tibble::tibble()
        },
        .package="phylogenize"
    )

    pz.options(out_dir="global-out", fdr_method="BH", error_to_file=FALSE)
    augmented <- augment_with_enrichments(core)

    expect_s3_class(augmented$enr_tbls, "tbl_df")
    expect_equal(seen$opts("out_dir"), "core-out")
    expect_equal(seen$opts("fdr_method"), "BY")
    expect_equal(pz.options("out_dir"), "global-out")
})
