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

    core <- phylogenize_core(
        do_enr=TRUE,
        error_to_file=FALSE
    )

    expect_null(core$list_signif)
    expect_null(core$enr_tbls)
    expect_false(seen$enrichment_called)
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
