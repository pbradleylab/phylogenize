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
