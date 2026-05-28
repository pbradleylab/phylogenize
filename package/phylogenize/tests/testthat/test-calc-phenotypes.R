test_that("calculate_phenotypes uses direct treemin override", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(which_phenotype="provided",
               treemin=10,
               error_to_file=FALSE)

    phenotype_file <- tempfile(fileext=".tsv")
    writeLines(c(
        "species\tphenotype",
        "s1\t1",
        "s2\t2",
        "s3\t3"
    ), phenotype_file)

    abd.meta <- list(
        mtx=Matrix::Matrix(
            matrix(
                c(1, 0, 0,
                  0, 1, 0,
                  0, 0, 1),
                nrow=3,
                byrow=TRUE,
                dimnames=list(paste0("s", 1:3), paste0("sample", 1:3))
            ),
            sparse=TRUE
        )
    )
    pz.db <- list(
        trees=list(TaxonA=ape::read.tree(text="(s1:1,s2:1,s3:1);")),
        gene.presence=list(
            TaxonA=Matrix::Matrix(
                matrix(
                    c(1, 1, 0),
                    nrow=1,
                    dimnames=list("gene1", paste0("s", 1:3))
                ),
                sparse=TRUE
            )
        )
    )

    phenotype_results <- calculate_phenotypes(
        abd.meta,
        pz.db,
        which_phenotype="provided",
        phenotype_file=phenotype_file,
        treemin=2,
        error_to_file=FALSE
    )

    expect_equal(phenotype_results$phenotype, c(s1=1, s2=2, s3=3))
})
