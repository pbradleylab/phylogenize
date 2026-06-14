test_that("resolved option objects remain isolated from later global changes", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    opts_a <- pz.resolve.options(
        out_dir="run-a",
        which_phenotype="prevalence",
        fdr_method="BH",
        error_to_file=FALSE
    )
    opts_b <- pz.resolve.options(
        out_dir="run-b",
        which_phenotype="specificity",
        fdr_method="BY",
        error_to_file=FALSE
    )

    pz.options(
        out_dir="global-run",
        which_phenotype="abundance",
        fdr_method="BH",
        error_to_file=FALSE
    )

    expect_equal(opts_a("out_dir"), "run-a")
    expect_equal(opts_b("out_dir"), "run-b")
    expect_equal(opts_a("which_phenotype"), "prevalence")
    expect_equal(opts_b("which_phenotype"), "specificity")
    expect_equal(phylogenize:::qvals(c(0.01, 0.02, 0.50), .opts=opts_b),
                 p.adjust(c(0.01, 0.02, 0.50), "BY"))
    expect_equal(pz.options("out_dir"), "global-run")
})

test_that("packaged database index loads without readr warnings", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(data_dir="", error_to_file=FALSE)

    expect_warning(
        expect_message(
            check_data_found(startup=FALSE),
            "Databases listed:"
        ),
        NA
    )
    expect_true(dir.exists(pz.options("data_dir")))
})

test_that("database index validation catches malformed files", {
    bad_index <- tempfile(fileext=".csv")
    writeLines(c(
        "genes,trees",
        "genes.rds,tree.rds"
    ), bad_index)

    expect_error(
        phylogenize:::read.internal.database.index(bad_index),
        "missing required column"
    )

    blank_index <- tempfile(fileext=".csv")
    writeLines(c(
        "database,genes,trees",
        ",genes.rds,tree.rds"
    ), blank_index)

    expect_error(
        phylogenize:::read.internal.database.index(blank_index),
        "blank database names"
    )
})
