test_that("install_data validates source files and extensions", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(error_to_file=FALSE)

    extdata <- file.path(tempdir(), "install-data-validation")
    dir.create(extdata, showWarnings=FALSE)

    expect_error(
        install_data(file.path(tempdir(), "missing.csv"), .extd_path=extdata),
        "Data archive not found"
    )

    unsupported <- tempfile(fileext=".txt")
    writeLines("not supported", unsupported)
    expect_error(
        install_data(unsupported, .extd_path=extdata),
        "Unsupported data file extension"
    )

    source_csv <- tempfile(fileext=".csv")
    writeLines("database,genes", source_csv)
    expect_error(
        install_data(source_csv, .extd_path=file.path(tempdir(), "missing-dir")),
        "extdata directory not found"
    )
})

test_that("install_data copies csv files and protects existing destinations", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(error_to_file=FALSE)

    extdata <- file.path(tempdir(), "install-data-copy")
    dir.create(extdata, showWarnings=FALSE)
    source_csv <- file.path(tempdir(), "test-database.csv")
    writeLines(c("database,genes", "db-one,genes.rds"), source_csv)

    installed <- install_data(source_csv, .extd_path=extdata)
    expect_equal(installed, file.path(extdata, basename(source_csv)))
    expect_equal(readLines(installed), c("database,genes", "db-one,genes.rds"))

    writeLines(c("database,genes", "db-two,genes.rds"), source_csv)
    expect_error(
        install_data(source_csv, .extd_path=extdata),
        "Destination file already exists"
    )
    expect_equal(readLines(installed), c("database,genes", "db-one,genes.rds"))

    overwritten <- install_data(source_csv, force=TRUE, .extd_path=extdata)
    expect_equal(overwritten, installed)
    expect_equal(readLines(installed), c("database,genes", "db-two,genes.rds"))
})
