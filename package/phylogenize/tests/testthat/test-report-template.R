test_that("report template creates output directories recursively", {
    report_path <- system.file("rmd", "phylogenize-report.Rmd", package="phylogenize")
    if (report_path == "") {
        report_path <- testthat::test_path("../../inst/rmd/phylogenize-report.Rmd")
    }

    report_text <- readLines(report_path)
    dir_create_lines <- grep("dir\\.create", report_text, value=TRUE)

    expect_true(all(grepl("recursive=TRUE", dir_create_lines, fixed=TRUE)))
    expect_true(all(grepl("showWarnings=FALSE", dir_create_lines, fixed=TRUE)))
})

test_that("report template uses valid core paths and guards missing associations", {
    report_path <- system.file("rmd", "phylogenize-report.Rmd", package="phylogenize")
    if (report_path == "") {
        report_path <- testthat::test_path("../../inst/rmd/phylogenize-report.Rmd")
    }

    report_text <- readLines(report_path)

    expect_false(any(grepl('core\\[\\["pz\\.db"\\]\\]\\[\\["list_pheno"\\]\\]',
                           report_text)))
    expect_true(any(grepl("is.null\\(signif_genes\\)", report_text)))
})

test_that("report template namespaces tidyverse helper calls", {
    report_path <- system.file("rmd", "phylogenize-report.Rmd", package="phylogenize")
    if (report_path == "") {
        report_path <- testthat::test_path("../../inst/rmd/phylogenize-report.Rmd")
    }

    report_text <- readLines(report_path)
    report_body <- paste(report_text, collapse="\n")

    expect_false(grepl("(?<!::)\\bfull_join\\(", report_body, perl=TRUE))
    expect_false(grepl("(?<!::)\\benframe\\(", report_body, perl=TRUE))
    expect_false(grepl("(?<!::)\\bmap\\(", report_body, perl=TRUE))
    expect_true(grepl("dplyr::full_join\\(", report_body))
    expect_true(grepl("tibble::enframe\\(", report_body))
    expect_true(grepl("purrr::map\\(", report_body))
})
