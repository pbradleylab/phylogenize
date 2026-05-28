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
