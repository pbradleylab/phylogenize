test_that("logging helpers infer local resolved options", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempdir()
    local_log <- file.path(tmp, "local-log-options.txt")
    global_log <- file.path(tmp, "global-log-options.txt")
    unlink(c(local_log, global_log))

    pz.options(
        out_dir=tmp,
        error_file=basename(global_log),
        error_to_file=FALSE,
        verbosity=3
    )

    log_with_local_opts <- function() {
        opts <- pz.resolve.options(
            out_dir=tmp,
            error_file=basename(local_log),
            error_to_file=TRUE,
            verbosity=0
        )
        pz.message("local message", level=1)
        suppressWarnings(pz.warning("local warning"))
    }

    expect_message(log_with_local_opts(), NA)
    expect_true(file.exists(local_log))
    expect_false(file.exists(global_log))
    expect_true(any(grepl("local message", readLines(local_log))))
    expect_true(any(grepl("local warning", readLines(local_log))))
})

test_that("message guard preserves file logging semantics", {
    no_output <- pz.resolve.options(
        error_to_file=FALSE,
        verbosity=0
    )
    file_output <- pz.resolve.options(
        error_to_file=TRUE,
        verbosity=0
    )

    expect_false(pz.should.message(level=1, .opts=no_output))
    expect_true(pz.should.message(level=1, .opts=file_output))
    expect_true(pz.should.message(level=0, .opts=no_output))
})
