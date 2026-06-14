test_that("run.vsearch passes strand option as separate arguments", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempdir()
    capture_file <- file.path(tmp, "vsearch-args.txt")
    fake_vsearch <- file.path(tmp, "fake-vsearch")
    writeLines(c(
        "#!/bin/sh",
        sprintf("printf '%%s\\n' \"$@\" > %s", shQuote(capture_file)),
        "exit 0"
    ), fake_vsearch)
    Sys.chmod(fake_vsearch, mode="0755")

    expect_true(run.vsearch(
        vsearch_dir=fake_vsearch,
        data_dir=tmp,
        vsearch_16sfile="db.fna",
        named_asv_file="input.fna",
        vsearch_cutoff=0.985,
        vsearch_outfile="hits.tsv",
        error_to_file=FALSE
    ))

    args <- readLines(capture_file)
    strand_idx <- match("--strand", args)

    expect_false("--strand both" %in% args)
    expect_false(is.na(strand_idx))
    expect_equal(args[strand_idx + 1], "both")
})
