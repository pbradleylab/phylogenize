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

test_that("16S helpers can use resolved options without global mutation", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempdir()
    local_asv_file <- file.path(tmp, "local-asvs.fna")
    local_hits_file <- file.path(tmp, "local-hits.tsv")
    global_asv_file <- file.path(tmp, "global-asvs.fna")
    global_hits_file <- file.path(tmp, "global-hits.tsv")

    opts <- pz.resolve.options(
        named_asv_file=local_asv_file,
        vsearch_outfile=local_hits_file,
        min_frac_16s=0.75,
        error_to_file=FALSE
    )
    pz.options(
        named_asv_file=global_asv_file,
        vsearch_outfile=global_hits_file,
        min_frac_16s=1,
        error_to_file=FALSE
    )

    mtx <- matrix(
        c(1, 2, 3, 4),
        nrow=2,
        dimnames=list(c("ACGTACGT", "TTTTCCCC"), c("sample1", "sample2"))
    )

    expect_true(prepare.vsearch.input(mtx, .opts=opts))
    expect_true(file.exists(local_asv_file))
    expect_false(file.exists(global_asv_file))

    writeLines(c(
        "Row1\tref;;genus;;species_a\t99",
        "Row2\tref;;genus;;species_b\t99"
    ), local_hits_file)
    writeLines("Row1\tref;;genus;;wrong_species\t99", global_hits_file)

    hits <- get.vsearch.results(.opts=opts)
    expect_equal(hits$targets, c("species_a", "species_b"))

    summed <- sum.nonunique.vsearch(hits, mtx, .opts=opts)
    expect_equal(rownames(summed), c("species_a", "species_b"))
    expect_equal(pz.options("named_asv_file"), global_asv_file)
    expect_equal(pz.options("vsearch_outfile"), global_hits_file)
    expect_equal(pz.options("min_frac_16s"), 1)
})
