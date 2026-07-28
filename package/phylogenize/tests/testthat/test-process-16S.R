test_that("run.vsearch passes strand option as separate arguments", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempfile("vsearch-")
    dir.create(tmp)
    capture_file <- file.path(tmp, "vsearch-args.txt")
    fake_bindir <- file.path(tmp, "bin")
    dir.create(fake_bindir)
    fake_vsearch <- file.path(fake_bindir, "vsearch")
    query_file <- file.path(tmp, "input.fna")
    db_file <- file.path(tmp, "db.fna")
    hits_file <- file.path(tmp, "hits.tsv")
    writeLines(c(">Row1", "ACGT"), query_file)
    writeLines(c(">ref;;genus;;species_a", "ACGT"), db_file)
    writeLines(c(
        "#!/bin/sh",
        sprintf("printf '%%s\\n' \"$@\" > %s", shQuote(capture_file)),
        "last_arg=''",
        "for arg in \"$@\"; do last_arg=\"$arg\"; done",
        "touch \"$last_arg\"",
        "exit 0"
    ), fake_vsearch)
    Sys.chmod(fake_vsearch, mode="0755")

    expect_true(run.vsearch(
        vsearch_dir=fake_bindir,
        data_dir=tmp,
        vsearch_16sfile="db.fna",
        named_asv_file=query_file,
        vsearch_cutoff=0.985,
        vsearch_outfile=hits_file,
        error_to_file=FALSE
    ))

    args <- readLines(capture_file)
    strand_idx <- match("--strand", args)

    expect_false("--strand both" %in% args)
    expect_false(is.na(strand_idx))
    expect_equal(args[strand_idx + 1], "both")
    expect_true(file.exists(hits_file))
})

test_that("run.appspam uses shared external command runner", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempfile("appspam-")
    dir.create(tmp)
    capture_file <- file.path(tmp, "appspam-args.txt")
    fake_appspam <- file.path(tmp, "appspam")
    aln_file <- file.path(tmp, "aln.fna")
    tree_file <- file.path(tmp, "tree.nwk")
    query_file <- file.path(tmp, "input.fna")
    jplace_file <- file.path(tmp, "out.jplace")
    writeLines(c(">ref", "ACGT"), aln_file)
    writeLines("(ref:1);", tree_file)
    writeLines(c(">Row1", "ACGT"), query_file)
    writeLines(c(
        "#!/bin/sh",
        sprintf("printf '%%s\\n' \"$@\" > %s", shQuote(capture_file)),
        "exit 0"
    ), fake_appspam)
    Sys.chmod(fake_appspam, mode="0755")

    expect_true(run.appspam(
        appspam_path=fake_appspam,
        aln_path_16s=aln_file,
        tree_path_16s=tree_file,
        named_asv_file=query_file,
        jplace_file=jplace_file,
        ncl=3,
        error_to_file=FALSE
    ))

    expect_equal(
        readLines(capture_file),
        c("-s", aln_file,
          "-t", tree_file,
          "-q", query_file,
          "--threads", "3",
          "-o", jplace_file)
    )

    writeLines(c("#!/bin/sh", "exit 7"), fake_appspam)
    Sys.chmod(fake_appspam, mode="0755")
    expect_error(
        run.appspam(
            appspam_path=fake_appspam,
            aln_path_16s=aln_file,
            tree_path_16s=tree_file,
            named_asv_file=query_file,
            jplace_file=jplace_file,
            ncl=3,
            error_to_file=FALSE
        ),
        "AppSpam failed with error code 7"
    )
})

test_that("run.appspam validates executable, inputs, and output options", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempfile("appspam-validation-")
    dir.create(tmp)
    fake_appspam <- file.path(tmp, "appspam")
    aln_file <- file.path(tmp, "aln.fna")
    tree_file <- file.path(tmp, "tree.nwk")
    query_file <- file.path(tmp, "input.fna")
    jplace_file <- file.path(tmp, "out.jplace")
    writeLines(c("#!/bin/sh", "exit 0"), fake_appspam)
    writeLines(c(">ref", "ACGT"), aln_file)
    writeLines("(ref:1);", tree_file)
    writeLines(c(">Row1", "ACGT"), query_file)

    expect_error(
        run.appspam(
            appspam_path=file.path(tmp, "missing-appspam"),
            aln_path_16s=aln_file,
            tree_path_16s=tree_file,
            named_asv_file=query_file,
            jplace_file=jplace_file,
            ncl=1,
            error_to_file=FALSE
        ),
        "AppSpam executable not found"
    )

    expect_error(
        run.appspam(
            appspam_path=fake_appspam,
            aln_path_16s=aln_file,
            tree_path_16s=tree_file,
            named_asv_file=query_file,
            jplace_file=jplace_file,
            ncl=1,
            error_to_file=FALSE
        ),
        "AppSpam executable is not runnable"
    )

    Sys.chmod(fake_appspam, mode="0755")

    expect_error(
        run.appspam(
            appspam_path=fake_appspam,
            aln_path_16s=file.path(tmp, "missing-aln.fna"),
            tree_path_16s=tree_file,
            named_asv_file=query_file,
            jplace_file=jplace_file,
            ncl=1,
            error_to_file=FALSE
        ),
        "16S alignment file not found"
    )

    expect_error(
        run.appspam(
            appspam_path=fake_appspam,
            aln_path_16s=aln_file,
            tree_path_16s=file.path(tmp, "missing-tree.nwk"),
            named_asv_file=query_file,
            jplace_file=jplace_file,
            ncl=1,
            error_to_file=FALSE
        ),
        "16S tree file not found"
    )

    expect_error(
        run.appspam(
            appspam_path=fake_appspam,
            aln_path_16s=aln_file,
            tree_path_16s=tree_file,
            named_asv_file=file.path(tmp, "missing-input.fna"),
            jplace_file=jplace_file,
            ncl=1,
            error_to_file=FALSE
        ),
        "16S query FASTA not found"
    )

    expect_error(
        run.appspam(
            appspam_path=fake_appspam,
            aln_path_16s=aln_file,
            tree_path_16s=tree_file,
            named_asv_file=query_file,
            jplace_file="",
            ncl=1,
            error_to_file=FALSE
        ),
        "jplace_file must be provided"
    )

    expect_error(
        run.appspam(
            appspam_path=fake_appspam,
            aln_path_16s=aln_file,
            tree_path_16s=tree_file,
            named_asv_file=query_file,
            jplace_file=file.path(tmp, "missing-dir", "out.jplace"),
            ncl=1,
            error_to_file=FALSE
        ),
        "AppSpam output directory not found"
    )

    expect_error(
        run.appspam(
            appspam_path=fake_appspam,
            aln_path_16s=aln_file,
            tree_path_16s=tree_file,
            named_asv_file=query_file,
            jplace_file=jplace_file,
            ncl=1.5,
            error_to_file=FALSE
        ),
        "ncl must be a positive integer"
    )
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

test_that("process.16s requires an explicit 16S mapping method", {
    abd.meta <- list(
        mtx=matrix(
            c(1, 0, 0, 1),
            nrow=2,
            dimnames=list(c("ACGTACGT", "TTTTCCCC"), c("sample1", "sample2"))
        ),
        metadata=data.frame(sample=c("sample1", "sample2"))
    )

    expect_error(
        process.16s(
            abd.meta,
            which_16s_method="",
            error_to_file=FALSE
        ),
        "type_16S=TRUE requires which_16s_method to be explicitly set"
    )
})

test_that("sum.nonunique.vsearch drops ambiguous ties and aggregates targets", {
    opts <- pz.resolve.options(
        min_frac_16s=0.5,
        error_to_file=FALSE
    )
    mtx <- Matrix::Matrix(
        matrix(
            c(1, 2, 3,
              4, 5, 6,
              7, 8, 9,
              10, 11, 12),
            nrow=4,
            byrow=TRUE,
            dimnames=list(paste0("asv", 1:4), paste0("sample", 1:3))
        ),
        sparse=TRUE
    )
    vsearch <- list(
        hits=c(1, 1, 2, 3, 3, 4),
        targets=c("target_b", "target_b", "target_a",
                  "target_c", "target_d", "target_a")
    )

    summed <- sum.nonunique.vsearch(vsearch, mtx, .opts=opts)

    expect_equal(rownames(summed), c("target_b", "target_a"))
    expect_equal(
        as.matrix(summed),
        rbind(
            target_b=c(1, 2, 3),
            target_a=c(14, 16, 18)
        ),
        ignore_attr="dimnames"
    )
})

test_that("run.vsearch validates executable and input files", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempfile("vsearch-validation-")
    dir.create(tmp)
    query_file <- file.path(tmp, "input.fna")
    db_file <- file.path(tmp, "db.fna")
    fake_vsearch <- file.path(tmp, "vsearch")
    writeLines(c(">Row1", "ACGT"), query_file)
    writeLines(c(">ref;;genus;;species_a", "ACGT"), db_file)
    writeLines(c("#!/bin/sh", "exit 0"), fake_vsearch)

    expect_error(
        run.vsearch(
            vsearch_dir=fake_vsearch,
            data_dir=tmp,
            vsearch_16sfile="db.fna",
            named_asv_file=query_file,
            vsearch_outfile=file.path(tmp, "hits.tsv"),
            error_to_file=FALSE
        ),
        "not runnable"
    )

    Sys.chmod(fake_vsearch, mode="0755")
    expect_error(
        run.vsearch(
            vsearch_dir=fake_vsearch,
            data_dir=tmp,
            vsearch_16sfile="db.fna",
            named_asv_file=file.path(tmp, "missing.fna"),
            vsearch_outfile=file.path(tmp, "hits.tsv"),
            error_to_file=FALSE
        ),
        "query FASTA not found"
    )

    expect_error(
        run.vsearch(
            vsearch_dir=fake_vsearch,
            data_dir=tmp,
            vsearch_16sfile="missing-db.fna",
            named_asv_file=query_file,
            vsearch_outfile=file.path(tmp, "hits.tsv"),
            error_to_file=FALSE
        ),
        "reference FASTA not found"
    )

    expect_error(
        run.vsearch(
            vsearch_dir=fake_vsearch,
            data_dir=tmp,
            vsearch_16sfile="db.fna",
            named_asv_file=query_file,
            vsearch_cutoff=1.5,
            vsearch_outfile=file.path(tmp, "hits.tsv"),
            error_to_file=FALSE
        ),
        "vsearch_cutoff must be"
    )
})

test_that("get.vsearch.results rejects malformed vsearch output", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempfile("vsearch-output-")
    dir.create(tmp)
    hits_file <- file.path(tmp, "hits.tsv")

    expect_error(
        get.vsearch.results(
            vsearch_outfile=hits_file,
            error_to_file=FALSE
        ),
        "output file not found"
    )

    writeLines("bad_query\tref;;genus;;species_a\t99", hits_file)
    expect_error(
        get.vsearch.results(
            vsearch_outfile=hits_file,
            error_to_file=FALSE
        ),
        "invalid query ID"
    )

    writeLines("Row1\tmalformed_target\t99", hits_file)
    expect_error(
        get.vsearch.results(
            vsearch_outfile=hits_file,
            error_to_file=FALSE
        ),
        "malformed target ID"
    )

    writeLines("Row1\tref;;genus;;species_a\tbad", hits_file)
    expect_error(
        get.vsearch.results(
            vsearch_outfile=hits_file,
            error_to_file=FALSE
        ),
        "nonnumeric percent identity"
    )
})
