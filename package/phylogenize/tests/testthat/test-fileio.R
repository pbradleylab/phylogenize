test_that("tabular import reads explicit abundance and metadata files", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempdir()
    abundance_file <- file.path(tmp, "abundance.tsv")
    metadata_file <- file.path(tmp, "metadata.tsv")

    writeLines(c(
        "taxon\ts1\ts2\ts3\ts4",
        "sp1\t1\t0\t2\t3",
        "sp2\t0\t4\t1\t5"
    ), abundance_file)
    writeLines(c(
        "sample\tdataset\tenv",
        "s1\td1\tA",
        "s2\td1\tA",
        "s3\td1\tB",
        "s4\td1\tB"
    ), metadata_file)

    abd.meta <- read.abd.metadata(
        input_format="tabular",
        abundance_file=abundance_file,
        metadata_file=metadata_file,
        sample_column="sample",
        dset_column="dataset",
        env_column="env",
        which_envir="A",
        which_phenotype="abundance",
        categorical=TRUE,
        error_to_file=FALSE
    )

    expect_equal(rownames(abd.meta$mtx), c("sp1", "sp2"))
    expect_equal(colnames(abd.meta$mtx), paste0("s", 1:4))
    expect_equal(as.numeric(abd.meta$mtx["sp1", ]), c(1, 0, 2, 3))
    expect_equal(as.character(abd.meta$metadata$sample), paste0("s", 1:4))
    expect_true(is.factor(abd.meta$metadata$env))
})

test_that("tabular import binarizes non-abundance phenotypes", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempdir()
    abundance_file <- file.path(tmp, "presence-abundance.tsv")
    metadata_file <- file.path(tmp, "presence-metadata.tsv")

    writeLines(c(
        "taxon\ts1\ts2\ts3\ts4",
        "sp1\t1\t0\t2\t3",
        "sp2\t0\t4\t0\t5"
    ), abundance_file)
    writeLines(c(
        "sample\tdataset\tenv",
        "s1\td1\tA",
        "s2\td1\tA",
        "s3\td1\tB",
        "s4\td1\tB"
    ), metadata_file)

    abd.meta <- read.abd.metadata(
        input_format="tabular",
        abundance_file=abundance_file,
        metadata_file=metadata_file,
        sample_column="sample",
        dset_column="dataset",
        env_column="env",
        which_envir="A",
        which_phenotype="prevalence",
        core_method="phylogenize",
        categorical=TRUE,
        error_to_file=FALSE
    )

    expect_s4_class(abd.meta$mtx, "Matrix")
    expect_equal(as.numeric(abd.meta$mtx["sp1", ]), c(1, 0, 1, 1))
    expect_equal(as.numeric(abd.meta$mtx["sp2", ]), c(0, 1, 0, 1))
})

test_that("tabular import rejects nonnumeric abundance values", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempdir()
    abundance_file <- file.path(tmp, "bad-abundance.tsv")
    metadata_file <- file.path(tmp, "bad-metadata.tsv")

    writeLines(c(
        "taxon\ts1\ts2\ts3\ts4",
        "sp1\t1\tbad\t2\t3",
        "sp2\t0\t4\t1\t5"
    ), abundance_file)
    writeLines(c(
        "sample\tdataset\tenv",
        "s1\td1\tA",
        "s2\td1\tA",
        "s3\td1\tB",
        "s4\td1\tB"
    ), metadata_file)

    expect_error(
        read.abd.metadata(
            input_format="tabular",
            abundance_file=abundance_file,
            metadata_file=metadata_file,
            sample_column="sample",
            dset_column="dataset",
            env_column="env",
            which_envir="A",
            which_phenotype="abundance",
            categorical=TRUE,
            error_to_file=FALSE
        ),
        "nonnumeric value"
    )
})

test_that("harmonization trims sample identifiers before matching", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    abd.meta <- list(
        mtx=matrix(
            c(1, 0, 2, 3,
              0, 4, 1, 5),
            nrow=2,
            byrow=TRUE,
            dimnames=list(c("sp1", "sp2"), c("s1 ", "s2", " s3", "s4"))
        ),
        metadata=tibble::tibble(
            sample=c("s1", " s2", "s3 ", "s4"),
            dataset=c("d1", "d1", "d1", "d1"),
            env=c("A", "A", "B", "B")
        )
    )

    harmonized <- harmonize.abd.meta(
        abd.meta,
        sample_column="sample",
        dset_column="dataset",
        env_column="env",
        which_phenotype="abundance",
        error_to_file=FALSE
    )

    expect_equal(colnames(harmonized$mtx), paste0("s", 1:4))
    expect_equal(harmonized$metadata$sample, paste0("s", 1:4))
    expect_equal(as.numeric(harmonized$mtx["sp1", ]), c(1, 0, 2, 3))
})

test_that("harmonization keeps metadata aligned to retained abundance columns", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    abd.meta <- list(
        mtx=matrix(
            c(1, 2, 3, 4, 0,
              5, 6, 7, 8, 0),
            nrow=2,
            byrow=TRUE,
            dimnames=list(c("sp1", "sp2"), paste0("s", 1:5))
        ),
        metadata=data.frame(
            sample=c("s3", "s1", "s5", "s4", "s2"),
            dataset=rep("d1", 5),
            env=c("B", "A", "A", "B", "A")
        )
    )

    harmonized <- harmonize.abd.meta(
        abd.meta,
        sample_column="sample",
        dset_column="dataset",
        env_column="env",
        which_phenotype="abundance",
        error_to_file=FALSE
    )

    expect_equal(harmonized$metadata$sample, colnames(harmonized$mtx))
    expect_equal(colnames(harmonized$mtx), c("s3", "s1", "s4", "s2"))
    expect_false("s5" %in% harmonized$metadata$sample)
})

test_that("BIOM import test documents optional external dependency", {
    skip("BIOM round-trip helper still depends on obsolete biom_dir/external biom tooling.")
})

test_that("16S mapping test documents optional external reference data", {
    skip("16S mapping requires external FASTA and vsearch/appspam data not distributed with tests.")
})
