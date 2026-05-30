test_that("filter.incomplete.metadata.samples drops incomplete metadata rows and matching matrix columns", {
    metadata <- data.frame(
        env=factor(c("case", "control", NA, "case")),
        dataset=factor(c("d1", "d1", "d2", NA)),
        row.names=paste0("sample", 1:4)
    )
    mtx <- Matrix::Matrix(
        matrix(
            seq_len(8),
            nrow=2,
            dimnames=list(c("taxon1", "taxon2"), rownames(metadata))
        ),
        sparse=TRUE
    )

    filtered <- phylogenize:::filter.incomplete.metadata.samples(
        metadata,
        mtx,
        c("env", "dataset")
    )

    expect_equal(rownames(filtered$metadata), c("sample1", "sample2"))
    expect_equal(colnames(filtered$mtx), c("sample1", "sample2"))
    expect_false(any(is.na(filtered$metadata$env)))
    expect_false(any(is.na(filtered$metadata$dataset)))
    expect_equal(levels(filtered$metadata$env), c("case", "control"))
    expect_equal(levels(filtered$metadata$dataset), "d1")
})
