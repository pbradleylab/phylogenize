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

test_that("matrix.plm serial loop matches previous row-wise apply behavior", {
    set.seed(42)
    taxa <- paste0("taxon", seq_len(8))
    extra_taxa <- paste0("extra", seq_len(4))
    all_taxa <- c(taxa, extra_taxa)
    genes <- paste0("gene", seq_len(12))
    tree <- ape::rtree(length(all_taxa), tip.label=all_taxa)
    pheno <- stats::rnorm(length(all_taxa))
    names(pheno) <- all_taxa
    gene_matrix <- matrix(
        stats::rbinom(length(genes) * length(all_taxa), 1, 0.4),
        nrow=length(genes),
        dimnames=list(genes, all_taxa)
    )
    opts <- pz.resolve.options(
        ncl=1,
        separate_process=FALSE,
        meas_err=FALSE,
        error_to_file=FALSE
    )
    old_result <- suppressWarnings(apply(
        gene_matrix[genes, taxa, drop=FALSE],
        1,
        lm.fx.pv,
        p=pheno,
        tr=tree,
        restrict=taxa,
        meas_err=opts("meas_err")
    ))
    new_result <- suppressWarnings(matrix.plm(
        tree,
        gene_matrix,
        pheno,
        method=lm.fx.pv,
        restrict.taxa=taxa,
        restrict.ff=genes,
        .opts=opts
    ))

    expect_identical(new_result, old_result)
})

test_that("matrix.plm preserves builtin and custom method inputs", {
    set.seed(43)
    taxa <- paste0("taxon", seq_len(8))
    extra_taxa <- paste0("extra", seq_len(4))
    all_taxa <- c(taxa, extra_taxa)
    genes <- paste0("gene", seq_len(12))
    tree <- ape::rtree(length(all_taxa), tip.label=all_taxa)
    pheno <- stats::rnorm(length(all_taxa))
    names(pheno) <- all_taxa
    gene_matrix <- matrix(
        stats::rbinom(length(genes) * length(all_taxa), 1, 0.4),
        nrow=length(genes),
        dimnames=list(genes, all_taxa)
    )
    opts <- pz.resolve.options(
        ncl=1,
        separate_process=FALSE,
        meas_err=FALSE,
        error_to_file=FALSE
    )

    old_phylolm <- suppressWarnings(apply(
        gene_matrix[genes, taxa, drop=FALSE],
        1,
        phylolm.fx.pv,
        p=pheno,
        tr=tree,
        restrict=taxa,
        meas_err=opts("meas_err")
    ))
    new_phylolm <- suppressWarnings(matrix.plm(
        tree,
        gene_matrix,
        pheno,
        method=phylolm.fx.pv,
        restrict.taxa=taxa,
        restrict.ff=genes,
        .opts=opts
    ))

    seen <- new.env(parent=emptyenv())
    seen$ok <- TRUE
    custom_method <- function(m, p, tr, restrict=NULL, meas_err=FALSE) {
        seen$ok <- seen$ok && identical(names(p), all_taxa) &&
            identical(restrict, taxa)
        c(Estimate=sum(m[restrict]),
          p.value=sum(p[restrict]),
          StdErr=0,
          df=length(restrict) - 2)
    }
    old_custom <- apply(
        gene_matrix[genes, taxa, drop=FALSE],
        1,
        custom_method,
        p=pheno,
        tr=tree,
        restrict=taxa,
        meas_err=opts("meas_err")
    )
    new_custom <- matrix.plm(
        tree,
        gene_matrix,
        pheno,
        method=custom_method,
        restrict.taxa=taxa,
        restrict.ff=genes,
        .opts=opts
    )

    expect_identical(new_phylolm, old_phylolm)
    expect_identical(new_custom, old_custom)
    expect_true(seen$ok)
})
