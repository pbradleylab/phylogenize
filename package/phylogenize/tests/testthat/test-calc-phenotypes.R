test_that("calculate_phenotypes uses direct treemin override", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(which_phenotype="provided",
               treemin=10,
               error_to_file=FALSE)

    phenotype_file <- tempfile(fileext=".tsv")
    writeLines(c(
        "species\tphenotype",
        "s1\t1",
        "s2\t2",
        "s3\t3"
    ), phenotype_file)

    abd.meta <- list(
        mtx=Matrix::Matrix(
            matrix(
                c(1, 0, 0,
                  0, 1, 0,
                  0, 0, 1),
                nrow=3,
                byrow=TRUE,
                dimnames=list(paste0("s", 1:3), paste0("sample", 1:3))
            ),
            sparse=TRUE
        )
    )
    pz.db <- list(
        trees=list(TaxonA=ape::read.tree(text="(s1:1,s2:1,s3:1);")),
        gene.presence=list(
            TaxonA=Matrix::Matrix(
                matrix(
                    c(1, 1, 0),
                    nrow=1,
                    dimnames=list("gene1", paste0("s", 1:3))
                ),
                sparse=TRUE
            )
        )
    )

    phenotype_results <- calculate_phenotypes(
        abd.meta,
        pz.db,
        which_phenotype="provided",
        phenotype_file=phenotype_file,
        treemin=2,
        error_to_file=FALSE
    )

    expect_equal(phenotype_results$phenotype, c(s1=1, s2=2, s3=3))
})

test_that("quantile_normalize preserves prevalence phenotypes with NULL phenoP", {
    phenotype_results <- list(
        phenotype=c(taxon1=2, taxon2=1, taxon3=3),
        phenoP=NULL
    )

    observed <- quantile_normalize(phenotype_results)
    expected <- quant_norm(phenotype_results$phenotype)

    old_normed <- quant_norm(c(phenotype_results$phenoP,
                               phenotype_results$phenotype))
    old_phenotype <- old_normed[-1]

    expect_equal(observed$phenotype, expected)
    expect_null(observed$phenoP)
    expect_equal(names(observed$phenotype), names(phenotype_results$phenotype))
    expect_false(identical(observed$phenotype, old_phenotype))
})

test_that("logit_auc_pheno indexes environments by sample ID", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    abd.meta <- list(
        mtx=matrix(
            c(100, 90, 10, 20,
              10, 20, 100, 90,
              35, 45, 55, 65),
            nrow=3,
            byrow=TRUE,
            dimnames=list(
                paste0("taxon", 1:3),
                paste0("s", 1:4)
            )
        ),
        metadata=data.frame(
            sample=c("s3", "s1", "s4", "s2"),
            dataset=rep("d1", 4),
            env=factor(c("B", "A", "B", "A"), levels=c("A", "B"))
        )
    )

    clr_mtx <- clr(abd.meta$mtx, pc=1)
    env.cols <- colnames(abd.meta$mtx) %in% c("s1", "s2")
    expected <- apply(clr_mtx, 1, function(x) {
        stat <- wilcox.test(x[env.cols], x[!env.cols])$statistic
        logit(stat / (sum(env.cols) * sum(!env.cols)))
    })

    old_env.rows <- abd.meta$metadata$env == "A"
    old_result <- apply(clr_mtx, 1, function(x) {
        stat <- wilcox.test(x[which(old_env.rows)],
                            x[which(!old_env.rows)])$statistic
        logit(stat / (sum(old_env.rows) * sum(!old_env.rows)))
    })

    observed <- logit_auc_pheno(
        abd.meta,
        sample_column="sample",
        dset_column="dataset",
        env_column="env",
        which_envir="A",
        error_to_file=FALSE
    )

    expect_equal(observed, expected)
    expect_false(isTRUE(all.equal(observed, old_result)))
})
