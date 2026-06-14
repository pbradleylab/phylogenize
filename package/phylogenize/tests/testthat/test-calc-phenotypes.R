calc_phenotype_fixture <- function() {
    list(
        abd.meta=list(
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
        ),
        pz.db=list(
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
    )
}

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

test_that("provided phenotype input is validated before use", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    fixture <- calc_phenotype_fixture()
    calculate_provided <- function(phenotype_file) {
        calculate_phenotypes(
            fixture$abd.meta,
            fixture$pz.db,
            which_phenotype="provided",
            phenotype_file=phenotype_file,
            treemin=2,
            error_to_file=FALSE
        )
    }

    expect_error(
        calculate_provided(file.path(tempdir(), "missing-phenotype.tsv")),
        "Phenotype file not found"
    )

    duplicate_file <- tempfile(fileext=".tsv")
    writeLines(c(
        "species\tphenotype",
        "s1\t1",
        " s1 \t2"
    ), duplicate_file)
    expect_error(
        calculate_provided(duplicate_file),
        "duplicate taxon IDs"
    )

    nonnumeric_file <- tempfile(fileext=".tsv")
    writeLines(c(
        "species\tphenotype",
        "s1\t1",
        "s2\tbad"
    ), nonnumeric_file)
    expect_error(
        calculate_provided(nonnumeric_file),
        "nonnumeric phenotype values"
    )

    no_overlap_file <- tempfile(fileext=".tsv")
    writeLines(c(
        "species\tphenotype",
        "x1\t1",
        "x2\t2"
    ), no_overlap_file)
    expect_error(
        calculate_provided(no_overlap_file),
        "No phenotype values matched database taxa"
    )
})

test_that("specificity prior files are resolved and validated", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    fixture <- calc_phenotype_fixture()
    calculate_specificity <- function(prior_file, working_dir=tempdir()) {
        calculate_phenotypes(
            fixture$abd.meta,
            fixture$pz.db,
            which_phenotype="specificity",
            prior_type="file",
            prior_file=prior_file,
            working_dir=working_dir,
            error_to_file=FALSE
        )
    }

    expect_error(
        calculate_specificity("missing-prior.tsv"),
        "Prior file not found"
    )

    tmp <- tempfile("prior-dir-")
    dir.create(tmp)

    missing_column_file <- file.path(tmp, "missing-column.tsv")
    writeLines(c(
        "env\tweight",
        "A\t0.5",
        "B\t0.5"
    ), missing_column_file)
    expect_error(
        calculate_specificity("missing-column.tsv", working_dir=tmp),
        "missing required column"
    )

    duplicate_env_file <- file.path(tmp, "duplicate-env.tsv")
    writeLines(c(
        "env\tprior",
        "A\t0.5",
        " A \t0.5"
    ), duplicate_env_file)
    expect_error(
        calculate_specificity("duplicate-env.tsv", working_dir=tmp),
        "duplicate environments"
    )

    nonnumeric_prior_file <- file.path(tmp, "nonnumeric-prior.tsv")
    writeLines(c(
        "env\tprior",
        "A\t0.5",
        "B\tbad"
    ), nonnumeric_prior_file)
    expect_error(
        calculate_specificity("nonnumeric-prior.tsv", working_dir=tmp),
        "nonnumeric prior values"
    )
})

test_that("specificity prior files must cover retained environments", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    abd.meta <- list(
        mtx=matrix(
            c(1, 0, 1, 0,
              0, 1, 0, 1),
            nrow=2,
            byrow=TRUE,
            dimnames=list(c("s1", "s2"), paste0("sample", 1:4))
        ),
        metadata=data.frame(
            sample=paste0("sample", 1:4),
            dataset=rep("d1", 4),
            env=factor(c("A", "A", "B", "B"))
        )
    )
    priors <- data.frame(env="A", prior=1)

    expect_error(
        calc.ess(
            abd.meta,
            pdata=priors,
            prior_type="file",
            env_column="env",
            dset_column="dataset",
            sample_column="sample",
            which_envir="A",
            error_to_file=FALSE
        ),
        "does not include retained environment"
    )
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

test_that("correl.clr returns finite values for perfect negative correlations", {
    trait <- c(-1, -0.5, 0.5, 1)
    baseline <- 1.5
    taxon1 <- exp(2 * trait) * baseline - 0.5
    taxon2 <- rep(1, length(trait))
    abd.meta <- list(
        mtx=matrix(
            c(taxon1, taxon2),
            nrow=2,
            byrow=TRUE,
            dimnames=list(c("taxon1", "taxon2"), paste0("s", 1:4))
        ),
        metadata=data.frame(
            sample=paste0("s", 1:4),
            dataset=rep("d1", 4),
            env=trait
        )
    )

    observed <- correl.clr(
        abd.meta,
        sample_column="sample",
        dset_column="dataset",
        env_column="env",
        error_to_file=FALSE
    )

    expect_true(all(is.finite(observed)))
    expect_lt(observed[["taxon2"]], 0)
})

test_that("differential abundance restores numeric taxon names without duplicates", {
    restored <- restore.diff.abund.taxon.names(
        sample_pheno=c(X123=1.5, taxonA=-0.5),
        sample_sd=c(X123=0.2, taxonA=0.3),
        original_names=c("123", "taxonA")
    )

    expect_equal(restored$old_name, c("123", "taxonA"))
    expect_equal(restored$new_name, c("X123", "taxonA"))
    expect_equal(restored$pheno, c(1.5, -0.5))
    expect_equal(restored$sd, c(0.2, 0.3))
    expect_equal(nrow(restored), 2)
})

test_that("differential abundance rejects ambiguous numeric taxon remapping", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(error_to_file=FALSE)

    expect_error(
        restore.diff.abund.taxon.names(
            sample_pheno=c(X123=1.5),
            sample_sd=c(X123=0.2),
            original_names=c("123", "X123")
        ),
        "ambiguous renamed value"
    )
})

test_that("nonparallel.results.generator returns empty results when all genes are filtered", {
    gene.matrix <- matrix(
        c(1, 1, 1,
          0, 0, 0),
        nrow=2,
        byrow=TRUE,
        dimnames=list(c("gene_present", "gene_absent"), paste0("s", 1:3))
    )
    tree <- ape::read.tree(text="(s1:1,s2:1,s3:1);")
    pheno <- c(s1=0, s2=1, s3=2)

    observed <- nonparallel.results.generator(
        gene.matrix=gene.matrix,
        tree=tree,
        taxa=paste0("s", 1:3),
        pheno=pheno,
        method=function(...) stop("method should not be called"),
        remove.low.variance=TRUE,
        use.for.loop=TRUE
    )

    expect_equal(dim(observed), c(4L, 0L))
    expect_equal(rownames(observed), c("Estimate", "p.value", "StdErr", "df"))
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
