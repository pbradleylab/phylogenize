test_that("make.results.matrix converts result matrices to long rows", {
    results <- list(
        TaxonA=matrix(
            c(1.5, -2.0,
              0.01, 0.20,
              0.10, 0.30,
              8.00, 8.00),
            nrow=4,
            byrow=TRUE,
            dimnames=list(c("Estimate", "p.value", "StdErr", "df"), c("g1", "g2"))
        )
    )

    result_tbl <- make.results.matrix(results)

    expect_s3_class(result_tbl, "tbl_df")
    expect_equal(names(result_tbl),
                 c("taxon", "gene", "effect.size", "p.value", "std.err", "df"))
    expect_equal(result_tbl$taxon, c("TaxonA", "TaxonA"))
    expect_equal(result_tbl$gene, c("g1", "g2"))
    expect_equal(as.numeric(result_tbl$effect.size), c(1.5, -2.0))
    expect_equal(as.numeric(result_tbl$p.value), c(0.01, 0.20))
})

test_that("make.results.matrix flattens repermulize-style result data frames", {
    results <- list(
        TaxonA=data.frame(
            g1=I(list(1.5, 0.01, NA_real_, NA_real_)),
            g2=I(list(-2.0, 0.20, NA_real_, NA_real_)),
            row.names=c("Estimate", "p.value", "StdErr", "df"),
            check.names=FALSE
        )
    )

    result_tbl <- make.results.matrix(results)

    expect_s3_class(result_tbl, "tbl_df")
    expect_false(any(vapply(result_tbl, is.list, logical(1))))
    expect_equal(result_tbl$taxon, c("TaxonA", "TaxonA"))
    expect_equal(result_tbl$gene, c("g1", "g2"))
    expect_equal(unname(result_tbl$effect.size), c(1.5, -2.0))
    expect_equal(unname(result_tbl$p.value), c(0.01, 0.20))
    expect_equal(unname(result_tbl$std.err), c(NA_real_, NA_real_))
    expect_equal(unname(result_tbl$df), c(NA_real_, NA_real_))
})

test_that("add.sig.descs annotates genes without tidyr warnings", {
    sigs <- list(
        TaxonA=c("g1", "g2"),
        TaxonB="g3",
        TaxonC=character(0)
    )
    gene_to_fxn <- tibble::tibble(
        gene=c("g1", "g2", "g3"),
        accession=c("K00001", "K00002", "K00003"),
        `function`=c("first gene", "second gene", "third gene")
    )

    expect_warning(
        out <- add.sig.descs(c("TaxonA", "TaxonB"), sigs, gene_to_fxn),
        NA
    )

    expect_equal(out$taxon, c("TaxonA", "TaxonA", "TaxonB"))
    expect_equal(out$gene, c("g1", "g2", "g3"))
    expect_equal(out$description, gene_to_fxn$`function`)
})

test_that("add.sig.descs handles taxa without significant genes", {
    sigs <- list(
        TaxonA=character(0),
        TaxonB=character(0)
    )
    gene_to_fxn <- tibble::tibble(
        gene=c("g1", "g2"),
        accession=c("K00001", "K00002"),
        `function`=c("first gene", "second gene")
    )

    out <- add.sig.descs(character(0), sigs, gene_to_fxn)

    expect_s3_class(out, "tbl_df")
    expect_equal(nrow(out), 0)
    expect_equal(names(out),
                 c("taxon", "gene", "accession", "description"))
})

test_that("qvals uses configured p.adjust methods", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    p <- c(0.01, 0.02, 0.50)

    pz.options(fdr_method="BH")
    expect_equal(phylogenize:::qvals(p), p.adjust(p, "BH"))

    pz.options(fdr_method="BY")
    expect_equal(phylogenize:::qvals(p), p.adjust(p, "BY"))

    pz.options(fdr_method="BH")
    expect_equal(phylogenize:::qvals(p, fdr_method="BY"), p.adjust(p, "BY"))
})

test_that("qvals can use resolved options without global mutation", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    p <- c(0.01, 0.02, 0.50)
    opts <- pz.resolve.options(fdr_method="BY")
    pz.options(fdr_method="BH")

    expect_equal(phylogenize:::qvals(p, .opts=opts), p.adjust(p, "BY"))
})

test_that("qvals validates methods and p-values", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    p <- c(g1=0.01, g2=NA, g3=0.50)
    observed <- phylogenize:::qvals(p, fdr_method="BH", error_to_file=FALSE)

    expect_equal(names(observed), names(p))
    expect_true(is.na(observed[["g2"]]))
    expect_equal(observed[c("g1", "g3")], p.adjust(p[c("g1", "g3")], "BH"))
    all_missing <- c(g1=NA_real_, g2=NA_real_)
    missing_observed <- phylogenize:::qvals(all_missing, fdr_method="BH",
                                            error_to_file=FALSE)
    expect_equal(names(missing_observed), names(all_missing))
    expect_true(all(is.na(missing_observed)))
    expect_equal(phylogenize:::qvals(numeric(0), fdr_method="BH"),
                 numeric(0))
    matrix_p <- matrix(
        0.01,
        nrow=1,
        dimnames=list("p.value", "g_matrix")
    )
    expect_equal(names(phylogenize:::qvals(matrix_p, fdr_method="BH")),
                 "g_matrix")

    expect_error(
        phylogenize:::qvals(c(0.01, 0.02), fdr_method="bonferroni",
                            error_to_file=FALSE),
        "Invalid fdr_method"
    )
    expect_warning(
        invalid_observed <- phylogenize:::qvals(
            c(g1=0.01, g2=Inf, g3=NaN, g4=0.20, g5=-0.1, g6=1.2),
            fdr_method="BH",
            error_to_file=FALSE
        ),
        "Marking 3 invalid p-value\\(s\\) as NA"
    )
    expect_equal(
        invalid_observed,
        c(g1=0.02, g2=NA_real_, g3=NA_real_,
          g4=0.20, g5=NA_real_, g6=NA_real_)
    )
})

test_that("FDR method snapshots preserve q-values, names, and missingness", {
    p <- c(
        g_strong=0.001,
        g_mid=0.020,
        g_border=0.049,
        g_null=0.200,
        g_missing=NA_real_,
        g_high=0.800
    )
    expected_bh <- c(
        g_strong=0.00500000,
        g_mid=0.05000000,
        g_border=0.0816666666666667,
        g_null=0.25000000,
        g_missing=NA_real_,
        g_high=0.80000000
    )
    expected_by <- c(
        g_strong=0.0114166666666667,
        g_mid=0.1141666666666667,
        g_border=0.1864722222222222,
        g_null=0.5708333333333333,
        g_missing=NA_real_,
        g_high=1.00000000
    )

    observed_bh <- phylogenize:::qvals(p, fdr_method="BH",
                                       error_to_file=FALSE)
    observed_by <- phylogenize:::qvals(p, fdr_method="BY",
                                       error_to_file=FALSE)
    expect_warning(
        observed_qvalue <- phylogenize:::qvals(p, fdr_method="qvalue",
                                               error_to_file=FALSE),
        "Trying lambda=0"
    )

    expect_named(observed_bh, names(p))
    expect_named(observed_by, names(p))
    expect_named(observed_qvalue, names(p))
    expect_equal(observed_bh, expected_bh, tolerance=1e-8)
    expect_equal(observed_by, expected_by, tolerance=1e-8)
    expect_equal(observed_qvalue, expected_bh, tolerance=1e-8)
    expect_true(is.na(observed_bh[["g_missing"]]))
    expect_true(is.na(observed_by[["g_missing"]]))
    expect_true(is.na(observed_qvalue[["g_missing"]]))
})

test_that("FDR method choice changes significant hits on fixed result objects", {
    results <- list(
        TaxonFDR=matrix(
            c(1.5, 1.2, -1.1, 0.8,
              0.001, 0.020, 0.049, 0.200,
              0.10, 0.10, 0.10, 0.10,
              20.0, 20.0, 20.0, 20.0),
            nrow=4,
            byrow=TRUE,
            dimnames=list(
                c("Estimate", "p.value", "StdErr", "df"),
                c("g_strong", "g_mid", "g_border", "g_null")
            )
        )
    )
    cuts <- c(strong=0.05)
    opts_bh <- pz.resolve.options(fdr_method="BH", error_to_file=FALSE)
    opts_by <- pz.resolve.options(fdr_method="BY", error_to_file=FALSE)

    sigs_bh <- make.sigs(results, cuts=cuts, min.fx=0, .opts=opts_bh)
    sigs_by <- make.sigs(results, cuts=cuts, min.fx=0, .opts=opts_by)
    signs <- make.signs(results)

    expect_named(signs$TaxonFDR,
                 c("g_strong", "g_mid", "g_border", "g_null"))
    expect_equal(sigs_bh$TaxonFDR$strong, c("g_strong", "g_mid"))
    expect_equal(sigs_by$TaxonFDR$strong, "g_strong")
    expect_equal(make.pos.sig(sigs_bh, signs)$TaxonFDR,
                 c("g_strong", "g_mid"))
    expect_equal(make.pos.sig(sigs_by, signs)$TaxonFDR, "g_strong")
})

test_that("calc.alpha.power compares named rejected tests to named truth sets", {
    pvs <- c(g_null=0.01, g_alt=0.02, g_miss=0.50, g_alt2=0.20)

    out <- calc.alpha.power(
        pvs,
        null=c("g_null", "g_miss"),
        alt=c("g_alt", "g_alt2"),
        alpha=0.05
    )

    expect_equal(unname(out["r"]), 0.5)
    expect_equal(unname(out["p"]), 0.5)
    expect_equal(unname(out["a"]), 0.5)

    filtered <- calc.alpha.power(
        pvs,
        null=c("g_null", "g_miss"),
        alt=c("g_alt", "g_alt2"),
        alpha=0.05,
        filter=c("g_alt", "g_miss")
    )

    expect_equal(unname(filtered["r"]), 0.25)
    expect_equal(unname(filtered["p"]), 1)
    expect_equal(unname(filtered["a"]), 0)
})

test_that("negative-only significant taxa are retained after thresholding", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(fdr_method="BH",
               min_fx=0,
               minimum=3,
               error_to_file=FALSE)

    results <- list(
        TaxonNeg=matrix(
            c(-1.5,
              0.001,
              0.10,
              8.00),
            nrow=4,
            dimnames=list(c("Estimate", "p.value", "StdErr", "df"), "g_neg")
        )
    )
    pz.db <- list(
        trees=list(TaxonNeg=ape::read.tree(text="(s1:1,s2:1,s3:1,s4:1,s5:1,s6:1);")),
        gene.presence=list(
            TaxonNeg=Matrix::Matrix(
                matrix(
                    c(1, 1, 1, 0, 0, 0),
                    nrow=1,
                    dimnames=list("g_neg", paste0("s", 1:6))
                ),
                sparse=TRUE
            )
        ),
        gene.to.fxn=tibble::tibble(
            gene="g_neg",
            accession="K00001",
            `function`="Mock negative gene"
        )
    )

    out <- get_signif_associated_genes(pz.db, results)

    expect_equal(out$phy.with.pos.sigs, character(0))
    expect_equal(out$phy.with.neg.sigs, "TaxonNeg")
    expect_equal(out$neg.sig.thresh$TaxonNeg, "g_neg")
    expect_equal(out$phy.with.sigs, "TaxonNeg")
})

test_that("significant gene thresholding uses binary gene_min_frac calls", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    opts <- pz.resolve.options(
        minimum=2,
        gene_min_frac=0.5,
        verbosity=0,
        error_to_file=FALSE
    )
    pz.db <- list(
        trees=list(TaxonA=ape::read.tree(text="(s1:1,s2:1,s3:1,s4:1);")),
        gene.presence=list(
            TaxonA=Matrix::Matrix(
                matrix(
                    c(1.0, 0.3, 0.3, 0.4,
                      1.0, 1.0, 0.0, 0.0),
                    nrow=2,
                    byrow=TRUE,
                    dimnames=list(c("g_fractional_singleton", "g_binary_ok"),
                                  paste0("s", 1:4))
                ),
                sparse=TRUE
            )
        )
    )
    pos.sig <- list(TaxonA=c("g_fractional_singleton", "g_binary_ok"))

    filtered <- threshold.pos.sigs(pz.db, "TaxonA", pos.sig, .opts=opts)

    expect_equal(filtered$TaxonA, "g_binary_ok")
})
