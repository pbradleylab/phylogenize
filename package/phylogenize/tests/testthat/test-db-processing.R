test_that("change.presence.tax.level splits sparse matrices by requested taxon", {
    tax <- tibble::tibble(
        cluster=paste0("s", 1:4),
        phylum="p1",
        family=c("f1", "f1", "f2", "f2")
    )
    binary <- list(
        p1=Matrix::Matrix(
            matrix(
                c(1, 0, 1, 0,
                  0, 1, 0, 1),
                nrow=2,
                byrow=TRUE,
                dimnames=list(c("gene1", "gene2"), paste0("s", 1:4))
            ),
            sparse=TRUE
        )
    )

    split_binary <- change.presence.tax.level(binary, "family", tax)

    expect_setequal(names(split_binary), c("f1", "f2"))
    expect_equal(colnames(split_binary$f1), c("s1", "s2"))
    expect_equal(colnames(split_binary$f2), c("s3", "s4"))
    expect_equal(rownames(split_binary$f1), c("gene1", "gene2"))
})

test_that("change.presence.tax.level skips missing taxonomy mappings", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(verbosity=0, error_to_file=FALSE)

    tax <- tibble::tibble(
        cluster=paste0("s", 1:2),
        phylum="p1",
        family=c("f1", "f1")
    )
    binary <- list(
        p1=Matrix::Matrix(
            matrix(
                c(1, 0),
                nrow=1,
                dimnames=list("gene1", paste0("s", 1:2))
            ),
            sparse=TRUE
        ),
        p2=Matrix::Matrix(
            matrix(
                c(1, 0),
                nrow=1,
                dimnames=list("gene1", paste0("s", 3:4))
            ),
            sparse=TRUE
        )
    )

    expect_warning(
        split_binary <- change.presence.tax.level(binary, "family", tax),
        "No data found for taxon classification p2"
    )

    expect_equal(names(split_binary), "f1")
    expect_equal(colnames(split_binary$f1), c("s1", "s2"))
})

test_that("change.tree.tax.level returns matching subtrees", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(verbosity=0, error_to_file=FALSE)

    tax <- tibble::tibble(
        cluster=paste0("s", 1:4),
        phylum="p1",
        family=c("f1", "f1", "f2", "f2")
    )
    tree <- list(p1=ape::read.tree(text="((s1:1,s2:1):1,(s3:1,s4:1):1);"))

    split_trees <- change.tree.tax.level(tree, "family", tax)

    expect_setequal(names(split_trees), c("f1", "f2"))
    expect_equal(split_trees$f1$tip.label, c("s1", "s2"))
    expect_equal(split_trees$f2$tip.label, c("s3", "s4"))
})

test_that("above_minimum_genes keeps genes with enough present and absent tips", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)
    pz.options(minimum=1, gene_min_frac=0.5, verbosity=0, error_to_file=FALSE)

    gene.presence <- list(
        TaxonA=Matrix::Matrix(
            matrix(
                c(1, 0, 1,
                  1, 1, 1,
                  0, 0, 0),
                nrow=3,
                byrow=TRUE,
                dimnames=list(c("g_keep", "g_all", "g_none"), c("s1", "s2", "s3"))
            ),
            sparse=TRUE
        )
    )
    trees <- list(TaxonA=ape::read.tree(text="(s1:1,s2:1,s3:1);"))

    filtered <- above_minimum_genes(gene.presence, trees)

    expect_equal(names(filtered), "TaxonA")
    expect_equal(rownames(filtered$TaxonA), "g_keep")
    expect_equal(colnames(filtered$TaxonA), c("s1", "s2", "s3"))
})

test_that("import.pz.db rejects missing database index", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    expect_error(
        import.pz.db(
            data_dir=tempdir(),
            db="missing-db",
            error_to_file=FALSE
        ),
        "Database index not found"
    )
})

test_that("import.pz.db validates database index columns", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempfile("db-index-")
    dir.create(tmp)
    writeLines(c(
        "database,genes,trees,taxonomy",
        "test-db,genes.rds,trees.rds,taxonomy.csv"
    ), file.path(tmp, "databases.csv"))

    expect_error(
        import.pz.db(
            data_dir=tmp,
            db="test-db",
            error_to_file=FALSE
        ),
        "missing required column"
    )
})

test_that("import.pz.db validates referenced database files before loading", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempfile("db-files-")
    dir.create(tmp)
    writeLines(c(
        "database,genes,trees,taxonomy,functions",
        "test-db,genes.rds,trees.rds,taxonomy.csv,functions.csv"
    ), file.path(tmp, "databases.csv"))

    expect_error(
        import.pz.db(
            data_dir=tmp,
            db="test-db",
            error_to_file=FALSE
        ),
        "Database file\\(s\\) not found"
    )
})

test_that("import.pz.db validates loaded taxonomy schema", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempfile("db-schema-")
    dir.create(tmp)
    writeLines(c(
        "database,genes,trees,taxonomy,functions",
        "test-db,genes.rds,trees.rds,taxonomy.csv,functions.csv"
    ), file.path(tmp, "databases.csv"))
    saveRDS(
        list(p1=Matrix::Matrix(
            matrix(
                c(1, 0),
                nrow=1,
                dimnames=list("gene1", c("s1", "s2"))
            ),
            sparse=TRUE
        )),
        file.path(tmp, "genes.rds")
    )
    saveRDS(
        list(p1=ape::read.tree(text="(s1:1,s2:1);")),
        file.path(tmp, "trees.rds")
    )
    writeLines(c(
        "cluster,species,genus,family,order,class,phylum",
        "s1,s1,g1,f1,o1,c1,p1",
        "s2,s2,g1,f1,o1,c1,p1"
    ), file.path(tmp, "taxonomy.csv"))
    writeLines(c(
        "node_head,accession,function",
        "gene1,acc1,fxn1"
    ), file.path(tmp, "functions.csv"))

    expect_error(
        import.pz.db(
            data_dir=tmp,
            db="test-db",
            error_to_file=FALSE
        ),
        "Taxonomy file is missing required column"
    )
})

test_that("import.pz.db validates loaded gene function schema", {
    old_opts <- pz.options()
    on.exit(do.call(pz.options, old_opts), add=TRUE)

    tmp <- tempfile("db-function-schema-")
    dir.create(tmp)
    writeLines(c(
        "database,genes,trees,taxonomy,functions",
        "test-db,genes.rds,trees.rds,taxonomy.csv,functions.csv"
    ), file.path(tmp, "databases.csv"))
    saveRDS(
        list(p1=Matrix::Matrix(
            matrix(
                c(1, 0),
                nrow=1,
                dimnames=list("gene1", c("s1", "s2"))
            ),
            sparse=TRUE
        )),
        file.path(tmp, "genes.rds")
    )
    saveRDS(
        list(p1=ape::read.tree(text="(s1:1,s2:1);")),
        file.path(tmp, "trees.rds")
    )
    writeLines(c(
        "cluster,species,genus,family,order,class,phylum,domain",
        "s1,s1,g1,f1,o1,c1,p1,d1",
        "s2,s2,g1,f1,o1,c1,p1,d1"
    ), file.path(tmp, "taxonomy.csv"))
    writeLines(c(
        "node_head,accession,description",
        "gene1,acc1,fxn1"
    ), file.path(tmp, "functions.csv"))

    expect_error(
        import.pz.db(
            data_dir=tmp,
            db="test-db",
            error_to_file=FALSE
        ),
        "Gene function file is missing required column"
    )
})
