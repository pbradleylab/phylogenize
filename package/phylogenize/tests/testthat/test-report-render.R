test_that("render_core_report renders a minimal partial-result report", {
    skip_if(system.file("rmd", "phylogenize-report.Rmd", package="phylogenize") == "",
            "report template is not available through system.file")

    out_dir <- file.path(tempdir(), "phylogenize-report-render-test")
    dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

    tree <- ape::read.tree(text="((s1:1,s2:1):1,s3:1);")
    core <- list(
        list_pheno=list(
            pz.db=list(
                trees=list(TaxonA=tree),
                gene.presence=list(
                    TaxonA=Matrix::Matrix(
                        matrix(
                            c(1, 0, 1),
                            nrow=1,
                            dimnames=list("gene1", paste0("s", 1:3))
                        ),
                        sparse=TRUE
                    )
                ),
                taxonomy=tibble::tibble(
                    cluster=paste0("s", 1:3),
                    species=paste0("Species ", 1:3),
                    family="TaxonA"
                ),
                species=list(TaxonA=paste0("s", 1:3)),
                ntaxa=1
            ),
            phenotype_results=list(
                phenotype=c(s1=logit(0.10), s2=logit(0.20), s3=logit(0.30)),
                phenoP=NULL,
                mapped.observed=paste0("s", 1:3)
            )
        ),
        list_signif=NULL,
        enr_tbls=NULL
    )

    output <- render_core_report(
        core,
        output_file="minimal-report.html",
        out_dir=out_dir,
        skip_graphs=TRUE,
        error_to_file=FALSE,
        taxon_level="family",
        reset_after=TRUE
    )

    expect_true(file.exists(output))
})
