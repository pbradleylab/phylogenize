#' @import Matrix
#' @importFrom ggtree %<+%
#' @importFrom settings clone_and_merge
NULL

#--- Main ways to invoke phylogenize: ---#

#' Run *phylogenize* start to finish, then render an interactive report on the
#' results.
#'
#' @param do_cache Turn on or off Rmarkdown's caching (default: TRUE).
#' @param reset_after Reset global options to package defaults after running?
#'   (default: TRUE)
#' @param do_enr Generate enrichment tables? (default: TRUE)
#' @param ... Other parameters to override defaults (see ?pz.options for a full
#'   list).
#' @export
phylogenize <- function(do_cache=TRUE,
                        reset_after=TRUE,
                        do_enr=TRUE,
                        ...) {
    
    opts <- pz.resolve.options(...)
    opts <- pz.resolve.options(
        .opts=opts,
        working_dir=normalizePath(getwd()),
        out_dir=normalizePath(opts("out_dir"), mustWork=FALSE)
    )
    poms_flag <- tolower(opts("core_method"))=="poms"
    if (!(dir.exists(opts("out_dir"))) &&
        !dir.create(opts("out_dir"), recursive=TRUE, showWarnings=FALSE)) {
        pz.error(paste0("Could not create output directory: ",
                        opts("out_dir")),
                 .opts=opts)
    }
    
    pz.message("Running the core of phylogenize...", .opts=opts)
    core <- phylogenize_core(do_POMS=poms_flag,
                             do_enr=do_enr,
                             ...,
                             .opts=opts)
    
    if (opts("rds_output_file") != "") {
        core_path <- file.path(opts("out_dir"),
                               opts("rds_output_file"))
        pz.message(paste0("Saving core results to ", core_path), .opts=opts)
        saveRDS(object=core, file=core_path)
    }
    
    pz.message("Generating report...", .opts=opts)
    report_path <- file.path(
        opts("out_dir"),
        basename(opts("output_file"))
    )
    render_core_report(core,
                       output_file=opts("output_file"),
                       do_cache=do_cache,
                       reset_after=reset_after,
                       .opts=opts)
    pz.message(paste0("Report written to ", report_path), .opts=opts)
    return(core)
}


#' Run just the core of the *phylogenize* pipeline start to finish, without
#' making all of the plots. Useful when incorporating phylogenize into a longer
#' analysis, or when you don't want to wait for everything to render.
#'
#' @details Note: this function uses package-wide options (see
#'   \code{?pz.options}), which can be overridden using the \code{...} argument.
#'
#' @param do_enr Run enrichment analysis. Can skip to save time (default: TRUE)
#' @param do_POMS to run POMS method. (default: FALSE)
#' @param force_return_data Return the input data and metadata, even if not
#'   using POMS (default: FALSE).
#' @param p.method Function that returns the effect size and p-value per gene
#'   (default: phylolm.fx.pv). Can be overridden to use (for example) a plain
#'   linear model instead, for the sake of comparison.
#' @param ... Parameters to override defaults.
#' @export
phylogenize_core <- function(
        do_enr=TRUE,
	do_POMS = FALSE,
        force_return_data=FALSE,
        p.method=phylogenize:::phylolm.fx.pv,
        ...,
        .opts=NULL
) {
    dots <- list(...)
    opts <- pz.resolve.options(..., .opts=.opts)
    do_POMS <- do_POMS || (tolower(opts('core_method')) == "poms")
    if (do_POMS) {
        dots[["core_method"]] <- "poms"
        opts <- pz.resolve.options(.opts=opts, core_method="poms")
    }
    pz.message("I. Generating phenotypes...")
    list_pheno <- do.call(data_to_phenotypes,
                          c(list(save_data = (!do_POMS || force_return_data)),
                            list(.opts=opts),
                            dots))
    pz.message(paste0(
        "  .....Phenotypes generated for ",
        length(list_pheno[["phenotype_results"]][["phenotype"]]),
        " taxon/taxa across ",
        list_pheno[["pz.db"]][["ntaxa"]],
        " testable group(s)"
    ))
    pz.message("II. Performing association tests...")
    list_signif <- do.call(get_all_associated_genes,
                            c(list(list_pheno=list_pheno,
                                   p.method=p.method,
                                   .opts=opts),
                              dots))
    if (!is.null(list_signif)) {
        pz.message(paste0(
            "  .....Association testing returned ",
            nrow(list_signif[["results.matrix"]]),
            " gene-level result(s)"
        ))
    }
    if (!do_enr) {
        return(list(list_pheno=list_pheno,
                    list_signif=list_signif,
                    options=opts))
    }
    if (is.null(list_signif)) {
        pz.warning("Association testing returned no results; skipping enrichment.")
        return(list(list_pheno=list_pheno,
                    list_signif=list_signif,
                    enr_tbls=NULL,
                    options=opts))
    }
    #list_signif <- readRDS("list_signif.rds")
    pz.message("III. Performing enrichment tests...")
    enr_tbls <- get_enrichment_tbls(list_signif[["signif"]],
                                    list_signif[["signs"]],
                                    list_pheno[["pz.db"]],
                                    list_signif[["results.matrix"]],
                                    export=TRUE,
                                    print_out=TRUE,
                                    .opts=opts,
                                    ...)
    pz.message(paste0(
        "  .....Enrichment testing returned ",
        ifelse(is.null(enr_tbls), 0, nrow(enr_tbls)),
        " filtered enrichment row(s)"
    ))
    pz.message("Done!")
    return(list(list_pheno=list_pheno,
                list_signif=list_signif,
                enr_tbls=enr_tbls,
                options=opts))
}

#' Take the output of `phylogenize_core` and generate a report.
#'
#' @details Note: this function uses package-wide options (see
#'   \code{?pz.options}), which can be overridden using the \code{...} argument.
#'
#' @param core Output from `phylogenize_core()` (a named list; see
#'   `?phylogenize_core`).
#' @param output_file Path giving what to name the resulting HTML file.
#' @param report_input Optionally specify which notebook to knit (useful for
#'   testing).
#' @param do_cache Turn on or off Rmarkdown's caching. (Default: TRUE)
#' @param reset_after Reset global options after running? (Default: TRUE)
#' @param ... Other options to be passed to `pz.options` that will override
#'   options in `core`.
#' @export
render_core_report <- function(core,
                               output_file="phylogenize-report.html",
                               report_input="phylogenize-report.Rmd",
                               do_cache=FALSE,
                               reset_after=TRUE,
                               verbose=FALSE,
                               ...,
                               .opts=NULL) {
    
    prev.options <- pz.options()
    if (reset_after) {
        on.exit(do.call(pz.options, prev.options), add=TRUE)
    }
    if (!("list_pheno" %in% names(core)) || !("list_signif" %in% names(core))) {
        pz.error(paste0(
	  "`core` does not look like Phylogenize2 output; did you pass ",
	  "a path to an RDS file instead of reading it with readRDS?"
	))
    }
    base_opts <- .opts
    if (is.null(base_opts) && "options" %in% names(core)) {
        base_opts <- core[["options"]]
    }
    opts <- pz.resolve.options(..., .opts=base_opts)
    opts <- pz.resolve.options(
        .opts=opts,
        working_dir=normalizePath(getwd()),
        out_dir=normalizePath(opts("out_dir"), mustWork=FALSE)
    )
    if (!dir.exists(opts('out_dir')) &&
        !dir.create(opts('out_dir'), recursive=TRUE, showWarnings=FALSE)) {
        pz.error(paste0("Could not create output directory: ",
                        opts("out_dir")),
                 .opts=opts)
    }
    p <- opts()
    if (verbose) {
        for (n in names(p)) {
            message(paste0(n, ": ", p[[n]]))
        }
    }
    file.copy(system.file("rmd",
                          report_input,
                          package="phylogenize"),
              opts("out_dir"),
              overwrite=TRUE)
    pz.message(paste0("Rendering report from ", report_input),
               level=1,
               .opts=opts)
    e <- environment()
    do.call(pz.options, opts())
    r <- rmarkdown::render(file.path(opts("out_dir"),
                                     report_input),
                           output_file=basename(output_file),
                           output_dir=opts("out_dir"),
                           output_options=list(
                               cache.path=opts("out_dir"),
                               cache=do_cache
                           ),
                           intermediates_dir=opts("out_dir"),
                           knit_root_dir=opts("out_dir"),
                           envir = e)
    r
}


#' Add enrichments after the fact to a phylogenize core object.
#' 
#' @details Note: this function uses package-wide options (see
#'   \code{?pz.options}), which can be overridden using the \code{...} argument.
#'
#' @param core The named list obtained from running `phylogenize_core()`.
#' @export
augment_with_enrichments <- function(core) {
    core[["enr_tbls"]] <- get_enrichment_tbls(
        core[["list_signif"]][["signif"]],
        core[["list_signif"]][["signs"]],
        core[["list_pheno"]][["pz.db"]],
        core[["list_signif"]][["results.matrix"]],
        export=TRUE,
        print_out=TRUE)
    core
}


#' Perform either phylogenetic regression or POMS.
#'
#' @param list_pheno A list, the result of data_to_phenotypes
#' @param abd.meta A list consisting of a species abundance matrix and the
#'     metadata (from read.abd.metadata or data_to_phenotypes).
#' @param ... Parameters to override defaults.
#' @export
get_all_associated_genes <- function(list_pheno,
    p.method=phylolm.fx.pv,
    ...,
    .opts=NULL) {
        opts <- pz.resolve.options(..., .opts=.opts)
        do_POMS <- (tolower(opts('core_method')) == "poms")
        spec_taxa <- opts('only_specific_taxa')
        if (!do_POMS) {
            pz.message("  A) Getting all associated genes")
            phenotype <- list_pheno$phenotype_results$phenotype
            taxaN <- names(which(pheno_nonzero_var(phenotype, list_pheno$pz.db$species)))
            pz.message(paste0("  .....All valid taxa: ", paste(taxaN, collapse = ", ")), level = 1)
            if (!is.null(spec_taxa)) { taxaN <- intersect(taxaN, spec_taxa )}
            pz.message(paste0("  .....Valid taxa (after filtering): ", paste(taxaN, collapse=", ")), level=1)
            if (length(taxaN) == 0) pz.error("Error: no taxa found. If you provided any, check that they are spelled correctly")
	    pz.message("  .....Running plm on valid taxa")
            pz.message(paste0("  ..........Testing ", length(taxaN), " taxon/taxa"))
            if (opts('ncl') > 1) {
                pz.message(paste0(
                    "  ..........Using ",
                    opts('ncl'),
                    " parallel worker(s)"
                ))
                results <- result.wrapper.plm(taxa=taxaN,
                    pheno=phenotype,
                    tree=list_pheno$pz.db$trees[taxaN],
                    clusters=list_pheno$pz.db$species[taxaN],
                    proteins=list_pheno$pz.db$gene.presence[taxaN],
                    method=p.method,
                    ncl=opts('ncl'),
                    .opts=opts)
            } else {
                results <- mapply(nonparallel.results.generator,
                    list_pheno$pz.db$gene.presence[taxaN],
                    list_pheno$pz.db$trees[taxaN],
                    list_pheno$pz.db$species[taxaN],
                    as.list(taxaN),
                        MoreArgs=list(pheno=phenotype,
                            method=p.method,
                            use.for.loop=FALSE,
                            .opts=opts),
                        SIMPLIFY=FALSE)
            }
        } else {
	    pz.message("  A) Getting all associated genes")
            taxaN <- names(list_pheno$pz.db$species)
            pz.message(paste0("  .....All valid taxa: ", paste(taxaN, collapse = ", ")), level = 2)
            if (!is.null(spec_taxa)) {
                taxaN <- intersect(taxaN, spec_taxa)
            }
            pz.message(paste0("  .....Valid taxa: ", paste(taxaN, collapse=", ")), level=1)
            if (length(taxaN) == 0) pz.error("Error: no taxa found. If you provided any, check that they are spelled correctly")
            pz.message(paste0("  ..........Running POMS for ", length(taxaN), " taxon/taxa"))
            results <- result.wrapper.plm(taxa=taxaN,
                pheno=NULL,
                tree=list_pheno$pz.db$trees[taxaN],
                clusters=list_pheno$pz.db$species[taxaN],
                proteins=list_pheno$pz.db$gene.presence[taxaN],
                method=p.method,
                poms=TRUE,
                abd.meta=list_pheno$abd.meta,
                .opts=opts,
                ...
            )
        }
    
    pz.message(paste0("  ..........Rows that remain: ", length(results)))
    pz.message("  ..........Trimming any results that didn't get dropped")
    # trim out any that didn't get dropped
    result_lens <- vapply(results, length, 1L)
    results <- results[names(which(na.omit(result_lens>0)))]

    pz.message(paste0("  ..........Rows that remain after trimming: ", length(results)))
    results <- results[!vapply(results, is.null, logical(1))]
    pz.message(paste0("  ..........Rows that remain after checking if not null: ", length(results)))
    pz.message("  B) Summarizing significant associations")
    return(get_signif_associated_genes(list_pheno$pz.db, results, ..., .opts=opts))
}

#' Process genes by significance threshold.
#'
#' @param pz.db A database (typically obtained with \code{import.pz.db}).
#' @param results Results object from \code{get_all_associated_genes}.
#' @param ... Parameters to override defaults.
#' @export
get_signif_associated_genes <- function(pz.db,
                                        results,
                                        ...,
                                        .opts=NULL) {
    pz.message("  ..........Processing genes by significance threshold")
    opts <- pz.resolve.options(..., .opts=.opts)
    if (length(results) == 0) {
        pz.message("  ..........No successful repermulize results to process.")
        return(NULL)
    }

    pz.message("  ..........Making sigs")
    signif <- make.sigs(results, ..., .opts=opts)
    pz.message("  ..........Making signs")
    signs <- make.signs(results)
    pz.message("  ..........Getting positive sigs")
    pos.sig <- nonequiv.pos.sig(results, min_fx=opts('min_fx'), .opts=opts)
    pz.message("  ..........Getting negative sigs")
    neg.sig <- nonequiv.pos.sig(results, min_fx=opts('min_fx'), dir=-1, .opts=opts)
    pz.message("  ..........Make results matrix")
    results.matrix <- make.results.matrix(results) %>%
        dplyr::filter(!is.na(p.value)) %>%
        dplyr::group_by(taxon) %>%
        dplyr::mutate(q.value = qvals(p.value, ..., .opts=opts)) %>%
        dplyr::arrange(taxon, q.value) %>%
        dplyr::ungroup()
    pz.message(paste0(
        "  ..........Results matrix contains ",
        nrow(results.matrix),
        " row(s)"
    ))
    phy.with.pos.sigs <- names(which(sapply(pos.sig, length) > 0))
    phy.with.neg.sigs <- names(which(sapply(neg.sig, length) > 0))
    pz.message("  ..........Add sig descriptions")
    pos.sig.descs <- add.sig.descs(phy.with.pos.sigs, pos.sig, pz.db$gene.to.fxn)
    pos.sig.thresh <- threshold.pos.sigs(
        pz.db,
        phy.with.pos.sigs,
        pos.sig,
        ...,
        .opts=opts
    )
    pos.sig.thresh.descs <- add.sig.descs(phy.with.pos.sigs,
                                          pos.sig.thresh,
                                          pz.db$gene.to.fxn)
    neg.sig.descs <- add.sig.descs(phy.with.neg.sigs, neg.sig, pz.db$gene.to.fxn)
    neg.sig.thresh <- threshold.pos.sigs(
        pz.db,
        phy.with.neg.sigs,
        neg.sig,
        ...,
        .opts=opts
    )
    neg.sig.thresh.descs <- add.sig.descs(
        phy.with.neg.sigs,
        neg.sig.thresh,
        pz.db$gene.to.fxn
    )
    pz.message("  ..........Recalculate sigs")
    # recalculate, since some of these may go away
    phy.with.pos.sigs <- names(which(sapply(pos.sig.thresh, length) > 0))
    phy.with.neg.sigs <- names(which(sapply(neg.sig.thresh, length) > 0))
    phy.with.sigs <- union(phy.with.pos.sigs, phy.with.neg.sigs)
    pz.message(paste0(
        "  ..........Taxa with positive significant genes: ",
        length(phy.with.pos.sigs)
    ))
    pz.message(paste0(
        "  ..........Taxa with negative significant genes: ",
        length(phy.with.neg.sigs)
    ))
    pz.message(paste0("  ..........Remaining results: ", length(results)))
    return(list(
        results = results, #1
        signif = signif, #2
        signs = signs, #3
        pos.sig = pos.sig, #4
        results.matrix = results.matrix, #5
        phy.with.sigs = phy.with.sigs, #6
        pos.sig.descs = pos.sig.descs, #7
        pos.sig.thresh = pos.sig.thresh, #8
        pos.sig.thresh.descs = pos.sig.thresh.descs, #9
        neg.sig = neg.sig, #10
        neg.sig.descs = neg.sig.descs, #11
        neg.sig.thresh = neg.sig.thresh, #12
        neg.sig.thresh.descs = neg.sig.thresh.descs, #13
        phy.with.pos.sigs = phy.with.pos.sigs, #14
        phy.with.neg.sigs = phy.with.neg.sigs #15
    ))
}


#' Get enrichment tables.
#'
#' @param signif Significant genes.
#' @param signs Signs of effect sizes.
#' @param pz.db A database (typically obtained with \code{import.pz.db}).
#' @param results.matrix Matrix of full results.
#' @param export Write enrichment tables to disk? (Default: FALSE)
#' @param kegg_pw_data Optional downloaded KEGG pathway data from
#'   \code{clusterProfiler::download_KEGG("ko", keggType="KEGG")}.
#' @param kegg_mod_data Optional downloaded KEGG module data from
#'   \code{clusterProfiler::download_KEGG("ko", keggType="MKEGG")}.
#' @param use_kegg_cache Reuse downloaded KEGG data from disk when available?
#'   (Default: TRUE)
#' @param kegg_cache_file File name for the KEGG cache, relative to
#'   \code{pz.options("out_dir")} unless an absolute path is supplied.
#' @param ... Parameters to override defaults.
#' @export
get_enrichment_tbls <- function(signif,
                                signs,
                                pz.db,
                                results.matrix,
                                export=FALSE,
                                kegg_pw_data=NULL,
                                kegg_mod_data=NULL,
                                use_kegg_cache=TRUE,
                                kegg_cache_file="kegg-cache.rds",
                                ...,
                                .opts=NULL) {
    opts <- pz.resolve.options(..., .opts=.opts)
    pretty.enr.tbl <- NULL
    enr.overlap <- NULL
    kegg_cache_path <- kegg_cache_file
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", kegg_cache_path)) {
        kegg_cache_path <- file.path(opts("out_dir"), kegg_cache_path)
    }
    if (use_kegg_cache && (is.null(kegg_pw_data) || is.null(kegg_mod_data)) &&
        file.exists(kegg_cache_path)) {
        pz.message(paste0("  .....Reading KEGG cache from ", kegg_cache_path))
        kegg_cache <- tryCatch(readRDS(kegg_cache_path),
                               error=function(e) {
                                   pz.warning(paste(
                                       "Could not read KEGG cache:",
                                       conditionMessage(e)
                                   ))
                                   NULL
                               })
        if (is.null(kegg_pw_data) && !is.null(kegg_cache[["KEGG"]])) {
            kegg_pw_data <- kegg_cache[["KEGG"]]
        }
        if (is.null(kegg_mod_data) && !is.null(kegg_cache[["MKEGG"]])) {
            kegg_mod_data <- kegg_cache[["MKEGG"]]
        }
    }
    if (is.null(kegg_pw_data)) {
        pz.message("  .....Downloading KEGG pathway annotations")
        kegg_pw_data <- clusterProfiler::download_KEGG("ko", keggType="KEGG")
    }
    if (is.null(kegg_mod_data)) {
        pz.message("  .....Downloading KEGG module annotations")
        kegg_mod_data <- clusterProfiler::download_KEGG("ko", keggType="MKEGG")
    }
    if (use_kegg_cache) {
        pz.message(paste0("  .....Writing KEGG cache to ", kegg_cache_path), level=2)
        tryCatch({
            cache_dir <- dirname(kegg_cache_path)
            if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive=TRUE)
            saveRDS(list(KEGG=kegg_pw_data, MKEGG=kegg_mod_data),
                    kegg_cache_path)
        }, error=function(e) {
            pz.warning(paste("Could not write KEGG cache:",
                             conditionMessage(e)))
        })
    }
    enrichment.tbl <- multi.kegg.enrich(signif,
                                        signs,
                                        pz.db$gene.to.fxn,
                                        kegg_pw = kegg_pw_data,
                                        kegg_mod = kegg_mod_data)
    pz.message(paste0(
        "  .....Raw enrichment rows: ",
        ifelse(is.null(enrichment.tbl), 0, nrow(enrichment.tbl))
    ))
    if (!is.null(enrichment.tbl)) {
        enrichment.tbl <- dplyr::filter(enrichment.tbl,
                                        qvalue <= 0.25, enr.estimate > 1)
        pz.message(paste0(
            "  .....Filtered enrichment rows: ",
            nrow(enrichment.tbl)
        ))
        pretty.enr.tbl <- dplyr::select(enrichment.tbl,
                                        taxon,
                                        cutoff,
                                        ID,
                                        Description,
                                        qvalue,
                                        enr.estimate) %>%
            dplyr::rename(Gene_significance=cutoff,
                          Taxon=taxon,
                          Description=Description,
                          q_value=qvalue,
                          Enrichment_odds_ratio=enr.estimate) %>%
            dplyr::arrange(factor(Gene_significance, levels=names(signif[[1]])),
                           Taxon,
                           q_value)
        if (export) {
            pz.message("  .....Writing enrichment table")
            write.csv(file=file.path(opts('out_dir'),
                                     "enr-table.csv"),
                      pretty.enr.tbl)
        }
    }
    if (!is.null(pretty.enr.tbl)) {
        accession_to_fxn <- pz.db$gene.to.fxn %>%
            dplyr::select(accession, gene, `function`) %>%
            dplyr::distinct()
        if (nrow(enrichment.tbl) > 0) {
            effect_lookup <- results.matrix %>%
                dplyr::select(taxon, gene, effect.size) %>%
                dplyr::rename(mapped_gene=gene,
                              effectsize=effect.size)
            enr.overlap <- dplyr::select(enrichment.tbl,
                                         taxon,
                                         cutoff,
                                         ID,
                                         Description,
                                         geneID) %>%
                tidyr::separate_longer_delim(geneID, delim="/") %>%
                dplyr::left_join(.,
                                 accession_to_fxn,
                                 by=c("geneID"="accession")) %>%
                dplyr::rename(gene=geneID,
                              mapped_gene=gene,
                              description=`function`) %>%
                dplyr::left_join(effect_lookup,
                                 by=c("taxon", "mapped_gene"))
            if (export) {
                pz.message("  .....Writing enrichment overlap tables")
                write.csv(file=file.path(opts('out_dir'),
                                         "enr-overlaps.csv"),
                          enr.overlap)
                write.csv(file=file.path(opts('out_dir'),
                                         "enr-overlaps-sorted.csv"),
                          dplyr::arrange(enr.overlap,
                                         taxon,
                                         factor(cutoff,
                                                levels=names(signif[[1]])),
                                         desc(effectsize)))
            }
        }
    }
    return(pretty.enr.tbl)
}  
