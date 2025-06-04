#' @import Matrix
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
    
    prev.options <- pz.options()
    do.call(pz.options, list(...))
    poms_flag <- tolower(pz.options("core_method"))=="poms"
    pz.options(working_dir=normalizePath(getwd()))
    pz.options(out_dir=normalizePath(pz.options("out_dir")))
    if (!(dir.exists(pz.options("out_dir")))) dir.create(pz.options("out_dir"))
    
    pz.message("Running the core of phylogenize...")
    core <- phylogenize_core(do_POMS=poms_flag,
                             do_enr=do_enr)
    
    if (pz.options("rds_output_file") != "") {
        core_path <- file.path(pz.options("out_dir"), 
                               pz.options("rds_output_file"))
        saveRDS(object=core, file=core_path)
    }
    
    pz.message("Generating report...")
    render_core_report(core,
                       output_file=pz.options("output_file"),
                       do_cache=do_cache,
                       reset_after=reset_after)
    if (reset_after) {
        do.call(pz.options, prev.options)
    }
    return(core)
}


#' Run just the core of the *phylogenize* pipeline start to finish, without
#' making all of the plots. Useful when incorporating phylogenize into a longer
#' analysis, or when you don't want to wait for everything to render.
#'
#' @details Note: this function uses package-wide options (see
#'   \code{?pz.options}), which can be overridden using the \code{...} argument.
#'
#' @param do_POMS Run the POMS algorithm instead of phylogenetic regression
#'   (default: FALSE).
#' @param do_enr Run enrichment analysis. Can skip to save time (default: TRUE)
#' @param force_return_data Return the input data and metadata, even if not
#'   using POMS (default: FALSE).
#' @param p.method Function that returns the effect size and p-value per gene
#'   (default: phylolm.fx.pv). Can be overridden to use (for example) a plain
#'   linear model instead, for the sake of comparison.
#' @param ... Parameters to override defaults.
#' @export
phylogenize_core <- function(
        do_POMS=FALSE,
        do_enr=TRUE,
        force_return_data=FALSE,
        p.method=phylogenize:::phylolm.fx.pv,
        ...
) {
    options <- clone_and_merge(PZ_OPTIONS, ...)
    list_pheno <- data_to_phenotypes(
        save_data = (!do_POMS || override_save_data),
        ...
    )
    list_signif <- get_all_associated_genes(
        list_pheno$pz.db,
        list_pheno$phenotype_results$phenotype,
        do_POMS,
        p.method,
        ...)
    if (!do_enr) {
        return(list(list_pheno=list_pheno,
                    list_signif=list_signif,
                    options=options))
    }
    enr_tbls <- get_enrichment_tbls(list_signif[["signif"]],
                                    list_signif[["signs"]],
                                    list_pheno[["pz.db"]],
                                    list_signif[["results.matrix"]],
                                    export=TRUE,
                                    print_out=TRUE,
                                    ...)
    return(list(list_pheno=list_pheno,
                list_signif=list_signif,
                enr_tbls=enr_tbls,
                options=options))
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
                               ...) {
    
    prev.options <- pz.options()
    if ("options" %in% names(core)) {
        do.call(pz.options, core[["options"]]())
    }
    do.call(pz.options, list(...))
    pz.options(working_dir=normalizePath(getwd()))
    pz.options(out_dir=normalizePath(pz.options("out_dir")))
    if (!dir.exists(pz.options('out_dir'))) {
        dir.create(pz.options('out_dir'))
    }
    p <- pz.options()
    if (verbose) {
        for (n in names(p)) {
            message(paste0(n, ": ", p[[n]]))
        }
    }
    file.copy(system.file("rmd",
                          report_input,
                          package="phylogenize"),
              pz.options("out_dir"),
              overwrite=TRUE)
    e <- environment()
    r <- rmarkdown::render(file.path(pz.options("out_dir"),
                                     report_input),
                           output_file=basename(output_file),
                           output_dir=pz.options("out_dir"),
                           output_options=list(
                               cache.path=pz.options("out_dir"),
                               cache=do_cache
                           ),
                           intermediates_dir=pz.options("out_dir"),
                           knit_root_dir=pz.options("out_dir"),
                           envir = e)
    if (reset_after) {
        do.call(pz.options, prev.options)
    }
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