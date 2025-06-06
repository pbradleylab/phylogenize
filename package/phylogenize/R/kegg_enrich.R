#' Given lists of significant genes (at different thresholds), effect sizes, and
#' gene set mappings, assemble a tbl of results using clusterProfiler.
#'
#' @param sigs List giving, per taxonomic group (outer) and per significance
#'   cutoff (inner), significant hits to test for enrichment.
#' @param signs List giving, per taxonomic group, signs of all gene effect
#'   sizes.
#' @param pid_to_ko Protein ID to KO mapping table (columns "node_head",
#'   "accession").
#' @param dirxn Count only genes with this effect sign as significant.
#' @param kegg_pw Optional downloaded data for KEGG pathways.
#' @param kegg_mod Optional downloaded data for KEGG modules.
#' @return A tbl giving Fisher's test p-values, q-values, effect sizes, and
#'   overlaps for all taxa and cutoffs.
#' @export
multi.kegg.enrich <- function(sigs, signs, pid_to_ko, dirxn=1,
                              kegg_pw=NULL, kegg_mod=NULL) {
  d <- sign(dirxn)
  if ((length(sigs) < 1) || (length(signs) < 1) || (nrow(pid_to_ko) < 1)) {
    pz.warning("Didn't find significant hits, signs, or mappings")
    return(c())
  }
  if (length(sigs[[1]]) < 1) {
    pz.warning("No significance cutoffs found")
    return(c())
  }
  tax_groups <- names(sigs)
  cutoffs <- names(sigs[[1]])
  background <- unique(pid_to_ko[["accession"]])
  enrichment_tbls <- purrr::pmap(
      list(s=sigs, sn=signs, tg=tax_groups),
      function(s, sn, tg) {
          all_cutoffs <- purrr::map2(s, names(s), function(sc, cn) {
              kegg.enrich.single(sc, sn, pid_to_ko, cn, d, tg,
                                 kegg_pw_data=kegg_pw,
                                 kegg_mod_data=kegg_mod)
          })
          dplyr::bind_rows(all_cutoffs)
      })
  out <- tryCatch({
      # Put together and add direction of enrichment
      dplyr::bind_rows(enrichment_tbls) %>% 
          tidyr::separate_wider_delim(cols = GeneRatio,
                                      delim = "/",
                                      names = c("GeneNum",
                                                "GeneDenom")) %>%
          tidyr::separate_wider_delim(cols = BgRatio,
                                      delim = "/",
                                      names = c("BgRatioNum",
                                                "BgRatioDenom")) %>%
          dplyr::mutate(enr.estimate =
                            (as.numeric(GeneNum) / as.numeric(GeneDenom)) /
                            (as.numeric(BgRatioNum) / as.numeric(BgRatioDenom)))
  }, error = function(e) {
    message("An error occurred: ", e$message)
    return(NULL)
  })
  return(out)
}

#' Convenience function to perform a clusterProfiler enrichment on downloaded KEGG pathway/module data.
#'
#' @param KOs KOs in significant set
#' @param background KOs in background set
#' @param qCut q-value cutoff
#' @param downloaded Data from `clusterProfiler::download_KEGG("ko", ...)`
#' @param ... Extra parameters passed to `clusterProfiler::enricher()`
#' @return See `clusterProfiler::enricher()`.
#' @export
enrich_downloaded_KEGG <- function(KOs, background, qCut, downloaded, ...) {
    clusterProfiler::enricher(
        gene = KOs,
        qvalueCutoff = qCut,
        universe = background,
        TERM2GENE = downloaded$KEGGPATHID2EXTID,
        TERM2NAME = downloaded$KEGGPATHID2NAME,
        ...
    )
}

#' Wraps a single enrichment call to clusterProfiler.
#'
#' @param sc String vector of significant genes.
#' @param sn Numeric vector of effect size signs.
#' @param p2k Protein ID to KO mapping table (columns "node_head", "accession").
#' @param cn String value to go in the "cutoff" column.
#' @param d Integer (1 or -1) giving effect size directions to keep.
#' @param tg String value to go in the "taxon" column.
#' @param dirxn Count only genes with this effect sign as significant.
#' @param kegg_pw_data Optionally provide already-downloaded KEGG KO pathway data.
#' @param kegg_mod_data Optionally provide already-downloaded KEGG KO module data.
#' @return A tbl giving Fisher's test p-values, q-values, effect sizes, and overlaps.
#' @export
kegg.enrich.single <- function(sc, sn, p2k, cn="test_cutoff",
                               d=1, tg="test_taxon", background=NULL,
                               kegg_pw_data=NULL,
                               kegg_mod_data=NULL) {
    genes <- intersect(sc, nw(sn == d))
    if (is.null(background)) { background <- unique(p2k[["accession"]]) }
    KOs <- unique(dplyr::filter(p2k, node_head %in% genes)[["accession"]])
    
    if (is.null(kegg_pw_data)) {
        pwy_enr <- clusterProfiler::enrichKEGG(gene = KOs,
                                               organism = "ko",
                                               universe = background,
                                               qvalueCutoff = 0.25)
    } else {
        pwy_enr <- enrich_downloaded_KEGG(KOs,
                                          background,
                                          qCut = 0.25,
                                          kegg_pw_data)
    }
    if (is.null(kegg_mod_data)) {
        mod_enr <- clusterProfiler::enrichMKEGG(gene = KOs,
                                                organism = "ko",
                                                universe = background,
                                                qvalueCutoff = 0.25)
    } else {
        mod_enr <- enrich_downloaded_KEGG(KOs,
                                          background,
                                          qCut = 0.25,
                                          kegg_mod_data)
    }
    if (is.null(pwy_enr) || is.null(mod_enr)) {
        return(NULL)
    } else if(!is.null(pwy_enr) && !is.null(mod_enr)) {
        total_enr_tbl <- dplyr::bind_rows(pwy_enr@result, mod_enr@result) %>%
            tibble::as_tibble() %>%
            dplyr::mutate(cutoff = cn, taxon = tg)
        return(total_enr_tbl)
    } else {
        pz.error(paste0(
            "Internal dataframe is malformed for kegg.enrich.single. ",
            "Please file a bug report."))
    }
}

