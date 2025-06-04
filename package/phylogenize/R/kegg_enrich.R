#' Given lists of significant genes (at different thresholds), effect sizes, and
#' gene set mappings, assemble a tbl of results using clusterProfiler.
#'
#' @param sigs List giving, per taxonomic group (outer) and per significance cutoff
#'     (inner), significant hits to test for enrichment.
#' @param signs List giving, per taxonomic group, signs of all gene effect sizes.
#' @param pid_to_ko Protein ID to KO mapping table (columns "node_head", "accession").
#' @param dirxn Count only genes with this effect sign as significant.
#' @return A tbl giving Fisher's test p-values, q-values, effect sizes, and overlaps for all taxa and cutoffs.
#' @export
multi.kegg.enrich <- function(sigs, signs, pid_to_ko, dirxn=1) {
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
  enrichment_tbls <- purrr::pmap(list(s=sigs, sn=signs, tg=tax_groups),
                                 function(s, sn, tg) {
                                   all_cutoffs <- purrr::map2(s, names(s), function(sc, cn) {
                                     kegg.enrich.single(sc, sn, pid_to_ko, cn, d, tg)
                                   })
                                   dplyr::bind_rows(all_cutoffs)
                                 })
  out <- tryCatch({
    dplyr::bind_rows(enrichment_tbls) %>%  # Put together and add direction of enrichment
      separate_wider_delim(cols = GeneRatio, delim = "/", names = c("GeneNum", "GeneDenom")) %>%
      separate_wider_delim(cols = BgRatio, delim = "/", names = c("BgRatioNum", "BgRatioDenom")) %>%
      mutate(enr.estimate =
             (as.numeric(GeneNum) / as.numeric(GeneDenom)) /
             (as.numeric(BgRatioNum) / as.numeric(BgRatioDenom)))
  }, error = function(e) {
    message("An error occurred: ", e$message)
    return(NULL)
  })
  return(out)
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
#' @return A tbl giving Fisher's test p-values, q-values, effect sizes, and overlaps.
#' @export
kegg.enrich.single <- function(sc, sn, p2k, cn="test_cutoff", d=1, tg="test_taxon", background=NULL) {
  genes <- intersect(sc, nw(sn == d))
  if (is.null(background)) { background <- unique(p2k[["accession"]]) }
  KOs <- unique(dplyr::filter(p2k, node_head %in% genes)[["accession"]])

  pwy_enr <- clusterProfiler::enrichKEGG(gene = KOs,
                                         organism = "ko",
                                         universe = background,
                                         qvalueCutoff = 0.25)
  mod_enr <- clusterProfiler::enrichMKEGG(gene = KOs,
                                          organism = "ko",
                                          universe = background,
                                          qvalueCutoff = 0.25)

  if (is.null(pwy_enr) || is.null(mod_enr)) {
	  return(NULL)
  } else if(!is.null(pwy_enr) && !is.null(mod_enr)) {
	  total_enr_tbl <- bind_rows(pwy_enr@result, mod_enr@result) %>%
  		  as_tibble() %>%
  		  mutate(cutoff = cn, taxon = tg)
	  return(total_enr_tbl)
  } else {
  	  pz.error("Internal dataframe is malformed for kegg.enrich.single. Please file a bug report.")
  }
}

