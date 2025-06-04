
# Functions for performing enrichment analysis. Note: deprecated in newest version, use kegg_enrich.R

#' Wrapper around \code{p.adjust('BY')}.
#'
#' @param x A vector of p-values.
#' @return A vector of BY FDR values.
#' @keywords internal
fdr.by <- function(x) { p.adjust(x, 'BY') }

#' Wrapper around \code{p.adjust('BH')}.
#'
#' @param x A vector of p-values.
#' @return A vector of BH FDR values.
#' @keywords internal
fdr.bh <- function(x) { p.adjust(x, 'BH') }


#' Wrapper around \code{fisher.test} that starts with string vectors instead of
#' a contingency table.
#'
#' @param list1 A vector of strings (e.g., names of significant gene hits).
#' @param list2 A vector of strings (e.g., names of genes in a pathway).
#' @param background A vector of strings that \code{list1} and \code{list2} were
#'     drawn from (e.g., names of all genes tested in an assay). Any strings in
#'     \code{list1} or \code{list2} not in this vector will be discarded.
#' @param alt The alternative hypothesis (see \code{?fisher.test}).
#' @return The results of \code{fisher.test} on a contingency table generated
#'     from \code{list1}, \code{list2}, and \code{background}, augmented with
#'     the additional field \code{overlap} which gives the intersection of
#'     \code{list1} and \code{list2}.
#' @export do.fisher
do.fisher <- function(list1, list2, background, alt="two.sided") {
    background <- unique(background)
    l1 <- unique(intersect(list1, background))
    l2 <- unique(intersect(list2, background))
    both <- intersect(l1, l2)
    only1 <- setdiff(l1, l2)
    only2 <- setdiff(l2, l1)
    neither <- setdiff(background, union(l1, l2))
    cont.table <- matrix(nr=2,
                         vapply(list(both, only1, only2, neither),
                                length,
                                1L))
    f <- fisher.test(cont.table, alternative = alt)
    f$overlap <- both
    f
}



#' Given lists of significant genes (at different thresholds), effect sizes, and
#' gene set mappings, assemble a tbl of results.
#'
#' @param sigs List giving, per taxonomic group (outer) and per significance cutoff
#'     (inner), significant hits to test for enrichment.
#' @param signs List giving, per taxonomic group, signs of all gene effect sizes.
#' @param mappings List of data.frames giving gene-to-gene-set mappings.
#' @param dirxn Count only genes with this effect sign as significant.
#' @return A tbl giving Fisher's test p-values, q-values, effect sizes, and overlaps.
#' @export
multi.enrich <- function(sigs, signs, mappings, dirxn=1) {
  d <- sign(dirxn)
  if ((length(sigs) < 1) || (length(signs) < 1) || (length(mappings) < 1)) {
    pz.warning("Didn't find significant hits, signs, or mappings")
    return(c())
  }
  if (length(sigs[[1]]) < 1) {
    pz.warning("No significance cutoffs found")
    return(c())
  }
  tbl.mappings <- tibble::enframe(purrr::map(mappings,
      ~dplyr::group_by_at(., colnames(.)[2]) %>%
        dplyr::nest() %>% dplyr::rename(term=1) %>%
        dplyr::mutate(data=purrr::map(data, unlist)))) %>%
        dplyr::rename(termset=name, terms=value)
     taxon <- names(sigs)
      cutoff <- names(sigs[[1]])
      map.bg <- Reduce(union, purrr::map(tbl.mappings$terms,
          ~Reduce(union, .$data)))
      full.table <- tidyr::crossing(taxon, cutoff, tbl.mappings) %>% unnest()
      full.table <- dplyr::mutate(full.table,
        enr=purrr::pmap(full.table, function(taxon,
            cutoff,
            termset,
            term,
            data) {
          s <- intersect(sigs[[taxon]][[cutoff]],
            nw(signs[[taxon]] == d))
          g <- names(which(!is.na(signs[[taxon]])))
          do.fisher(data,
            s,
            intersect(map.bg, g))
          }))
      full.table <- dplyr::mutate(full.table,
        enr.pval=purrr::map_dbl(enr, ~.$p.value),
        enr.estimate=purrr::map_dbl(enr, ~.$estimate),
        enr.overlap=purrr::map(enr, ~.$overlap))
      full.table <- full.table %>%
        dplyr::group_by(taxon, cutoff, termset) %>%
        dplyr::mutate(enr.qval=qvals(pv1(enr.pval))) %>%
        dplyr::ungroup()
}

#' Get q-values for results in tbl format.
#'
#' @param results A tbl of results with columns for "taxon" and "p.value".
#' @param method A function: method for obtaining q-values from p-values.
#' @param ... Additional parameters to pass to \code{method}.
#' @return A tbl with an additional "q.value" column.
#' @export tbl.result.qvs
tbl.result.qvs <- function(results, method=qvals, ...) {
    results <- results %>%
        dplyr::group_by(taxon) %>%
        dplyr::nest() %>%
        dplyr::mutate(q.value=purrr::map(data, ~ {
            method(.x$p.value)
        })) %>%
        dplyr::unnest
    results
}

#' An alternative way to get the enrichment table, using a tbl of results
#' instead of separate inputs for significant results, effect signs, etc.
#'
#' @param results A tidyverse tbl of results with at least the following
#'     columns: "taxon", "gene", "effect.size", and "q.value" (if only
#'     "p.value" is present, first apply function \code{tbl.result.qvs}).
#' @param mappings List of data.frames giving gene-to-gene-set mappings.
#' @param dirxn Count only genes with this effect sign as significant.
#' @return A tbl giving Fisher's test p-values, q-values, effect sizes, and overlaps.
#' @export
alt.multi.enrich <- function(results, mappings, dirxn=1,
                             qcuts=c(strong=0.05, med=0.1, weak=0.25),
                             future=FALSE) {
    d <- sign(dirxn)
    if (length(mappings) < 1) {
        pz.warning("Didn't find significant hits, signs, or mappings")
        return(c())
    }
    tbl.mappings <- tibble::enframe(
        purrr::map(mappings, ~group_by_at(., colnames(.)[2]) %>%
                       dplyr::nest() %>%
                       dplyr::rename(term=1) %>%
                       dplyr::mutate(data=purrr::map(data, unlist)))) %>%
        dplyr::rename(termset=name, terms=value)
    taxa <- unique(results$taxon)
    cutoff <- names(qcuts)
    map.bg <- Reduce(union, purrr::map(tbl.mappings$terms,
                                       ~Reduce(union, .$data)))
    full.table <- tidyr::crossing(
        taxon=taxa, cutoff, tbl.mappings) %>% dplyr::unnest()
    map_fxn <- purrr::pmap
    full.table <- dplyr::mutate(full.table,
                         enr=map_fxn(full.table,
                                     indiv.enr,
                                     results=results,
                                     qcuts=qcuts,
                                     dirxn=dirxn,
                                     bg=map.bg))
    full.table <- dplyr::mutate(
        full.table,
        enr.pval=purrr::map_dbl(enr, ~.$p.value),
        enr.estimate=purrr::map_dbl(enr, ~.$estimate),
        enr.overlap=purrr::map(enr, ~.$overlap))
    full.table <- full.table %>%
        dplyr::group_by(taxon, cutoff, termset) %>%
        dplyr::mutate(enr.qval=qvals(pv1(enr.pval))) %>%
        dplyr::ungroup()
}

#' Internal function to perform an individual enrichment.
#' @keywords internal
indiv.enr <- function(taxon, cutoff, termset, term, data, results,
                      qcuts=c(strong=0.05), dirxn=1, bg=NULL) {
    s <- dplyr::filter(results,
                       q.value <= qcuts[cutoff],
                       taxon == !!taxon,
                       sign(effect.size) == dirxn) %>%
        dplyr::select(gene) %>%
        unlist()
    g <- dplyr::filter(results,
                       taxon == !!taxon,
                       !is.na(effect.size)) %>%
        dplyr::select(gene) %>%
        unlist()
    if (is.null(bg)) {
        bg <- g
    } else {
        bg <- intersect(bg, g)
    }
    do.fisher(data,
              s,
              bg)
}

#' Fix p-values that are above 1.
#'
#' Sometimes, p-values from the Fisher test can apparently be slightly larger
#' than one for some reason; this works around that problem.
#'
#' @param x A vector of p-values.
#' @return The same vector with all p-values above 1 changed to exactly 1.
#' @export
pv1 <- function(x) { x[which(x > 1)] <- 1; return(x) }

