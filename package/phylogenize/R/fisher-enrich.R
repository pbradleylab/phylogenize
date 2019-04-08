# Functions for performing enrichment analysis

#' Wrapper around \code{qvalue} that extracts only q-values. If there is an
#' error in estimating q-values, will automatically fall back to a
#' Benjamini-Hochberg-style correction (by setting lambda to zero), finally
#' returning a vector of NAs if this still does not work.
#'
#' @param x A vector of p-values.
#' @return A vector of q-values.
qvals <- function(x, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    tryCatch(qvalue::qvalue(x, fdr=T, lambda=seq(0.001, 0.95, 0.005))$qvalues,
             error=function(e) {
                 pz.warning("Falling back to BH...", ...)
                 tryCatch(qvalue::qvalue(x, fdr=T, lambda=0)$qvalues,
                          error=function(e) {
                              pz.warning(e, ...)
                              q <- rep(NA, length(x))
                              names(q) <- names(x)
                              q
                          })
             })
}

#' Wrapper around \code{p.adjust('BY')}.
#'
#' @param x A vector of p-values.
#' @return A vector of BY FDR values.
fdr.by <- function(x) { p.adjust(x, 'BY') }

#' Wrapper around \code{p.adjust('BH')}.
#'
#' @param x A vector of p-values.
#' @return A vector of BH FDR values.
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
#' @param sigs List giving, per phylum (outer) and per significance cutoff
#'     (inner), significant hits to test for enrichment.
#' @param signs List giving, per phylum, signs of all gene effect sizes.
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
    tbl.mappings <- enframe(map(mappings, ~group_by_at(., colnames(.)[2]) %>%
                                 nest() %>% rename(term=1) %>%
                                     mutate(data=map(data, unlist)))) %>%
        rename(termset=name, terms=value)
    phylum <- names(sigs)
    cutoff <- names(sigs[[1]])
    map.bg <- Reduce(union, map(tbl.mappings$terms,
                            ~Reduce(union, .$data)))
    full.table <- crossing(phylum, cutoff, tbl.mappings) %>% unnest()
    full.table <- mutate(full.table,
                         enr=pmap(full.table, function(phylum,
                                                       cutoff,
                                                       termset,
                                                       term,
                                                       data) {
                             s <- intersect(sigs[[phylum]][[cutoff]],
                                            nw(signs[[phylum]] == d))
                             g <- names(which(!is.na(signs[[phylum]])))
                             do.fisher(data,
                                       s,
                                       intersect(map.bg, g))
                         }))
    full.table <- mutate(full.table,
                         enr.pval=map_dbl(enr, ~.$p.value),
                         enr.estimate=map_dbl(enr, ~.$estimate),
                         enr.overlap=map(enr, ~.$overlap))
    full.table <- full.table %>%
        group_by(phylum, cutoff, termset) %>%
        mutate(enr.qval=qvals(pv1(enr.pval))) %>%
        ungroup
}


#' Fix p-values that are above 1.
#'
#' Sometimes, p-values from the Fisher test can apparently be slightly larger
#' than one for some reason; this works around that problem.
#'
#' @param x A vector of p-values.
#' @return The same vector with all p-values above 1 changed to exactly 1.
pv1 <- function(x) { x[which(x > 1)] <- 1; return(x) }

