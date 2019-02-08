#' Outermost wrapper to perform enrichments at different levels of significance.
#'
#' @param sigs List of lists of string vectors, giving which genes were
#'     significant at different significance levels in different phyla.
#'     Outermost list is per-phylum, inner list is per-significance level.
#' @param signs Named list of numeric vectors, giving the signs of the effect
#'     sizes for all genes.
#' @param results Named list of result matrices, one per phylum (list names).
#'     First row is effect size, second row is p-value, and genes are in
#'     columns.
#' @param mappings Named list of mappings, one per level of the mapping
#'     hierarchy. For each, first column gives gene names, second column gives
#'     annotations.
#' @param dirxn An integer. Consider only positive (1) or negative (-1) effect
#'     sizes.
#' @param exclude A character vector optionally giving genes to exclude from the
#'     analysis.
#' @return A list of \code{dirxn.enrich} results, one per level of significance
#'     in \code{sigs}.
#' @export
perform.enrichments <- function(sigs,
                                signs,
                                results,
                                mappings,
                                dirxn = 1,
                                exclude = NULL,
                                ...) {
    lapply.across.names(names(sigs[[1]]), function(x) {
        dirxn.enrich(sigs,
                     signs,
                     results,
                     siglevel=x,
                     exclude=exclude,
                     fdr.method=fdr.bh,
                     dirxn=dirxn,
                     mappings=mappings,
                     ...)
    })
}

#' Wrapper to perform enrichments and calculate q-values.
#'
#' @param sig Vector of names of significant gene hits (strings).
#' @param tested Vector of names of all genes tested (strings).
#' @param bg Optional vector of names of genes to intersect with both \code{sig}
#'     and \code{tested}.
#' @param mapping Mapping of gene names (column 1) to annotations (column 2).
#' @param descs Optional one-column data frame with row names equal to
#'     annotation labels and values equal to descriptions of each annotation
#'     label.
#' @param ms Integer; minimum size of enrichment set to consider.
#' @param cutoff Double; anything with an FDR less than/equal to this value will
#'     be returned.
#' @param fdr.method A function to convert p-values into FDRs.
#' @return A list, with components:
#'   \item{full}{The full results of \code{enrich}.}
#'   \item{qv}{q-values or adjusted p-values for each enrichment set.}
#'   \item{dep.enr}{The names of all enrichment sets that were significantly
#'     enriched or depleted.}
#'   \item{enr}{The names of all enrichment sets that were significantly
#'     enriched.}
#'   \item{table}{A table giving the names of all significantly enriched sets in
#'     the first column and their q/p-values in the second column.}
#' @export
enr.wrapper <- function(sig,
                        tested,
                        bg = NULL,
                        mapping,
                        descs = NULL,
                        ms = 3,
                        cutoff = 0.25,
                        fdr.method = qvals,
                        ...) {
  if (is.null(bg)) { bg = tested }
  bg.present <- intersect(tested, bg)
  enriched.and.depleted.list <- enrich(sig,
                                       mapping,
                                       bg.present,
                                       minsize = ms,
                                       alt = "two.sided",
                                       ...)
  enr.dep.qvals <- enriched.and.depleted.list$enrichments %>%
      pv1 %>%
      fdr.method
  enriched.and.depleted <- nw(enr.dep.qvals <= cutoff)
  enriched <- intersect(enriched.and.depleted,
                        nw(enriched.and.depleted.list$odds.ratios > 1))
  if (!is.null(descs)) {
    table <- cbind(enriched, enr.dep.qvals[enriched], descs[enriched,])
  } else {
    table <- cbind(enriched, enr.dep.qvals[enriched])
  }
  t.order <- order(enr.dep.qvals[enriched])
  return(list(full = enriched.and.depleted.list,
              qv = enr.dep.qvals,
              dep.enr = enriched.and.depleted,
              enr = enriched,
              table = table[t.order,]))
}

#' Wrapper around \code{enr.wrapper} to get pathway enrichments of genes having
#' either positive or negative effect sizes.
#'
#' @param sigs List of lists of string vectors, giving which genes were
#'     significant at different significance levels in different phyla.
#'     Outermost list is per-phylum, inner list is per-significance level.
#' @param signs Named list of numeric vectors, giving the signs of the effect
#'     sizes for all genes.
#' @param results Named list of result matrices, one per phylum (list names).
#'     First row is effect size, second row is p-value, and genes are in
#'     columns.
#' @param siglevel String; gives which of the significance levels in \code{sigs}
#'     to take as significant.
#' @param mappings Named list of mappings, one per level of the mapping
#'     hierarchy. For each, first column gives gene names, second column gives
#'     annotations.
#' @param dirxn An integer. Consider only positive (1) or negative (-1) effect
#'     sizes.
#' @param exclude A character vector optionally giving genes to exclude from the
#'     analysis.
#' @return A list of lists (per mapping, per phylum) of \code{enr.wrapper}
#'     results.
#' @export
dirxn.enrich <- function(sigs,
                         signs,
                         results,
                         siglevel = "strong",
                         exclude = NULL,
                         dirxn = 1,
                         mappings,
                         ...) {
    lapply(names(mappings), function(l) {
                                        # dirxn = 1 is POSITIVE
        message(l)
        enrichments <- lapply.across.names(names(sigs), function(x) {
            enr.wrapper(intersect(sigs[[x]][[siglevel]],
                                  nw((dirxn * signs[[x]]) > 0)),
                        tested = setdiff(names(na.omit(results[[x]][1,])),
                                         exclude[[x]]),
                        mapping = mappings[[l]], ...)
        })
        return(enrichments)
    })
}

# Functions for performing enrichment analysis

#' Wrapper around \code{qvalue} that extracts only q-values.
#'
#' @param x A vector of p-values.
#' @return A vector of q-values.
qvals <- function(x) {
  qvalue(x,
         fdr = T,
         lambda = seq(0.001, 0.95, 0.005)
  )$qvalues

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
#'     from \code{list1}, \code{list2}, and \code{background}.
do.fisher <- function(list1, list2, background, alt) {
    background <- unique(background)
    l1 <- unique(intersect(list1, background))
    l2 <- unique(intersect(list2, background))
    both <- intersect(l1, l2) %>% length
    only1 <- setdiff(l1, l2) %>% length
    only2 <- setdiff(l2, l1) %>% length
    neither <- setdiff(background, union(l1, l2)) %>% length
    cont.table <- matrix(nr=2, c(both, only1, only2, neither))
    return(fisher.test(cont.table, alternative = alt))
}

#' Test enrichment of a gene set with various pathways, using Fisher's exact
#' test.
#'
#' This function wraps \code{do.fisher}, and also yields the genes that
#' overlapped, which can be useful for downstream analyses.
#'
#' @param hits A string vector of significant hits.
#' @param mapping Mapping of gene names (column 1) to annotations (column 2).
#' @param background A vector of strings giving all potential hits (e.g. genes
#'     tested in an assay).
#' @param minsize Minimum pathway size to test.
#' @param alt Alternative for Fisher's test (see \code{?fisher.test})
#' @param verbose Boolean; display a progress bar?
#' @return A list:
#'   \item{enrichments}{Unadjusted p-values of enrichments, per pathway.}
#'   \item{odds.ratios}{Odds ratios of enrichments, per pathway.}
#'   \item{overlaps}{List of character vectors of hits that overlapped with each
#'     pathway.}
#'   \item{sets}{List of lists giving each entry in the contingency table.}
#' @export
enrich <- function(hits,
                   mapping,
                   background,
                   minsize = 1,
                   alt = "greater",
                   verbose = T) {
  m.hits <- intersect(hits, background)
  possible <- intersect(mapping[,1], background)
  pways.total <- unique(mapping[which(mapping[,1] %in% background), 2])
  pw.lengths <- sapply(pways.total, function(y) {
      (mapping[,2] == y) %>%
          which %>%
          (function(z) length(intersect(mapping[z, 1], background)))})
  pways <- nw(pw.lengths >= minsize)
  if (length(pways) == 0) {
      return(list(enrichments = NULL,
                  odds.ratios = NULL,
                  overlaps = NULL,
                  sets = NULL)) }
  enrichments <- vector("numeric", length(pways))
  odds.ratios <- vector("numeric", length(pways))
  overlaps <- list()
  hits.total <- length(m.hits)
  progress.bit <- ceiling(length(pways)/20)
  sets <- list()
  if (verbose) cat("|----|----|----|----|\n")
  for (np in 1:length(pways)) {
    if (verbose) {
      if ( (np %% progress.bit) == 0 ) { cat (".") }
    }
    p <- pways[np]
    names(enrichments)[np] <- p
    in.pway <- intersect(background, mapping[which(mapping[,2] == p), 1])
    ft <- do.fisher(in.pway, m.hits, background, alt)
    pway.nohits <- length(setdiff(in.pway, m.hits))
    pway.hits <- length(intersect(in.pway, m.hits))
    nopway.hits <- length(setdiff(m.hits, in.pway))
    nopway.nohits <- length(setdiff(background, union(m.hits, in.pway)))
    cont.table <- (matrix(nr=2,
                          c(pway.hits,
                            nopway.hits,
                            pway.nohits,
                            nopway.nohits)))
    enrichments[np] <- ft$p.value
    odds.ratios[np] <- ft$estimate
    overlaps[[p]] <- intersect(in.pway, m.hits)
    sets[[p]] <- list(PH = pway.hits,
                      PNH = pway.nohits,
                      NPH = nopway.hits,
                      NPNH = nopway.nohits)
  }
  names(odds.ratios) <- pways
  if(verbose) cat(".\n")
  return(list(enrichments=enrichments,
              odds.ratios = odds.ratios,
              overlaps = overlaps,
              sets = sets))
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

#' Produce a table of significant overlaps.
#'
#' Given results of an enrichment and results of an original test, return a
#' table showing which genes contributed to each enrichment and what their
#' original (i.e., not enrichment) q-values were.
#'
#' @param enr Results of \code{enrich}.
#' @param results Result matrix for the original genewise test. First row is
#'     effect size, second row is p-value, and genes are in columns.
#' @return A list of data frames, one per pathway, each containing:
#'   \item{genes}{Gene names that overlapped with a particular pathway.}
#'   \item{descs}{Description of genes.}
#'   \item{estimate}{Per-gene effect sizes.}
#'   \item{qval}{Per-gene q-values.}
#' @export
signif.overlaps <- function(enr, result) {
  if (!is.null(enr$table %>% dim)) {
    overlap.list <- enr$full$overlaps[enr$table %>% rownames]
  } else if ("enriched" %in% names(enr$table)) {
    overlap.list <- enr$full$overlaps[enr$table["enriched"]]
  } else {
    overlap.list <- list()
  }
  x <- lapply(names(overlap.list), function(n) {
    g <- overlap.list[[n]]
    if (length(g) > 0) {
      cbind(genes = overlap.list[[n]],
        descs = gene.annot(overlap.list[[n]]),
        estimate = result[1, g],
        qval = qvals(result[2, ])[g])
    } else {
      cbind(NA, NA, NA, NA)
    }
  })
  names(x) <- names(overlap.list)
  x
}
