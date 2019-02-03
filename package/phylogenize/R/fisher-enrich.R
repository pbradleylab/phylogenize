
perform.enrichments <- function(sigs, signs, results, mapping, dirxn = 1, exclude = NULL, ...) {
  lapply.across.names(c("strong", "med", "weak"), function(x) {
    dirxn.enrich(sigs, signs, results, siglevel = x,
      exclude = exclude, fdr.method = fdr.bh, dirxn = dirxn,
      mapping = mapping, ...)
  })
}

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

dirxn.enrich <- function(sigs,
                         signs,
                         results,
                         siglevel = "strong",
                         exclude = NULL,
                         dirxn = 1,
                         mapping,
                         ...) {
    lapply(c(l1 = "level1", l2 = "level2", l3 = "level3"), function(l) {
        # dirxn = 1 is POSITIVE
        message(l)
        enrichments <- lapply.across.names(names(sigs), function(x) {
            enr.wrapper(intersect(sigs[[x]][[siglevel]],
                                  nw((dirxn * signs[[x]]) > 0)),
                        tested = setdiff(names(na.omit(results[[x]][1,])),
                                         exclude[[x]]),
                        mapping = mapping[[l]], ...)
        })
        return(enrichments)
    })
}


# Functions for performing enrichment analysis

qvals <- function(x) {
  qvalue(x,
         fdr = T,
         lambda = seq(0.001, 0.95, 0.005)
  )$qvalues
}
fdr.by <- function(x) { p.adjust(x, 'BY') }
fdr.bh <- function(x) { p.adjust(x, 'BH') }

do.fisher <- function(list1, list2, background, alt) {
  l1 <- intersect(list1, background)
  l2 <- intersect(list2, background)
  both <- intersect(l1, l2) %>% length
  only1 <- setdiff(l1, l2) %>% length
  only2 <- setdiff(l2, l1) %>% length
  neither <- setdiff(background, union(l1, l2)) %>% length
  cont.table <- matrix(nr=2, c(both, only1, only2, neither))
  return(fisher.test(cont.table, alternative = alt))
}

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
                  odds.ratios= NULL,
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



# Sometimes p-values from the Fisher test can be 1+epsilon for some reason;
# this fixes that problem

pv1 <- function(x) { x[which(x > 1)] <- 1; return(x) }

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
