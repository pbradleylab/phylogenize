# significance testing

make.pos.sig <- function(sigs, signs, cut = "strong") {
  lapply.across.names(names(sigs), function(x) {
    intersect(sigs[[x]][[cut]], nw(signs[[x]] > 0))
  })
}

qvals <- function(x) qvalue(na.omit(x),
                            pfdr = T,
                            pi0.method = "bootstrap")$qvalues

fdr.bh <- function(x) p.adjust(x, "BH")
fdr.by <- function(x) p.adjust(x, "BY")

make.sigs <- function(results,
                      cuts = c(strong = 0.05, med = 0.1, weak = 0.25),
                      method = qvals,
                      exclude = NULL) {
  lapply.across.names(names(results), function(x) { 
    lapply(cuts, function(cut) {
      if (!is.null(exclude)) {
      valid <- setdiff(colnames(results[[x]]), 
        exclude[[x]])
      } else { valid <- colnames(results[[x]]) }
      tested <- na.omit(results[[x]][2, valid])
      tryCatch(nw(method(tested) <= cut), error = function(e) character(0))
      })
  })
}

make.signs <- function(results) {
  lapply(results, function(r) {
    sign(r[1, ]) %>% na.omit
  })
}

calc.alpha.power <- function(pvs, null, alt, alpha = 0.05, filter = NULL) {
  reject = which(pvs <= alpha)
  if (!is.null(filter)) {
    alt <- intersect(alt, filter)
    null <- intersect(null, filter)
    reject <- intersect(reject, filter)
  }   # allows you to look at e.g. just predicted - fx
  pwr = length(intersect(reject, alt)) / length(alt)
  a = length(intersect(reject, null)) / length(null)
  return(c(r = length(reject) / length(pvs), p = pwr, a = a))
}


results.report <- function(results, sigs, signs) {
  pos.sigs <- lapply(c(strong="strong",
                       med="med",
                       weak="weak"),
                     function(level) {
                       make.pos.sig(sigs, signs, level)
                     })
  sapply(pos.sigs, function(x) {
    sapply(x, function(y) length(y))
  })
}

get.top.N <- function(p,
                      sigs,
                      signs,
                      results,
                      level = "strong",
                      exclude = NULL,
                      N = 25,
                      total.n.cutoff = 0,
                      genomes.per.protein = NULL) {
  sig.up <- intersect(sigs[[p]][[level]], nw(signs[[p]] > 0))
  if (!is.null(genomes.per.protein)) {
    gn <- names(which(genomes.per.protein[[p]] >= total.n.cutoff))
    sig.up <- intersect(sig.up, gn)
  }
  sig.pv <- results[[p]][2, sig.up]
  setdiff(names(sort(sig.pv, dec = F)), exclude)[1:N]
}
