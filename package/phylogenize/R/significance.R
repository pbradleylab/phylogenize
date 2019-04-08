#' Get vectors of significant genes with positive effect sizes.
#'
#' @param sigs Output of \code{make.sigs}
#' @param signs Output of \code{make.signs}
#' @param cut String giving named significance level to use.
#' @return List (per phylum) of string vectors of positive significant hits.
#' @export
make.pos.sig <- function(sigs, signs, cut = "strong") {
    lapply.across.names(names(sigs), function(x) {
        intersect(sigs[[x]][[cut]], nw(signs[[x]] > 0))
    })
}

#' Get vectors of significant genes from result tables.
#'
#' @param results List of result matrices with two rows (effect size and
#'     p-value) and one column per gene tested.
#' @param cuts Named numeric vector giving different significance cutoffs.
#' @param method Function that will be used to adjust raw p-values in
#'     \code{results}.
#' @param exclude String vector of genes to exclude (optional).
#' @return List (per phylum) of string vectors of significant hits.
#' @export
make.sigs <- function(results,
                      cuts = c(strong = 0.05,
                               med = 0.1,
                               weak = 0.25),
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

#' Get effect sizes of genes from result tables.
#'
#' @param results List of result matrices with two rows (effect size and
#'     p-value) and one column per gene tested.
#' @return List (per phylum) of numeric vectors of signs of hits.
#' @export
make.signs <- function(results) {
  lapply(results, function(r) {
    sign(r[1, ]) %>% na.omit
  })
}

#' Calculate the FPR and (1 - FNR) for results of a set of tests.
#'
#' @param pvs A named vector of p-values, one per test.
#' @param null A vector of strings giving the tests in \code{pvs} for which the
#'     null was true.
#' @param alt A vector of strings giving the tests in \code{pvs} for which the
#'     alternative hypothesis was true.
#' @param filter Optional vector of strings giving the tests to which the
#'     analysis should be restricted.
#' @return A numeric vector:
#'   \item{r}{Proportion of p-values where the null was rejected.}
#'   \item{p}{Power (1 - FNR)}
#'   \item{a}{Alpha (FPR)}
#' @export
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

#' Create a summary table giving how many tests were significant.
#'
#' @param results List of result matrices, one per phylum.
#' @param sigs The output of \code{make.sigs}.
#' @param signs The output of \code{make.signs}.
#' @return A table with the number of positive significant results per phylum at
#'     each significance level in \code{sigs}.
#' @export
results.report <- function(results, sigs, signs) {
    if (length(sigs) < 1) { error("sigs must have at least one element") }
    pos.sigs <- lapply(names(sigs[[1]]),
                       function(level) {
                           make.pos.sig(sigs, signs, level)
                       })
    sapply(pos.sigs, function(x) {
        sapply(x, function(y) length(y))
    })
}

#' Return the significant hits with the N smallest p-values.
#'
#' @param p A phylum
#' @param sigs The output of \code{make.sigs}.
#' @param signs The output of \code{make.signs}.
#' @param results List of result matrices, one per phylum.
#' @param N Integer; how many hits to return.
#' @param level Significance level (must be in \code{sigs[[1]]}).
#' @param exclude Optional: exclude these genes from any list.
#' @param genomes.per.protein Optional: list (one per phylum) of named numeric
#'     vectors giving the number of genomes that each protein was found in.
#' @param total.n.cutoff Optional: if \code{genomes.per.protein} provided, only
#'     return hits found in at least this many genomes.
#' @return A named vector of N significant hits in descending order of
#'     significance.
#' @export
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
