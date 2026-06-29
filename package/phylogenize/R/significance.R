#' Get vectors of significant genes with positive effect sizes.
#'
#' @param sigs Output of \code{make.sigs}
#' @param signs Output of \code{make.signs}
#' @param cut String giving named significance level to use.
#' @return List (per taxon) of string vectors of positive significant hits.
#' @export
make.pos.sig <- function(sigs, signs, cut = "strong") {
    lapply.across.names(names(sigs), function(x) {
        intersect(sigs[[x]][[cut]], nw(signs[[x]] > 0))
    })
}

#' Get vectors of significant genes with positive (or negative) effect sizes.
#'
#' @param results List of result matrices (4 x N).
#' @param method Method for performing multiple test adjustment.
#' @param qcut_sig Desired q-value cutoff for significance test.
#' @param qcut_eq Desired q-value cutoff for equivalence test (N.B.:
#'     significantly equivalent hits will be *excluded*).
#' @param min_fx Minimum effect size for equivalence test.
#' @param exclude Optional list of character vectors of genes to exclude.
#' @param dir Direction of test (1 is positive, -1 is negative; default: 1).
#' @return List of character vectors of hits that were significantly different
#'     from zero, had positive effect sizes, and not significantly equivalent to
#'     a minimum effect size.
#' @export
nonequiv.pos.sig <- function(results,
                             method=qvals,
                             qcut_sig=0.05,
                             qcut_eq=0.05,
                             min_fx=0.25,
                             exclude=NULL,
                             dir=1,
                             ...,
                             .opts=NULL) {
    opts <- pz.resolve.options(..., .opts=.opts)
    call_method <- function(vals) {
        tryCatch(
            method(vals, .opts=opts),
            error=function(e) {
                if (grepl("unused argument.*\\.opts", conditionMessage(e))) {
                    return(method(vals))
                }
                stop(e)
            }
        )
    }
    dir <- sign(dir)
    if (dir==0) { stop("dir must not be zero") }
    if (is.null(exclude)) exclude <- lapply(results, function(.) NULL)
    sig_fxn <- function(r, ex) {
        valid <- colnames(r)[which(!is.na(unlist(r[1, , drop = FALSE])))]
        if (!is.null(ex)) {
            valid <- setdiff(valid, ex)
        }
        #tryCatch(
        #    {
                tested <- unlist(r[2, valid, drop = FALSE])
                if (is.null(names(tested))) names(tested) <- valid
                tested_names <- names(tested)
                tested <- as.numeric(tested)
                names(tested) <- tested_names
                tested <- stats::na.omit(tested)
                qv <- call_method(tested)
                fx_sig <- nw(qv <= qcut_sig)
                neq_sig <- valid
                if (min_fx > 0) {
                    if (nrow(r) >= 4) {
                        neq <- apply(
                            unlist(r[, valid, drop = FALSE]),
                            2,
                            function(x) equiv_test(x[1], x[3], x[4], min_fx)
                        )
                        neq_qv <- call_method(neq)
                        neq_sig <- nw(neq_qv > qcut_eq)
                    } else {
                        pz.warning(
                            "Cannot perform non-equivalence test with permulation/permutration"
                        )
                    }
                }
                fx_sizes <- unlist(r[1, valid, drop=FALSE])
                if (is.null(names(fx_sizes))) names(fx_sizes) <- valid
                fx_names <- names(fx_sizes)
                fx_sizes <- as.numeric(fx_sizes)
                names(fx_sizes) <- fx_names
                which_pos <- nw((dir * fx_sizes) > 0)
                Reduce(intersect, list(fx_sig, neq_sig, which_pos))
          #  },
          #  error = function(e) character(0)
        #)
    }
    mapply(sig_fxn, results, exclude, SIMPLIFY=FALSE)
}

#' Equivalence test based on two one-sided tests.
#'
#' @param fx Estimated effect size.
#' @param se Estimated standard error of the effect size.
#' @param df Degrees of freedom.
#' @param min_fx Minimum effect size we care about.
#' @return A p-value; reject non-equivalence if below alpha.
#' @keywords internal
equiv_test <- function(fx, se, df, min_fx=0.25) {
    # test 1: H0 is fx >= min_fx, HA is fx < min_fx
    # test 2: H0 is fx <= -min_fx, HA is fx > -min_fx
    t_stat1 <- -(fx - min_fx) / se # more positive as fx << min_fx
    t_stat2 <- -(-min_fx - fx) / se # more positive as fx >> -min_fx
    pv1 <- pt(t_stat1, df, lower.tail=FALSE)
    pv2 <- pt(t_stat2, df, lower.tail=FALSE)
    max(pv1, pv2)
}

#' Get vectors of significant genes from result tables.
#'
#' @param results List of result matrices with two rows (effect size and
#'     p-value) and one column per gene tested.
#' @param cuts Named numeric vector giving different significance cutoffs.
#' @param method Function that will be used to adjust raw p-values in
#'     \code{results}.
#' @param exclude String vector of genes to exclude (optional).
#' @param min.fx Minimum effect size for calling something significant.
#' @param ... Extra parameters to be passed to `method`.
#' @return List (per taxon) of string vectors of significant hits.
#' @export
make.sigs <- function(
    results,
    cuts = c(strong = 0.05, med = 0.1, weak = 0.25),
    method = qvals,
    exclude = NULL,
    min.fx = 0,
    ...,
    .opts=NULL
) {
    opts <- pz.resolve.options(..., .opts=.opts)
    call_method <- function(vals) {
        tryCatch(
            method(vals, ..., .opts=opts),
            error=function(e) {
                if (grepl("unused argument.*\\.opts", conditionMessage(e))) {
                    return(method(vals, ...))
                }
                stop(e)
            }
        )
    }
    sig_fxn <- function(x, cut) {
        if (!is.null(exclude)) {
            valid <- setdiff(colnames(results[[x]]), exclude[[x]])
        } else {
            valid <- colnames(results[[x]])
        }
        # should work with repermulize output
        all_tested <- na.omit(
                purrr::map_dbl(results[[x]][2, valid, drop = FALSE], ~.x) |>
                    setNames(valid)
        )
        above.min <- nw(
            purrr::map_lgl(
                results[[x]][1, valid, drop = FALSE],
                ~ abs(.x) >= min.fx
            ) |>
                setNames(valid)
        )
        tryCatch(
            intersect(nw(call_method(all_tested) <= cut), above.min),
            error = function(e) character(0)
        )
    }
    lapply.across.names(names(results), function(x) {
        lapply(cuts, \(cut) sig_fxn(x, cut))
    })
}

#' Get effect sizes of genes from result tables.
#'
#' @param results List of result matrices with two rows (effect size and
#'     p-value) and one column per gene tested.
#' @return List (per taxon) of numeric vectors of signs of hits.
#' @export
make.signs <- function(results) {
    lapply(results, function(r) {
        fx <- as.numeric(r[1, , drop=FALSE])
        names(fx) <- colnames(r)
        stats::na.omit(sign(fx))
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
    reject <- which(pvs <= alpha)
    if (!is.null(names(pvs))) {
        reject <- names(reject)
    }
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
#' @param results List of result matrices, one per taxon.
#' @param sigs The output of \code{make.sigs}.
#' @param signs The output of \code{make.signs}.
#' @return A table with the number of positive significant results per taxon at
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
#' @param p A taxon
#' @param sigs The output of \code{make.sigs}.
#' @param signs The output of \code{make.signs}.
#' @param results List of result matrices, one per taxon.
#' @param N Integer; how many hits to return.
#' @param level Significance level (must be in \code{sigs[[1]]}).
#' @param exclude Optional: exclude these genes from any list.
#' @param genomes.per.protein Optional: list (one per taxon) of named numeric
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
    sig.pv <- unlist(results[[p]][2, sig.up, drop=FALSE])
    setdiff(names(sort(sig.pv, dec = F)), exclude)[1:N]
}

#' Wrapper around \code{qvalue} and \code{p.adjust} that extracts only q-values.
#' If there is an error in `qvalue()` with default parameters, this function
#' will automatically fall back to a Benjamini-Hochberg-style correction (by
#' setting lambda to zero), using p.adjust if that still doesn't work, and
#' finally returning a vector of NAs if all else fails.
#'
#' @param x A vector of p-values.
#' @param ... Extra parameters to override defaults, especially `fdr_method`
#'   (which can be "BH", "BY", or "qvalue).
#' @return A vector of q-values.
#' @keywords internal
qvals <- function(x, ..., .opts=NULL) {
    opts <- pz.resolve.options(..., .opts=.opts)
    method <- tolower(opts("fdr_method"))
    valid_methods <- c("bh", "by", "qvalue")
    if (!(method %in% valid_methods)) {
        pz.error(paste0(
            "Invalid fdr_method: ",
            opts("fdr_method"),
            ". Expected one of: BH, BY, qvalue"
        ), .opts=opts)
    }
    if (length(x) == 0) {
        return(x)
    }
    x_unlisted <- unlist(x)
    x_numeric <- suppressWarnings(as.numeric(x_unlisted))
    x_names <- names(x_unlisted)
    if (is.null(x_names) && !is.null(dim(x)) && ncol(x) == length(x_numeric)) {
        x_names <- colnames(x)
    }
    names(x_numeric) <- x_names
    nonmissing <- !is.na(x_numeric)
    invalid <- nonmissing & (!is.finite(x_numeric) |
                             x_numeric < 0 |
                             x_numeric > 1)
    if (any(invalid)) {
        pz.warning(
            paste0(
                "Marking ",
                sum(invalid),
                " invalid p-value(s) as NA before FDR correction"
            ),
            .opts=opts
        )
        x_numeric[invalid] <- NA_real_
        nonmissing <- !is.na(x_numeric)
    }
    q <- rep(NA_real_, length(x_numeric))
    names(q) <- names(x_numeric)
    if (!any(nonmissing)) {
        return(q)
    }
    x_valid <- x_numeric[nonmissing]
    if (method=="bh") {
        q[nonmissing] <- p.adjust(x_valid, 'BH')
        return(q)
    }
    if (method=="by") {
        q[nonmissing] <- p.adjust(x_valid, 'BY')
        return(q)
    }
    if (method=="qvalue") {
        q[nonmissing] <- tryCatch(qvalue::qvalue(x_valid,
                                     fdr=T,
                                     lambda=seq(0.001, 0.95, 0.005))$qvalues,
                      error=function(e) {
                          pz.warning("Trying lambda=0...", .opts=opts)
                          tryCatch(qvalue::qvalue(x_valid,
                                                  fdr=T,
                                                  lambda=0)$qvalues,
                                   error=function(e) {
                                       pz.warning(e, .opts=opts)
                                       pz.warning("Falling back to BH", .opts=opts)
                                       p.adjust(x_valid, 'BH')
                                   })
                      })
        return(q)
    }
}


#' Convert results into a long (vs. wide) format.
#'
#' @param results Output of result.wrapper.plm.
#' @return A single data frame with entries from \code{results}.
#' @export
make.results.matrix <- function(results) {
    result_names <- names(results)
    n_genes <- vapply(results, ncol, integer(1))
    n_rows <- sum(n_genes)

    genes <- character(n_rows)
    effect_sizes <- numeric(n_rows)
    p_values <- numeric(n_rows)
    std_errs <- numeric(n_rows)
    dfs <- numeric(n_rows)

    result_row <- function(result, row) {
        unname(as.numeric(unlist(result[row, , drop=FALSE],
                                 recursive=TRUE,
                                 use.names=FALSE)))
    }

    row_start <- 1L
    for (rn in result_names) {
        n_gene <- ncol(results[[rn]])
        if (n_gene == 0) {
            next
        }
        rows <- row_start:(row_start + n_gene - 1L)
        genes[rows] <- colnames(results[[rn]])
        effect_sizes[rows] <- result_row(results[[rn]], 1)
        p_values[rows] <- result_row(results[[rn]], 2)
        std_errs[rows] <- result_row(results[[rn]], 3)
        dfs[rows] <- result_row(results[[rn]], 4)
        row_start <- row_start + n_gene
    }

    names(effect_sizes) <- genes
    names(p_values) <- genes
    names(std_errs) <- genes
    names(dfs) <- genes

    tibble::tibble(
        taxon=rep(result_names, n_genes),
        gene=genes,
        effect.size=effect_sizes,
        p.value=p_values,
        std.err=std_errs,
        df=dfs
    )
}
