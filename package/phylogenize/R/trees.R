#' Fudge trees with unresolved polytomies.
#'
#' \code{fix.tree} converts polytomies to dichotomies with very small branch
#' lengths.
#'
#' @param len Zero-length branches will be replaced with this fraction of the
#'     maximum node height.
#' @keywords internal
fix.tree <- function(phy, len=1e-6) {
  phy <- ape::multi2di(phy)
  ## from Liam Revell's blog, June 23 2015:
  ## http://blog.phytools.org/2015/06/update-to-rerootingmethod-for-ancestral.html
  phy$edge.length[phy$edge.length==0] <- max(phytools::nodeHeights(phy)) * len
  phy
}

#' Keep some tips on a tree.
#'
#' \code{keep.tips} keeps only the set of specified tips in a tree.
#'
#' @param tree A \code{phylo} object.
#' @param keep A character vector of tip labels. Any tip not in this vector will
#'     be dropped.
#' @export keep.tips
keep.tips <- function(tree, keep) {
  ape::drop.tip(tree, setdiff(tree$tip.label, keep))
}

#' Get tip-to-root distance.
#'
#' \code{tipToRoot} returns all tip-to-root distances.
#'
#' @param phy A \code{phylo} object.
#' @return A vector of root-to-tip distances.
#' @keywords internal
tipToRoot <- function(phy) phy %>% ape::vcv.phylo %>% diag

#' Get tip-to-tip distances.
#'
#' \code{tree.to.dist} returns all tip-to-tip distances in distance matrix
#' format.
#'
#' @param tree A \code{phylo} object.
#' @return A \code{dist} object representing tip-to-tip distances.
#' @keywords internal
tree.to.dist <- function(tree) {
  ((1 - (ape::vcv(tree, TRUE))) / 2) %>% as.dist
}

#' Plot continuous trait on a tree.
#'
#' \code{gg.cont.tree} paints a continuous trait along a tree.
#'
#' @param phy A \code{phylo} object.
#' @param ctrait A named numeric vector assigning trait values to tree tips.
#' @param cAnc Calculated ancestry for continuous trait; if this is NULL, it is calculated.
#' @param model Model for calculating ancestry (see phytools::fastAnc).
#' @param cLimits Scale bar limits for plotting continuous trait.
#' @param n Character vector giving which nodes to display; if NULL, defaults to intersection of \code{phy$tip.label} and \code{names(ctrait)}.
#' @param reduced.phy Dichotomous tree with only the nodes in \code{n} represented; if NULL, this is calculated.
#' @param colors Named character vector with at least "low.col" and "high.col" and optionally "mid.col" defined, giving colors to use for plotting.
#' @param plot Whether to plot the tree object or just return it.
#' @param restrict Character vector giving which continuous trait values to plot; if NULL, all are used.
#' @param cName String giving the title of the plot.
#' @param reverse Mirror the resulting plotted tree.
#' @param ladderize Ladderize the plotted tree.
#' @param ... Additional arguments passed to ggplot2.
#' @return A ggtree plot of a continuous trait plotted along a tree.
#' @export
gg.cont.tree <- function(phy,
                         ctrait,
                         cAnc = NULL,
                         model = "ARD",
                         cLimits = logit(c(0.025, 0.1)),
                         n = NULL,
                         reduced.phy = NULL,
                         colors = c(
                           low.col = "slateblue",
                           mid.col = "black",
                           high.col = "orange2"),
                         plot = T,
                         restrict = NULL,
                         cName = "prevalence",
                         reverse = F,
                         ladderize = T,
                         ...) {
  if (!is.null(restrict)) {
    ctrait <- ctrait[intersect(names(ctrait), restrict)]
  }
  if (is.null(n)) { n <- intersect(phy$tip.label, names(ctrait)) }
  if (is.null(reduced.phy)) { reduced.phy <- fix.tree(keep.tips(phy, n)) }
  pz.message("getting continuous trait ancestry")
  if (is.null(reduced.phy)) {
      pz.error("all tips were dropped")
  }
  if (length(reduced.phy$tip.label) < 2) {
      pz.error("need at least two tips for a tree")
  }
  if (is.null(cAnc)) cAnc <- phytools::fastAnc(reduced.phy, ctrait[n])
  # concatenate tip values and node values
  cDisplay <- truncated(c(ctrait[n], cAnc), cLimits)
  # make color scales
  if ("mid.col" %in% names(colors)) {
      cColors <- ggplot2::scale_color_gradient2(low = colors["low.col"],
                                                high = colors["high.col"],
                                                mid = colors["mid.col"],
                                                midpoint = mean(cLimits),
                                                guide = "colorbar",
                                                name = cName)
  } else {
      cColors <- ggplot2::scale_color_gradient(low = colors["low.col"],
                                               high = colors["high.col"],
                                               guide = "colorbar",
                                               name = cName)
  }
  # plot trees
  if (reverse) {
      ctree <- ggtree::ggtree(reduced.phy,
                              ladderize = ladderize,
                              ggplot2::aes_string(color = cDisplay), ...) +
          cColors + ggplot2::scale_x_reverse() + ggplot2::theme(legend.position = "bottom")
  } else {
      ctree <- ggtree::ggtree(reduced.phy,
                              ladderize = ladderize,
                              ggplot2::aes_string(color = cDisplay), ...) +
          cColors + ggplot2::theme(legend.position = "bottom")
  }
  if (plot) print(ctree)
  return(list(tree = ctree,
              cAnc = cAnc,
              rphy = reduced.phy,
              n = n,
              cols = colors,
              lims = cLimits,
              disp = cDisplay))
}


#' Remove clades above a certain size with no observed tips.
#'
#' @param phy A phylo object.
#' @param observed A character vector giving all the observed tips. May contain
#'     tips not in the tree, which will be ignored.
#' @param pct Double; a cutoff. For example, \code{pct}=0.1 means that any
#'     unobserved clade containing at least 10% of the total number of tips in
#'     the tree will be dropped.
#' @return A new, trimmed phylo object.
#' @keywords internal
rem_unobs_clades <- function(phy, observed, pct=0.1) {
    ntips <- length(phy$tip.label)
    nodes <- (1:(phy$Nnode)) + ntips
    tiplist <- lapply(nodes, function(n) {
        clade <- ape::extract.clade(phy, n)
        clade.tips <- clade$tip.label
        nc <- length(clade.tips)
        if (((nc / ntips) >= pct) && (!any(clade.tips %in% observed))) {
            return(clade.tips)
        } else {
            return(character(0))
        }
    })
    to.drop <- Reduce(union, tiplist)
    ape::drop.tip(phy, to.drop)
}
