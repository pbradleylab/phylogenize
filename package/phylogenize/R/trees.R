#--- Functions to manipulate trees ---#

#' Fudge trees with unresolved polytomies.
#'
#' \code{fix.tree} converts polytomies to dichotomies with very small branch
#' lengths.
#'
#' @param len Zero-length branches will be replaced with this fraction of the
#'     maximum node height.
#' @keywords internal
fix.tree <- function(phy, len=1e-6) {
  if (!is.null(phy) && inherits(phy, "phylo")) {
      phy <- ape::multi2di(phy)
      if (any(phy$edge.length == 0)) {
          max_height <- max(castor::get_all_distances_to_root(phy))
          phy$edge.length[phy$edge.length == 0] <- max_height * len
      }
  }
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
