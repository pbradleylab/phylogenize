# Tree functions

fix.tree <- function(phy) {
  phy <- multi2di(phy)
  # from Liam Reveil's blog June 23 2015
  phy$edge.length[phy$edge.length==0] <- max(nodeHeights(phy)) * 1e-6
  phy
}

keep.tips <- function(tree, keeplist) {
  drop.tip(tree, setdiff(tree$tip.label, keeplist))
}

tipToRoot <- function(phy) phy %>% vcv.phylo %>% diag

tree.to.dist <- function(tree) {
  ((1 - (vcv(tree, TRUE))) / 2) %>% as.dist
}

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
  message("getting continuous trait ancestry")
  if (is.null(cAnc)) cAnc <- fastAnc(reduced.phy, ctrait[n])
  # concatenate tip values and node values
  cDisplay <- truncated(c(ctrait[n], cAnc), cLimits)
  # make color scales
  if ("mid.col" %in% names(colors)) {
    message("found mid.col")
    cColors <- scale_color_gradient2(low = colors["low.col"],
                                     high = colors["high.col"],
                                     mid = colors["mid.col"],
                                     midpoint = mean(cLimits),
                                     guide = "colorbar",
                                     name = cName)
  } else {
    cColors <- scale_color_gradient(low = colors["low.col"],
                                    high = colors["high.col"],
                                    guide = "colorbar",
                                    name = cName)
  }
  # plot trees
  if (reverse) {
    ctree <- ggtree(reduced.phy,
                    ladderize = ladderize,
                    aes_string(color = cDisplay), ...) +
      cColors + scale_x_reverse() + theme(legend.position = "bottom")
  } else {
    ctree <- ggtree(reduced.phy,
                    ladderize = ladderize,
                    aes_string(color = cDisplay), ...) +
      cColors + theme(legend.position = "bottom")
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
