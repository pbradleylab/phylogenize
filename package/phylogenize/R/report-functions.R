#--- Plotting ---#

#' Get phenotype plotting scales.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence"
#'   or "specificity" or "correlation").}
#'   \item{prev_color_low}{String. When graphing prevalence on a tree, this
#'   color is the lowest value.} \item{prev_color_high}{String. When graphing
#'   prevalence on a tree, this color is the highest value.}
#'   \item{spec_color_high}{String. When graphing specificity on a tree, this
#'   color is the lowest value (most anti-specific).}
#'   \item{spec_color_med}{String. When graphing specificity on a tree, this
#'   color denotes the prior (no association).} \item{spec_color_high}{String.
#'   When graphing specificity on a tree, this color is the highest value (most
#'   specific).}
#' }
#' @param phenotype A named vector giving the phenotype for each taxon ID.
#' @param trees A list of tree objects.
#' @param phenoP An optional value giving the prior probability for the
#'     environment of interest.
#' @return A list of overall limits (\code{limits}), taxon-specific limits
#'     (\code{phy.limits}), a color scale (\code{colors}), and the zero point
#'     (\code{zero}).
#' @export
get.pheno.plotting.scales <- function(phenotype, trees, phenoP=0, ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    if (opts('which_phenotype') == 'prevalence') {
        get.pheno.plotting.scales.prevalence(phenotype,
                                             trees,
                                             phenoP,
                                             ...)
    } else if (opts('which_phenotype') %in% c('specificity', 'abundance', 'provided')) {
        get.pheno.plotting.scales.specificity(phenotype,
                                              trees,
                                              phenoP,
                                              ...)
    }
}


#' Get phenotype plotting scales (prevalence-specific).
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{prev_color_low}{String. When graphing prevalence on a tree, this
#'   color is the lowest value.}
#'   \item{prev_color_high}{String. When graphing prevalence on a tree, this
#'   color is the highest value.}
#' }
#' @param phenotype A named vector giving the phenotype for each taxon ID.
#' @param trees A list of tree objects.
#' @param phenoP An optional value giving the prior probability for the
#'     environment of interest.
#' @return A list of overall limits (\code{limits}), taxon-specific limits
#'     (\code{phy.limits}), a color scale (\code{colors}), and the zero point
#'     (\code{zero}).
#' @keywords internal
get.pheno.plotting.scales.prevalence <- function(phenotype,
                                                 trees,
                                                 phenoP=0,
                                                 ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    phenoLimits <- quantile(unique(phenotype), c(0.2, 0.8))
    phenoLimitsTaxon <- lapply(trees, function(tr) {      
        phi <- phenotype[intersect(names(phenotype), tr$tip.label)] %>%
            na.omit
	if (length(unique(phi)) > 1) {
		quantile(unique(phi), c(0.2, 0.8))
        } else {
                return(NULL)
        }
    })
    #phenoLimitsTaxon[!sapply(phenoLimitsTaxon, is.null)]
    phenoColors <- c(low.col=opts('prev_color_low'),
                     high.col=opts('prev_color_high'))
    return(list(limits=phenoLimits,
                phy.limits=phenoLimitsTaxon,
                colors=phenoColors,
                zero=0))
}


#' Get phenotype plotting scales (specificity or abundance).
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{spec_color_high}{String. When graphing specificity on a tree, this
#'   color is the lowest value (most anti-specific).}
#'   \item{spec_color_med}{String. When graphing specificity on a tree, this
#'   color denotes the prior (no association).} \item{spec_color_high}{String.
#'   When graphing specificity on a tree, this color is the highest value (most
#'   specific).}
#' }
#' @param phenotype A named vector giving the phenotype for each taxon ID.
#' @param trees A list of tree objects.
#' @param phenoP An optional value giving the prior probability for the
#'     environment of interest.
#' @return A list of overall limits (\code{limits}), taxon-specific limits
#'     (\code{phy.limits}), a color scale (\code{colors}), and the zero point
#'     (\code{zero}).
#' @keywords internal
get.pheno.plotting.scales.specificity <- function(phenotype,
                                                  trees,
                                                  phenoP=0,
                                                  ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    phenoLimits <- c(phenoP - 1 * sd(phenotype),
                     phenoP + 1 * sd(phenotype))
    phenoLimitsTaxon <- lapply(trees, function(tr) {
        phi <- phenotype[intersect(names(phenotype), tr$tip.label)] %>%
            na.omit %>%
            unique
        c(phenoP - (1 * sd(phi)), phenoP + (1 * sd(phi)))
    })
    phenoColors <- c(low.col=opts('spec_color_low'),
                     mid.col=opts('spec_color_mid'),
                     high.col=opts('spec_color_high'))
    return(list(limits=phenoLimits,
                phy.limits=phenoLimitsTaxon,
                colors=phenoColors,
                zero=phenoP))
}

#' Plot a phenotype along a list of trees.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence"
#'   or "specificity").}
#' }
#'
#' @param phenotype A named vector with the phenotype values for each taxon.
#' @param trees A list of trees.
#' @param taxonomy A dataframe with the taxonomix information.
#' @param scale A list returned from \code{get.pheno.plotting.scales}.
#' @return A list of ggtree objects in which the phenotype has been plotted
#'     across each tree in \code{trees}.
#' @export plot.phenotype.trees
plot.phenotype.trees <- function(phenotype,
                                 trees,
                                 taxonomy,
                                 scale,
                                 ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    if (any(!(names(trees) %in% names(scale$phy.limits)))) {
        pz.error("taxon-specific limits must be calculated for every tree")
    }
    
    plotted.pheno.trees <- lapply(names(trees), function(tn) {
        gg.cont.tree(trees[[tn]],
                     phenotype,
                     taxonomy,
                     cLimits=scale$phy.limits[[tn]],
                     colors=scale$colors,
                     cName=tn,
                     plot=FALSE)#,
    })
    if(!is.null(plotted.pheno.trees)) {
        names(plotted.pheno.trees) <- names(trees)
        plotted.pheno.trees <- plotted.pheno.trees[
            vapply(plotted.pheno.trees, is.list, TRUE)
        ]
        if (length(plotted.pheno.trees) == 0) {
            pz.warning("No trees were plotted: too few taxa for any taxon?")
        }
        return(plotted.pheno.trees)
    }
}
    

#' Plot distributions of a phenotype across taxa.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence"
#'   or "specificity").}
#' }
#' @param phenotype A named vector with the phenotype values for each species.
#' @param pz.db A database containing a \code{taxonomy} and \code{trees}.
#' @return A ggplot object with the phenotype distribution plotted per taxon.
#' @export plot.pheno.distributions
plot.pheno.distributions <- function(phenotype,
                                     pz.db,
                                     ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    kept.species <- Reduce(c, lapply(pz.db$trees, function(x) x$tip.label))
    pheno.taxon <- pz.db$taxonomy[match(names(phenotype),
                                        pz.db$taxonomy$cluster),
                                   pz.options("taxon_level"),
                                   drop=TRUE]
    pheno.characteristics <- data.frame(pheno=phenotype,
                                        taxon=pheno.taxon,
                                        cluster=names(phenotype))
    sub.pheno <- subset(pheno.characteristics, cluster %in% kept.species)
    taxon_uniq <- unique(sub.pheno$taxon)
    # Get color pallet
    color_palette <- scales::hue_pal()(length(taxon_uniq))

    # Remove any rows with NA as they were unable to be found by by that
    # taxonomic classification in the codebase
    sub.pheno <- sub.pheno[!(sub.pheno$taxon == "" | is.na(sub.pheno$taxon)), ]

    distros <- sub.pheno %>%
	    dplyr::group_by(taxon) %>%
	    tidyr::nest() %>%
	    dplyr::mutate(n_datapoints = purrr::map_dbl(data, nrow)) %>%
	    dplyr::filter(n_datapoints >= 3) %>%
	    dplyr::arrange(-n_datapoints) %>%
        dplyr::mutate(plots = purrr::map2(data, taxon, ~ {
            p <- ggplot2::ggplot(.x, ggplot2::aes(pheno, fill=.y)) +
                ggplot2::geom_density() +  
                ggplot2::xlab(opts('which_phenotype')) +
                ggplot2::ggtitle(.y) +
                ggplot2::theme_minimal() + 
                ggplot2::scale_fill_manual(values = setNames(color_palette,
                                                             taxon_uniq)) + 
                ggplot2::guides(fill = "none")
            return(p)
        })) %>%
        dplyr::ungroup() %>%
        dplyr::select(plots) %>%
        tibble::deframe()
    
    # Group the ggplots into sets of 10 with them already being ordered from the
    # groups with the most individuals to the least for node tips in the kept
    # taxon
    combine_plots <- function(plots, ncol = 2) {
	    patchwork::wrap_plots(plots, ncol = ncol)
    }
    # This number can be adjusted for how many plots per facet.
    # Anything above 6 starts to squish the axis however.
    grouped_plots <- split(distros, ceiling(seq_along(distros) / 6)) 
    distros <- lapply(grouped_plots, combine_plots)

    return(distros)
}


#' Edit a list of plotted trees to add fancy highlight labels.
#'
#' This function adds fancy SVG highlight labels to ggtree objects and then
#' plots them. If there's an error, it will fall back to a regular plot.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence"
#'   or "specificity").}
#' }
#'
#' @param plotted.pheno.trees A named list of ggtree plots (per taxon).
#' @param phenotype Phenotype to plot/label.
#' @param label Label to give to the phenotype.
#' @param stroke.scale How thick to make the highlight.
#' @param units A string appended to each label, used to give units of
#'   phenotype.
#' @return list containing plot.labeled.phenotype.trees
#' @export
plot.labeled.phenotype.trees <- function(plotted.pheno.trees,
                                         phenotype,
                                         label='prevalence',
                                         stroke.scale=0.3,
                                         units='%') {
    if (is.null(plotted.pheno.trees)) {
        pz.message("warning: no trees found")
        return(NULL)
    }
    if (length(plotted.pheno.trees) == 0) {
        pz.message("warning: no trees found")
        return(NULL)
    }

    plots <- list()
    taxon_names <- names(plotted.pheno.trees)
    for (tree in 1:length(plotted.pheno.trees)) {
	    plotted_tree <- plotted.pheno.trees[[tree]]
	    fn <- knitr::fig_path('svg', number = tree)
	    name <- taxon_names[[tree]]
	    tryCatch(
		     plots[[name]] <- interactive.plot(plotted_tree, fn, name),
		     error = function(e) {
			    pz.message(e)
			    plots[[name]] <- non.interactive.plot(plotted_tree, fn, name)
	    })
    }
    return(plots)
}

#' Make a hybrid tree-heatmap plot showing the taxon distribution of significant
#' hits.
#'
#' This function wraps single.cluster.plot, running it in a separate process.
#' This is because there can be problems with memory leaks. For relevant global
#' options, see the documentation for `pz.options`.
#'
#' @param gene.presence Gene presence/absence matrix.
#' @param sig.genes Character vector of the significant genes.
#' @param tree A tree object.
#' @param plotted.tree A ggtree plot of \code{tree}.
#' @param taxon Name of the taxon represented by \code{tree}
#' @param verbose Whether to report debugging information (boolean).
#' @export do.clust.plot
do.clust.plot <- function(gene.presence,
                          sig.genes,
                          tree,
                          plotted.tree,
                          taxon,
                          verbose=FALSE,
                          ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    # Run these on a separate process to avoid memory leak
    cl <- parallel::makeCluster(1)
    if (verbose) message("exporting data...")
    parallel::clusterExport(cl,
                  c("gene.presence",
                    "sig.genes",
                    "tree",
                    "plotted.tree",
                    "taxon",
                    "verbose",
                    "PZ_OPTIONS"),
                  envir=environment())
    if (verbose) message("importing source...")
    # parallel::clusterCall(cl, library, phylogenize)
    cluster.load.pkg(cl, opts("devel"), opts("devel_pkgdir"))
    if (verbose) message("performing call...")
    tmpL <- parallel::clusterCall(cl,
                        single.cluster.plot,
                        gene.presence,
                        sig.genes,
                        tree,
                        plotted.tree,
                        taxon,
                        verbose,
                        ...)
    print(tmpL[[1]])
    rm(tmpL)
    parallel::stopCluster(cl)
    gc()
    return(NULL) # Avoid wasting memory since we never touch these
}

#' Make a hybrid tree-heatmap plot showing the taxon distribution of significant
#' hits.
#'
#' Some particularly relevant global options are:
#' \describe{
#'   \item{which_phenotype}{String. Which phenotype to calculate ("prevalence"
#'   or "specificity").}
#'   \item{gene_color_absent}{String. When graphing gene presence/absence,
#'   this color indicates absence.}
#'   \item{gene_color_present}{String. When graphing gene presence/absence, this
#'   color indicates presence.}
#' }
#'
#' @param gene.presence Gene presence/absence matrix.
#' @param sig.genes Character vector of the significant genes.
#' @param tree A tree object.
#' @param plotted.tree A ggtree plot of \code{tree}.
#' @param taxon Name of the taxon represented by \code{tree}
#' @param verbose Whether to report debugging information (boolean).
#' @return A faceted ggplot object.
#' @export
single.cluster.plot <- function(gene.presence,
                                sig.genes,
                                tree,
                                plotted.tree,
                                taxon,
                                verbose=FALSE,
                                ...) {
    opts <- clone_and_merge(PZ_OPTIONS, ...)
    sig.bin <- gene.presence[intersect(rownames(gene.presence),
                                       sig.genes), , drop=FALSE]
    if (is.null(dim(sig.bin))) {
        sig.bin <- as.matrix(sig.bin) %>% t
        rownames(sig.bin) <- sig.genes
    }
    names(dimnames(sig.bin)) <- c("gene", "id")
    col.df <- data.frame(id=names(plotted.tree$disp),
                         disp=plotted.tree$disp)
    p <- ggtree::ggtree(tree, ladderize = TRUE) %<+%
        col.df +
        ggtree::geom_tippoint(ggplot2::aes(color=disp, shape='.', guide="none"))
    if (opts('which_phenotype') == "prevalence") {
        p <- p + ggplot2::scale_color_gradient(
            low=plotted.tree$cols["low.col"],
            high=plotted.tree$cols["high.col"])
    } else if (opts('which_phenotype') %in% c("specificity",
                                              "correlation",
                                              "abundance",
                                              "provided")) {
        p <- p + ggplot2::scale_color_gradient2(
            low=plotted.tree$cols["low.col"],
            mid=plotted.tree$cols["mid.col"],
            high=plotted.tree$cols["high.col"],
            midpoint=mean(plotted.tree$lims))
    }
    
    if (length(sig.genes) == 0) { return(p) }
    if (length(sig.genes) > 1) {
        clust <- hclust(dist(sig.bin, method = "binary"))
        sig.ord <- sparseMelt(t(sig.bin)[, clust$order, drop=FALSE])
    } else {
        sig.ord <- sparseMelt(t(sig.bin))
        sig.ord <- sig.ord[order(sig.ord[, 3]), , drop=FALSE]
    }
    tmp <- ggtree::facet_plot(p,
                              panel=paste0('heatmap: ', taxon),
                              data=sig.ord,
                              geom=ggplot2::geom_tile,
                              mapping=ggplot2::aes(
                                  x = as.numeric(as.factor(gene)),
                                  fill = as.numeric(as.factor(value)))) +
        ggplot2::scale_fill_gradient(low = opts('gene_color_absent'),
                                     high = opts('gene_color_present'),
                                     na.value = opts('gene_color_absent')) +
        ggplot2::labs(color=opts("which_phenotype"), fill="gene presence") +
        ggplot2::scale_shape(guide="none")
    tmp
}


#' Make a pretty enrichment table.
#'
#' @param enr.table Input enrichment table.
#' @return Mutated enrichment table with better-labeled columns and significance
#'   coloring.
#' @export
output.enr.table <- function(enr.table) {
    enr.table %>%
        dplyr::mutate(
            Gene_significance=capwords(Gene_significance),
            Taxon=capwords(Taxon)
        ) %>%
        dplyr::mutate(
            q_value=kableExtra::cell_spec(
                prettyNum(q_value, digits=2),
                "html",
                background=kable.recolor(
                    -log10(q_value),
                    colors=c("#FFFFFF", "#FF8888"),
                    limits=c(0, 10)))) %>%
        dplyr::mutate(
            Enrichment_odds_ratio=kableExtra::cell_spec(
                prettyNum(Enrichment_odds_ratio, digits=3),
                "html",
                background=kable.recolor(
                    Enrichment_odds_ratio,
                    colors=c("#FFFFFF", "#FFFF44"),
                    limits=c(1, 10)))) %>%
        knitr::kable("html", escape=FALSE, align="l") %>%
        kableExtra::kable_styling(c("striped", "condensed"))
}

                                     
#' Helper function for gg.cont.tree
#'
#' \code{gg.cont.tree} generates the right internals for later use in plotting
#' interactive trees.
#'
#' @param taxonomy A dataframe with the taxonomic information.
#' @param reduced.phy A reduced phylogenetic tree (phylo)
#' @param ctree a graphed phylogenetic tree object
#' @return a graphed phylogenetic tree object with edited labels and color
#'   column
change_tree_plot_internals <- function(taxonomy, reduced.phy, ctree) {
    # Map the colors to the tree object directly to use ggplotly later for
    # interactive trees
    tax_filt <- taxonomy[taxonomy$cluster %in% ctree$data$label,
                         c("species", "cluster")]
    ctree_data <- merge(ctree$data,
                        tax_filt,
                        by.x = "label",
                        by.y = "species",
                        all.x = TRUE)
    tax_colors <- data.frame(label = as.character(names(ctree$mapping$colour)),
                             color = ctree$mapping$colour
    )
    ctree_data <- merge(ctree_data, tax_colors, by = "label", all.x = TRUE)
    ctree_data$color <- ifelse(is.na(ctree_data$color), 0, ctree_data$color)
    
    colnames(taxonomy)[colnames(taxonomy) == "cluster"] <- "label"
    ctree_data <- merge(ctree_data, taxonomy[c("label", "species")],
                        by = "label",
                        all.x = TRUE)
    ctree_data$label <- ctree_data$species
    
    ctree_data <- merge(ctree$data,
                        ctree_data[c("node", "species", "color")],
                        by = "node",
                        all.x = TRUE)
    ctree_data$label <- ctree_data$species
    ctree$data <- ctree_data
    
    ctree$data$rounded <- signif(ctree$data$color, digits = 2)
    ctree$data$label <- ifelse(
        is.na(ctree$data$label),
        ctree$data$label,
        paste0(ctree$data$label, " (phenotype: ", ctree$data$rounded, ")")
    )
    
    #Change color section in tree to avoid error
    swap <- tax_filt %>%
        dplyr::mutate(cluster = as.character(cluster)) %>%
        dplyr::select(cluster, species) 
    
    color_head <- data.frame(cluster = names(ctree$mapping$colour),
                             stringsAsFactors = FALSE)
    color_head <- color_head %>%
        dplyr::left_join(swap, by = "cluster")
    
    color_head$final_name <- ifelse(is.na(color_head$species),
                                    color_head$cluster,
                                    color_head$species)
    
    # Only assign names when lengths match
    names(ctree$mapping$colour) <- color_head$final_name
    
    return(ctree)
}


#' Plot continuous trait on a tree.
#'
#' \code{gg.cont.tree} paints a continuous trait along a tree.
#'
#' @param phy A \code{phylo} object.
#' @param ctrait A named numeric vector assigning trait values to tree tips.
#' @param taxonomy A dataframe with the taxonomic information.
#' @param cAnc Calculated ancestry for continuous trait; if this is NULL, it is
#'   calculated using castor::get_ancestral_nodes.
#' @param cLimits Scale bar limits for plotting continuous trait.
#' @param n Character vector giving which nodes to display; if NULL, defaults to
#'   intersection of \code{phy$tip.label} and \code{names(ctrait)}.
#' @param reduced.phy Dichotomous tree with only the nodes in \code{n}
#'   represented; if NULL, this is calculated.
#' @param colors Named character vector with at least "low.col" and "high.col"
#'   and optionally "mid.col" defined, giving colors to use for plotting.
#' @param plot Whether to plot the tree object or just return it.
#' @param restrict Character vector giving which continuous trait values to
#'   plot; if NULL, all are used.
#' @param cName String giving the title of the plot.
#' @param reverse Mirror the resulting plotted tree.
#' @param ladderize Ladderize the plotted tree.
#' @param ... Additional arguments passed to ggplot2.
#' @return A ggtree plot of a continuous trait plotted along a tree.
#' @export
gg.cont.tree <- function(phy,
                         ctrait,
                         taxonomy,
                         cAnc = NULL,
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
    
    kept_tips <- keep.tips(phy, n)
    
    if(!is.null(kept_tips)) {
        if(is.null(reduced.phy)) {reduced.phy <- fix.tree(kept_tips)}
        pz.message("getting continuous trait ancestry")
        if (is.null(reduced.phy)) {
            pz.warning(paste0("all tips were dropped from ", cName))
            return(NA)
        }
        if (length(reduced.phy$tip.label) < 2) {
            pz.warning(paste0(
                "need at least two tips for a tree; skipping ", cName))
            return(NA)
        }
        if (is.null(cAnc)) {
            #cAnc <- phytools::fastAnc(reduced.phy, ctrait[n])
            cAnc <- castor::asr_independent_contrasts(
                reduced.phy,
                ctrait[n][reduced.phy[["tip.label"]]]
            )[["ancestral_states"]]
        }
        cDisplay <- truncated(c(unlist(ctrait[n]),
                                unlist(cAnc)),
                              unlist(cLimits))
        
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
                                    ggplot2::aes_string(color = cDisplay),
                                    ...) +
                cColors + ggplot2::scale_x_reverse() +
                ggplot2::theme(legend.position = "bottom")
        } else {
            ctree <- ggtree::ggtree(reduced.phy,
                                    ladderize = ladderize,
                                    ggplot2::aes_string(color = cDisplay),
                                    ...) +
                cColors + ggplot2::theme(legend.position = "bottom")
        }
        
        # Make the plots so that they can be graphed interactively or non-interactively later on
        ctree <- change_tree_plot_internals(taxonomy, reduced.phy, ctree)
        
        if (plot) {ctree}
        
        return(list(tree = ctree,
                    cAnc = cAnc,
                    rphy = reduced.phy,
                    n = n,
                    cols = colors,
                    lims = cLimits,
                    disp = cDisplay))
    } else {
        pz.message(paste0(cName, " has been removed from the phylogenetic ",
                          "trees: No phenotype found associated"))
        return(NA)
    }
}


#' Create interactive tree diagrams.
#'
#' This function replaces a very old and brittle SVG hack with ggplotly!
#'
#' @param tree.obj A ggtree representation of a tree.
#' @param file A filename where the final SVG output will be written.
#' @param name String. the name of the taxon being added. Used for title.
#' @export
interactive.plot <- function(tree.obj, file, name) {
    tree <- non.interactive.plot(tree.obj, file, name)
    interactive_tree <- plotly::ggplotly(tree, tooltip = "text")
    return(interactive_tree)
}


#' A fall-back plotting option for when \code{hack.tree.labels} fails, designed
#' to produce the same kind of output.
#'
#' @param tree.obj A ggtree object.
#' @param file File to which an SVG representation of this tree object will be
#'   written.
#' @param name String. the name of the taxon being added. Used for title.
#' @export
non.interactive.plot <- function(tree.obj, file, name) {
    warning(paste0("replotting to: ", file))
    
    valid_labels <- subset(tree.obj$tree$data, !is.na(label))
    low_color <- tree.obj$cols["low.col"]
    high_color <- tree.obj$cols["high.col"]
    
    tree <- ggtree::ggtree(ape::as.phylo(tree.obj$tree)) +
        ggtree::geom_point(data = valid_labels,
                   ggplot2::aes(text = label, color = color)) +
        ggtree::geom_tiplab(data = valid_labels,
                    ggplot2::aes(color = color)) +
        ggplot2::ggtitle(name) +
        ggplot2::labs(color="phenotype")
    
    # Write to an svg
    svg <- svglite::xmlSVG(print(tree), standalone = TRUE)
    writeLines(as.character(svg), file)
    
    return(tree)
}
