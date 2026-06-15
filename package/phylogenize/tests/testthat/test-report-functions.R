test_that("non.interactive.plot adds an explicit visible branch layer", {
    tr <- ape::read.tree(text="(s1:1,s2:1);")
    tree.obj <- list(
        tree=list(data=data.frame(label=c("s1", "s2"), color=c(1, 2))),
        rphy=tr,
        cols=c(low.col="black", high.col="orange2"),
        lims=c(0, 2)
    )

    p <- suppressWarnings(non.interactive.plot(tree.obj))
    layer_colors <- vapply(p$layers, function(layer) {
        color <- layer$aes_params[["colour"]]
        if (is.null(color)) color <- layer$aes_params[["color"]]
        if (is.null(color)) return(NA_character_)
        color
    }, character(1))

    expect_true("grey35" %in% layer_colors)
})

test_that("change_tree_plot_internals maps cluster tip labels to species labels", {
    taxonomy <- data.frame(
        cluster=c("cluster1", "cluster2"),
        species=c("Species one", "Species two")
    )
    ctree <- list(
        data=data.frame(
            node=c(1, 2, 3),
            label=c("cluster1", "cluster2", NA)
        ),
        mapping=list(colour=c(cluster1=1.25, cluster2=-0.5))
    )

    updated <- change_tree_plot_internals(taxonomy, NULL, ctree)

    tip_labels <- updated$data$label[updated$data$node %in% c(1, 2)]
    expect_true(any(grepl("Species one", tip_labels, fixed=TRUE)))
    expect_true(any(grepl("Species two", tip_labels, fixed=TRUE)))
    expect_false(any(grepl("cluster1", tip_labels, fixed=TRUE)))
    expect_equal(names(updated$mapping$colour), c("Species one", "Species two"))
})

test_that("tree.branch.segments extracts plotly-safe tree coordinates", {
    tree <- list(data=data.frame(
        node = c(3, 1, 2),
        parent = c(3, 3, 3),
        x = c(0, 1, 1),
        y = c(1.5, 1, 2)
    ))

    segments <- tree.branch.segments(tree)

    expect_equal(nrow(segments), 3)
    expect_true(any(segments$x == 0 & segments$xend == 1 &
                    segments$y == 1 & segments$yend == 1))
    expect_true(any(segments$x == 0 & segments$xend == 1 &
                    segments$y == 2 & segments$yend == 2))
    expect_true(any(segments$x == 0 & segments$xend == 0 &
                    segments$y == 1 & segments$yend == 2))
})

test_that("add.plotly.tree.branches appends branch traces to interactive trees", {
    tree <- list(data=data.frame(
        node = c(3, 1, 2),
        parent = c(3, 3, 3),
        x = c(0, 1, 1),
        y = c(1.5, 1, 2)
    ))
    p <- plotly::plot_ly()

    p <- add.plotly.tree.branches(p, tree)

    expect_s3_class(p, "plotly")
    expect_true(length(p$x$attrs) > 0 || length(p$x$data) > 0)
})
