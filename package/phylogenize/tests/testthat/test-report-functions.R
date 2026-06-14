test_that("non.interactive.plot adds an explicit visible branch layer", {
    tr <- ape::read.tree(text="(s1:1,s2:1);")
    tree.obj <- list(
        tree=list(data=data.frame(label=c("s1", "s2"), color=c(1, 2))),
        rphy=tr,
        cols=c(low.col="black", high.col="orange2"),
        lims=c(0, 2)
    )

    p <- non.interactive.plot(tree.obj)
    layer_colors <- vapply(p$layers, function(layer) {
        color <- layer$aes_params[["colour"]]
        if (is.null(color)) color <- layer$aes_params[["color"]]
        if (is.null(color)) return(NA_character_)
        color
    }, character(1))

    expect_true("grey35" %in% layer_colors)
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
