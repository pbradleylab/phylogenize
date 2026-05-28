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
