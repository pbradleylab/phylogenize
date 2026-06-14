test_that("fix.tree repairs missing and invalid branch lengths", {
    no_lengths <- ape::read.tree(text="(s1,s2,(s3,s4));")
    fixed_missing <- fix.tree(no_lengths)

    expect_true(inherits(fixed_missing, "phylo"))
    expect_length(fixed_missing$edge.length, nrow(fixed_missing$edge))
    expect_true(all(is.finite(fixed_missing$edge.length)))
    expect_true(all(fixed_missing$edge.length > 0))

    bad_lengths <- ape::read.tree(text="((s1:1,s2:1):1,s3:1);")
    bad_lengths$edge.length <- c(Inf, 0, -1, NA)
    fixed_bad <- fix.tree(bad_lengths)

    expect_true(all(is.finite(fixed_bad$edge.length)))
    expect_true(all(fixed_bad$edge.length > 0))
})

test_that("fix.tree leaves non-tree inputs unchanged", {
    expect_null(fix.tree(NULL))
    expect_equal(fix.tree("not a tree"), "not a tree")
})
