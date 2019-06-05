### IO TESTS

context("IO")
library(phylogenize)

test.abd.meta <- generate.fake.abd.meta(n.samples=100,
                                        n.taxa=1000,
                                        n.envs=3,
                                        n.dsets=1,
                                        n.reads=1e6,
                                        env.fx.sizes=c(0,1,-1,0.5,-0.5),
                                        dset.fx.sizes=c(0,0.1,0.5),
                                        env.frac.affected=0.3,
                                        dset.frac.affected=0.1)

test_that("biom import matches original table", {
    write.test.biom(test.abd.meta,
                    biom_file='test.biom',
                    in_dir=tempdir(),
                    overwrite=TRUE)
    compare.abd.meta <- read.abd.metadata(input_format='biom',
                                          biom_file='test.biom',
                                          in_dir=tempdir())
    expect_equal(compare.abd.meta$mtx,
                 test.abd.meta$mtx)
    expect_equal(compare.abd.meta$metadata$sample,
                 as.character(test.abd.meta$metadata$sample))
    expect_equal(compare.abd.meta$metadata[[pz.options('env_column')]],
                 test.abd.meta$metadata[[pz.options('env_column')]])
    expect_equal(compare.abd.meta$metadata[[pz.options('dset_column')]],
                 test.abd.meta$metadata[[pz.options('dset_column')]])
})

test_that("tabular import matches original table", {
    write.test.tabular(test.abd.meta,
                       in_dir=tempdir(),
                       overwrite=TRUE)
    compare.abd.meta <- read.abd.metadata(input_format='tabular',
                                          in_dir=tempdir())
    expect_equal(compare.abd.meta$mtx,
                 test.abd.meta$mtx)
    expect_equal(compare.abd.meta$metadata$sample,
                 as.character(test.abd.meta$metadata$sample))
    expect_equal(compare.abd.meta$metadata[[pz.options('env_column')]],
                 test.abd.meta$metadata[[pz.options('env_column')]])
    expect_equal(compare.abd.meta$metadata[[pz.options('dset_column')]],
                 test.abd.meta$metadata[[pz.options('dset_column')]])
})

# These BURST assignments should be unique because they have been clustered at
# 90 percent similarity with any duplicate species removed

test.abd.meta.16s <- generate.fake.abd.meta(n.samples=100,
                                            n.taxa=750,
                                            n.envs=3,
                                            n.dsets=1,
                                            n.reads=1e6,
                                            env.fx.sizes=c(0,1,-1,0.5,-0.5),
                                            dset.fx.sizes=c(0,0.1,0.5),
                                            env.frac.affected=0.3,
                                            dset.frac.affected=0.1,
                                            make.16s=TRUE,
                                            data_dir='../../data',
                                            burst_16sfile=
                                                paste('16s_centroids',
                                                      '90_filt500_nodups.fa',
                                                      sep='_'),
                                            tag.length=200)

test_that("burst mapping works appropriately", {
    processed.test.16s <- process.16s(abd.meta=test.abd.meta.16s,
                                      data_dir=system.file(package="phylogenize",
                                                           'data'),
                                      burst_dir='/home/pbradz/bin/',
                                      burst_16sfile=
                                          '16s_centroids_90_filt500_nodups.fa')
    expect_equal(processed.test.16s$n,
                 test.abd.meta.16s$n)
    expect_equal(sort(processed.test.16s$mtx %>% rownames),
                 sort(test.abd.meta.16s$n))
    testmtx <- Matrix::Matrix(test.abd.meta.16s$mtx > 0)
    rownames(testmtx) <- test.abd.meta.16s$n
    expect_equal(testmtx[rownames(processed.test.16s$mtx), ],
                Matrix::Matrix(processed.test.16s$mtx > 0))
})
