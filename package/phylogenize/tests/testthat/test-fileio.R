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

test.abd.meta.16s <- generate.fake.abd.meta(n.samples=100,
                                            n.taxa=1000,
                                            n.envs=3,
                                            n.dsets=1,
                                            n.reads=1e6,
                                            env.fx.sizes=c(0,1,-1,0.5,-0.5),
                                            dset.fx.sizes=c(0,0.1,0.5),
                                            env.frac.affected=0.3,
                                            dset.frac.affected=0.1,
                                            make.16s=TRUE,
                                            data_dir='../../data')

s
