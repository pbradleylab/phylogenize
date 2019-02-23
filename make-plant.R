### 2018 Jul 09
rmarkdown::render("phylogenize-report.Rmd",
  output_file = "./emp/emp-plant-specificity-retest.html",
  params = list(out_dir = "./emp/output-plant-specificity-2/",
    in_dir = ".",
    type = "16S",
    db_version = "midas_v1.2",
    which_phenotype = "specificity",
    which_envir = "Plant",
    abundance_file = "",
    metadata_file = "",
    biom_file = "emp_deblur_freeliving.biom",
    input_format = "biom",
		burst_dir = "./",
    ncl = 10))

devtools::load_all("package/phylogenize")
phylogenize::render.report.alt(
  out_dir = "/home/pbradz/projects/phylogenize/jordan/jlp16s-3/",
  in_dir = "/home/pbradz/projects/phylogenize",
  output_file ="/home/pbradz/projects/phylogenize/jordan/jlp16s-4.html",
  type = "16S",
  db_version = "midas_v1.2",
  which_phenotype = "prevalence",
  which_envir = "Diet",
  abundance_file = "j16s-abundance-reduced.tab",
  metadata_file = "metadata-diet-date-interaction.tab",
  biom_file = "",
  input_format = "tabular",
  burst_dir = "/home/pbradz/projects/phylogenize/",
  data_dir = "/home/pbradz/projects/phylogenize/data/",
  devel = TRUE,
  devel_pkgdir = "/home/pbradz/projects/phylogenize/package/phylogenize/",
  ncl = 10)

library(phylogenize)


devtools::load_all("package/phylogenize")
setwd("~/projects/phylogenize")
dir.create("/home/pbradz/projects/phylogenize/emp-new/output-plant-specificity-retest2/")
phylogenize::render.report.alt(
  output_file = "/home/pbradz/projects/phylogenize/emp-new/emp-plant-specificity-retest2.html",
  out_dir = "/home/pbradz/projects/phylogenize/emp-new/output-plant-specificity-retest2/",
  in_dir = "/home/pbradz/projects/phylogenize",
  type = "16S",
  db_version = "midas_v1.2",
  which_phenotype = "specificity",
  which_envir = "Plant",
  abundance_file = "",
  metadata_file = "",
  biom_file = "emp_deblur_freeliving.biom",
  input_format = "biom",
  burst_dir = "/home/pbradz/bin/",
  ncl = 10,
  data_dir = "/home/pbradz/projects/phylogenize/data/",
  devel = TRUE,
  devel_pkgdir = "/home/pbradz/projects/phylogenize/package/phylogenize/",
  use_rmd_params = FALSE)

###TESTING###

pz.options(
                 out_file = "/home/pbradz/projects/phylogenize/emp-new/emp-plant-specificity-retest.html",
                 out_dir = "/home/pbradz/projects/phylogenize/emp-new/output-plant-specificity-retest/",
                 in_dir = "/home/pbradz/projects/phylogenize",
                 type = "16S",
                 db_version = "midas_v1.2",
                 which_phenotype = "specificity",
                 which_envir = "Plant",
                 abundance_file = "",
                 metadata_file = "",
                 biom_file = "emp_deblur_freeliving.biom",
                 input_format = "biom",
                 burst_dir = "/home/pbradz/bin/",
                 ncl = 1,
                 data_dir = "/home/pbradz/projects/phylogenize/data/",
                 burst_dir = "/home/pbradz/bin/",
                 use_rmd_params = FALSE)




pz.options(
  output_file = "/home/pbradz/projects/phylogenize/emp/emp-plant-specificity-retest.html",
  out_dir = "/home/pbradz/projects/phylogenize/emp/output-plant-specificity-retest/",
  in_dir = "/home/pbradz/projects/phylogenize",
  type = "16S",
  db_version = "midas_v1.2",
  which_phenotype = "specificity",
  which_envir = "Plant",
  abundance_file = "",
  metadata_file = "",
  biom_file = "emp_deblur_freeliving.biom",
  input_format = "biom",
  burst_dir = "/home/pbradz/bin/",
  ncl = 10,
  use_rmd_params = FALSE)

devtools::document("package/phylogenize")
devtools::load_all("package/phylogenize")
setwd("~/projects/phylogenize")
pz.options(
    out_dir = "/home/pbradz/projects/phylogenize/testing-out/",
    in_dir = "/home/pbradz/projects/phylogenize/testing-in/",
    type = "16S",
    db_version = "midas_v1.2",
    which_phenotype = "prevalence",
    which_envir = "env1",
    abundance_file = "",
    metadata_file = "",
    data_dir = "/home/pbradz/projects/phylogenize/data",
    biom_file = "test.biom",
    input_format = "biom",
    biom_dir = "/usr/local/bin/",
    burst_dir = "/home/pbradz/bin/",
    burst_16sfile = '16s_centroids_90_filt500_nodups.fa',
    ncl = 10,
    use_rmd_params = FALSE)
generate.test.pzdb(db_version="midas_v1.2")
test.data <- generate.fake.abd.meta(n.samples=50, n.taxa=150, make.16s=TRUE)
dir.create("./testing-in")
dir.create("./testing-out")
write.test.biom(test.data, overwrite=TRUE)





devtools::load_all("package/phylogenize")
setwd("~/projects/phylogenize")
phylogenize::render.report.alt(
  output_file = "/home/pbradz/projects/phylogenize/testing-out-2/testing-prev-output.html",
  in_dir = "/home/pbradz/projects/phylogenize/testing-in/",
  out_dir = "/home/pbradz/projects/phylogenize/testing-out-2/",
  type = "16S",
  db_version = "midas_v1.2",
#  db_version = "test",
  which_phenotype = "prevalence",
  which_envir = "env1",
  abundance_file = "",
  metadata_file = "",
  data_dir = "/home/pbradz/projects/phylogenize/data/",
  biom_file = "test.biom",
  input_format = "biom",
  burst_dir = "/home/pbradz/bin/",
  burst_16sfile = '16s_centroids_90_filt500_nodups.fa',
  ncl = 10,
  use_rmd_params = FALSE,
  devel = TRUE,
  devel_pkgdir = file.path(getwd(), "package/phylogenize"),
  pryr = FALSE)

pz.options(
    out_dir = "/home/pbradz/projects/phylogenize/testing-out-3/",
    in_dir = "/home/pbradz/projects/phylogenize/testing-in-3/",
    type = "16S",
    db_version = "midas_v1.2",
    which_phenotype = "specificity",
    ncl = 10,
    biom_dir = "/home/pbradz/.local/bin/",
    burst_dir = "/home/pbradz/bin/",
    use_rmd_params = FALSE)
dir.create("./testing-in-3")
dir.create("./testing-out-3")
test.data <- generate.fake.abd.meta(n.samples=50, n.taxa=150, make.16s=TRUE, n.dsets=1)
write.test.biom(test.data, overwrite=TRUE)

devtools::load_all("package/phylogenize")
setwd("~/projects/phylogenize")
phylogenize::render.report.alt(
  output_file = "/home/pbradz/projects/phylogenize/testing-out-3/testing-spec-output.html",
  in_dir = "/home/pbradz/projects/phylogenize/testing-in-3/",
  out_dir = "/home/pbradz/projects/phylogenize/testing-out-3/",
  type = "16S",
  db_version = "midas_v1.2",
  which_phenotype = "prevalence",
  which_envir = "env1",
  abundance_file = "",
  metadata_file = "",
  data_dir = "/home/pbradz/projects/phylogenize/data/",
  biom_file = "test.biom",
  input_format = "biom",
  burst_dir = "/home/pbradz/bin/",
  burst_16sfile = '16s_centroids_90_filt500_nodups.fa',
  ncl = 10,
  use_rmd_params = FALSE,
  devel = TRUE,
  devel_pkgdir = file.path(getwd(), "package/phylogenize"),
  pryr = FALSE)
