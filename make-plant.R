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

rmarkdown::render("phylogenize-report.Rmd",
  output_file = "jordan-late-prevalence-16s.html",
  params = list(out_dir = "/home/pbradz/projects/phylogenize/jordan/jlp16s/",
    in_dir = "/home/pbradz/projects/phylogenize",
    source_dir = "/home/pbradz/projects/phylogenize",
    type = "16S",
    db_version = "midas_v1.2",
    which_phenotype = "prevalence",
    which_envir = "Diet",
    abundance_file = "j16s-abundance-reduced.tab",
    metadata_file = "metadata-diet-date-interaction.tab",
    biom_file = "",
    input_format = "tabular",
		burst_dir = "/home/pbradz/projects/phylogenize/",
    ncl = 1))

library(phylogenize)


devtools::load_all("package/phylogenize")
setwd("~/projects/phylogenize")
phylogenize::render.report.alt(
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
    biom_file = "test.biom",
    input_format = "biom",
    burst_dir = "/home/pbradz/bin/",
    burst_16sfile = '16s_centroids_90_filt500_nodups.fa',
    ncl = 10,
    use_rmd_params = FALSE)
test.data <- generate.fake.abd.meta(n.samples=50, n.taxa=150, make.16s=TRUE)
write.test.biom(test.data, overwrite=TRUE)
phylogenize::render.report.alt(
  output_file = "/home/pbradz/projects/phylogenize/testing-out/testing-prev-output.html",
  in_dir = "/home/pbradz/projects/phylogenize/testing-in/",
  type = "16S",
  db_version = "midas_v1.2",
  which_phenotype = "prevalence",
  which_envir = "env1",
  abundance_file = "",
  metadata_file = "",
  biom_file = "test.biom",
  input_format = "biom",
  burst_dir = "/home/pbradz/bin/",
  burst_16sfile = '16s_centroids_90_filt500_nodups.fa',
  ncl = 10,
  use_rmd_params = FALSE)
