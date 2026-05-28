# Phylogenize2 (v2.0.0-alpha)

Phylogenize2 allows users to link microbial genes to environments, accounting for phylogeny. More specifically, given community composition data Phylogenize2 links patterns of microbes in a given environment to genes in those microbes' pangenomes, while taking into account an important potential confounder: the phylogenetic relationships between microbes. We allow several different patterns to be calculated, including prevalence, specificity, and differential abundance (using [MaAsLin2](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009442) or [ANCOM-BC2](https://www.nature.com/articles/s41467-020-17041-7)). By default, we use phylogenetic regression, but we also allow users to apply the alternative method [POMS](https://academic.oup.com/bioinformatics/article/38/22/5055/6731923). The method is described in a forthcoming preprint (Kananen et al., in preparation).

In addition, we highly recommend using v2.0.0-alpha or later. We no longer support the use of earlier versions, as there are significant improvements since v0.91.

## Installing Phylogenize2

The easiest way to install all the dependencies needed is by using mamba or conda. We recommend using [miniforge3](https://github.com/conda-forge/miniforge). Please make sure you are using miniforge v3-23.3.1-0 or later. Miniforge3 is available for MacOS, Linux, and Windows OS. Phylogenize is not tested on Windows (proceed with caution); for all future examples, unless otherwise stated, we are assuming you are using Linux.

To install miniforge, run `wget -c https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh` and then in a terminal type `bash Miniforge3-Linux-x86_64.sh`. You'll need to run through the prompts and then give it a download location if the default is not to your liking. Finally, you can let miniforge initialize itself if you want conda to always be in your "base" when you open the terminal. Otherwise, if you say `no` then you will have to manually source the executable for conda which can be done similarly as `source /your/path/to/miniforge3/bin/activate`.

### Now you are ready to start installing the dependencies.

Create a new environment by typing `conda create -n phylogenize` and `conda activate phylogenize`. Then you can install phylogenize by running `conda install bioconda::phylogenize`. For any future analysis, all you have to do is activate this environment to have the dependencies run.

#### Locally - Command line and Rstudio (MacOS/Linux)

Please note, we assume in these instructions you are working off of base-r and NOT Rstudio. We describe at the bottom of this section how to use Rstudio while still installing the dependencies with mamba.

#### Install with mamba - configuration file

You can make a conda environment using the supplied yaml file and not worry about installing any dependencies. Run `conda env create -f environment.yml` and then `conda activate phylogenize`. Open base-r and then type `devtools::install_github("biocore/phylogenize")`.

#### Install with conda - no configuration file

1.  Make sure you have R installed. You can verify if you type `R --version`. If you don't you can get the latest version [here](https://www.r-project.org/) or install it using conda [here](https://anaconda.org/r/r).
2.  Create a new environment in conda by running `conda create -n phylogenize`
3.  Activate your new environment with `conda activate phylogenize`
4.  Install the dependencies with the bioconda and conda-forge channels as shown below

```         
mamba install -y bioconda::phylogenize
```

4.  Open R and then run `library("phylogenize")`.

#### Locally - Rstudio

After creating a `phylogenize` environment with conda using `conda create -n phylogenize` and installing phylogenize `conda install bioconda::phylogenize`, to use Rstudio run `conda install r::rstudio`. Then you can activate it by typing `rstudio` in your terminal. This will launch an Rstudio IDE. There, if you haven't already, you can run followed by `library("phylogenize")`.

### Installing Phylogenize2 package for use on AWS

We recommend you install Phylogenize2 in a conda environment as above. However, because the default Amazon images are meant for headless operation they are missing some tools to deal with fonts that Phylogenize2 uses to generate its plots. You can install those as follows:

```         
conda install -c conda-forge xorg-libxt
sudo apt install zlib
sudo apt show zlib1g
sudo apt install fontconfig
```

## Selecting a database

We have several premade databases that you can select from depending on what is expected to match your host's system. If you are unsure what database to use, then we recommend using GTDB as the default.

| Environment        | Version | Database | Number of families | Number of species |
|--------------------|---------|----------|--------------------|-------------------|
| chicken gut        | v1.0.1  | MGnify   | 142                | 1007              |
| cow rumen          | v1.0.1  | MGnify   | 121                | 1914              |
| honeybee gut       | v1.0.1  | MGnify   | 31                 | 131               |
| human gut          | v2.0.2  | MGnify   | 215                | 3445              |
| human oral         | v1.0.1  | MGnify   | 52                 | 260               |
| human vaginal      | v1.0    | MGnify   | 52                 | 189               |
| marine eukaryotes  | vbeta   | MGnify   | 250                | 250               |
| marine             | v2.0    | MGnify   | 1192               | 7408              |
| mouse gut          | v1.0    | MGnify   | 136                | 1639              |
| non model fish gut | v2.0    | MGnify   | 60                 | 87                |
| pig gut            | v1.0    | MGnify   | 138                | 800               |
| sheep rumen        | v1.0    | MGnify   | 117                | 2122              |
| zebrafish fecal    | v1.0    | MGnify   | 41                 | 24                |
| global             |         | GlobDB   |                    |                   |

All databases have been been matched against the UniRef50, FesNov, and UHGP databases, and any remaining protein sequences have been clustered *de novo*. Functional annotations have been obtained using [anvi'o](https://peerj.com/articles/1319/) and [KEGG](https://www.genome.jp/kegg/pathway.html) KOfams as described in Kananen et al., 2025.

Databases can be downloaded manually and decompressed from our Zenodo page [here](), or they can be downloaded and decompressed using Phylogenize2's `phylogenize::download.zenodo.db("your/html/link/here.zip")`. The default if no database is available is GTDB. If using a custom database, then all the database files must be placed into a directory called `package/inst/extdata/`.

## Preparing your data

If you are using shotgun metagenomes, you will need to first quantify species abundances. The species definitions and names must match the database you plan to use. We recommend using Kraken2 with Bracken, as there are Kraken2 databases for every MGnify database. (Make sure that the version numbers match!) For example:

 - Human gut v2.0.2 Kraken2 database: https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/kraken2_db_uhgg_v2.0.2/
 - Mouse gut Kraken2 database: https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/mouse-gut/v1.0/kraken2_db_mouse-gut_v1.0/
 - Marine Kraken2 database: https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/marine/v2.0/kraken2_db_marine_v2.0/

An example workflow for UHGG written in Snakemake can be seen under `shotgun_kraken2_example`.

Finally, you will want to make the taxon names from Bracken match the IDs in Phylogenize2, and also to merge any technical replicates for the same biological sample (as these will lead to overconfident predictions). There is a script to perform this under `shotgun_kraken2_example` called `parse-bracken.R`. You can run this script as follows:

```
Rscript parse_bracken.R -t [path to taxonomy file] -i [path to bracken output files] -o [path to output tab-separated file] -m [path to metadata file]
```

The last option (`-m`) is optional, but allows you to provide a tab-separated file with "sample" and "run" columns that will merge any runs belonging to the same sample. The taxonomy file provided should be the one in the Phylogenize2 database that you are using (e.g. `mouse-gut-taxonomy.csv`). (If you are having trouble finding the path where a database was installed, try looking under the directory where Phylogenize2 was installed, which you should be able to see by running `system.file(package="phylogenize")` in R.)

## Running Phylogenize2

Congratulations! Phylogenize2 should now be installed.

### Running Phylogenize2 using the R interface

Most users should start with the `phylogenize()` function. It runs the full workflow: reading the abundance and metadata files, calculating microbial phenotypes, testing genes for association with those phenotypes, running enrichment tests, saving an RDS object with the complete results, and rendering an HTML report.

Before running, decide:

1. Which database matches your taxa, for example `human-gut`, `mouse-gut`, `marine`, or `gtdb`.
2. Which phenotype you want to test:
   - `prevalence`: whether taxa are present in `which_envir`.
   - `specificity`: whether taxa are specific to `which_envir` compared with other environments.
   - `abundance`: differential abundance across groups or along a continuous variable.
   - `provided`: a phenotype table you calculated elsewhere.
3. Which taxonomic level to run, usually `family` for a first pass.
4. Which metadata columns contain sample IDs, environment/group labels, and dataset labels.

The abundance table should have taxa as rows and samples as columns. The first column should contain taxon IDs matching the selected database. The metadata table should have one row per sample and must include the columns named by `sample_column`, `env_column`, and `dset_column`. If all samples come from one study or batch, set `single_dset=TRUE` and you do not need a real dataset column.

Compared to some R packages, passing options to Phylogenize2 works a little differently under the hood. `phylogenize()` and related functions read global options that can be set with `pz.options()` or overridden directly as extra arguments. For a one-off run, passing options directly to `phylogenize()` is usually clearest. To see every available option, run `?pz.options` in R.

#### Minimal tabular run

This is a good first run for shotgun data that has already been mapped to species IDs in the selected Phylogenize2 database:

```r
library(phylogenize)

results <- phylogenize(
  abundance_file = "data/abundance.tsv",
  metadata_file = "data/metadata.tsv",
  input_format = "tabular",
  db = "human-gut",
  taxon_level = "family",
  which_phenotype = "prevalence",
  which_envir = "case",
  sample_column = "sample",
  env_column = "status",
  single_dset = TRUE,
  out_dir = "output/human_gut_prevalence",
  output_file = "phylogenize-report.html",
  rds_output_file = "core_output.rds",
  ncl = 4
)
```

This writes:

- `output/human_gut_prevalence/phylogenize-report.html`: the interactive report.
- `output/human_gut_prevalence/core_output.rds`: the full result object, which can be reused later.
- `output/human_gut_prevalence/errmsg.txt`: progress messages and warnings, when file logging is enabled.
- Enrichment CSV files such as `enr-table.csv` and `enr-overlaps.csv`, when enrichment results are available.

#### Differential abundance run

Use `which_phenotype="abundance"` when you want gene associations with differential abundance estimates rather than simple presence/prevalence. The `which_envir` value should match the case, treatment, or focal group in your metadata:

```r
results <- phylogenize(
  abundance_file = "data/cirrhosis-abundance.tsv",
  metadata_file = "data/cirrhosis-metadata.tsv",
  input_format = "tabular",
  db = "human-gut",
  taxon_level = "family",
  which_phenotype = "abundance",
  diff_abund_method = "ANCOMBC2",
  which_envir = "case",
  sample_column = "sampleid",
  env_column = "disease_status",
  dset_column = "study",
  out_dir = "output/cirrhosis_abundance_family",
  output_file = "cirrhosis-fam-abd.html",
  rds_output_file = "cirrhosis-fam-abd.rds",
  ncl = 4
)
```

Use `diff_abund_method="Maaslin2"` instead if that method is better suited to your analysis.

#### BIOM input

If your abundance and sample metadata are in one BIOM file, switch `input_format` and provide `biom_file`:

```r
results <- phylogenize(
  biom_file = "data/table-with-metadata.biom",
  input_format = "biom",
  db = "gtdb",
  taxon_level = "family",
  which_phenotype = "prevalence",
  which_envir = "soil",
  sample_column = "sample",
  env_column = "habitat",
  single_dset = TRUE,
  out_dir = "output/biom_run"
)
```

If the BIOM file contains only the abundance matrix and your metadata are in a separate TSV file, also set `separate_metadata=TRUE` and `metadata_file="data/metadata.tsv"`.

#### Re-rendering a report

For long analyses, it can be useful to separate computation from report rendering. `phylogenize_core()` runs the analysis and returns the full result object. `render_core_report()` creates a report from an existing result object.

```r
core <- phylogenize_core(
  abundance_file = "data/abundance.tsv",
  metadata_file = "data/metadata.tsv",
  input_format = "tabular",
  db = "human-gut",
  taxon_level = "family",
  which_phenotype = "prevalence",
  which_envir = "case",
  sample_column = "sample",
  env_column = "status",
  single_dset = TRUE,
  ncl = 4
)

saveRDS(core, "output/human_gut_prevalence/core_output.rds")

core <- readRDS("output/human_gut_prevalence/core_output.rds")
render_core_report(
  core,
  output_file = "phylogenize-report.html",
  out_dir = "output/human_gut_prevalence"
)
```

#### Common options

| Option | Typical value | Description |
|----|----|----|
| `abundance_file` | `"data/abundance.tsv"` | Taxon-by-sample abundance table for tabular input. |
| `metadata_file` | `"data/metadata.tsv"` | Sample metadata table. |
| `biom_file` | `"data/table.biom"` | BIOM file, if using `input_format="biom"`. |
| `input_format` | `"tabular"` | Either `"tabular"` or `"biom"`. |
| `db` | `"human-gut"` | Database to use. Choose one matching your taxon IDs. |
| `data_dir` | package extdata directory | Directory containing `databases.csv` and database files. Set this for custom database locations. |
| `taxon_level` | `"family"` | Taxonomic level to test: `"phylum"`, `"class"`, `"order"`, `"family"`, or `"genus"`. |
| `which_phenotype` | `"prevalence"` | Phenotype to calculate: `"prevalence"`, `"specificity"`, `"abundance"`, or `"provided"`. |
| `which_envir` | `"case"` | Focal environment/group. Must match a value in `env_column`. |
| `sample_column` | `"sample"` | Metadata column containing sample IDs. |
| `env_column` | `"status"` | Metadata column containing environment, group, or numeric phenotype values. |
| `dset_column` | `"study"` | Metadata column containing study/batch labels. |
| `single_dset` | `TRUE` | Use when all samples are from one dataset and no dataset column is needed. |
| `diff_abund_method` | `"ANCOMBC2"` | Differential abundance method for `which_phenotype="abundance"`. |
| `ncl` | `4` | Number of worker processes. Increase when running on a machine or cluster node with more cores. |
| `out_dir` | `"output/my_run"` | Output directory. |
| `output_file` | `"phylogenize-report.html"` | Report file name. |
| `rds_output_file` | `"core_output.rds"` | Saved RDS result file name. Set to `""` to disable. |
| `verbosity` | `1` | Progress-message detail. Increase to `2` or `3` for more diagnostics. |

#### Troubleshooting first runs

- If no samples match between metadata and the abundance matrix, check `sample_column`, abundance table column names, and whitespace in sample IDs.
- If no taxa are retained, confirm that `db` matches the taxon IDs in your abundance table.
- If an environment is not found, check that `which_envir` exactly matches a value in `env_column`.
- If you have only one dataset, set `single_dset=TRUE`.
- If report rendering is slow or memory-intensive, set `skip_graphs=TRUE` for a lighter report.

## Acknowledgements

-   Principal investigator: [Patrick H. Bradley](https://bradleylab.science)
-   Development: Kathryn Kananen, Nia Tran, [Patrick H. Bradley](https://bradleylab.science)
-   Funding:
    -   Startup funds from The Ohio State University
    -   National Institutes of Health, NIGMS R35GM151155

## Contact

If you have questions or comments, please contact [support\@phylogenize.org](mailto:support@phylogenize.org). If Phylogenize2 is giving you an error, please also feel free to file a bug using our [issue tracker](https://bitbucket.org/pbradz/phylogenize/issues?status=new&status=open). Thanks for your feedback!
