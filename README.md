# Phylogenize2 (v2.0.0-alpha)

Phylogenize2 allows users to link microbial genes to environments, accounting for phylogeny. More specifically, given community composition data Phylogenize2 links patterns of microbes in a given environment to genes in those microbes' pangenomes, while taking into account an important potential confounder: the phylogenetic relationships between microbes. We allow several different patterns to be calculated, including prevalence, specificity, and differential abundance (using [MaAsLin2](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009442) or [ANCOM-BC2](https://www.nature.com/articles/s41467-020-17041-7)). By default, we use phylogenetic regression, but we also allow users to apply the alternative method [POMS](https://academic.oup.com/bioinformatics/article/38/22/5055/6731923). The method is described in a forthcoming preprint (Kananen et al., in preparation).

In addition, we highly recommend using v2.0.0-alpha or later. We no longer support the use of earlier versions, as there are significant improvements since v0.91.

## Installing Phylogenize2

The easiest way to install all the dependencies needed is by using mamba or conda. We recommend using [miniforge3](https://github.com/conda-forge/miniforge). Please make sure you are using miniforge v3-23.3.1-0 or later. Miniforge3 is available for MacOS, Linux, and Windows OS. Phylogenize is not tested on Windows (proceed with caution); for all future examples, unless otherwise stated, we are assuming you are using Linux.

To install miniforge, run `wget -c https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh` and then in a terminal type `bash Miniforge3-Linux-x86_64.sh`. You'll need to run through the prompts and then give it a download location if the default is not to your liking. Finally, you can let miniforge initialize itself if you want conda to always be in your "base" when you open the terminal. Otherwise, if you say `no` then you will have to manually source the executable for conda which can be done similarly as `source /your/path/to/miniforge3/bin/activate`.

You can also install Phylogenize2 using [Pixi](https://pixi.sh) with the provided `pixi.toml` file.

### Now you are ready to start installing the dependencies.

Create a new environment by typing `conda create -n phylogenize` and `conda activate phylogenize`. Then you can install phylogenize by running `conda install bioconda::phylogenize`. For any future analysis, all you have to do is activate this environment to have the dependencies run.

To use Pixi, first use `git clone https://github.com/pbradleylab/phylogenize`, enter the `phylogenize` directory, then type `pixi install` to download all of the dependencies. You can then use Phylogenize2 from within an R session that you start by typing `pixi run R` from within the `phylogenize` directory.

#### Locally - Command line and Rstudio (MacOS/Linux)

Please note, we assume in these instructions you are working off of command line R and NOT Rstudio. We describe at the bottom of this section how to use Rstudio while still installing the dependencies with mamba.

#### Install with mamba - configuration file

You can make a conda environment using the supplied yaml file and not worry about installing any dependencies. Run `conda env create -f environment.yml` and then `conda activate phylogenize`. Open R from the terminal, and then type `devtools::install_github("biocore/phylogenize")`.

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

Currently, two databases can be used with Phylogenize2:

| Name | Environment        | Version | Database | Number of families | Number of species |
|------|---------------|---------------|---------------|---------------|---------------|
| uhgp | human gut          | v1.0  | MGnify   | 202                | 4542              |
| gtdb | mixed environment  | v202    | GTDB     | 3003               | 43058             |

These databases can be downloaded from our Zenodo page [here](https://zenodo.org/communities/bradley_phylogenize), then installed using Phylogenize2's `phylogenize::install_data("path/to/downloaded/file")`. Note that you will also have to download and install the `databases.csv` file in addition to one or both databases.

The default if no database is available is GTDB. If using a custom database, then all the database files must be placed into a directory called `package/inst/extdata/`.

We are currently preparing more biome-specific databases that may be better suited for particular environments (in progress):

| Environment        | Version | Database | Number of families | Number of species |
|---------------|---------------|---------------|---------------|---------------|
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
| mixed environment  | v202    | GTDB     | 3003               | 43058             |

All of the above databases will have been been matched against the UniRef50, FesNov, and UHGP databases, and any remaining protein sequences clustered *de novo*. Functional annotations have been obtained using [anvi'o](https://peerj.com/articles/1319/) and [KEGG](https://www.genome.jp/kegg/pathway.html) KOfams as described in Kananen et al., 2025.

## Preparing your data

If you are using shotgun metagenomes, you will need to first quantify species. The species definitions and names must match the database you plan to use. We recommend using Kraken2/Bracken with one of the following databases:

 - UHGG v1.0 Kraken2 database: https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgg_kraken2-db/
 - GTDB v202 Kraken2 database: http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release202/kraken2/
    * Thanks to Nick Youngblut who generated this database using [Struo2](https://github.com/leylabmpi/Struo2).
  
An example workflow for UHGG written in Snakemake can be seen under `shotgun_kraken2_example`. (Note that the names in the UHGG v1.0 database does not exactly match the database used in Phylogenize2, so they get processed further in `read-bracken.R`. We will make this easier in a future release.)

Finally, we also recommend that you merge any technical replicates at this point, as leaving in multiple measurements per experimental unit will lead to overconfident predictions. `read-bracken.R` has an example of how to do this using run info downloaded from the SRA (provided as an example).

## Running Phylogenize2

Congratulations! Phylogenize2 should now be installed.

### Running Phylogenize2 using the R interface

The main function in Phylogenize2 is called `phylogenize`. The parameters that you are the most likely to use are as follows:

| Option | Default | Description |
|--------------------------|------------------------|----------------------|
| in_dir | "." | String. Path to input directory (i.e., where to look for input files. |
| out_dir | "output" | String. Path to output directory. |
| abundance_file | "test-abundance.tab" | String. Name of abundance tabular file. |
| metadata_file | "test-metadata.tab" | String. Name of metadata tabular file. |
| biom_file | "test.biom" | String. Name of BIOM abundance-and-metadata file, if using BIOM instead of tabular data. |
| input_format | "tabular" | String. Whether to look for tabular or BIOM-formatted data ("tabular" or "biom"). |
| ncl | 1 | Integer. Number of cores to use for parallel computation. |
| type_16S | FALSE | Boolean. Set to true if your species names are 16S ASV sequences, instead of species IDs from your database of interest. |
| db | "uhgp" | String. Gives the database to use. Some options are "uhgp" and "gtdb"; see above for others. |
| env_column | "env" | String. Name of column in metadata file containing the environment annotations. |
| dset_column | "dataset" | String. Name of column in metadata file containing the dataset annotations. |
| sample_column | "sample_id" | Name of column in metadata file containing the sample IDs. |
| single_dset | FALSE | Boolean. If true, will assume that all samples come from a single dataset called `"dset1"` no matter what, if anything, is in `dset_column`. |
| diff_abund_method | "maaslin2" | String. Which tool to use to give differential abundance estimates ("Maaslin2" or "ANCOMBC2"; case insensitive). |
| which_phenotype | "prevalence" | String. Which phenotype to calculate ("prevalence", "abundance", "specificity", or "provided"). |
| taxon_level | "family" | String. Run analyses for each of these taxonomic units (can be "phylum", "class", "order", "family", or "genus"; "family" is recommended). |
| which_envir | "Stool" | String. Environment in which to calculate prevalence or specificity. Must match annotations in metadata. |

Compared to some R packages, passing options to Phylogenize2 works a little differently under the hood. Instead of having its own parameters, `phylogenize` and other Phylogenize2 functions look for global options that can either be set using the function `pz.options` or overridden as extra arguments. This allows you to set parameters once and then work with the Phylogenize2 functions without retyping them, and therefore makes the code easier to read. To see the full list of parameters that can be overridden, see `?pz.options`.

Here is an example invocation:

```         
library(phylogenize)
cirrhosis_family_abundance <- phylogenize(
  output_file="cirrhosis-fam-abd.html",
  output_rds_file="cirrhosis-fam-abd.rds",
  out_dir=file.path("output", "cirrhosis_uhgp_abd_family"),
  db="uhgp",
  taxon_level="family",
  type_16S=FALSE,
  which_phenotype="abundance",
  diff_abund_method="maaslin2",
  which_envir="case",
  abundance_file="test_data/cirr/cirrhosis-abundance.tab",
  metadata_file="test_data/cirr/cirrhosis-metadata.tab", 
  input_format="tabular",
  sample_column="sampleid",
  ncl=4)
```

This invocation will run Phylogenize2 with four cores, using Maaslin2 to get differential abundance of microbes between cases and controls, and using the UHGP human gut database. It will then output the report to `output/cirrhosis_uhgp_abd_family/cirrhosis-fam-abd.html` and will also generate a so-called RDS object under `output/cirrhosis_uhgp_abd_family/cirrhosis-fam-abd.rds` that contains the full output generated by Phylogenize2, so that you can later re-generate just the report if desired.

You can also run just the analysis part of Phylogenize2 using the function `phylogenize_core()`, or just render a new report from an existing analysis run of Phylogenize2 using `render_core_report()`. (Note that `phylogenize_core()` does not save a RDS file of its results by default, but you can save it with `saveRDS`.) The above call would be equivalent to:

```         
cirrhosis_family_abundance <- phylogenize_core(
  db="uhgp",
  taxon_level="family",
  type_16S=FALSE,
  which_phenotype="abundance",
  diff_abund_method="maaslin2",
  which_envir="case",
  abundance_file="test_data/cirr/cirrhosis-abundance.tab",
  metadata_file="test_data/cirr/cirrhosis-metadata.tab", 
  input_format="tabular",
  sample_column="sampleid",
  ncl=4)
  
saveRDS(cirrhosis_family_abundance, 
  output_rds_file="output/cirrhosis_uhgp_abd_family/cirrhosis-fam-abd.rds")
  
# To load this output back into memory after writing to disk:
# cirrhosis_family_abundance <- readRDS("output/cirrhosis_uhgp_abd_family/cirrhosis-fam-abd.rds")

render_core_report(
  cirrhosis_family_abundance,
  output_file="cirrhosis-fam-abd.html",
  out_dir=file.path("output", "cirrhosis_uhgp_abd_family"))
```

## Acknowledgements

-   Principal investigator: [Patrick H. Bradley](https://bradleylab.science)
-   Development: Kathryn Kananen, Nia Tran, [Patrick H. Bradley](https://bradleylab.science)
-   Funding:
    -   Startup funds from The Ohio State University
    -   National Institutes of Health, NIGMS R35GM151155

## Contact

If you have questions or comments, please contact [support\@phylogenize.org](mailto:support@phylogenize.org). If Phylogenize2 is giving you an error, please also feel free to file a bug using our [issue tracker](https://bitbucket.org/pbradz/phylogenize/issues?status=new&status=open). Thanks for your feedback!
