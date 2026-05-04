# Phylogenize2 (v2.0.1)

Phylogenize2 allows users to link microbial genes to environments, accounting for phylogeny. More specifically, given community composition data Phylogenize2 links patterns of microbes in a given environment to genes in those microbes' pangenomes, while taking into account an important potential confounder: the phylogenetic relationships between microbes. We allow several different patterns to be calculated, including prevalence, specificity, and differential abundance (using [MaAsLin2](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009442) or [ANCOM-BC2](https://www.nature.com/articles/s41467-020-17041-7)). By default, we use phylogenetic regression, but we also allow users to apply the alternative method [POMS](https://academic.oup.com/bioinformatics/article/38/22/5055/6731923). The method is described in a forthcoming preprint (Kananen et al., in preparation).

In addition, we highly recommend using v2.0.1 or later. We no longer support the use of earlier versions, as there are significant improvements since v0.91.

## Installing Phylogenize2

The easiest way to install all the dependencies needed is by using mamba or conda. We recommend using [miniforge3](https://github.com/conda-forge/miniforge). Please make sure you are using miniforge v3-23.3.1-0 or later. Miniforge3 is available for MacOS, Linux, and Windows OS. Phylogenize is not tested on Windows (proceed with caution); for all future examples, unless otherwise stated, we are assuming you are using Linux.

To install miniforge, run `wget -c https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh` and then in a terminal type `bash Miniforge3-Linux-x86_64.sh`. You'll need to run through the prompts and then give it a download location if the default is not to your liking. Finally, you can let miniforge initialize itself if you want conda to always be in your "base" when you open the terminal. Otherwise, if you say `no` then you will have to manually source the executable for conda which can be done similarly as `source /your/path/to/miniforge3/bin/activate`.

### Now you are ready to start installing the dependencies.

Create a new environment by typing `conda create -n phylogenize` and `conda activate phylogenize`. Then you can install phylogenize by running `conda install bioconda::phylogenize`. For any future analysis, all you have to do is activate this environment to have the dependencies run.

#### Locally - Command line and Rstudio (MacOS/Linux)

Please note, we assume in these instructions you are working off of base-r and NOT Rstudio. We describe at the bottom of this section how to use Rstudio while still installing the dependencies with mamba.

#### Install with mamba/conda - no configuration file

1.  Make sure you have R installed. You can verify if you type `R --version`. If you don't you can get the latest version [here](https://www.r-project.org/) or install it using conda [here](https://anaconda.org/r/r).
2.  Create a new environment in mamba/conda by running `conda create -n phylogenize`
3.  Activate your new environment with `conda activate phylogenize`
4.  Install the dependencies with the bioconda and conda-forge channels as shown below

```         
mamba install -y bioconda::phylogenize
```
If you are running phylogenize2 and plan to use abundance phenotype calculations. The conda version comes with ancombc2 preinstalled. To use maaslin2, you will have to install that separately. 

Additionally, you should install these packages to ensure a smooth workflow for abundance runs:
```
mamba install -c bioconda \
  bioconductor-mia \
  bioconductor-phyloseq \
  bioconductor-microbiom
```

4.  Open R and then run `library("phylogenize")`. You should be all set to run phylogenize!

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

| Environment        | Version | Database | Number of families | Number of species | Archaea Included? | Zenodo                                    |
|--------------------|---------|----------|--------------------|-------------------|-------------------|-------------------------------------------|
| barley rhizosphere | v2.0    | MGnify   | 34                 | 66                |                   |[here](https://zenodo.org/records/19094275)|
| chicken gut        | v1.0.1  | MGnify   | 142                | 1007              |                   |[here](https://zenodo.org/records/18706394)|
| cow rumen          | v1.0.1  | MGnify   | 121                | 1914              |                   |[here](https://zenodo.org/records/19075302)|
| multiple           | v226    | GlobDB   | 10906              | 306261            | True              |[here](https://zenodo.org/records/19607691)|
| honeybee gut       | v1.0.1  | MGnify   | 31                 | 131               |                   |[here](https://zenodo.org/records/18706413)|
| human gut          | v2.0.2  | MGnify   | 215                | 3445              |                   |[here](https://zenodo.org/records/18706603)|
| human oral         | v1.0.1  | MGnify   | 52                 | 260               |                   |[here](https://zenodo.org/records/18706578)|
| human skin         | v1.0    | MGnify   | 86                 | 552               |                   |[here](https://zenodo.org/records/18706597)|
| human vaginal      | v1.0    | MGnify   | 52                 | 189               |                   |[here](https://zenodo.org/records/18706478)|
| maize rhizosphere  | v1.0    | MGnify   | 153                | 268               |                   |[here](https://zenodo.org/records/18706584)|
| marine             | v2.0    | MGnify   | 1192               | 7408              | True              |[here](https://zenodo.org/records/18706676)|
| marine sediment    | v1.0    | MGnify   | 1571               | 4362              | True              |[here](https://zenodo.org/records/19154534)|
| mouse gut          | v1.0    | MGnify   | 136                | 1639              |                   |[here](https://zenodo.org/records/18706664)|
| non model fish gut | v2.0    | MGnify   | 60                 | 87                |                   |[here](https://zenodo.org/records/18706644)|
| pig gut            | v1.0    | MGnify   | 138                | 800               | True              |[here](https://zenodo.org/records/18706629)|
| sheep rumen        | v1.0    | MGnify   | 117                | 2122              |                   |[here](https://zenodo.org/records/18706500)|
| soil               | v1.0    | MGnify   | 1353               | 9122              | True              |[here](https://zenodo.org/records/18707064)|
| tomato rhizosphere | v1.0    | MGnify   | 153                | 268               | True              |[here](https://zenodo.org/records/18706452)|
| zebrafish fecal    | v1.0    | MGnify   | 41                 | 24                |                   |[here](https://zenodo.org/records/18706621)|

### GlobDB v226 Special Note
[GlobDB](https://globdb.org/) is a dereplicated database from multiple sources that are processed by Speth et al, 2025 (1). The project includes 14 genome consolidated resources: GTDB, mOTU, SPIRE, BCRBG, GEM, 13 MGnify Biome Mag catalogs, GOMC, SMAG, TPMC, cFMD, MRGM, HRGM2, sheep and goat gut microbiome compendium, genome catalog of anammox microbiotas, and GFS.

For phylogenize, all databases have been been matched against the UniRef50, FesNov, and UHGP databases, and any remaining protein sequences have been clustered *de novo*. Functional annotations have been obtained using [anvi'o](https://peerj.com/articles/1319/) and [KEGG](https://www.genome.jp/kegg/pathway.html) KOfams as described in Kananen et al., 2025.

Databases can be downloaded manually and decompressed from our Zenodo pages in the table above. All the database files must be placed into a directory called `package/inst/extdata/`. Older database versions can also be located on the Zenodo in the phylogenize community. 

### Making Your Own Database
We recommend using MGnify's v3.0.0 pipeline [here](https://github.com/EBI-Metagenomics/genomes-catalogue-pipeline/releases/tag/v3.0.0) for processing raw files into workable databases. If the files follow standard MGnify format, then they will work in our custom workflow. After you have run their pipeline - a custom database can be generated using our snakemake workflow [here](https://github.com/pbradleylab/phylogenize-db-prep). 

## Preparing your data

If you are using shotgun metagenomes, you will need to first quantify species abundances. The species definitions and names must match the database you plan to use. We recommend using Kraken2 with Bracken, as there are Kraken2 databases for every MGnify database. (Make sure that the version numbers match!) For example:

 - Human gut v2.0.2 Kraken2 database: https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/kraken2_db_uhgg_v2.0.2/
 - Mouse gut Kraken2 database: https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/mouse-gut/v1.0/kraken2_db_mouse-gut_v1.0/
 - Marine Kraken2 database: https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/marine/v2.0/kraken2_db_marine_v2.0/

To use GlobDB, you will need to run taxonomic assignment using Sylph using GlobDB's pregenerated Sylph database found [here](https://fileshare.lisc.univie.ac.at/globdb/globdb_r226/taxonomic_profiling/)

Finally, you will want to make the taxon names from Bracken match the IDs in Phylogenize2. You can check this by seeing if the sampleid column's values match the values in the selected databases `cluster` column in the taxonomy file (i.e mouse-gut-taxonomy.csv). Additionally, you may wish to merge any technical replicates for the same biological sample (as these will lead to overconfident predictions). There is a script to perform this under `shotgun_kraken2_example` called `parse-bracken.R`. You can run this script as follows:

```
Rscript parse_bracken.R -t [path to taxonomy file] -i [path to bracken output files] -o [path to output tab-separated file] -m [path to metadata file]
```

The last option (`-m`) is optional, but allows you to provide a tab-separated file with "sample" and "run" columns that will merge any runs belonging to the same sample. The taxonomy file provided should be the one in the Phylogenize2 database that you are using (e.g. `mouse-gut-taxonomy.csv`). (If you are having trouble finding the path where a database was installed, try looking under the directory where Phylogenize2 was installed, which you should be able to see by running `system.file(package="phylogenize")` in R.)

## Running Phylogenize2

Congratulations! Phylogenize2 should now be installed.

### Running Phylogenize2 using the R interface

The main function in Phylogenize2 is called `phylogenize`. The parameters that you are the most likely to use are as follows:

| Option | Default | Description |
|----|----|----|
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

## Citations
1. Daan R Speth, Nick Pullen, Samuel T N Aroney, Benjamin L Coltman, Jay Osvatic, Ben J Woodcroft, Thomas Rattei, Michael Wagner, GlobDB: a comprehensive species-dereplicated microbial genome resource, Bioinformatics Advances, Volume 5, Issue 1, 2025, vbaf280, https://doi.org/10.1093/bioadv/vbaf280
