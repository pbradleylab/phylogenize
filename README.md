# Phylogenize2 (v2.0.0-alpha)

Phylogenize2 allows users to link microbial genes to environments, accounting for phylogeny. More specifically, given community composition data Phylogenize2 links patterns of microbes in a given environment to genes in those microbes' pangenomes, while taking into account an important potential confounder: the phylogenetic relationships between microbes. We allow several different patterns to be calculated, including prevalence, specificity, and differential abundance (using [MaAsLin2](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009442) or [ANCOM-BC2](https://www.nature.com/articles/s41467-020-17041-7)). By default, we use phylogenetic regression, but we also allow users to apply the alternative method [POMS](https://academic.oup.com/bioinformatics/article/38/22/5055/6731923). The method is described in a forthcoming preprint (Kananen et al., in preparation).

In addition, we highly recommend using v2.0.0-alpha or later. We no longer support the use of earlier versions, as there are significant improvements since v0.91.

## Installing Phylogenize2

The easiest way to install all the dependencies needed is by using mamba or conda. We recommend using [miniforge3](https://github.com/conda-forge/miniforge). Please make sure you are using miniforge v3-23.3.1-0 or later. Miniforge3 is available for MacOS, Linux, and Windows OS. Phylogenize is not tested for an Windows OS, so proceed at your own paril. For all future examples, unless otherwise stated, we are assuming you are using Linux.

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
|---------------|--------------|--------------|---------------|--------------|
| chicken gut        | v1.0.1  | MGnify   | 142                | 1007              |
| cow rumen          | v1.0.1  | MGnify   | 121                | 1914              |
| honeybee gut       | v1.0.1  | MGnify   | 31                 | 131               |
| human gut          | v2.0.2  | MGnify   | 203                | 4543              |
| human oral         | v1.0.1  | MGnify   | 52                 | 260               |
| human vaginal      | v1.0    | MGnify   | 52                 | 189               |
| marine eukaryotes  | vbeta   | MGnify   | 250                | 250               |
| marine             | v2.0    | MGnify   | 1192               | 7408              |
| mouse gut          | v1.0    | MGnify   | 136                | 1639              |
| non model fish gut | v2.0    | MGnify   | 60                 | 87                |
| pig gut            | v1.0    | MGnify   | 138                | 800               |
| sheep rumen        | v1.0    | MGnify   | 117                | 2122              |
| zebrafish fecal    | v1.0    | MGnify   | 41                 | 24                |
| mixed environment  | v214    | GTDB     | 3003               | 43058             |

All databases have been been matched against the UniRef50, FesNov, and UHGP databases, and any remaining protein sequences have been clustered *de novo*. Functional annotations have been obtained using [anvi'o](https://peerj.com/articles/1319/) and [KEGG](https://www.genome.jp/kegg/pathway.html) KOfams as described in Kananen et al., 2025.

Databases can be downloaded manually and decompressed from our Zenodo page [here](), or they can be downloaded and decompressed using Phylogenize2's `phylogenize::download_zenodo_db("your/html/link/here.zip")`

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
| db | "uhgp" | String. Gives the database to use. Some options are "uhgp" and "gtdb"; see below. |
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

This invocation will run Phylogenize2 with four cores, using Maaslin2 to get differential abundance of microbes between cases and controls, and using the UHGP human gut database. It will then output the report to `output/cirrhosis_uhgp_abd_family/cirrhosis-fam-abd.html` and will also generate an RDS object under `output/cirrhosis_uhgp_abd_family/cirrhosis-fam-abd.rds`, so that you could re-generate just the report if desired.

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
render_core_report(
  cirrhosis_family_abundance,
  output_file="cirrhosis-fam-abd.html",
  out_dir=file.path("output", "cirrhosis_uhgp_abd_family"))
```

### Running Phylogenize2 locally with the web interface

Phylogenize2 also can be used with its own graphical user interface and job manager, which can be accessed with a web browser. We provide an installation of *Phylogenize2* at [<https://www.phylogenize.org>], but you can also run this interface locally. (The next part of the guide assumes a \*-nix environment like Ubuntu or OS X, but you may also be able to run this on Windows using the Windows Subsystem for Linux or Cygwin.)

The Phylogenize2 web interface is written as a Flask WSGI app. If you are just going to be running it on your own computer or within a trusted intranet, you can probably use the built-in Flask server. (If you are concerned about security and will be allowing potentially untrusted users to use the server, or if the built-in Flask server is inadequate for any other reason, we recommend using a production web server like Apache2. The Flask documentation has [some information](http://flask.pocoo.org/docs/1.0/deploying/mod_wsgi/) on how to get started hosting a WSGI app on Apache2 using `mod_wsgi`.)

We recommend installing Flask using a virtual environment as per the instructions [here](http://flask.pocoo.org/docs/1.0/installation/). Once you have activated the virtual environment, before running the app but after cloning the repository with `git clone`, you will need to launch the Beanstalk-based queueing system. From the repository root, run:

```         
nohup beanstalkd -l 127.0.0.1 -p 14711 &
nohup python3 phylogenize_app/worker.py &
```

To keep these jobs running if the terminal is closed, you may want to `disown` these jobs or alternatively, start them in a `tmux` or `screen` session.

To allow multiple simultaneous jobs, you will need to edit `worker.py` and change MaxJobs as appropriate. It may be worth starting with a single job, particularly if each job runs multi-threaded, to make sure memory use stays reasonable.

Next, you will need to start the server as follows:

```         
 FLASK_APP=phylogenize_app flask run
```

The application should then be accessible from <http://localhost:5000>.

## Acknowledgements

-   Principal investigator: [Patrick H. Bradley](https://bradleylab.science)
-   Development: Kathryn Kananen, Nia Tran, [Patrick H. Bradley](https://bradleylab.science)
-   Funding:
    -   Startup funds from The Ohio State University
    -   National Institutes of Health, NIGMS R35GM151155

## Contact

If you have questions or comments, please contact [[support\@phylogenize.org](mailto:support@phylogenize.org){.email}]. If Phylogenize2 is giving you an error, please also feel free to file a bug using our [issue tracker](https://bitbucket.org/pbradz/phylogenize/issues?status=new&status=open). Thanks for your feedback!
