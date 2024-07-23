# phylogenize (v0.94 beta)

*phylogenize* is a tool that allows users to link microbial genes to environments, accounting for phylogeny. More specifically, given community composition data, *phylogenize* links genes in microbial genomes to either microbial prevalence in, or specificity for, a given environment, while also taking into account an important potential confounder: the phylogenetic relationships between microbes. *phylogenize* comes with [web](https://www.phylogenize.org), QIIME 2, and R interfaces.

The method is described fully in [Bradley, Nayfach, and Pollard (2018)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006242).

## Installing phylogenize
The easiest way to install all the dependencies needed is by using mamba or conda. We recommend using mamba's maintained [miniforge](https://github.com/conda-forge/miniforge). Miniforge is available for MacOS, Linux, and Windows. For all future examples, unless otherwise stated, we are assuming you are using Linux. 

To install miniforge, run `wget -c https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh` and then `bash Miniforge3-Linux-x86_64.sh`. You'll need to run through the prompts and then give it a download location if the default is not to your liking. Finally, you can let mamba initialize itself if you want mamba to always be in your "base" when you open the terminal. Otherwise, if you say `no` then you will have to manually source the executable for mamba which can be done similarly as `source /your/path/to/miniforge3/bin/activate`.

Now you are ready to start installing the dependencies.

Important::: We have replaced the BURST package with Vsearch. If you are having issues with BURST, please switch over to vsearch! 

### Locally - Command line and Rstudio (MacOS/Linux)
Please note, we assume in these instructions you are working off of base-r and NOT Rstudio. We describe at the bottom of this section how to use Rstudio while still installing the dependencies with mamba.

#### Install with mamba - configuration file
You can make a conda evnironment using the supplied yaml file and not worry about installing any dependencies. Run `mamba env create -f environment.yml
` and then `mamba activate phylogenize`. Open base-r (look below for how to use Rstudio) and then type `devtools::install_bitbucket('pbradz/phylogenize/package/phylogenize')` followed by `library("phylogenize")` and then `phylogenize::install.data.figshare()`.
   
#### Install with mamba - no configuration file
1. Make sure you have R installed. You can verify if you type `R --version`. If you don't you can get the latest version [here](https://www.r-project.org/) or install it using mamba [here](https://anaconda.org/r/r).
   
   *P.S use this website to look for any packages you need to install. Conda is the older version of mamba and the commands are the same. For R the command is like so `mamba install -c r r`*
2. Create a new environment in mamba by running `mamba create -n phylogenize`
3. Activate your new environment with `mamba activate phylogenize`
4. Install the dependencies with the bioconda and conda-forge channels as shown below
```
mamba install -y -c bioconda \
	bioconductor-qvalue \
	bioconductor-ggtree \
	bioconductor-biomformat \
	vsearch
```
```
mamba install -y -c conda-forge \
	r-devtools \
	r-ragg \
	r-phylolm \
	r-phangorn	
```
5. Now you can install phylogenize by running either `R -e "devtools::install_bitbucket('pbradz/phylogenize/package/phylogenize')"` or by opening an R session and then running `devtools::install_bitbucket('pbradz/phylogenize/package/phylogenize')`.
6. Run `library("phylogenize")` in your R session.
7. Then download the necessary databases with `phylogenize::install.data.figshare()`
##### Locally - Rstudio
After creating a `phylogenize` environment in mamba, to install and use Rstudio run `mamba install -c r rstudio`. Then you can activate it by typing `rstudio`. This will launch an Rstudio IDE. There, if you haven't already, you can run `devtools::install_bitbucket('pbradz/phylogenize/package/phylogenize')` followed by `library("phylogenize")` and then `phylogenize::install.data.figshare()`.

### QIIME 2
QIIME 2 runs in a conda environment, meaning it has its own installation of R and related packages. To run phylogenize with QIIME 2, you will need to install phylogenize within the QIIME 2 conda environment, then install the [q2-phylogenize plugin](https://bitbucket.org/pbradz/q2-phylogenize). The instructions are similar for installing in any other conda environment, except you won't need the plugin.

1. First, switch to the correct environment. For QIIME 2, this is accomplished with `source activate qiime2-2019.4` (see [here](https://docs.qiime2.org/2019.4/install/native/#activate-the-conda-environment)). Note that you may need to replace "2019.4" with the most recent version of QIIME2, e.g., "2020.6".

```
options(unzip="internal")
Sys.setenv(TAR="/bin/tar")   # replace with path to tar in your installation, if necessary
```

These options only stick around for your current session, so re-run them if you quit and come back to the installation process. You shouldn't need to run them once phylogenize is installed, though.

Next try:

```
ackages("BiocManager")
BiocManager::install(c("qvalue","biomformat","ggtree"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
devtools::install_bitbucket("pbradz/phylogenize/package/phylogenize")
phylogenize::install.data.figshare()    # nb: this may take a while to download
```

To avoid causing problems elsewhere in QIIME2, we recommend *not* updating any other packages if prompted.

*Note* (updated 8/4/2020): If you get an error message about the package "mnormt" and/or "phylolm" not being available for R 3.5.1, try installing a specific version with the command `remotes::install_version("mnormt", "1.5-5"); remotes::install_version("phylolm", "2.6")`. Then retry the above commands starting from the `devtools::install_bitbucket` line.

*Note* (updated 8/10/2020): If you get errors installing the package "igraph" that relate to "-lgomp", first install libgomp from the command line (outside of R, but in your QIIME2 environment) with `conda install libgomp`. Because of some quirks in where conda expects libraries to be, you may then also need to manually put the `libgomp` library in your QIIME2 environment. You can do that with a soft-link, e.g., `ln -s ~/miniconda3/lib/libgomp.so ~/miniconda3/envs/qiime2-2020.6/x86_64-conda-linux-gnu/lib/`. Substitute your anaconda directory for "~/miniconda3" and your QIIME2 environment name for "qiime2-2020.6" as needed.

*Note*: If you are having trouble with the `install_bitbucket` command on Windows or within the QIIME2 VM on Windows, you can try instead [downloading](https://bitbucket.org/pbradz/phylogenize/downloads/) the repository and unzipping it, then running:

```
devtools::install(pkg="[...]/package/phylogenize")
```
... where `[...]` should be replaced with the path where you unzipped the repository. Then run `install.data.figshare()` as above.

Finally, install the plugin. From the UNIX command line (i.e., not in R -- you can use the command `quit()` or press Ctrl+D to leave the R environment; don't worry about "saving" since the packages will remain installed):

```
git clone https://bitbucket.org/pbradz/q2-phylogenize
cd q2-phylogenize
python setup.py build
python setup.py install
```

Further information about how to use q2-phylogenize can be found on its [git repository](https://bitbucket.org/pbradz/q2-phylogenize). If you are just using conda and not QIIME 2, you can just proceed to the section entitled "Running *phylogenize* using the R interface."

### Installing *phylogenize* package for use on AWS

We recommend you install *phylogenize* in a conda environment as above. However, because the default Amazon images are meant for headless operation they are missing some tools to deal with fonts that *phylogenize* uses to generate its plots. You can install those as follows:

```
conda install -c conda-forge xorg-libxt
sudo apt install zlib
sudo apt show zlib1g
sudo apt install fontconfig
```

### Re-installing *phylogenize* after an upgrade

To **reinstall** *phylogenize*, the following should be all that's necessary:

```
devtools::install_bitbucket("pbradz/phylogenize/package/phylogenize")
phylogenize::install.data.figshare()
```

## Running *phylogenize*

Congratulations! *phylogenize* should now be installed.

### Running *phylogenize* using the R interface

The main function in *phylogenize* is `render.report`. The parameters that you are the most likely to use are as follows:

| Option            | Default         |   Description  |
|-------------------|-----------------|----------------|
| in_dir | "." | String. Path to input directory (i.e., where to look for input files. |
| out_dir | "output" | String. Path to output directory. |
| abundance_file | "test-abundance.tab" | String. Name of abundance tabular file. |
| metadata_file | "test-metadata.tab" | String. Name of metadata tabular file. |
| biom_file | "test.biom" | String. Name of BIOM abundance-and-metadata file. |
| input_format | "tabular" | String. Whether to look for tabular or BIOM-formatted data ("tabular" or "biom"). |
| ncl | 1 | Integer. Number of cores to use for parallel computation. |
| type | "midas" |  Type of data to use, either "midas" (shotgun) or "16S" (amplicon). |
| env_column | "env" | String. Name of column in metadata file containing the environment annotations. |
| dset_column | "dataset" | String. Name of column in metadata file containing the dataset annotations. |
| sample_column | "sample_id" | Name of column in metadata file containing the sample IDs. |
| single_dset | FALSE | Boolean. If true, will assume that all samples come from a single dataset called `"dset1"` no matter what, if anything, is in `dset_column`. |
| db_version | "midas_v1.2" | String. Which version of the MIDAS database to use ("midas_v1.2" or "midas_v1.0"). |
| which_phenotype | "prevalence" | String. Which phenotype to calculate ("prevalence" or "specificity"). |
| which_envir | "Stool" | String. Environment in which to calculate prevalence or specificity. Must match annotations in metadata. |

Compared to some R packages, passing options to *phylogenize* works a little differently under the hood. Instead of having its own parameters, `render.report` and other *phylogenize* functions look for global options that can either be set using the function `pz.options` or overridden as extra arguments. This allows you to set parameters once and then work with the *phylogenize* functions without retyping them, and therefore makes the code easier to read. To see the full list of parameters that can be overridden, see `?pz.options`.

Here is an example invocation of `render.report`: 

~~~~
library(phylogenize)
render.report(
    output_file="16S-results.html",
    in_dir="/home/kananen13/workspace/bradleyLab/tools/phylogenize/phylogenize/hmp/",
    out_dir=file.path("hmp", "16S-results"),
    type="midas",
    db_version="midas_v1.0",
    which_phenotype="prevalence",
    which_envir="Stool",
    abundance_file="hmp-shotgun-bodysite.tab",
    metadata_file="hmp-shotgun-bodysite-metadata.tab", 
    input_format="tabular",
    ncl=10)
~~~~

This invocation will generate a report under "./hmp/16S-results" called "16S-results.html".

Note that for now it is necessary to call `phylogenize::set_data_internal()` if you don't explicitly load the package with `library(phylogenize)`; that function is automatically triggered when the package is loaded with `library()`.


### Running *phylogenize* locally with the web interface

*phylogenize* also can be used with its own graphical user interface and job manager, which can be accessed with a web browser. We provide an installation of *phylogenize* at [https://www.phylogenize.org], but you can also run this interface locally. (The next part of the guide assumes a \*-nix environment like Ubuntu or OS X, but you may also be able to run this on Windows using the Windows Subsystem for Linux or Cygwin.)

The *phylogenize* web interface is written as a Flask WSGI app. If you are just going to be running it on your own computer or within a trusted intranet, you can probably use the built-in Flask server. (If you are concerned about security and will be allowing potentially untrusted users to use the server, or if the built-in Flask server is inadequate for any other reason, we recommend using a production web server like Apache2. The Flask documentation has [some information](http://flask.pocoo.org/docs/1.0/deploying/mod_wsgi/) on how to get started hosting a WSGI app on Apache2 using `mod_wsgi`.)

We recommend installing Flask using a virtual environment as per the instructions [here](http://flask.pocoo.org/docs/1.0/installation/). Once you have activated the virtual environment, before running the app but after cloning the repository with `git clone`, you will need to launch the Beanstalk-based queueing system. From the repository root, run:

    nohup beanstalkd -l 127.0.0.1 -p 14711 &
    nohup python3 phylogenize_app/worker.py &

To keep these jobs running if the terminal is closed, you may want to `disown` these jobs or alternatively, start them in a `tmux` or `screen` session.

To allow multiple simultaneous jobs, you will need to edit `worker.py` and change MaxJobs as appropriate. It may be worth starting with a single job, particularly if each job runs multi-threaded, to make sure memory use stays reasonable.

Next, you will need to start the server as follows:

     FLASK_APP=phylogenize_app flask run

The application should then be accessible from http://localhost:5000.

## Acknowledgements

 * Project lead and repository maintainer: [Patrick H. Bradley](http://docpollard.org/people/patrick-j-h-bradley/)
 * Principal investigator: [Katherine S. Pollard](http://docpollard.org)
 * Testing, troubleshooting, and debugging AWS/Mac installation: [Chunyu Zhao](https://github.com/zhaoc1)
 * More testing: [Jordan Bisanz](https://github.com/jbisanz) and other members of the [Turnbaugh lab](https://turnbaughlab.ucsf.edu/)
 * Funding: 
     - National Science Foundation (DMS-106[]9303, DMS-156[]3159)
     - Gordon & Betty Moore Foundation (#3300)

## Contact

If you have questions or comments, please contact [support@phylogenize.org]. If *phylogenize* is giving you an error, please also feel free to file a bug using our [issue tracker](https://bitbucket.org/pbradz/phylogenize/issues?status=new&status=open). Thanks for your feedback!
 
