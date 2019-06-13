[TOC]

# phylogenize (v0.92 beta)

*phylogenize* is a tool that allows users to link microbial genes to environments, accounting for phylogeny. More specifically, given community composition data, *phylogenize* links genes in microbial genomes to either microbial prevalence in, or specificity for, a given environment, while also taking into account an important potential confounder: the phylogenetic relationships between microbes. *phylogenize* comes with [web](https://www.phylogenize.org), QIIME2, and R interfaces.

The method is described fully in [Bradley, Nayfach, and Pollard (2018)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006242).

## Installing the *phylogenize* package and its dependencies

The "core" of *phylogenize* is an R package, so to run *phylogenize* locally, you will need to first install that package from this repository. There are separate instructions depending on whether you are installing *phylogenize* for use 1. with R or the web interface on a local machine, 2. with QIIME2 or in a conda environment, or 3. on AWS.

### Installing BURST

BURST is needed for all 16S analyses using *phylogenize*. BURST is a high-speed pairwise aligner that *phylogenize* uses to map 16S amplicon sequence variants back to a database of genomes.

You can download the binaries from the [BURST Github repository](github.com/knights-lab/BURST). By default *phylogenize* expects these binaries to be copied into the directory `/usr/local/bin`, but you can override this (see below).

**Note for Mac**: Binaries for the latest version of BURST are available only on request, but *phylogenize* should work with an [earlier version of BURST](https://github.com/knights-lab/BURST/releases/tag/v0.99.4a) as long as the binary is renamed `burst12` and copied to `/usr/local/bin` or your preferred path (again, see below for how to override this).

### Installing *phylogenize* package for use with R or the web interface on a local machine

Because *phylogenize* is an analysis and visualization pipeline, it has more dependencies than the average package, so this will require a little patience. `devtools` is capable of tracking down most of these dependencies automatically, but there are a few packages from BioConductor that it will not be able to find. It's a good idea to install those first:

~~~~
install.packages("BiocManager")
BiocManager::install(c("qvalue","biomformat","ggtree"))
install.packages("devtools")
devtools::install_bitbucket("pbradz/phylogenize/package/phylogenize")
~~~~

(Note that you need to tell R to look in a specific subdirectory of this repository -- i.e., `package/phylogenize` -- and not the root.)

***Important note***: If you are using R >= 3.6.0, you will need to install the latest development versions of `phytools` and `treeio`, which have an important [bug fix](https://github.com/liamrevell/phytools/issues/47) (thanks to Liam Revell and Guangchuang Yu for fixing this so quickly):

```
devtools::install_github("liamrevell/phytools")
devtools::install_github("GuangchuangYu/treeio")
```

Finally, you will need to download and install the data files that *phylogenize* needs to run. That has been automated so that all you should need to do is run the following function:

```
library(phylogenize)
phylogenize::install.data.figshare()
```


### Installing *phylogenize* package for use with QIIME2 or in another conda environment

QIIME2 runs in a conda environment, meaning it has its own installation of R and related packages. To run *phylogenize* with QIIME2, you will need to install *phylogenize* within the QIIME2 conda environment, then install the [q2-phylogenize plugin](https://bitbucket.org/pbradz/q2-phylogenize). The instructions are similar for installing in any other conda environment, except of course you won't need the plugin.

First, switch to the correct environment. For QIIME2, this is accomplished with `source activate qiime2-2019.4` (see [here](https://docs.qiime2.org/2019.4/install/native/#activate-the-conda-environment)).

Installing *phylogenize* within conda is a little tricky. You will need to manually install a few libraries and packages that are either not included, or difficult to install from source. From the UNIX command line:

```
conda install libcurl
conda install r-devtools
conda install -c bioconda bioconductor-rhdf5lib
conda install -c conda-forge r-magick
conda install -c r r-git2r
conda install -c r r-shiny
```

Next, run R within the same environment and install the *phylogenize* library. However, you will probably need to work around a [known issue in conda](https://github.com/r-lib/devtools/issues/1722) before calling `devtools::install_bitbucket`. The following should work (from within R):

```
options(unzip="internal")
Sys.setenv(TAR="/bin/tar")   # replace with path to tar in your installation, if necessary
install.packages("BiocManager")
BiocManager::install(c("qvalue","biomformat","ggtree"))
devtools::install_bitbucket("pbradz/phylogenize/package/phylogenize")
phylogenize::install.data.figshare()    # nb: this may take a while to download
```

Finally, install the plugin. From the UNIX command line (i.e., not in R):

```
git clone bitbucket.org/pbradz/q2-phylogenize
cd q2-phylogenize
python setup.py build
python setup.py install
```

Further information about how to use q2-phylogenize can be found on its [git repository](https://bitbucket.org/pbradz/q2-phylogenize). If you are just using conda and not QIIME2, you can just proceed to the section entitled "Running *phylogenize* using the R interface."

### Installing *phylogenize* package for use on AWS

We recommend you install *phylogenize* in a conda environment as above. However, because the default Amazon images are meant for headless operation they are missing some tools to deal with fonts that *phylogenize* uses to generate its plots. You can install those as follows:

```
conda install -c conda-forge xorg-libxt
sudo apt install zlib
sudo apt show zlib1g
sudo apt install fontconfig
```

### Special instructions for Mac users

If you see the error `Fontconfig error: unable to match font pattern`, you may need to install or reinstall some packages using [Homebrew](https://brew.sh/). From the command line, after installing Homebrew:

```
brew update
brew upgrade
brew reinstall cairo
```

Then, open R and in a fresh session, run:

```
devtools::install_github("davidgohel/gdtools") 
devtools::install_github("davidgohel/rvg")
devtools::install_github("davidgohel/ggiraph")
```

### Re-installing *phylogenize* after an upgrade

To **reinstall** *phylogenize*, the following should be all that's necessary:

```
devtools::install_bitbucket("pbradz/phylogenize/package/phylogenize")
phylogenize::install.data.figshare()
```

## Running *phylogenize*

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
| burst_dir | "/usr/local/bin" | String. Path to the BURST binaries. |

Compared to some R packages, passing options to *phylogenize* works a little differently under the hood. Instead of having its own parameters, `render.report` and other *phylogenize* functions look for global options that can either be set using the function `pz.options` or overridden as extra arguments. This allows you to set parameters once and then work with the *phylogenize* functions without retyping them, and therefore makes the code easier to read. To see the full list of parameters that can be overridden, see `?pz.options`.

Here is an example invocation of `render.report`: 

~~~~
library(phylogenize)
render.report(
    output_file="16S-results.html",
    in_dir=hmp_dir,
    out_dir=file.path("hmp", "16S-results"),
    type="16S",
    db_version="midas_v1.2",
    which_phenotype="prevalence",
    which_envir="Stool",
    abundance_file="hmp-16s-dada2-full.tab",
    metadata_file="hmp-16s-phylogenize-metadata-full.tab",
    input_format="tabular",
    burst_dir="/home/pbradz/bin/",
    ncl=10)
~~~~

This invocation will generate a report under "./hmp/16S-results" called "16S-results.html".

Note that for now it is necessary to call `phylogenize::set_data_internal()` if you don't explicitly load the package with `library(phylogenize)`; that function is automatically triggered when the package is loaded with `library()`. If you get an error `ERROR: invalid input file.` after BURST is invoked, this is probably the reason why.


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
 * Testing and troubleshooting: [Chunyu Zhao](https://github.com/zhaoc1)
 * More testing: [Jordan Bisanz](https://github.com/jbisanz) and other members of the [Turnbaugh lab](https://turnbaughlab.ucsf.edu/)
 * Funding: 
   * National Science Foundation [DMS-1069303, DMS-1563159]
   * Gordon & Betty Moore Foundation [#3300]
   
## Contact

If you have questions or comments, please contact [support@phylogenize.org]. If *phylogenize* is giving you an error, please also feel free to file a bug using our [issue tracker](https://bitbucket.org/pbradz/phylogenize/issues?status=new&status=open). Thanks for your feedback!
 
