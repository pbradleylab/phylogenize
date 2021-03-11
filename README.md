[TOC]

# phylogenize (v0.94 beta)

*phylogenize* is a tool that allows users to link microbial genes to environments, accounting for phylogeny. More specifically, given community composition data, *phylogenize* links genes in microbial genomes to either microbial prevalence in, or specificity for, a given environment, while also taking into account an important potential confounder: the phylogenetic relationships between microbes. *phylogenize* comes with [web](https://www.phylogenize.org), QIIME 2, and R interfaces.

The method is described fully in [Bradley, Nayfach, and Pollard (2018)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006242).

## Installing the *phylogenize* package and its dependencies

The "core" of *phylogenize* is an R package, so to run *phylogenize* locally, you will need to first install that package from this repository. There are separate instructions depending on whether you are installing *phylogenize* for use 1. with R or the web interface on a local machine, 2. with QIIME 2 or in a conda environment, or 3. on AWS.

### Installing BURST and/or vsearch

(Update: 5/4/2020) Either BURST or, as of phylogenize 0.94, vsearch is needed for all 16S analyses using *phylogenize*. BURST is a high-speed pairwise aligner that *phylogenize* uses to map 16S amplicon sequence variants back to a database of genomes.

You can download the binaries from the [BURST Github repository](github.com/knights-lab/BURST) or the [vsearch repository](https://github.com/torognes/vsearch). By default *phylogenize* expects these binaries to be copied into the directory `/usr/local/bin`, but you can override this (see below).

**Note for Mac** (Update: 5/4/2020): On a Mac, using vsearch may be more straightforward than using BURST. Mac binaries for the latest version of BURST are available only on request, but *phylogenize* should work with an [earlier version of BURST](https://github.com/knights-lab/BURST/releases/tag/v0.99.4a) as long as the binary is renamed `burst12` and copied to `/usr/local/bin` or your preferred path (again, see below for how to override this). If that version of BURST still doesn't work, follow the instructions for older computers.

**Note for older computers** (Update: 5/4/2020): Some of the optimizations that BURST uses are only available on newer architectures. You may be able to run an older version of BURST called [EMBALMER](https://github.com/knights-lab/BURST/releases): pick a version tagged as "buzzard," which means it is compiled to be slower but compatible with more machines. If you use such a version of BURST, instead of renaming it, specify the original name by setting `burst_bin` when calling *phylogenize* (see below). This way, *phylogenize* won't use command-line options that only work in later BURST releases.

### Installing *phylogenize* package for use with R or the web interface on a local machine

Because *phylogenize* is an analysis and visualization pipeline, it has more dependencies than the average package, so this will require a little patience. `devtools` is capable of tracking down most of these dependencies automatically, but there are a few packages from BioConductor that it will not be able to find. It's a good idea to install those first:

~~~~
install.packages("BiocManager")
BiocManager::install(c("qvalue","biomformat","ggtree"))
install.packages("devtools")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
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


### Installing *phylogenize* package for use with QIIME 2 or in another conda environment

QIIME 2 runs in a conda environment, meaning it has its own installation of R and related packages. To run *phylogenize* with QIIME 2, you will need to install *phylogenize* within the QIIME 2 conda environment, then install the [q2-phylogenize plugin](https://bitbucket.org/pbradz/q2-phylogenize). The instructions are similar for installing in any other conda environment, except of course you won't need the plugin.

First, switch to the correct environment. For QIIME 2, this is accomplished with `source activate qiime2-2019.4` (see [here](https://docs.qiime2.org/2019.4/install/native/#activate-the-conda-environment)). Note that you may need to replace "2019.4" with the most recent version of QIIME2, e.g., "2020.6".

Installing *phylogenize* within conda is a little tricky. You will need to manually install a few libraries and packages that are either not included, or difficult to install from source. From the UNIX command line:

```
conda install -c conda-forge -c bioconda -c r libcurl r-devtools bioconductor-rhdf5lib r-magick r-git2r r-shiny
```

Next, run R within the same environment and install the *phylogenize* library. However, you will probably need to work around a [known issue in conda](https://github.com/r-lib/devtools/issues/1722) before calling `devtools::install_bitbucket`. The following should work (from within R).

*Note (added 8/4/2020)*: On recent versions of OS X, you will probably need to change "/bin/tar" to "/usr/bin/tar" when calling `Sys.setenv` in this code block.

```
options(unzip="internal")
Sys.setenv(TAR="/bin/tar")   # replace with path to tar in your installation, if necessary
```

These options only stick around for your current session, so re-run them if you quit and come back to the installation process. You shouldn't need to run them once phylogenize is installed, though.

Next try:

```
install.packages("BiocManager")
BiocManager::install(c("qvalue","biomformat","ggtree"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
devtools::install_bitbucket("pbradz/phylogenize/package/phylogenize")
phylogenize::install.data.figshare()    # nb: this may take a while to download
```

To avoid causing problems elsewhere in QIIME2, I recommend *not* updating any other packages if prompted.

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
| burst_dir | "/usr/local/bin" | String. Path to the BURST or vsearch binaries. |
| burst_bin | "burst12" | String. Name of the BURST or vsearch binary. |

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

## FAQ

### I get an error stating `BURST failed with error code 127`.

This indicates that the BURST binary couldn't be found. Check that you have the correct binary for your operating system and architecture, that it is installed to wherever `burst_dir` is set to point, and that you can run the BURST binary by itself.

### I get a different BURST error.

Update (5/4/2020): This often means that BURST is not compatible with your hardware or operating system. First, make sure you have the right version for your operating system (OS X, Linux, or Windows). If you do, try following the instructions for older computers above. As an alternative, as of phylogenize 0.94 you can now use vsearch instead of BURST: you may have more luck going that route.

### I get an error stating `Error in plotted.pheno.trees[[pn]] : subscript out of bounds`.

Update (11/4/2019): in the latest version of *phylogenize* you should get a more informative error message; please file a bug if this still happens!

This means that you have too few taxa mapped to MIDAS IDs in every phylum *phylogenize* tests. "Too few" here is defined by the options `treemin` (minimum taxa per phylum as an absolute number, default: 5) and `pctmin` (minimum taxa per phylum as a fraction of all taxa in that phylum, default: 0.01, range: 0 to 1.0). You can try lowering the cutoffs, but they are in place because with so few observations, even if you were to get gene hits with *phylogenize*, the results probably would not be that meaningful.

This situation can happen if you, for example, are using a "testing" dataset with very few ASVs, if your read depth was extremely low and few ASVs were detected, or if you are sequencing a community where almost nothing maps to the MIDAS database. If you are using a testing dataset, using the full dataset should solve the problem (if you really want to have something fast to run, you could try only keeping ASVs that map to a single phylum). Assuming you started with lots of ASVs, taking a look at the intermediate BURST output file `output_assignments.txt` should tell you how many of them were successfully mapped to MIDAS IDs: if there are only a few entries with a percent identity above 98.5\% (default), lack of reference genomes is likely your problem.

### When I try to install phylogenize in conda, it spins on "solving environment" for hours without finishing.

There are a couple of things to try:
 * If your conda/QIIME2 installation is not brand new, try removing your existing installation and installing miniconda3 again from scratch. This has helped me before, possibly because there was stuff installed in the 'base' conda environment that was conflicting with requirements for phylogenize.
 * Try using the [mamba](https://github.com/TheSnakePit/mamba) environment solver, instead of the one built into conda.
 * Try using [strict channel priority](https://www.anaconda.com/blog/understanding-and-improving-condas-performance).

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
 
