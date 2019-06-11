# phylogenize (v0.91 beta)


## Running *phylogenize* locally

To run *phylogenize* locally, you will need to first install the R package from within this repository. Because *phylogenize* is an analysis and visualization pipeline, it has more dependencies than the average package, so this will require a little patience. `devtools` is capable of tracking down most of these dependencies automatically, but there are a few packages from BioConductor that it will not be able to find. It's a good idea to install those first:

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

You will also need a copy of the BURST binaries from [their Github](github.com/knights-lab/BURST). BURST is a high-speed pairwise aligner that *phylogenize* uses to map 16S amplicon sequence variants back to a database of genomes. By default *phylogenize* expects these binaries to be in `/usr/local/bin`, but you can override this (see below).

Finally, you will need to download and install the data files that *phylogenize* needs to run. That has been automated so that all you should need to do is run the following function:

```
library(phylogenize)
phylogenize::install.data.figshare()
```

You may use this package by itself, using the web interface, or using the [QIIME2 interface](https://bitbucket.org/pbradz/q2-phylogenize).

To **reinstall** *phylogenize*, the following should be all that's necessary:

```
devtools::install_bitbucket("pbradz/phylogenize/package/phylogenize")
phylogenize::install.data.figshare()
```


### Running *phylogenize* locally in R

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

### Running *phylogenize* locally with QIIME2

To run *phylogenize* with QIIME2, you will need to install *phylogenize* within the QIIME2 conda environment, then install the [q2-phylogenize plugin](https://bitbucket.org/pbradz/q2-phylogenize).

First, switch to the correct environment using `source activate qiime2-2019.4` (see [here](https://docs.qiime2.org/2019.4/install/native/#activate-the-conda-environment)).

Installing *phylogenize* within conda is a little tricky. You will need to manually install a few libraries and packages that are either not included, or difficult to install from source. From the UNIX command line:

```
conda install libcurl
conda install r-devtools
conda install -c bioconda bioconductor-rhdf5lib
conda install -c conda-forge r-magick
```

Next, run R within the same environment and install the *phylogenize* library. However, you will need to work around a [known issue in conda](https://github.com/r-lib/devtools/issues/1722) before calling `devtools::install_bitbucket`. The following should work (from within R):

```
options(unzip="internal")
Sys.setenv(TAR="/bin/tar")   # replace with path to tar in your installation, if necessary
install.packages("BiocManager")
BiocManager::install(c("qvalue","biomformat","ggtree"))
devtools::install_bitbucket("pbradz/phylogenize/package/phylogenize")
phylogenize::install.data.figshare()    # nb: this may take a while to download
```

Finally, install the plugin. From the UNIX command line (i.e., not R):

```
git clone bitbucket.org/pbradz/q2-phylogenize
cd q2-phylogenize
python setup.py build
python setup.py install
```

Further information about how to use q2-phylogenize can be found on its [git repository](https://bitbucket.org/pbradz/q2-phylogenize).


### Running *phylogenize* locally with the web interface

This part of the guide is written assuming you are in a \*nix environment like Ubuntu or OS X. Windows users may be able to accomplish the command-line steps using Cygwin.

*phylogenize* is written as a Flask WSGI app. If you are going to be running it on your own computer, you can probably use the built-in Flask server. (If you are concerned about security and will be allowing potentially untrusted users to use the server, or if the built-in Flask server is inadequate for any other reason, we recommend using a production web server like Apache2. The Flask documentation has [some information](http://flask.pocoo.org/docs/1.0/deploying/mod_wsgi/) on how to get started hosting a WSGI app on Apache2 using `mod_wsgi`.)

We recommend installing Flask using a virtual environment as per the instructions [here](http://flask.pocoo.org/docs/1.0/installation/). Once you have activated the virtual environment, before running the app, you will need to launch the Beanstalk-based queueing system:

    nohup beanstalkd -l 127.0.0.1 -p 14711 &
    nohup python3 phylogenize_app/worker.py &

To keep these jobs running if the terminal is closed, you may want to `disown` these jobs or alternatively, start them in a `tmux` or `screen` session.

To allow multiple simultaneous jobs, you will need to edit `worker.py` and change MaxJobs as appropriate. It may be worth starting with a single simultaneous job, particularly if each job runs multi-threaded, to make sure memory use is appropriate.

Next, you will need to start the server as follows:

     FLASK_APP=phylogenize_app flask run

The application should then be accessible from http://localhost:5000.
