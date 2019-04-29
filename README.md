# phylogenize (v0.9 beta)


## Running locally without the web server

To run *phylogenize* locally without the web server, you will need to first install the R package:

`devtools::install_bitbucket("pbradz/phylogenize/package/phylogenize")`

You may use this package directly, or together with the QIIME2 interface (see: [https://bitbucket.org/pbradz/q2-phylogenize]).



To run *phylogenize* locally without the web server, you will need to render the RMarkdown notebook "phylogenize-report.Rmd," overriding the values in the header. An example script is provided in `example-job.R`. Here is a full list of fields that can be overridden and their descriptions:

| Field            | Description |
|------------------|----------------------------------------------------------------------------------------------------------------------------------|
| ncl              | Maximum number of threads to run simultaneously. If you are short on memory, keep this at 1 to run phylogenize single-threaded. |
| type             | Can be "midas" (shotgun data) or "16S" (16S amplicon data). |
| out_dir          | Directory for storing output. Will be created if it doesn't exist. |
| in_dir           | Directory where input files are stored. |
| abundance_file   | If providing tabular data, this is the tab-delimited data matrix. |
| metadata_file    | If providing tabular data, this is the tab-delimited metadata matrix. |
| biom_file        | If providing BIOM-formatted data, this is the .biom file containing both data and sample metadata. |
| input_format     | Can be "tabular" or "biom" |
| env_column       | Which column of the metadata matrix contains environment annotations? Defaults to "env" |
| dset_column      | Which column of the metadata matrix contains dataset annotations? Defaults to "dataset" |
| phenotype_file   | Optional phenotype tabular file (allows you to skip calculating it) |
| db_version       | Can be "midas_v1.0" or "midas_v1.2". Defaults to "midas_v1.2" |
| which_phenotype  | Can be "prevalence" or "specificity." Defaults to "prevalence" |
| which_envir      | Which environment should be chosen to calculate the phenotype? Needs to match an annotation in the sample metadata. |
| prior_type       | For calculating specificity, what should the prior probability of the environment be? Can be "uninformative" or "provided" (if prior_file is provided). Defaults to "uninformative" |
| prior_file       | Tab-delimited file giving priors per environment. |
| minimum          | Lower bound for the number of species a gene must appear in, in order to be counted as positive and significant. Defaults to 3 |
| treemin          | Lower bound for the number of observed taxa per phylum-level tree. Phyla with fewer observed taxa will be skipped. Defaults to 5 |
| assume_below_LOD | Can be TRUE or FALSE. If true, taxa with zero observations will be assumed to be absent everywhere; if false, taxa with zero observations will be removed. Defaults to TRUE |
| burst_dir        | Location of BURST binaries. Defaults to "`./bin`" |

You will probably need to install dependencies first (run `install-dependencies.R`, plus there may be some system-wide dependencies like libxml and libcairo). You will also need to download a version of BURST (see [github.com/knights-lab/BURST](https://github.com/knights-lab/BURST)).

## Running locally with the web server

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
