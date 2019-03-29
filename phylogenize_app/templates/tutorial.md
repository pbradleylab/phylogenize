# Tutorial

# What is *phylogenize*?

*phylogenize* is a web tool that allows users to link microbial genes to either microbial prevalence in, or specificity for, a given environment, while taking into account the phylogenetic relationship between microbes. The method is described in [Bradley, Nayfach, and Pollard (2018)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006242).

# How do I use *phylogenize*?

To use *phylogenize*, you need the following things:

  * **Your data.** This should be a matrix where:
	  * each row is a species ID (shotgun data) or a sequence variant (16S data)
		* each column is a sample ID of your choosing
		* each entry is any of the following:
		  * a number of reads
			* a relative abundance
			* 0 for absence and 1 for presence
	* **Your sample metadata.** This should be a matrix mapping sample IDs to environments and/or datasets, where:
		* each row corresponds to a sample in your data
		* there are at least the following three columns:
		  * *sample*: this column contains the sample IDs, matching the sample IDs in your data
			* *env*: this column contains the environment labels corresponding to each sample
			* *dataset*: this column contains the batch or dataset to which each sample belongs
			  * Note: if there is only one dataset or environment, it is fine to have all of these entries be the same.
				* Note: right now, you can only calculate specificity if there is just one dataset.
	* **Information about your data**. This tells *phylogenize* how your data were processed (shotgun vs. 16S) and what version of the database to use (most people can leave this set to the default).
	* **A choice of phenotype.** This is the variable you are linking to gene presence. It can be one of:
	  * prevalence: how frequently is a microbe observed in a given environment?
		* specificity: how specific is a microbe for a given environment compared to all others?
		  * Note: to calculate specificity, there must be more than one environment present.
	* **A choice of environment.** Prevalence and specificity are defined with respect to a given environment, like "soil" or "stool" or "marine."
	  * Note: this entry needs to match the *env* column in your sample metadata above.

![Example of form](phylogenize-form-example.png)


## What types of data does *phylogenize* accept?

*phylogenize* requires taxonomic abundances. These can either be derived from shotgun data or 16S data.

### Shotgun data

Shotgun data should be processed with [MIDAS](https://github.com/snayfach/MIDAS), with database version either v1.2 (recommended) or v1.0 (deprecated). You do not need to run the full MIDAS pipeline; you will just need to do metagenomic species profiling as described [here](https://github.com/snayfach/MIDAS/blob/master/docs/species.md). The output file you need should be called `species_profile.txt`. 

### 16S data

16S data should be processed using a denoising algorithm. Two common options are [Deblur](https://github.com/biocore/deblur) and [DADA2](https://benjjneb.github.io/dada2/). Both of these algorithms will result in abundances for "amplicon sequence variants" (ASVs), which correspond to the distinct DNA sequences that are obtained by denoising.

*phylogenize* works by mapping these ASVs back to [MIDAS](https://www.github.com/snayfach/MIDAS) genome clusters using the fast aligner [BURST](https://github.com/knights-lab/BURST). The default threshold for matching a genome cluster is 98.5% sequence identity, since we found that this gives a reasonable sensitivity-specificity threshold. Ambiguous matches are discarded, and reads/abundances mapping to the same genome cluster are summed together. (In practice, though, we have never actually seen ambiguous matches.)

## Gathering your data and metadata

Your data and metadata should be in either *tabular* or *BIOM* format. Tabular format is simpler, but takes up more space. The maximum file upload is 35M, so if your data are larger, you will need to convert them to BIOM format (see below) or, alternatively, to run *phylogenize* locally.

You can upload either two tab-delimited files (data and metadata) or one BIOM file (data plus metadata). Right now you cannot, for example, provide separate tabular metadata for a BIOM-formatted file.

### Tabular format data

If you choose the tabular format, your data and metadata will be in two separate tab-delimited files. They should look something like this:

![Top corner of data](data-corner.png)

![Top corner of metadata](metadata-corner.png)

### BIOM format data

[BIOM format](http://biom-format.org/) data is more compact for tabular data and is more suitable for larger datasets. This is because 1. BIOM represents matrices that are mostly zeroes (i.e., "sparse") more efficiently, and most taxonomic matrices are sparse, and 2. BIOM files can use a binary representation (HDF5) that takes up less space than plain text.

A BIOM file can contain multiple tables, so you will only need to upload one file. This one file needs to have the following tables:

  * Observation matrix or "OTU table"
	  * Despite the name, rows should actually be either MIDAS species IDs or amplicon sequence variants (ASVs), not OTUs
	* Sample metadata matrix
	  * Make sure the sample metadata includes columns labeled *env* and *dataset*, just as with tabular data.

To convert your tabular data into a BIOM file, [install the `biom` utility](http://biom-format.org/index.html). Then run the conversion as follows (adapted from the biom-format.org documentation):

    biom convert -i your_data.tab -o your_data_and_metadata.biom --to-hdf5 --table-type="OTU table" --sample-metadata-fp metadata.tab

Replace `your_data.tab` and `metadata.tab` with the names of your tab-delimited data files and, optionally, `your_data_and_metadata` with the name of your dataset.

For more help using the `biom` utility, try looking at the pages entitled ["Converting between file formats"](http://biom-format.org/documentation/biom_conversion.html) and ["Adding sample and observation metadata to biom files"](http://biom-format.org/documentation/adding_metadata.html). These documents also cover how to rename columns to match what *phylogenize* accepts, if you already have a BIOM file, for example, because you've run QIIME.

## Information about your data

"Data type" is relatively self-explanatory: is your data 16S amplicon sequencing data, or shotgun metagenomics data processed with MIDAS?

For most people, the MIDAS database version will be the latest, v1.2. We also provide the option to use a previous version of the database, since species are named differently in MIDAS v1.2 and v1.0, but for most people the default (1.2) will be appropriate.

## Choosing a phenotype

*phylogenize* can associate genes with either microbial *prevalence* or microbial *specificity* for a given environment. Prevalence gives how often a microbe is found in an environment. Specificity compares this prevalence to prevalences in all other provided environments. The following toy example gives an intuition for what this means:

![Cartoon of prevalence and specificity calculations](prevalence-specificity.png)

Microbe 1 is detected (black boxes) often in both samples from healthy (A-G) and sick (H-N) individuals. Microbe 2 is detected more frequently in healthy than sick samples, and microbe 3 is detected relatively seldom in either. This means that microbes 1 and 2 will have *high prevalence* in the "healthy" environment, while microbe 3 will have *low prevalence*. However, microbe 2 is the only microbe to have *high specificity* for "healthy" over "sick," because of the difference in prevalence between these environments.

Prevalence will tend to capture both cosmopolitan microbes and ones that are particularly well-adapted to a given environment. Specificity is useful for contrasting two sites (like gut vs. skin) or states (healthy vs. sick).

The exact details of how prevalence and specificity are calculated can be read in our [paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006242). Basically, we perform additional steps to re-weight datasets equally (for prevalence), to make sure that small numbers of observations don't skew the results (for specificity), and to transform these quantities into something approximately normally-distributed (both).

Note: it is possible with low numbers of samples that no microbes will be specific to a given environment, in which case *phylogenize* will return an error because there is not enough signal to compute specificity.

Note: while you do not have to rarefy your data, you should make sure that read depth does not systematically differ across environments. If it does, then rarefying may be a good idea because some microbes will simply appear more by chance in deeper-sequenced data.

## Choosing an environment

Simply enter the name of the environment for which you are interested in calculating a phenotype. This should match an entry in the *env* column of the sample metadata exactly. At least two samples per environment are necessary.

## Results pages

Your results will be delivered at a URL that should look like: `www.phylogenize.org/results/<result_id>` with a long string replacing `<result_id>`. This is where you should eventually be able to pick up your results. We do not make the result IDs searchable, so it is a good idea to bookmark or otherwise note this URL.

Anyone with this URL will be able to see your results, so only share it with members of your team and collaborators. (Also, even though we don't make these URLs publicly searchable, you still should not upload any data with personally-identifiable information, sensitive personal information, or any other information that needs to remain 100% completely secure.)

The results have two components. The first is an HTML report summarizing and visualizing the data. The second is a .tgz file (compressed archive) that should contain the full set of results, including:

  * `phenotype.tab`: Tab-delimited file giving the calculated phenotype values for every microbe.
	* `pos-sig-thresholded.csv`: Comma-delimited file giving all FIGfam gene families significantly positively associated with the phenotype, within a given phylum. This file also contains descriptions of the significant genes.
	* `all-results.csv`: Comma-delimited file giving the effect size and uncorrected p-value for the association of every FIGfam gene family in every phylum with the calculated phenotype. This can be useful if, for example, you want to apply your own threshold for significance or a different p-value correction method, or if you want to further filter your results by effect size.
	* `enr-table.csv`: FDR-corrected SEED subsystem enrichments for the significantly positively associated FIGfam gene families. These are pre-computed for three gene-wise levels of significance: "strong" (0.05), "med" (0.1), and "weak" (0.25). SEED subsystems significant at an FDR of 25% are returned.
	* `enr-overlaps.csv` and `enr-overlaps-sorted.csv`: These files give the individual significant FIGfams underlying any SEED subsystem enrichments, which can aid interpretation.
	* `progress.txt` and `stderr.txt`: These are the outputs of *phylogenize*; if your job finished without an error message, you probably don't need to refer to them.

To make *phylogenize* available for other users, results will be deleted approximately 7 days after completion, so it is a good idea to save the report and results to disk before then.

# Troubleshooting

  * *phylogenize* appears to be stuck a little over halfway through.
	  * This is probably normal; the progress bar is approximate and this is likely to be the step where associations are actually calculated. If you check back in several minutes, you should see more subtle signs of progress under "monitor warnings/errors."
	* *phylogenize* doesn't appear to be running my job: the bar is stuck at 0%.
		* Because *phylogenize* takes a lot of memory to run, only one job can run at a time. If you reload the page and look under "monitor warnings/errors" you should see information about how many other jobs are in the queue. If you don't see this information and the job is still stuck at 0% for a prolonged period of time, contact support (see About/Contact in the navigation bar).
	* *phylogenize* gave me an error that I don't understand.
	  * Check that your data tables satisfy the above requirements (in particular, are your metadata columns named correctly and do your sample IDs match up?).
		* If your job ran out of memory, try converting it to BIOM format (see above).
		* If there's still a problem, e-mail the webmaster (see About/Contact in the navigation bar) and attach the output that appears when you click "monitor warnings/errors".

# How do I run *phylogenize* locally?

If you are getting an "out of memory" error, you have a lot of jobs, or you have a more specific use case that is not covered by the server, we recommend that you download *phylogenize* from our code repository at https://www.bitbucket.com/pbradz/phylogenize. The core of *phylogenize* is an RMarkdown notebook, `phylogenize-report.Rmd`, which can be run either with or without the web app (written in Python). The web app provides a user interface and a scheduling system. 

# How should I cite *phylogenize*?

See About/Contact for more information.


