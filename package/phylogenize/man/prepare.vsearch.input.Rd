% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process-16S.R
\name{prepare.vsearch.input}
\alias{prepare.vsearch.input}
\title{Prepare input file for alignment}
\usage{
prepare.vsearch.input(mtx, ...)
}
\arguments{
\item{mtx}{A presence/absence or abundance matrix, with row names equal to
amplicon sequence variant DNA sequences.}
}
\description{
\code{prepare.vsearch.input} outputs a FASTA file of the sequences in the
input 16S data for analysis using vsearch or vsearch.
}
\details{
Some particularly relevant global options are:
\describe{
  \item{vsearch_infile}{String. File name of the sequences written to disk
  and then read into vsearch/vsearch.}
}
}
