% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/significance.R
\name{get.top.N}
\alias{get.top.N}
\title{Return the significant hits with the N smallest p-values.}
\usage{
get.top.N(
  p,
  sigs,
  signs,
  results,
  level = "strong",
  exclude = NULL,
  N = 25,
  total.n.cutoff = 0,
  genomes.per.protein = NULL
)
}
\arguments{
\item{p}{A taxon}

\item{sigs}{The output of \code{make.sigs}.}

\item{signs}{The output of \code{make.signs}.}

\item{results}{List of result matrices, one per taxon.}

\item{level}{Significance level (must be in \code{sigs[[1]]}).}

\item{exclude}{Optional: exclude these genes from any list.}

\item{N}{Integer; how many hits to return.}

\item{total.n.cutoff}{Optional: if \code{genomes.per.protein} provided, only
return hits found in at least this many genomes.}

\item{genomes.per.protein}{Optional: list (one per taxon) of named numeric
vectors giving the number of genomes that each protein was found in.}
}
\value{
A named vector of N significant hits in descending order of
    significance.
}
\description{
Return the significant hits with the N smallest p-values.
}
