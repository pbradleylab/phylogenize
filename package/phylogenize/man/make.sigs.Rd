% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/significance.R
\name{make.sigs}
\alias{make.sigs}
\title{Get vectors of significant genes from result tables.}
\usage{
make.sigs(
  results,
  cuts = c(strong = 0.05, med = 0.1, weak = 0.25),
  method = qvals,
  exclude = NULL,
  min.fx = 0,
  ...
)
}
\arguments{
\item{results}{List of result matrices with two rows (effect size and
p-value) and one column per gene tested.}

\item{cuts}{Named numeric vector giving different significance cutoffs.}

\item{method}{Function that will be used to adjust raw p-values in
\code{results}.}

\item{exclude}{String vector of genes to exclude (optional).}

\item{min.fx}{Minimum effect size for calling something significant.}

\item{...}{Extra parameters to be passed to `method`.}
}
\value{
List (per taxon) of string vectors of significant hits.
}
\description{
Get vectors of significant genes from result tables.
}
