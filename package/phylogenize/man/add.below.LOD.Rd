% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc-phenotypes.R
\name{add.below.LOD}
\alias{add.below.LOD}
\title{Add in taxa that were not observed, assuming this means they were
zero-prevalence.}
\usage{
add.below.LOD(pz.db, abd.meta, ...)
}
\arguments{
\item{pz.db}{A database (typically obtained with \code{import.pz.db}).}

\item{abd.meta}{A list consisting of a taxon abundance matrix and the
metadata.}
}
\value{
An updated version of \code{abd.meta}.
}
\description{
Add in taxa that were not observed, assuming this means they were
zero-prevalence.
}
