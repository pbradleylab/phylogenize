% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dummy.R
\name{build.fx.matrix}
\alias{build.fx.matrix}
\title{Assemble an effect-size matrix from a set of distinct effects.}
\usage{
build.fx.matrix(samples, species, fx, wspecies, wsamples)
}
\arguments{
\item{samples}{Vector of strings representing the sample names.}

\item{species}{Vector of strings representing species names.}

\item{fx}{Numeric vector of sizes of effects to be added.}

\item{wspecies}{List of vectors of strings giving the subset of species to be
affected in each distinct effect.}

\item{wsamples}{List of vectors of strings giving which samples to be
affected, one for every distinct effect.}
}
\value{
An effect size matrix.
}
\description{
Assemble an effect-size matrix from a set of distinct effects.
}
\keyword{internal}
