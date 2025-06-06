% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plm-functions.R
\name{rndpos.fx.pv}
\alias{rndpos.fx.pv}
\title{Wrapper to return random effect sizes and p-values with a bias; useful for
testing.}
\usage{
rndpos.fx.pv(
  m,
  p,
  tr,
  coefname = "m",
  restrict = NULL,
  pos.pct = 0.1,
  pos.fx = 2,
  pos.beta = c(shape1 = 1, shape2 = 20),
  meas_err = FALSE,
  se = 0.1
)
}
\arguments{
\item{m}{Named numeric vector of gene presence/absences per taxon.}

\item{p}{Named numeric vector of phenotype values per taxon.}

\item{tr}{Phylogeny relating taxa (class \code{"phylo"}).}

\item{coefname}{Which coefficient from the phylolm to return (meaningless
here)?}

\item{restrict}{If not NULL, a character vector of taxa to consider
(meaningless here).}

\item{pos.pct}{Fraction of times to return a positive effect.}

\item{pos.fx}{Shift in mean for positive effects, in standard deviations.}

\item{pos.beta}{Length-2 numeric vector of beta parameters for true effect
p-values.}
}
\value{
Length-2 numeric vector with names \code{"Estimate"} and
    \code{"p.value"}, with distributions N(0,1) and U(0,1), respectively.
}
\description{
Wrapper to return random effect sizes and p-values with a bias; useful for
testing.
}
\keyword{internal}
