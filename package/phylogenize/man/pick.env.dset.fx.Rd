% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dummy.R
\name{pick.env.dset.fx}
\alias{pick.env.dset.fx}
\title{Function for picking species to simulate.}
\usage{
pick.env.dset.fx(
  n.species,
  env.n.affected,
  dset.n.affected,
  divided,
  env.mean = 0,
  env.sd = 0,
  dset.mean = 0,
  dset.sd = 0
)
}
\arguments{
\item{n.species}{Number of species to simulate}

\item{env.n.affected}{Number of species to be affected by environment effects.}

\item{dset.n.affected}{Number of species to be affected by dataset effects.}

\item{divided}{Output of \code{divide.samples}}
}
\value{
A list:
  \item{mtx}{An effect size matrix.}
  \item{fx}{The chosen, randomly generated effect sizes (see
    \code{pick.env.dset.fx}.}
}
\description{
Function for picking species to simulate.
}
\keyword{internal}
