% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plm-functions.R
\name{calc.ess}
\alias{calc.ess}
\title{Main function to calculate environmental specificity scores.}
\usage{
calc.ess(abd.meta, pdata = NULL, b.optim = NULL, ...)
}
\arguments{
\item{abd.meta}{A list giving an abundance matrix and metadata.}

\item{pdata}{Named numeric vector giving priors per environment.}

\item{b.optim}{If not NULL, use this value for the regularization parameter
$b$, otherwise optimize it.}
}
\value{
A list with the following components:
\describe{
  \item{b.optim}{Shrinkage parameter b obtained through optimization}
  \item{ess}{Logit-transformed, shrunken estimates of specificity.}
  \item{regularized}{Non-transformed regularized values.}
  \item{priors}{Values of priors.}
  \item{phenoP}{Prior for environment of interest.}
}
}
\description{
Some particularly relevant global options are:
\describe{
  \item{env_column}{String. Name of column in metadata file containing the
  environment annotations.}
  \item{dset_column}{String. Name of column in metadata file containing the
  dataset annotations.}
  \item{which_envir}{String. Environment in which to calculate prevalence or
  specificity. Must match annotations in metadata.}
  \item{prior_type}{String. What type of prior to use ("uninformative" or
  "file").}
}
}
