% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report-functions.R
\name{non.interactive.plot}
\alias{non.interactive.plot}
\title{A fall-back plotting option for when \code{hack.tree.labels} fails, designed
to produce the same kind of output.}
\usage{
non.interactive.plot(tree.obj, file, name)
}
\arguments{
\item{tree.obj}{A ggtree object.}

\item{file}{File to which an SVG representation of this tree object will be
written.}

\item{name}{String. the name of the taxon being added. Used for title.}
}
\description{
A fall-back plotting option for when \code{hack.tree.labels} fails, designed
to produce the same kind of output.
}
