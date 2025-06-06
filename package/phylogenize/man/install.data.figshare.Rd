% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load-process.R
\name{install.data.figshare}
\alias{install.data.figshare}
\title{Download data from figshare (or provide it locally) and un-gzip it into the
package directory so that it can be imported.}
\usage{
install.data.figshare(
  data_path = NULL,
  figshare_url = paste0("https://ndownloader.figshare.com/files/",
    "43692576?private_link=987aeecdfebd2da02302")
)
}
\arguments{
\item{data_path}{Optional: provide a path to a local file containing a
compressed .tar archive of data. Must extract to the subdirectory
\code{extdata/}.}

\item{figshare_url}{Optional: override the URL from which to obtain the data.}
}
\description{
Download data from figshare (or provide it locally) and un-gzip it into the
package directory so that it can be imported.
}
