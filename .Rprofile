library(utils)

# Automatically select the cran mirror that is closest to the user
local({r <- getOption("repos")
       r["CRAN"] <- "https://cloud.r-project.org" 
       options(repos=r)
})

# Install renv for package management
if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
}
renv::install("arrow", prompt=FALSE)

# Install phylogenize
if (!requireNamespace("phylogenize", quietly = TRUE)) {
  renv::install("pbradleylab/phylogenize/package/phylogenize@dev")
}
library(phylogenize)
phylogenize::install.data.figshare()
