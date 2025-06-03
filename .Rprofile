# Assume that R and pandoc are already installed.
library(utils)

# Automatically select the cran mirror that is closest to the user
local({r <- getOption("repos")
       r["CRAN"] <- "https://cloud.r-project.org" 
       options(repos=r)
})

# Install renv for package management
if (!requireNamespace("renv", quietly = TRUE)) {    
    install.packages("renv", quietly = TRUE, type = "source")
}
# Install the dependancies that rely on biocmanager now
if (!requireNamespace("qvalue", quietly = TRUE)) {
    renv::install("bioc::qvalue", prompt=FALSE)
}
if (!requireNamespace("ggtree", quietly = TRUE)) {
    renv::install("bioc::ggtree", prompt=FALSE)
}
if (!requireNamespace("biomformat", quietly = TRUE)) {
    renv::install("bioc::biomformat", prompt=FALSE)
}

# Install the dev branch of phylogenize
if (!requireNamespace("phylogenize", quietly = TRUE)) {
    renv::install("pbradleylab/phylogenize/package/phylogenize@dev", prompt=FALSE)
}

if (!requireNamespace("Maaslin2", quietly = TRUE)) {
    renv::install("bioc::Maaslin2", prompt=FALSE)
}

#phylogenize::install.data.figshare()

