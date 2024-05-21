# Returns the last element of a vector
last_elem <- function(x) x[length(x)]

# Calculates q-values for p-values within groups specified by a given variable
# in a table. Takes a table (tbl) and optional arguments for grouping (.nestby),
# p-value column name (.pvcol), and q-value column name (.qvcol). Then the table
# is grouped by the specified variable, data is nested within each group, 
# and q-values are calcaulted for each. Finally, the dataframe is unnested and
# the original structure is returned.
make_qvs <- function(tbl, .nestby = "phylum", .pvcol = "p.value", .qvcol = "q.value") {
    tbl %>%
        group_by({{ .nestby }}) %>%
        nest() %>%
        mutate(data = map(data, ~ mutate(., {{ .qvcol }} := phylogenize:::qvals({{ .pvcol }}, error_to_file = FALSE)))) %>%
        unnest()
}

# Calculates equivalent p-values for effect size, standard error, and degrees
# of freedom within a table. The function takes a table (tbl) and an optional
# argument for the minimum effect size (mfx). It creates a new column (equiv.pv)
# using the phylogenize package's equivalent test function, which calculates
# p-values based on effect size, standard error, and degrees of freedom.
make_equivs <- function(tbl, mfx=0.5) {
    mutate(tbl,
           equiv.pv = pmap_dbl(list(effect.size, std.err, df),
                               phylogenize:::equiv_test,
                               min_fx=mfx))
}

# Processes data by calculating equivalent p-values, q-values, and adjusted q-values.
# The function takes a table (tbl) and applies a series of operations: first, it
# calculates equivalent p-values using make_equivs, then it calculates q-values
# using make_qvs, and finally, it calculates adjusted q-values using make_qvs
# again with the equivalent p-values and corresponding q-values
process_data <- function(tbl) {
    tbl %>% 
        make_equivs %>% 
        make_qvs %>%
        make_qvs(., .pvcol="equiv.pv", .qvcol="equiv.qv")
}

# Retrieves the alpha parameter from a phylogenetic tree for a given phylum and
# gene. The function takes the phylum name (phylum), gene name (gene),
# phylogenetic tree database (db), and an optional argument for the lab parameter
# (lab) which sets the upper limit for log(alpha). It extracts the phylogenetic
# tree corresponding to the specified phylum from the provided database, selects
# the tips common to the tree and gene presence data, constructs a data frame, and
# performs a phylogenetic linear regression using the phyloglm function from the
# ape package. It returns the estimated alpha parameter if successful, otherwise
# raises a warning and returns NA.
get_alpha <- function(phylum, gene, db=pz.db, lab=8) {
    p <- db$trees[[phylum]]
    tips <- intersect(colnames(db$gene.presence[[phylum]]), p$tip.label)
    df <- data.frame(g=db$gene.presence[[phylum]][gene, tips])
    if (var(df$g) == 0) { 
        warning(paste0("zero-variance: ", gene))
        return(NA)
    }
    tryCatch({
        phyloglm(g ~ 1,
                data=df,
                phy=keep.tip(p, tips),
                log.alpha.bound=lab)$alpha
    }, error=function(e) { 
        warning(e) 
        NA 
    })
}

# Sets the class attribute of an object to "phylo". This is used to mark an object
# as a phylogenetic tree.
ap_tree <- function(x) {
    attr(x, "class") <- "phylo"
    x
}
