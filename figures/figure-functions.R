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
        mutate(data = map(data, ~ mutate(.,{{ .qvcol }} := phylogenize:::qvals({{ .pvcol }}, error_to_file = FALSE)))) %>%
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

# Unzip 16S or shotgun data if necessary
unzip_file <- function(path, name) {
    if (!file.exists(file.path(path, name))) {
        if (file.exists(file.path(path,name))) {
            message(paste("unzipping data: ", name))
            system2("xz", paste0("-d ",
                        file.path(path, paste0(name".xz"))),
                        "-k")
        } else {
            stop(paste(name, "file not found")
        }
    }
}

# Makes a phylogenize rendered report for various associations. Note, this does
# not contain all of the options in the program but those that are used for 
# this report rendering as done in figures.R
generate_report <- function(association_type, output_file_suffix = NULL, linearize = FALSE) {
    # Initialize variables used by majority associations
    out_dir_suffix = NULL; abundance_file = ""; metadata_file = ""
    db_version = "midas_v1.2"
    # Add in linear variables suffix and if should be linear
    if (association_type == "hmp16s-linear" || association_type == "emp-linear") {
        output_file_suffix="-linear"
        linearize = TRUE
    } else if (association_type == "hmpshotgun") {
        db_version = "midas_v1.0"
    }

    param_table <- list(
        "hmp16s" = list(
            abundance_file = "hmp-16s-dada2-full.tab",
            metadata_file = "hmp-16s-phylogenize-metadata-full.tab"),
        "hmp16s-linear" = list(
            abundance_file = "hmp-16s-dada2-full.tab",
            metadata_file = "hmp-16s-phylogenize-metadata-full.tab"),
        "hmpshotgun" = list(
            abundance_file = "hmp-shotgun-bodysite.tab",
            metadata_file = "hmp-shotgun-bodysite-metadata.tab"))

    output_file <- switch(association_type,
                            paste0(file.path(hmp_dir, "16S-results"), output_file_suffix, ".html"),
                            paste0(file.path(hmp_dir, "16S-linear-results"), output_file_suffix, ".html"),
                            paste0(file.path(hmp_dir, "shotgun-results"), output_file_suffix, ".html"),
                            paste0(file.path(emp_dir, "emp-plant-rhizosphere"), output_file_suffix, ".html"),
                            paste0(file.path(emp_dir, "emp-plant-rhizosphere-linear"), output_file_suffix, ".html"))

    out_dir <- switch(association_type,
                        file.path(hmp_dir, "16S-results", output_file_suffix),
                        file.path(hmp_dir, "16S-linear-output", output_file_suffix),
                        file.path(hmp_dir, "shotgun-output", output_file_suffix),
                        file.path(emp_dir, "plant-rhizosphere-phylo", output_file_suffix),
                        file.path(emp_dir, "plant-rhizosphere-linear", output_file_suffix))

    params <- c(
        output_file = output_file,
        out_dir = out_dir,
        in_dir = ifelse(association_type %in% c("hmp16s", "hmp16s-linear", "hmpshotgun"), hmp_dir, emp_dir),
        type = ifelse(association_type %in% c("hmp16s", "hmp16s-linear"), "16S", "midas"),
        db_version = db_version,
        which_phenotype = ifelse(association_type %in% c("hmp16s", "hmp16s-linear", "hmpshotgun"), "prevalence", "specificity"),
        which_envir = ifelse(association_type %in% c("hmp16s", "hmp16s-linear", "hmpshotgun"), "Stool", "Plant rhizosphere"),
        abundance_file = ifelse(association_type %in% c("hmp16s", "hmp16s-linear", "hmpshotgun"), param_table[[association_type]]$abundance_file, abundance_file),
        metadata_file = ifelse(association_type %in% c("hmp16s", "hmp16s-linear", "hmpshotgun"), param_table[[association_type]]$metadata_file, metadata_file),
        data_dir = system.file(package="phylogenize", "extdata"),
        input_format = "tabular",
        vsearch_dir = VSEARCH_DIR,
        ncl = NCL,
        meas_err = TRUE,
        pryr = FALSE,
        linearize = linearize,
        single_dset = association_type %in% c("emp", "emp-linear"),
        env_column = ifelse(association_type %in% c("emp", "emp-linear"),"empo_3",NA),
        biom_file = ifelse(association_type %in% c("emp", "emp-linear"),"emp_deblur_orig_metadata.biom",NA),
        use_rmd_params = FALSE
    )
    return(params)
}