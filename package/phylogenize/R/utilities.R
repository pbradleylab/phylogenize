# utilities

# helper functions

#' Test whether a value is between two other values (non-inclusive).
#'
#' @param x Value(s) to test (numeric vector).
#' @param y Numeric vector of length 2, giving minimum and maximum values of \code{x}.
#' @export
`%btwn%` <- function(x, y) { (x > min(y)) & (x < max(y)) }

#' Test whether a value is between two other values (inclusive).
#'
#' @param x Value(s) to test (numeric vector).
#' @param y Numeric vector of length 2, giving minimum and maximum values of \code{x}.
#' @export
`%btwn.inc%` <- function(x, y) { (x >= min(y)) & (x <= max(y)) }

#' Intersect two vectors.
#'
#' @param x First vector.
#' @param y Second vector.
#' @export
`%intr%` <- function(x, y) intersect(x,y)

#' Assign names to a vector.
#'
#' @param x Vector to have names assigned (any type).
#' @param y Character vector giving new names for values in \code{x}.
#' @export
`%withnames%` <- function(x, y) { names(x) <- y; x }

#' Abbreviation for \code{names(which(x))}.
#'
#' @param x Boolean vector.
#' @export
nw <- function(x) { names(which(x)) }

#' Find the minimum value of a vector that is still greater than zero.
#'
#' @param x Numeric vector.
#' @export
min.nonzero <- function(x) min(x[x > 0])

#' Count the number of instances of every unique value of a vector.
#'
#' @param x Vector to be counted.
#' @export
count.each <- function(x, na.rm = FALSE) {
    u <- unique(x)
    simplify2array(lapply.across.names(u, function(y)
        sum(x == y, na.rm = na.rm)))
}

#' Abbreviation for \code{grep} with \code{value=TRUE}.
#'
#' @param x Pattern to find.
#' @param y Values to search for pattern \code{x}.
#' @export
grepv <- function(x, y, ...) grep(x, y, ..., value = TRUE)

#' Apply a function to a vector of names, with the returned list having those
#' names.
#'
#' @param X Vector of names.
#' @param FUN A function to apply to the names in \code{X} (typically using them
#'     as list indices).
#' @return A list of the results of applying \code{FUN} to \code{X}, with
#'     \code{X} as the list elements' names.
#' @export
lapply.across.names <- function(X, FUN, ...) {
    r <- lapply(X, FUN, ...)
    names(r) <- X
    r
}

#' Apply a function to a vector of names, with the returned list having those
#' names (progress bar).
#'
#' @param X Vector of names.
#' @param FUN A function to apply to the names in \code{X} (typically using them
#'     as list indices).
#' @return A list of the results of applying \code{FUN} to \code{X}, with
#'     \code{X} as the list elements' names.
#' @export
pblapply.across.names <- function(X, FUN, ...) {
    r <- pblapply(X, FUN, ...)
    names(r) <- X
    r
}

#' Merge two vectors/matrices by (row) names, returning a matrix with row names.
#'
#' @param x First vector or matrix.
#' @param y Second vector or matrix.
#' @return A merged matrix with row names.
#' @export
small.merge <- function(x, y, ...) {
    z <- merge(x, y, by = 0)
    r <- z[, -1]
    rownames(r) <- z[, 1]
    r
}

#' Calculate geometric mean of a vector of values.
#'
#' @param x Numeric vector of values.
#' @export
geommean <- function(x) exp(sum(log(x)) / length(x))

#' Calculate logit of a value or vector of values.
#'
#' @param x Numeric value, or numeric vector of numeric values.
#' @export
logit <- function(x) (log(x / (1 - x)))

#' Calculate inverse-logit of a value or vector of values.
#'
#' @param x Numeric value, or numeric vector of numeric values.
#' @export
logistic <- function(x) exp(x) / (1 + exp(x))

#' Truncate a vector with a particular lower and upper limit.
#'
#' @param x Vector to truncate.
#' @param lim Two-element numeric vector giving the lower and upper limits.
#' @return A vector where elements below (or above) the lower (or upper) limit
#'     have been replaced with that limit.
#' @export
truncated <- function(x, lim = c(logit(0.001), logit(0.25))) {
    y <- x
    y[which(x < lim[1])] <- lim[1]
    y[which(x > lim[2])] <- lim[2]
    y
}

#' Wrapper around fastread that preserves row names.
#'
#' @param location Path to file to be read.
#' @param cn Whether to check that the row names are valid, non-duplicated R row names.
#' @return A data matrix with rownames equal to the first column of the input
#'     file and colnames equal to the first row.
#' @export
fastread <- function(location, cn = TRUE) {
    # rownames are useful
    master <- data.frame(data.table::fread(location, header = T),
                         check.names = cn)
    rn <- master[,1]
    rest <- master[,-1]
    rownames(rest) <- rn
    return(data.matrix(rest))
}

# Parallelization

#' Apply a function over a margin of a matrix in parallel.
#'
#' @param X A matrix.
#' @param MARGIN A number specifying whether to apply over rows (1) or columns (2).
#' @param FUN A function to be applied.
#' @param mc.cores Number of cores to use.
#' @param simplify Whether to simplify the results using \code{simplify2array}.
#' @export
mcapply <- function(X, MARGIN, FUN, mc.cores = 10, simplify = TRUE, ...) {
    if (MARGIN == 1) {
        mlist <- lapply(seq_len(nrow(X)), function(i) X[i, ])
        names(mlist) <- rownames(X)
    } else if (MARGIN == 2) {
        mlist <- lapply(seq_len(ncol(X)), function(i) X[, i])
        names(mlist) <- colnames(X)
    } else {
        stop("invalid MARGIN value")
    }
    r <- mclapply(mlist, FUN, mc.cores = mc.cores, ...)
    if (simplify) simplify2array(r) else r
}

#' Apply a function over a margin of a matrix in parallel (with progress bar).
#'
#' @param X A matrix.
#' @param MARGIN A number specifying whether to apply over rows (1) or columns
#'     (2).
#' @param FUN A function to be applied.
#' @param mc.cores Number of cores to use.
#' @param simplify Whether to simplify the results using \code{simplify2array}.
#' @export
pbmcapply <- function(X, MARGIN, FUN, mc.cores = 10, simplify = TRUE, ...) {
    if (MARGIN == 1) {
        mlist <- lapply(seq_len(nrow(X)), function(i) X[i, ])
        names(mlist) <- rownames(X)
    } else if (MARGIN == 2) {
        mlist <- lapply(seq_len(ncol(X)), function(i) X[, i])
        names(mlist) <- colnames(X)
    } else {
        stop("invalid MARGIN value")
    }
    r <- pbmclapply(mlist, FUN, mc.cores = mc.cores, ...)
    if (simplify) simplify2array(r) else r
}

# Annotation (FIGfams and taxa)

#' Annotate genes using a gene-to-function table.
#'
#' @param x A gene (string) or vector of genes (strings).
#' @param gene.to.fxn A data frame with at least "gene" and "function" as columns.
#' @return A character vector of gene functions, with names equal to \code{x}.
#' @export
gene.annot <- function(x, gene.to.fxn) {
    gf <- data.frame(gene.to.fxn[match(x, gene.to.fxn$gene), ])
    gf$"function" %withnames% gf$"gene"
}


#' Annotate taxa using a taxonomy table.
#'
#' @param tns A vector of taxon IDs.
#' @param taxonomy A data frame with at least "cluster" and "species" columns;
#'     "cluster" is used to match the identifiers in \code{tns}.
#' @return A character vector of species names.
#' @export
tax.annot <- function(tns, taxonomy) {
    taxonomy$species[match(tns, taxonomy$cluster)]
}


#' Calculate the maximum depth of a list of lists.
#'
#' \code{list.depth} calculates the maximum depth of a list of lists. Code is
#' modified from original by user "flodel" on stackoverflow: see
#' https://stackoverflow.com/questions/13432863/determine-level-of-nesting-in-r.
#'
#' @param this A list, possibly containing other lists.
#' @export
list.depth <- function(this) {
    if ((is.list(this)) && (length(this) == 0)) { return(0L) }
    ifelse(is.list(this), 1L + max(sapply(this, list.depth)), 0L)
}

#' Make sure that a list of lists has even depth all the way down.
#'
#' @param this A list, possibly containing other lists.
#' @export
list.even <- function(this) {
    if (!is.list(this)) return(TRUE)
    d <- vapply(this, list.depth, FUN.VALUE=1L)
    this.even <- all(d == d[1])
    return(this.even & all(vapply(this, list.even, TRUE)))
}

#' Find the first list element at a desired depth.
#'
#' \code{first.element.at.depth} returns the first list element if \code{n} is
#' 1; otherwise, the function recursively dives into nested lists in the first
#' list element until it reaches level \code{n}.
#'
#' @param l A list.
#' @param n A number giving a particular depth of nesting.
#' @return The first element of a list
#' @export
first.element.at.depth <- function(l, n) {
    if (n == 1) { l[[1]] } else { (first.element.at.depth(l[[1]], n - 1)) }
}

#' Transform nested lists into a data frame.
#'
#' \code{annotate.nested} takes a nested list (i.e., a list of lists, each of
#' which lists may also be a list of lists), and transforms it into a tabular
#' structure. This is useful for visualizing and processing the results of
#' nested \code{lapply} statements, such as those that might be used when
#' performing a grid search across multiple parameters.
#'
#' @param nested A nested list.
#' @param summarize If NULL, indicates that no summarizing should be done on the
#'     most basic list elements; if a function or closure, this function or
#'     closure will be applied to the most basic list elements (for example,
#'     \code{mean} would take the mean of vectors at the "bottom" of the nested
#'     list), and the return value are incorporated into the final table instead
#'     of the originals.
#' @param n An optional vector of values to bind to the result table.
#' @param n.names An optional vector of names for the columns in the final table.
#' @param stop.at Stop at this depth and summarize.
#' @return A 2D data frame representation of the (summarized) values in the
#'     nested lists.
annotate.nested <- function(nested,
                            summarize = NULL,
                            n = NULL,
                            n.names = NULL,
                            stop.at = 0) {
    if (is.null(nested)) {
       return(c())
    }
    nestedness <- list.depth(nested)
    if (nestedness == stop.at) {
        if (is.null(summarize)) {
            if (is.matrix(nested)) {
                values <- data.frame(nested)
            } else {
                values <- data.frame(value = nested)
            }
            rownames(values) <- names(nested)
        } else {
            val.raw <- summarize(nested)
            if (class(val.raw) != "numeric") {
                if ((!is.null(nrow(val.raw))) || (nrow(val.raw) > 0)) {
                    rownames(val.raw) <- NULL
                }
            }
            if (!is.matrix(val.raw)) {
                values <- data.frame(value = val.raw)
            } else {
                values <- data.frame(val.raw)
            }
        }
        if (is.null(nrow(values))) {
            n.mtx <- matrix(rep(n, (length(values))),
                            nc = length(n),
                            byrow = TRUE)
        } else {
            n.mtx <- matrix(rep(n, nrow(values)), nc = length(n), byrow = TRUE)
        }
        if (!is.null(n.names)) {
            if (length(n.names) > ncol(n.mtx)) {
                pz.warning(paste0(
                    "more names provided than columns; just using first ",
                    ncol(n.mtx)))
            }
            colnames(n.mtx) <- n.names[1:ncol(n.mtx)]
        }
        if (is.matrix(nested)) {
            final.df <- cbind(names = rownames(values), n.mtx, values)
        } else {
            final.df <- cbind(names = rownames(values), n.mtx, value = values)
        }
        rownames(final.df) <- NULL
        return(final.df)
    } else {
        if (is.null(names(nested))) { names(nested) <- 1:length(nested) }
        Reduce(rbind,
               lapply.across.names(names(nested), function(x)
                   annotate.nested(nested[[x]],
                                   summarize = summarize,
                                   n = c(n, x),
                                   n.names = n.names,
                                   stop.at = stop.at)))
    }

}

#' ``Zip'' a list of data together with its names.
#'
#' This function takes a list and constructs a new structured list with fields
#' "name" and "data". This is typically called by \code{zipLapply} or
#' \code{zipSapply}. The purpose is to allow the user to iterate over a list
#' using lapply- or sapply-like syntax while still having access to the names of
#' the list elements.
#'
#' @param iterover List to iterate over.
#' @return A list with the same number of elements as \code{iterover}. Each list
#'     element is itself a list with the slots "name" and "data". For the ith
#'     element of the returned list, "name" contains \code{names(iterover)[i]}
#'     and "data" contains \code{iterover[[i]]}.
zipData <- function(iterover) {
    n <- names(iterover)
    names(n) <- n
    return(lapply(n, function(x) {
        list(name = n[x], data = iterover[[n[x]]])
    }))
}

#' A wrapper for lapply that allows access to list element names.
#'
#' @param iterover List to iterate over. This list will be transformed by
#'     \code{zipData}.
#' @param fxn A function to apply to each element \code{x} of
#'     \code{zipData(iterover)}, where the name is accessible as \code{x$name}
#'     and the original list element is accessible as \code{x$data}.
#' @export
zipLapply <- function(iterover, fxn, ...) {
    lapply(zipData(iterover), fxn, ...)
}

#' A wrapper for lapply that allows access to list element names, with simplification at the end to an array.
#'
#' @param iterover List to iterate over. This list will be transformed by \code{zipData}.
#' @param fxn A function to apply to each element \code{x} of
#'     \code{zipData(iterover)}, where the name is accessible as \code{x$name}
#'     and the original list element is accessible as \code{x$data}.
#' @export
zipSapply <- function(iterover, fxn, ...) {
    simplify2array(zipLapply(iterover, fxn), ...)
}

#' Function to obtain HTML colors for a particular value.
#'
#' @param x A value or vector of values to colorize.
#' @param direction If 1, use the color scale as given; if -1, reverse it.
#' @param na_color Set NA elements to this color.
#' @param scale_from Instead of the minimum and maximum of \code{x}, scale
#'     values from this minimum and maximum (see \code{?rescale}).
#' @param colors A vector of two strings giving the low and high color,
#'     respectively.
#' @param limits A vector of two numbers giving the minimum and maximum value
#'     outside which values will be represented by the bottom or top of the
#'     color scale, respectively.
#' @return An HTML color.
#' @export
kable.recolor <- function(x,
                          direction = 1,
                          option = "D",
                          na_color = "#BBBBBB",
                          scale_from = NULL,
                          colors = c("#000000","#FFFFFF"),
                          limits = c(-Inf, Inf)
                          )
{
    if (any(na.omit(x) < limits[1])) {
        x[x < limits[1]] <- limits[1]
    }
    if (any(na.omit(x) > limits[2])) {
        x[x < limits[2]] <- limits[2]
    }
    if (is.null(scale_from)) {
        x <- round(scales::rescale(x, c(1, 256)))
    }
    else {
        x <- round(scales::rescale(x, to = c(1, 256), from = scale_from))
    }
    if (direction == -1) { x <- (257 - x) }
    cmap <- grDevices::colorRampPalette(colors)(256)
    color_code <- cmap[x]
    color_code[is.na(color_code)] <- na_color
    return(color_code)
}

#' Capitalize the first letter of a word/vector of words.
#'
#' @param s Word or vector of words.
#' @param strict If TRUE, also force non-initial letters to be lowercase.
#' @export
capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s, 1, 1)),
    {s <- substring(s, 2); if(strict) tolower(s) else s},
    sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
###

# Super gross XML hack follows to make SVGs of trees interactive

#' XML hack to make interactive tree diagrams.
#'
#' This hack is very ugly but works most of the time. However, it is a good idea
#' to wrap it in a tryCatch so that you can fall back to a less flashy
#' implementation, because it relies on editing a poorly-annotated SVG file as
#' if it were an XML document.
#'
#' @param tree.obj A ggtree representation of a tree.
#' @param file A filename where the final SVG output will be written.
#' @param stroke.scale Multiplier of stroke width in dendrogram.
#' @param pheno A vector with names corresponding to the tips of the tree and
#'     values corresponding to the phenotype value at that tip.
#' @param pheno.name The name of the phenotype being calculated (e.g.
#'     "prevalence").
#' @param native.tooltip Instead of using mouseover, use SVG tooltips (less
#'     powerful).
#' @param units Postfix for the values in \code{pheno} (e.g. "%" for
#'     percentages).
#' @export
hack.tree.labels <- function(tree.obj,
                             file,
                             stroke.scale = 0.7,
                             pheno = NULL,
                             pheno.name = NULL,
                             native.tooltip = FALSE,
                             units = "",
                             ...) {
    tip.labels <- with(tree.obj$data, label[isTip])
    xml <- svglite::xmlSVG(print(tree.obj), standalone = TRUE, ...)
    new.style.text <- " \n .faketip:hover ~ .realtip { \n stroke-width: 5; \n opacity: 1; \n  } \n .faketip:hover ~ .specieslabel { \n opacity: 1; \n } \n "
    style.nodes <- xml2::xml_find_all(xml, "//*[local-name()='style']")
    xml2::xml_set_text(style.nodes[1], new.style.text)
    xml.text <- xml2::xml_find_all(xml, "//*[local-name()='text']")
    xml.text.contents <- sapply(xml.text, xml2::xml_text)
    xml.label.indices <- which(xml.text.contents %in% tip.labels)
    xml.label.heights <- sapply(xml.text, function(x) {
        xml2::xml_attrs(x)["y"]
    })
    xml.label.pair <- cbind(label = xml.text.contents,
                            y = xml.label.heights)[xml.label.indices, ]
    ordered.labels <- xml.label.pair[order(xml.label.pair[, "y"] %>% as.numeric),
                                     "label"]
    xml.lines <- xml2::xml_find_all(xml, "//*[local-name()='line']")
    xml.line.props <- sapply(xml.lines, function(x) xml2::xml_attrs(x))
    xml.x2 <- xml.line.props["x2", ] %>% as.numeric
    xml.y2 <- xml.line.props["y2", ] %>% as.numeric
    # terminus <- max(xml.x2)
    uniq.x2 <- unique(xml.x2)
    count.x2 <- sapply(uniq.x2, function(x) sum(xml.x2 == x))
    terminus <- uniq.x2[which(count.x2 == length(tip.labels))]
    xml.y2.sorted <- sort(xml.y2[xml.x2 == terminus])
    # skootch over, remove tip, add title
    for (x in xml.text) {
        label <- xml2::xml_text(x)
        if (label %in% ordered.labels) {
            xml2::xml_set_attr(x, "x", as.character(terminus))
            xml2::xml_set_text(x, " ")
            if (native.tooltip) {
                xml2::xml_add_sibling(x, xml2::read_xml(paste0("<title>",
                                                               label,
                                                               "</title>")))
            }
        }
    }
    for (l in xml.lines) {
        l.attr <- xml2::xml_attrs(l)
        l.y2 <- as.numeric(l.attr["y2"])
        l.x2 <- as.numeric(l.attr["x2"])
        l.x1 <- as.numeric(l.attr["x1"])
        # This step is necessary because otherwise mouseover won't work
        style <- xml2::xml_attrs(l)["style"]
        s.parsed <- style.parse(style)
        for (n in 1:length(s.parsed)) {
            xml2::xml_set_attr(l, names(s.parsed)[n], s.parsed[n])
        }
        xml2::xml_set_attr(l, "style", "")
        if ("stroke-width" %in% names(s.parsed)) {
            xml2::xml_set_attr(l,
                         "stroke-width",
                         (as.numeric(s.parsed["stroke-width"]) * stroke.scale) %>%
                         as.character)
        }
        if (!("stroke" %in% names(s.parsed))) {
            xml2::xml_set_attr(l, "stroke", "#000000")
        }
        if (l.x2 == terminus) {
            label <- ordered.labels[which(xml.y2.sorted == l.y2)]
            xml2::xml_set_attr(l, "id", label)
            xml2::xml_set_attr(l, "class", "realtip")
            new.group <- xml2::read_xml("<g class=\"tip\"> </g>")
            l2 <- xml2::xml_add_child(new.group, l)
            xml2::xml_add_child(new.group, l)
            xml2::xml_set_attr(l2, "opacity", "0")
            xml2::xml_set_attr(l2, "pointer-events", "all")
            xml2::xml_set_attr(l2, "stroke-width", 5)
            xml2::xml_set_attr(l2, "class", "faketip")
            xml2::xml_set_attr(l2, "x1",
                         as.character(
                             l.x1 - 500
                         ))
            xml2::xml_set_attr(l2, "x2",
                         as.character(
                             terminus + 500
                         ))
            extra.info <-  ""
            if (!is.null(pheno)) {
                if (is.null(pheno.name)) pheno.name <- "phenotype"
                if (label %in% names(pheno)) {
                    phi <- format(pheno[label], digits = 3)
                } else {
                    phi <- "NA"
                }
                extra.info <- paste0("(", pheno.name, " = ", phi, units,  ")")
            }
            xml2::xml_add_child(new.group, xml2::read_xml(paste0(
              "<text x=\"",
              l.x2 + 5,
              "\" y = \"",
              l.y2 + 3,
              "\" opacity=\"0\" pointer-events=\"all\"" ,
              " style=\"font-family: Arial; font-size: 10px;",
              " fill: ",
              xml2::xml_attr(l, "stroke"),
              ";",
              "\" class=\"specieslabel\"> ",
              label,
              " ",
              extra.info,
              " ",
              "</text>")))
            if (native.tooltip) {
                xml2::xml_add_child(new.group, xml2::read_xml(paste0("<title>",
                                                                     label,
                                                                     "</title>")))
            }
            xml2::xml_replace(l, new.group)
        }
    }
    xml2::write_xml(x = xml, file)
}


#' Helper function to parse SVG styles.
#'
#' @param str A style string to parse.
#' @return A named vector of style attributes.
style.parse <- function(str) {
    semi.split <- strsplit(str, ";") %>% sapply(., trimws)
    if (is.null(dim(semi.split))) {
        c.split <- trimws(strsplit(semi.split, ":")[[1]])
        c.output <- c.split[2]
        names(c.output) <- c.split[1]
    } else {
        c.split <- apply(semi.split, 1, function(x)
            strsplit(x, ":")[[1]]) %>% trimws
        c.output <- c.split[2, ]
        names(c.output) <- c.split[1, ]
    }
    return(c.output)
}

#' A fall-back plotting option for when \code{hack.tree.labels} fails, designed
#' to produce the same kind of output.
#'
#' @param tree.obj A ggtree object.
#' @param file File to which an SVG representation of this tree object will be
#'     written.
#' @export
non.interactive.plot <- function(tree.obj, file) {
    warning(paste0("replotting to: ", file))
    non.int <- svglite::xmlSVG(print(tree.obj), standalone = TRUE)
    xml2::write_xml(x = non.int, file)
}

#' A wrapper around \code{apply} and \code{parApply} that allows them to be
#' called with a single syntax.
#'
#' @param mtx Matrix to apply a function over.
#' @param margin Margin of matrix to iterate over.
#' @param fun Function to apply over matrix.
#' @param cl If NULL, \code{apply} will be called; otherwise, should be the
#'     object returned by \code{makeCluster}.
#' @export
maybeParApply <- function(mtx, margin, fun, cl=NULL, ...) {
    if (!is.null(cl)) {
        parallel::parApply(cl, mtx, margin, fun, ...)
    } else {
        apply(mtx, margin, fun, ...)
    }
}

