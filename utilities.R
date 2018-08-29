# utilities

# helper functions

`%btwn%` <- function(x, y) { (x > min(y)) & (x < max(y)) }
`%btwn.inc%` <- function(x, y) { (x >= min(y)) & (x <= max(y)) }
`%intr%` <- function(x, y) intersect(x,y)
`%withnames%` <- function(x, y) { names(x) <- y; x }

nw <- function(x) { names(which(x)) }
min.nonzero <- function(x) min(x[x > 0])

count.each <- function(x, na.rm = FALSE) {
  u <- unique(x)
  simplify2array(lapply.across.names(u, function(y)
    sum(x == y, na.rm = na.rm)))
}

grepv <- function(x, y, ...) grep(x, y, ..., value = TRUE)

lapply.across.names <- function(X, FUN, ...) {
  r <- lapply(X, FUN, ...)
  names(r) <- X
  r
}
pblapply.across.names <- function(X, FUN, ...) {
  r <- pblapply(X, FUN, ...)
  names(r) <- X
  r
}

small.merge <- function(x, y, ...) {
  z <- merge(x, y, by = 0)
  r <- z[, -1]
  rownames(r) <- z[, 1]
  r
}

logit <- function(x) (log(x / (1 - x)))
logistic <- function(x) exp(x) / (1 + exp(x))

truncated <- function(x, lim = c(logit(0.001), logit(0.25))) {
  # truncate a vector by values
  y <- x
  y[which(x < lim[1])] <- lim[1]
  y[which(x > lim[2])] <- lim[2]
  y
}

fastread <- function(location, cn = TRUE) {
  # rownames are useful
  master <- data.frame(fread(location, header = T), check.names = cn)
  rn <- master[,1]
  rest <- master[,-1]
  rownames(rest) <- rn
  return(data.matrix(rest))
}
# Parallelization

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

gene.annot <- function(x) {
  merge(data.table(x),
        gene.to.fxn,
        by.x = colnames(data.table(x))[1],
        by.y = "gene",
        all = T)[x]$"function" %withnames% x
}


tax.annot <- function(tns, taxonomy) {
  taxonomy$species[match(tns, taxonomy$cluster)]
}

# from flodel @ stackoverflow
list.depth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, list.depth)), 0L)
first.element.at.depth <- function(l, n) {
  if (n == 1) { l[[1]] } else { (first.element.at.depth(l[[1]], n - 1)) }
}

annotate.nested <- function(nested,
  summarize = NULL,
  n = NULL,
  n.names = NULL,
  stop.at = 0, 
  summarize.values = 1) {

  nestedness <- list.depth(nested)
  if (nestedness == stop.at) {
    if (is.null(summarize)) {
      values <- data.frame(value = nested)
      rownames(values) <- names(nested)  
    } else {
      val.raw <- summarize(nested)
      if ((!is.null(nrow(val.raw))) || (nrow(val.raw) > 0)) {
        rownames(val.raw) <- NULL
      }
      values <- data.frame(value = val.raw)
    }
    if (is.null(nrow(values))) {
      n.mtx <- matrix(rep(n, (length(values) / summarize.values)),
        nc = length(n),
        byrow = TRUE)
    } else {
      n.mtx <- matrix(rep(n, nrow(values)), nc = length(n), byrow = TRUE)
    }
    if (!is.null(n.names)) { colnames(n.mtx) <- n.names }
    final.df <- cbind(names = rownames(values), n.mtx, value = values)
    rownames(final.df) <- NULL
    return(final.df)
  } else {
    if (is.null(names(nested))) { names(nested) <- 1:length(nested) }
    Reduce(rbind, lapply.across.names(names(nested),
      function(x) annotate.nested(nested[[x]],
        summarize = summarize,
        n = c(n, x),
        n.names = n.names,
        stop.at = stop.at)))
  }

}

zipData <- function(iterover) {
  n <- names(iterover)
  names(n) <- n
  return(lapply(n, function(x) {
    list(name = n[x], data = iterover[[n[x]]])
  }))
}

zipLapply <- function(iterover, fxn) {
  lapply(zipData(iterover), fxn)
}

zipSapply <- function(iterover, fxn) {
  simplify2array(zipLapply(iterover, fxn))
}
kable.recolor <- function(x, direction = 1, option = "D",
  na_color = "#BBBBBB", scale_from = NULL, colors = c("#000000","#FFFFFF"),
  limits = c(-Inf, Inf))
{
  if (any(na.omit(x) < limits[1])) {
    x[x < limits[1]] <- limits[1]
  }
  if (any(na.omit(x) > limits[2])) {
    x[x < limits[2]] <- limits[2] 
  }
  if (is.null(scale_from)) {
    x <- round(rescale(x, c(1, 256)))
  }
  else {
    x <- round(rescale(x, to = c(1, 256), from = scale_from))
  }
  if (direction == -1) { x <- (257 - x) }
  cmap <- colorRampPalette(colors)(256)
  color_code <- cmap[x]
  color_code[is.na(color_code)] <- na_color
  return(color_code)
}
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
    {s <- substring(s, 2); if(strict) tolower(s) else s},
    sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
###

# Super gross XML hack follows to make SVGs of trees interactive

hack.tree.labels <- function(tree.obj,
  file,
  stroke.scale = 0.7,
  pheno = NULL,
  pheno.name = NULL,
  native.tooltip = FALSE,
  units = "",
  ...) {

  tip.labels <- with(tree.obj$data, label[isTip])
  xml <- xmlSVG(print(tree.obj), standalone = TRUE, ...)
  new.style.text <- " \n .faketip:hover ~ .realtip { \n stroke-width: 5; \n opacity: 1; \n  } \n .faketip:hover ~ .specieslabel { \n opacity: 1; \n } \n " 
  style.nodes <- xml_find_all(xml, "//*[local-name()='style']")
  xml_set_text(style.nodes[1], new.style.text)
  xml.text <- xml_find_all(xml, "//*[local-name()='text']")
  xml.text.contents <- sapply(xml.text, xml_text)
  xml.label.indices <- which(xml.text.contents %in% tip.labels)
  xml.label.heights <- sapply(xml.text, function(x) {
    xml_attrs(x)["y"]
  })
  xml.label.pair <- cbind(label = xml.text.contents,
    y = xml.label.heights)[xml.label.indices, ]
  ordered.labels <- xml.label.pair[order(xml.label.pair[, "y"] %>% as.numeric), "label"]
  xml.lines <- xml_find_all(xml, "//*[local-name()='line']")
  xml.line.props <- sapply(xml.lines, function(x) xml_attrs(x))
  xml.x2 <- xml.line.props["x2", ] %>% as.numeric
  xml.y2 <- xml.line.props["y2", ] %>% as.numeric
  # terminus <- max(xml.x2)
  uniq.x2 <- unique(xml.x2)
  count.x2 <- sapply(uniq.x2, function(x) sum(xml.x2 == x))
  terminus <- uniq.x2[which(count.x2 == length(tip.labels))]
  xml.y2.sorted <- sort(xml.y2[xml.x2 == terminus])
  # skootch over, remove tip, add title
  for (x in xml.text) {
    label <- xml_text(x)
    if (label %in% ordered.labels) {
      xml_set_attr(x, "x", as.character(terminus))
      xml_set_text(x, " ")
      if (native.tooltip) {
        xml_add_sibling(x, read_xml(paste0("<title>",
              label,
              "</title>")))
      }
    }
  }
  for (l in xml.lines) {
    l.attr <- xml_attrs(l)
    l.y2 <- as.numeric(l.attr["y2"])
    l.x2 <- as.numeric(l.attr["x2"])
    l.x1 <- as.numeric(l.attr["x1"])
    # This step is necessary because otherwise mouseover won't work
    style <- xml_attrs(l)["style"]
    s.parsed <- style.parse(style)
    for (n in 1:length(s.parsed)) {
      xml_set_attr(l, names(s.parsed)[n], s.parsed[n])
    }
    xml_set_attr(l, "style", "")
    if ("stroke-width" %in% names(s.parsed)) {
      xml_set_attr(l, "stroke-width", (as.numeric(s.parsed["stroke-width"]) * stroke.scale) %>% as.character)
    }
    if (!("stroke" %in% names(s.parsed))) {
      xml_set_attr(l, "stroke", "#000000")
    }
    if (l.x2 == terminus) {
      label <- ordered.labels[which(xml.y2.sorted == l.y2)]
      xml_set_attr(l, "id", label)
      xml_set_attr(l, "class", "realtip")
      new.group <- read_xml("<g class=\"tip\"> </g>")
      l2 <- xml_add_child(new.group, l)
      xml_add_child(new.group, l)
      xml_set_attr(l2, "opacity", "0")
      xml_set_attr(l2, "pointer-events", "all")
      xml_set_attr(l2, "stroke-width", 5)
      xml_set_attr(l2, "class", "faketip")
      xml_set_attr(l2, "x1",
        as.character(
          l.x1 - 100
          ))
      xml_set_attr(l2, "x2",
        as.character(
          terminus + 300
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
      xml_add_child(new.group, read_xml(paste0(
            "<text x=\"",
            l.x2 + 5,
            "\" y = \"",
            l.y2 + 3, 
            "\" opacity=\"0\" pointer-events=\"all\" style=\"font-family: Arial; font-size: 10px;",
            " fill: ",
            xml_attr(l, "stroke"),
            ";",
            "\" class=\"specieslabel\"> ",
            label, 
            " ",
            extra.info,
            " ",
            "</text>")))
      if (native.tooltip) {
        xml_add_child(new.group, read_xml(paste0("<title>",
              label,
              "</title>")))
      }
      xml_replace(l, new.group)
    }
  }
  write_xml(x = xml, file) 
}

style.parse <- function(str) {
  semi.split <- strsplit(str, ";") %>% sapply(., trimws)
  if (is.null(dim(semi.split))) {
    c.split <- trimws(strsplit(semi.split, ":")[[1]])
    c.output <- c.split[2]
    names(c.output) <- c.split[1]
  } else {
    c.split <- apply(semi.split, 1, function(x) strsplit(x, ":")[[1]]) %>% trimws
    c.output <- c.split[2, ]
    names(c.output) <- c.split[1, ]
  }
  return(c.output)
}
non.interactive.plot <- function(tree.obj, file) {
  warning(paste0("replotting to: ", file))
  non.int <- xmlSVG(print(tree.obj), standalone = TRUE)
  write_xml(x = non.int, file)
}

generic.make.tables <- function(enr, depth = 3) {
  annotate.nested(enr,
                  stop.at = list.depth(enr) - depth,
                  summarize = function(x) {
                    if (is.null(nrow(x$table))) {
                      names(x$table) <- c("enriched", "V2")
                      data.frame(t(x$table))
                    } else if (nrow(x$table) > 0) { 
                      x$table 
                    } else {
                      rbind(x$table, c(enriched = NA, V2 = NA))
                    }}
  )
}

maybeParApply <- function(mtx, margin, fun, cl=NULL, ...) {
  if (!is.null(cl)) {
    parApply(cl, mtx, margin, fun, ...)
  } else {
    apply(mtx, margin, fun, ...)
  }
}


sparseMelt <- function(mtx) {
  mtxT <- as(mtx, "TsparseMatrix")
  df <- data.frame(row=mtxT@Dimnames[[1]][mtxT@i + 1],
    col=mtxT@Dimnames[[2]][mtxT@j + 1],
    value=mtxT@x)
  if (!is.null(names(dimnames(mtxT)))) {
    colnames(df)[1:2] <- names(dimnames(mtxT))
  }
  df
}
