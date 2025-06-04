
#--- Functions rewritten from clinUtils to lighten dependencies ---#

# See documentation for clinUtils. Incorporated here under the MIT license

#' See clinUtils::knitPrintListPlots.
#'
#' @export
knitPrintListPlots <- function(
        plotsList,
        generalLabel = "plotsList",
        type = c("ggplot2", "plotly"), ...) 
{
    type <- match.arg(type)
    knitPrintListObjects(xList = plotsList, generalLabel = generalLabel, 
                         printObject = (type == "ggplot2"), ...)
}

#' See clinUtils::knitPrintListPlots.
#'
#' @export
knitPrintListObjects <- function(
        xList,
        generalLabel = "objectsList",
        labels = paste0(generalLabel, seq_along(xList)),
        titles = NULL,
        titleLevel = 2,
        printObject = FALSE, 
        ...) 
{
    labels <- paste0("'", labels, "'")
    argsChunk <- list(...)
    if (length(argsChunk) > 0) {
        idxArgsChunkCharac <- which(sapply(argsChunk, is.character))
        if (length(idxArgsChunkCharac) > 0) {
            argsChunk[idxArgsChunkCharac] <- lapply(
                argsChunk[idxArgsChunkCharac], 
                function(x) {
                    paste0("'", x, "'")
                })
        }
        argsChunkTxt <- paste0(names(argsChunk), "=", "{{", names(argsChunk), 
                               "}}")
    }
    else argsChunkTxt <- NULL
    chunkTemplate <- paste0(
        "```{r {{label}}",
        if (!is.null(argsChunkTxt)) paste0(", ", toString(argsChunkTxt)),
        ", results = 'asis', echo = FALSE}\n", 
        if (!is.null(titles)) 
            paste0("cat(\"\\n\", paste(rep(\"#\", titleLevel), ",
                   "collapse = \"\"), \" {{title}}\\n\", sep = \"\")\n"), 
        if (printObject) 
            "print(", "xList[[{{i}}]]", if (printObject) 
                ")", "\n", "```\n")
    argsKnitExpand <- c(list(FUN = knitr::knit_expand,
                             text = chunkTemplate, 
                             i = seq_along(xList),
                             label = labels),
                        if (!is.null(titles)) list(title = titles), 
                        argsChunk)
    chunkTxt <- do.call(mapply, argsKnitExpand)
    cat(knitr::knit(text = paste(chunkTxt, collapse = "\n"), quiet = TRUE))
}

#' See clinUtils::getClinDT.
#' 
#' @export
getClinDT <- function (data,
                       nonVisibleVar = NULL,
                       nonVisible = NULL,
                       percVar = NULL,
                       barVar = NULL,
                       barColorThr = NULL,
                       barRange = NULL,
                       filter = "top",
                       searchBox = FALSE,
                       pageLength,
                       fixedColumns = NULL,
                       columnsWidth = NULL,
                       options = list(),
                       expandVar = NULL,
                       expandIdx = NULL,
                       escape = TRUE,
                       rowGroup = NULL,
                       rowGroupVar = NULL,
                       vAlign = "top",
                       callback = NULL,
                       buttons = getClinDTButtons(),
                       scrollX = TRUE,
                       file = NULL,
                       verbose = TRUE,
                       ...)
{
    
    extraArgs <- list(...)
    isSharedData <- inherits(x = data, what = "SharedData")
    dataContent <- if (isSharedData) {
        data$origData()
    } else data
    if (inherits(dataContent, "tbl_df")) {
        dataContent <- as.data.frame(dataContent)
    }
    if (!inherits(dataContent, c("data.frame", "matrix"))) 
        stop("'data' should be a data.frame, a matrix, a tibble or a SharedData object.")
    colnames <- extraArgs$colnames
    if (!is.null(colnames)) {
        colnames <- colnames[colnames %in% colnames(dataContent)]
        if (length(colnames) == 0) {
            colnames <- NULL
            warning("'colnames' doesn't contain labels for columns in data. ", 
                    "Are you sure you have specified it correctly (c([newName] = [oldName], ...)?")
        }
        extraArgs$colnames <- colnames
    }
    if (!is.null(nonVisible)) 
        warning("'nonVisible' is deprecated, please use: 'nonVisibleVar' instead.")
    nonVisibleVar <- checkVarInData(var = nonVisibleVar, data = dataContent, 
                                    label = "non-visible")
    if (!is.null(nonVisibleVar)) {
        if (!is.null(nonVisible)) 
            warning("'nonVisible' or 'nonVisibleVar' should be specified, 'nonVisibleVar' is used")
        nonVisible <- match(nonVisibleVar, colnames(dataContent)) - 
            1
    }
    if (missing(pageLength)) {
        pageLength <- ifelse(nrow(dataContent) <= 10, Inf, 10)
    }
    if (!is.null(rowGroup)) {
        warning("'rowGroup' is deprecated, please use: 'rowGroupVar' instead.")
        rowGroupVar <- rowGroup
    }
    rowGroupVar <- checkVarInData(var = rowGroupVar, data = dataContent, 
                                  label = "row group")
    if (!is.null(rowGroupVar)) {
        rowGroup <- match(rowGroupVar, colnames(dataContent)) - 
            1
        if (length(rowGroup) == 0) 
            rowGroup <- NULL
    }
    else rowGroup <- NULL
    if (is.logical(escape)) {
        if (length(escape) != 1) {
            stop("If escape is logical, it should be of length 1.")
        }
        else {
            if (escape) {
                escape <- seq(from = 1, to = ncol(dataContent))
            }
            else {
                escape <- numeric()
            }
        }
    }
    else if (is.numeric(escape)) {
        idxEscNotInData <- escape[!abs(escape) %in% seq_len(ncol(dataContent))]
        if (length(idxEscNotInData) > 0) {
            stop("'Escape' contains columns not in data: ", toString(idxEscNotInData), 
                 ".")
        }
        if (any(escape < 0)) {
            if (!all(escape < 0)) 
                stop("If 'escape' contains negative elements, they should all be negative.")
            escape <- setdiff(seq(from = 1, to = ncol(dataContent)), 
                              -escape)
        }
    }
    if (!is.null(rowGroup)) 
        nonVisible <- union(nonVisible, rowGroup)
    idxControl <- NULL
    expandVar <- checkVarInData(var = expandVar, data = dataContent, 
                                label = "expandable")
    isExpandIdxWrong <- !is.null(expandIdx) && ((!is.matrix(expandIdx)) || 
                                                    !all(c("row", "col") %in% colnames(expandIdx)))
    if (isExpandIdxWrong) {
        stop("'expandIdx' should be a matrix with columns: ", 
             "'row' and 'col'.")
    }
    if (!is.null(expandVar) | !is.null(expandIdx)) {
        if (!is.null(expandIdx)) {
            idxExpandVar <- unique(expandIdx[, "col"])
            for (iCol in seq_along(idxExpandVar)) {
                idxCol <- idxExpandVar[iCol]
                idxColNew <- idxCol + iCol - 1
                expandIdxCol <- expandIdx[which(expandIdx[, "col"] %in% 
                                                    idxCol), , drop = FALSE]
                expandIdxCol[, "col"] <- idxColNew
                expandRow <- rep(NA_character_, nrow(dataContent))
                expandRow[expandIdxCol[, "row"]] <- dataContent[expandIdxCol]
                dataContent[expandIdxCol] <- "&oplus;"
                idxBefore <- seq_len(idxColNew)
                idxAfter <- setdiff(seq_len(ncol(dataContent)), 
                                    idxBefore)
                dataContent <- cbind(dataContent[, idxBefore, 
                                                 drop = FALSE], expandRow = expandRow, dataContent[, 
                                                                                                   idxAfter, drop = FALSE])
            }
            newIdxForExpandVar <- idxExpandVar + seq_along(idxExpandVar) - 
                1
            getCol <- function(x) {
                x
            }
            body(getCol) <- bquote({
                xNew <- sapply(x, function(xI) {
                    idxDiff <- xI - .(idxExpandVar)
                    idxDiff <- idxDiff[idxDiff > 0]
                    ifelse(length(idxDiff) > 0, xI + which.min(idxDiff), 
                           xI)
                })
                return(xNew)
            })
            getColFormatStyle <- function(x) {
                x
            }
            body(getColFormatStyle) <- bquote(.(getCol)(x) - 
                                                  1)
            idxControl <- getCol(idxExpandVar) - 1
            escapeExpand <- getCol(idxExpandVar)
            nonVisibleExpand <- getCol(idxExpandVar)
            expandJS <- paste0("' + d[iCol + 1]+ '")
            callback <- JS(paste0("\n\t\t\t\t\ttable.column(1).nodes().to$().css({cursor: 'pointer'});\n\t\t\t\t\tvar format = function(d, iCol) {\n\t\t\t\t\t\treturn '<div>", 
                                  expandJS, "</div>';\n\t\t\t\t\t};\n\t\t\t\t\ttable.on('click', 'td.details-control', function() {\n\t\t\t\t\t\tvar td = $(this), row = table.row(td.closest('tr')), iCol = td[0]._DT_CellIndex['column'];\n\t\t\t\t\t\tif (row.child.isShown()) {\n\t\t\t\t\t\t\trow.child.hide();\n\t\t\t\t\t\t\ttd.html('&oplus;');\n\t\t\t\t\t\t} else {\n\t\t\t\t\t\t\toldVal = format(row.data(), iCol-1);\n\t\t\t\t\t\t\tif(oldVal === '<div>&oplus;</div>'){\n\t\t\t\t\t\t\t\trow.child(format(row.data(), iCol)).show();\n\t\t\t\t\t\t\t\ttd.html('&CircleMinus;');\n\t\t\t\t\t\t\t}\n\t\t\t\t\t\t}\n\t\t\t\t\t});"), 
                           callback)
        }
        else if (!is.null(expandVar)) {
            idxExpandVar <- which(colnames(dataContent) %in% 
                                      expandVar)
            getCol <- function(x) return(x + 1)
            getColFormatStyle <- function(x) return(x)
            dataContent <- cbind(` ` = "&oplus;", dataContent)
            idxControl <- 0
            escapeExpand <- 1
            nonVisibleExpand <- idxExpandVar
            expandJS <- paste(sapply(idxExpandVar, function(i) {
                labelI <- colnames(dataContent)[getCol(i)]
                if (!is.null(colnames)) {
                    labelCNI <- names(colnames)[match(labelI, colnames)]
                    if (!is.na(labelCNI)) 
                        labelI <- labelCNI
                }
                paste0(labelI, ": ' + d[", i, "] + '")
            }), collapse = "<br>")
            callback <- JS(paste0("\n\t\t\t\t\ttable.column(1).nodes().to$().css({cursor: 'pointer'});\n\t\t\t\t\tvar format = function(d) {\n\t\t\t\t\t\treturn '<div>", 
                                  expandJS, "</div>';\n\t\t\t\t\t};\n\t\t\t\t\ttable.on('click', 'td.details-control', function() {\n\t\t\t\t\t\tvar td = $(this), row = table.row(td.closest('tr'));\n\t\t\t\t\t\tif (row.child.isShown()) {\n\t\t\t\t\t\t\trow.child.hide();\n\t\t\t\t\t\t\ttd.html('&oplus;');\n\t\t\t\t\t\t} else {\n\t\t\t\t\t\t\trow.child(format(row.data())).show();\n\t\t\t\t\t\t\ttd.html('&CircleMinus;');\n\t\t\t\t\t\t}\n\t\t\t\t\t});"), 
                           callback)
        }
        escape <- setdiff(getCol(escape), escapeExpand)
        nonVisible <- union(getCol(nonVisible), nonVisibleExpand)
    }
    else {
        getColFormatStyle <- getCol <- function(x) return(x)
        callback <- callback
    }
    if (any(nonVisible >= ncol(dataContent))) 
        stop(paste("'nonVisible' should contain indices of columns within data (< ncol(data)).", 
                   "Are you sure you are using Javascript indexing", 
                   "(0 for first column, 1 for second column and so on)?"))
    if (!is.null(options$columnDefs)) {
        options$columnDefs <- sapply(options$columnDefs, function(x) {
            if (is.list(x) && "targets" %in% names(x)) {
                x[["targets"]] <- getCol(x[["targets"]])
            }
            x
        }, simplify = FALSE)
    }
    columnDefs <- c(options$columnDefs, if (!is.null(columnsWidth)) {
        list({
            columnsWidths <- rep(columnsWidth, length.out = ncol(dataContent))
            lapply(seq_along(columnsWidths), function(i) list(targets = getCol(i), 
                                                              columnsWidth = columnsWidths[i]))
        })
    }, if (!is.null(nonVisible)) list(list(targets = nonVisible, 
                                           visible = FALSE, className = "noVis")), if (!is.null(idxControl)) columnDefs <- list(list(orderable = FALSE, 
                                                                                                                                     className = "details-control", targets = idxControl)))
    isOptionAvailable <- function(options, label) {
        isOptionAvailable <- !label %in% names(options)
        if (!isOptionAvailable & verbose) {
            message("The", sQuote(label), " specified in 'options' overwrites the default.")
        }
        return(isOptionAvailable)
    }
    if (isOptionAvailable(options, "dom")) {
        domDefault <- paste0(if (length(buttons) > 0) 
            "B", if (pageLength < Inf) 
                "l", if (searchBox) 
                    "f", "rt", if (pageLength < Inf) 
                        "ip")
        options[["dom"]] <- domDefault
    }
    if (!is.null(fixedColumns)) {
        idx <- which(names(fixedColumns) %in% c("leftColumns", 
                                                "rightColumns"))
        if (length(idx) > 0) 
            fixedColumns[idx] <- sapply(fixedColumns[idx], getCol, 
                                        simplify = FALSE)
        if (isOptionAvailable(options, "fixedColumns")) {
            options[["fixedColumns"]] <- fixedColumns
        }
    }
    if (isOptionAvailable(options, "fixedHeader")) {
        options[["fixedHeader"]] <- if (is.null(fixedColumns)) 
            TRUE
        else FALSE
    }
    if (isOptionAvailable(options, "buttons")) {
        options[["buttons"]] <- buttons
    }
    if (isOptionAvailable(options, "searching")) {
        options[["searching"]] <- TRUE
    }
    if (isOptionAvailable(options, "scrollX")) {
        options[["scrollX"]] <- scrollX
    }
    if (isOptionAvailable(options, "autoWidth")) {
        options[["autoWidth"]] <- (!is.null(columnsWidth))
    }
    if (isOptionAvailable(options, "pageLength")) {
        options[["pageLength"]] <- ifelse(pageLength == Inf, 
                                          nrow(dataContent), pageLength)
    }
    if (length(rowGroup) > 0 && isOptionAvailable(options, "rowGroup")) {
        rowGroup <- getCol(rowGroup)
        options[["rowGroup"]] <- list(dataSrc = rowGroup)
        columnDefs <- c(columnDefs, list(list(targets = rowGroup, 
                                              className = "rowGroup")))
    }
    if (length(columnDefs) > 0) {
        options[["columnDefs"]] <- columnDefs
    }
    if (length(options) == 0) 
        options <- NULL
    extensions <- c(if (!is.null(rowGroup)) "RowGroup", if (length(buttons) > 
                                                            0) "Buttons", if (!is.null(fixedColumns)) c("FixedColumns", 
                                                                                                        "Scroller"), if (is.null(fixedColumns)) "FixedHeader")
    dataDT <- if (isSharedData) {
        if (nrow(dataContent) != length(data$key())) 
            stop("Key vector is of different length than the number of records in the data.")
        keySD <- data$.__enclos_env__$private$.key
        crosstalk::SharedData$new(data = dataContent, key = keySD, group = data$groupName())
    }
    else dataContent
    argsDT <- list(data = dataDT, rownames = FALSE, filter = filter, 
                   extensions = extensions, options = options, escape = escape)
    if (!is.null(callback)) 
        argsDT <- c(argsDT, list(callback = callback))
    extraArgsSpec <- intersect(names(extraArgs), names(argsDT))
    if (length(extraArgsSpec) > 0) {
        warning(paste("Extra parameter(s)", toString(sQuote(extraArgsSpec)), 
                      "are ignored because some internal defaults are set for these parameters."))
        extraArgs <- extraArgs[setdiff(names(extraArgs), extraArgsSpec)]
    }
    argsDT <- c(argsDT, extraArgs)
    tableDT <- do.call(DT::datatable, argsDT)
    if (!is.null(percVar)) 
        tableDT <- DT::formatPercentage(tableDT, columns = percVar, 
                                        digits = 2)
    tableDT <- formatDTBarVar(tableDT = tableDT, data = dataContent, 
                              barVar = barVar, barColorThr = barColorThr, barRange = barRange, 
                              getCol = getColFormatStyle)
    if (!is.null(vAlign)) {
        tableDT <- tableDT %>%
            DT::formatStyle(columns = seq_len(ncol(dataContent)), 
                            `vertical-align` = vAlign)
    }
    if (!is.null(file)) {
        if (file_ext(file) != "html") 
            stop("'file' should be of extension 'html'.")
        wdInit <- getwd()
        on.exit(setwd(wdInit))
        setwd(dirname(file))
        htmlwidgets::saveWidget(widget = tableDT, file = basename(file))
    }
    return(tableDT)
}


formatDTBarVar <- function (tableDT, data, barVar = NULL, barColorThr = NULL, barRange = NULL, 
                            getCol = function(x) x) 
{
    barVar <- checkVarInData(var = barVar, data = data, label = "bar")
    if (!is.null(barVar)) {
        barVarNotNum <- barVar[!sapply(data[, barVar, drop = FALSE], 
                                       is.numeric)]
        if (length(barVarNotNum) > 0) {
            warning(paste(toString(barVarNotNum), "variable(s)", 
                          "not represented as bar because they are not numeric."))
            barVar <- setdiff(barVar, barVarNotNum)
        }
        getElFromList <- function(param, var) {
            if (!is.null(param)) {
                if (!is.null(names(param))) {
                    if (var %in% names(param)) {
                        param[[var]]
                    }
                }
                else param
            }
        }
        for (var in barVar) {
            idxVar <- getCol(match(var, colnames(data)))
            barColorThrVar <- getElFromList(param = barColorThr, 
                                            var = var)
            barRangeVar <- getElFromList(param = barRange, var = var)
            if (is.null(barRangeVar)) 
                barRangeVar <- range(as.numeric(data[, var]), 
                                     na.rm = TRUE)
            barRangeVar[1] <- barRangeVar[1] - diff(barRangeVar) * 
                0.01
            barColor <- if (!is.null(barColorThrVar)) {
                styleInterval(cuts = barColorThrVar, values = viridis(length(barColorThrVar) + 
                                                                          1))
            }
            else "black"
            barBg <- styleColorBar(data = barRangeVar, color = "green")
            tableDT <- tableDT %>% formatStyle(columns = idxVar, 
                                               color = barColor, background = barBg)
        }
    }
    return(tableDT)
}

checkVarInData <- function (var, data, label) 
{
    varNotInData <- setdiff(var, colnames(data))
    if (length(varNotInData) > 0) 
        warning(paste(label, "variable(s):", sQuote(toString(varNotInData)), 
                      "not used because not available in the data."), call. = FALSE)
    var <- intersect(var, colnames(data))
    if (length(var) == 0) 
        var <- NULL
    return(var)
}

getClinDTButtons <- function(type = c("copy", "csv", "excel", "pdf", "print"),
                             typeExtra = NULL,
                             opts = NULL) {
    type <- unique(c(type, typeExtra))
    type <- match.arg(type, choices = c("copy", "csv", "excel", 
                                        "pdf", "print", "colvis"),
                      several.ok = TRUE)
    getExportButton <- function(typeBtn, ...) {
        if (typeBtn %in% type) {
            c(list(extend = typeBtn, ...,
                   exportOptions = list(columns = list(".rowGroup", 
                                                       ":visible"))),
              opts[[typeBtn]], ...)
        }
    }
    buttons <- list(getExportButton(typeBtn = "copy"),
                    getExportButton(typeBtn = "csv"), 
                    getExportButton(typeBtn = "excel"),
                    getExportButton(typeBtn = "pdf"), 
                    getExportButton(typeBtn = "print"),
                    if ("colvis" %in% type) c(
                        list(extend = "colvis",
                             columns = ":not(.noVis)", 
                             text = "Show/hide columns"),
                        opts[["colvis"]]))
    buttons <- buttons[!sapply(buttons, is.null)]
    return(buttons)
}