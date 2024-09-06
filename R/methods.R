###############################################################
# Setting validity and show methods.

.check_inputs <- function(anchor1, anchor2, nameBait, regions,same.length=TRUE) {
  if (!all(is.finite(anchor1)) || !all(is.finite(anchor2))) {
    return("all anchor indices must be finite integers")
  }
  if (!all(anchor1 >= 1L) || !all(anchor2 >= 1L)) {
    return('all anchor indices must be positive integers')
  }
  nregs1 <- length(nameBait)
  nregs2 <- length(regions)
  if (!all(anchor1 <= nregs2) || !all(anchor2 <= nregs2)) {
    return("all anchor indices must refer to entries in 'regions'")
  }
  if (same.length && length(nameBait)!=length(anchor2)) {
    return("first and second anchor vectors have different lengths")
  }
  return(TRUE)
}

setValidity2("linkSet", function(object) {
  if (is.unsorted(regions(object))) { # Don't move into .check_inputs, as resorting comes after checking validity in various methods.
    return("'regions' should be sorted")
  }
  msg <- .check_inputs(anchor1(object), anchor2(object), bait(object),regions(object))
  if (is.character(msg)) { return(msg) }

  ### Length of anchors versus object is automatically checked by 'parallel_slot_names.'

  if (!is.null(names(object))) {
    if (length(names(object))!=length(object)) {
      stop("'NAMES' must be NULL or have length equal to that of the object")
    }
  }

  return(TRUE)
})

#' @export
setMethod("show", "linkSet", function(object){
  showLinkSet(object, margin="  ", print.seqinfo=TRUE, print.classinfo=TRUE)
})

showLinkSet <- function(x, margin="", print.seqinfo=FALSE, print.classinfo=FALSE) {
  lx <- length(x)
  nr <- lx
  nc <- .safeNMcols(x)
  cat(class(x), " object with ",
      lx, " ", ifelse(lx == 1L, "interaction", "interactions"), " and ",
      nc, " metadata ", ifelse(nc == 1L, "column", "columns"),
      ":\n", sep="")
  out <- makePrettyMatrixForCompactPrinting(x, .makeNakedMatFromGInteractions)

  # Ripped from GenomicRanges:::showGenomicRanges (with some mods).
  if (print.classinfo) {
    .COL2CLASS <- c(bait = "character", "   "="", seqnames_oe="Rle", ranges_oe="IRanges")
    extraColumnNames <- GenomicRanges:::extraColumnSlotNames(x)
    .COL2CLASS <- c(.COL2CLASS, getSlots(class(x))[extraColumnNames])
    classinfo <- makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
    classinfo[,"   "] <- ""
    stopifnot(identical(colnames(classinfo), colnames(out)))
    out <- rbind(classinfo, out)
  }

  if (nrow(out) != 0L) {
    rownames(out) <- paste0(margin, rownames(out))
  }
  print(out, quote=FALSE, right=TRUE, max=length(out))
  if (print.seqinfo) {
    cat(margin, "-------\n", sep="")
    ncr <-  .safeNMcols(regions(x))
    cat(margin, "regions: ", nr, " ranges and ", ncr, " metadata ", ifelse(ncr==1L, "column", "columns"), "\n", sep="")
    cat(margin, "seqinfo: ", summary(seqinfo(x)), "\n", sep="")
  }
}


.safeNMcols <- function(x) {
  #return column number safely
  nc <- ncol(mcols(x))
  if (is.null(nc)) { nc <- 0L }
  return(nc)
}

.makeNakedMatFromGInteractions <- function(x) {
  lx <- length(x)
  nc <- .safeNMcols(x)
  ans <- cbind(.pasteAnchor(anchors(x, type="bait"), append="bait"),
               "   "=rep.int("---", lx),
               .pasteAnchor(anchors(x, type="oe"), append="oe"))
  if (nc > 0L) {
    tmp <- do.call(data.frame, c(lapply(mcols(x), showAsCell), list(check.names=FALSE)))
    ans <- cbind(ans, `|`=rep.int("|", lx), as.matrix(tmp))
  }
  ans
}


.pasteAnchor <- function(x, append) {
  if(is.character(x)){
    out <- as.matrix(x)
    colnames(out) <- "bait"
  } else{
    out <- cbind(as.character(seqnames(x)), showAsCell(ranges(x)))
    colnames(out) <- paste0(c("seqnames", "ranges"),"_", append)
  }
  out
}
###############################################################
# Constructors

.resort_regions <- function(anchor1, anchor2, regions) {
  if (is.unsorted(regions)) {
    o <- order(regions)
    new.pos <- seq_along(o)
    new.pos[o] <- new.pos
    anchor1 <- new.pos[anchor1]
    anchor2 <- new.pos[anchor2]
    regions <- regions[o]
  }
  return(list(anchor1=anchor1, anchor2=anchor2, regions=regions))
}


.new_LK <- function(anchor1, anchor2, nameBait, regions, metadata) {
  elementMetadata <- make_zero_col_DFrame(length(anchor1))

  # Checking odds and ends.
  anchor1 <- as.integer(anchor1)
  anchor2 <- as.integer(anchor2)
  if (is.null(nameBait)){
    nameBait <- paste(gr1)
  }
  msg <- .check_inputs(anchor1, anchor2, nameBait, regions)
  if (is.character(msg)) { stop(msg) }

  #out <- .resort_regions(anchor1, anchor2, regions)
  # anchor1 <- out$anchor1
  # anchor2 <- out$anchor2
  # regions <- out$regions

  cls <- "linkSet"

  new(cls,
      anchor1=anchor1,
      anchor2=anchor2,
      nameBait=nameBait,
      regions=regions,
      elementMetadata=elementMetadata,
      metadata=as.list(metadata))
}

#' @export
setMethod("linkSet", c("character", "GRanges","character_Or_missing"),
          function(anchor1, anchor2, specificCol,metadata=list(),  ...) {
            mcol2 <- mcols(anchor2)
            mcols(anchor2) <- NULL
            colnames(mcol2) <- sprintf("anchor2.%s", colnames(mcol2))
            extraCols <- DataFrame(...)
            if (ncol(extraCols) == 0L) {
              extraCols <- make_zero_col_DFrame(length(anchor1))
            }
            mcolBind <- cbind(extraCols, mcol2)
            nameBait <- anchor1
            anchor1 <- seq_along(nameBait)
            regions <- anchor2
            anchor2 <- seq_along(regions)
            out <- .new_LK(anchor1=anchor1, anchor2=anchor2,
                           nameBait=nameBait,regions=regions,
                           metadata=metadata)
            mcols(out) <- mcolBind
            out
          }
)

.collate_GRanges <- function(...) {
  incoming <- list(...)
  obj.dex <- rep(factor(seq_along(incoming)), lengths(incoming))
  combined <- do.call(c, incoming)

  # Sorting and re-indexing.
  o <- order(combined)
  refdex <- integer(length(o))
  refdex[o] <- seq_along(combined)
  combined <- combined[o]

  # Removing duplicates and re-indexing.
  is.first <- !duplicated(combined)
  new.pos <- cumsum(is.first)
  combined <- combined[is.first]
  refdex <- new.pos[refdex]
  return(list(indices=split(refdex, obj.dex), ranges=combined))
}

#' @export
setMethod("linkSet", c("GRanges", "GRanges","character_Or_missing"),
          function(anchor1, anchor2, specificCol,metadata=list(),  ...) {
            # Stripping metadata and putting it somewhere else.
            mcol1 <- mcols(anchor1)
            mcols(anchor1) <- NULL
            colnames(mcol1) <- sprintf("anchor1.%s", colnames(mcol1))
            mcol2 <- mcols(anchor2)
            mcols(anchor2) <- NULL
            colnames(mcol2) <- sprintf("anchor2.%s", colnames(mcol2))

            # Additional Interaction-specific metadata
            extraCols <- DataFrame(...)
            if (ncol(extraCols) == 0L) {
              extraCols <- make_zero_col_DFrame(length(anchor1))
            }

            mcolBind <- cbind(extraCols, mcol1, mcol2)

            if (!missing(specificCol)){
              specificColName <- paste0("anchor1.",specificCol)
              if (specificColName %in% colnames(mcolBind)){
                nameBait <- mcolBind[specificColName]
              } else{
                warning(paste0("Can't find ", specificCol, "in metadata............"))
                nameBait <- paste(anchor1)
              }
            } else{
              nameBait <- paste(anchor1)
            }


            collated <- .collate_GRanges(anchor1, anchor2)
            regions <- collated$ranges
            anchor1 <- collated$indices[[1]]
            anchor2 <- collated$indices[[2]]

            # regionBait <- anchor1
            # anchor1 <- seq_along(regionBait)
            # regionOE <- anchor2
            # anchor2 <- seq_along(regionOE)

            out <- .new_LK(anchor1=anchor1, anchor2=anchor2,
                           nameBait=nameBait,
                           region= regions,
                           metadata=metadata)
            mcols(out) <- mcolBind
            out
          }
)

