###############################################################
# Setting validity and show methods.

.check_inputs <- function(anchor1, anchor2, nameBait, regionBait, regionOE,same.length=TRUE) {
  if (!all(is.finite(anchor1)) || !all(is.finite(anchor2))) { 
    return("all anchor indices must be finite integers")
  }
  if (!all(anchor1 >= 1L) || !all(anchor2 >= 1L)) {
    return('all anchor indices must be positive integers')
  } 
  nregs1 <- length(nameBait)
  nregs2 <- length(regionOE)
  if ( !all(anchor1 <= nregs1) || !all(anchor2 <= nregs2)) {
    return("all anchor indices must refer to entries in 'regions'")
  } 
  if (same.length && length(anchor1)!=length(anchor2)) { 
    return("first and second anchor vectors have different lengths")
  }
  return(TRUE)
}

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
    ncr <- .safeNMcols(c(regions(x)[[1]],regions(x)[[2]]))
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

.new_LK <- function(anchor1, anchor2, nameBait, regionBait, regionOE, metadata) {
  elementMetadata <- make_zero_col_DFrame(length(anchor1))
  
  # Checking odds and ends.
  anchor1 <- as.integer(anchor1)
  anchor2 <- as.integer(anchor2)
  if (is.null(nameBait)){
    nameBait <- paste(gr1)
  }
  msg <- .check_inputs(anchor1, anchor2, nameBait, regionBait,regionOE)
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
      regionBait=regionBait,
      regionOE = regionOE,
      elementMetadata=elementMetadata,
      metadata=as.list(metadata))
}

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
            regionOE <- anchor2
            anchor2 <- seq_along(regionOE)
            regionBait <- GRanges()
            out <- .new_LK(anchor1=anchor1, anchor2=anchor2,  
                           nameBait=nameBait,regionBait=regionBait,
                           regionOE=regionOE, metadata=metadata)
            mcols(out) <- mcolBind
            out
            
          }
)
          


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
            
            regionBait <- anchor1
            anchor1 <- seq_along(regionBait)
            regionOE <- anchor2
            anchor2 <- seq_along(regionOE)
            
            out <- .new_LK(anchor1=anchor1, anchor2=anchor2,  
                           nameBait=nameBait,
                           regionBait = regionBait,
                           regionOE=regionOE, metadata=metadata)
            mcols(out) <- mcolBind
            out
          }
)

