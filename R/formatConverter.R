
#' Convert GInteractions to linkSet
#' @aliases Convert
#' @description Convert other data formats to linkSet. Currently supported: GInteractions, data.frame.
#' 
#' @param x A GInteractions object
#' @param specificCol A character string specifying the column to use for bait naming
#' @param ... Additional arguments (not used)
#' 
#' @rdname Convert
#' 
#' @return A linkSet object
#' @export
setMethod("Convert", signature(x = "GInteractions"), function(x, baitCol = NULL, ...) {
  anchor1 <- x@anchor1
  anchor2 <- x@anchor2
  regions <- x@regions
  metadata <- as.data.frame(x@elementMetadata)

  if (!is.null(baitCol) && !baitCol %in% colnames(metadata)) {
    warning("baitCol not found in metadata, using first regions as bait")
    baitCol <- NULL
  }

  if (is.null(baitCol)) {
    anchor1Region <- regions[anchor1]
    nameBait <- paste0(anchor1Region)
  } else {
    nameBait <- metadata[[baitCol]]
  }
  ls <- linkSet(
    anchor1 = regions[anchor1],
    anchor2 = regions[anchor2],
    specificCol =  nameBait
  )
  mcols(ls) <- metadata
  return(ls)
})

#' Convert GenomicInteractions to linkSet
#' @param x A GenomicInteractions object
#' @param baitCol A character string specifying the column to use for bait naming (optional)
#' 
#' @rdname Convert
#'
#' @importFrom GenomicInteractions GenomicInteractions
#' @return A linkSet object
#' @export
setMethod("Convert", signature(x = "GenomicInteractions"), function(x, baitCol = NULL, ...) {
  # Extract anchor1 and anchor2 GRanges
  anchor1 <- x@regions[x@anchor1]
  anchor2 <- x@regions[x@anchor2]
  
  # Extract metadata
  metadata <- as.data.frame(x@elementMetadata)
  
  # Create specificCol (using 'InteractionID' if available, otherwise use ranges)
  if ("InteractionID" %in% colnames(metadata)) {
    specificCol <- metadata$InteractionID
  } else {
    specificCol <- paste0(anchor1)
  }

  if (is.null(baitCol)) {
    nameBait <- paste0(anchor1)
  } else if (baitCol %in% colnames(metadata)) {
    nameBait <- metadata[[baitCol]]
  } else {
    warning("baitCol not found in metadata, using first regions as bait")
    nameBait <- paste0(anchor1)
  }

  # Create linkSet object
  ls <- linkSet(
    anchor1 = anchor1,
    anchor2 = anchor2,
    specificCol = nameBait
  )
  
  # Add metadata to linkSet
  mcols(ls) <- metadata
  
  return(ls)
})

.convert_to_grange <- function(intervals) {
  # convert "chr1.816066.816566" or "chr1:816066-816566" to grange format
  parts <- strsplit(intervals, "[.:\\-]")
  
  # Check if all parts have exactly 3 elements
  if (!all(sapply(parts, length) == 3)) {
    stop("Please input peak format like chr1.816066.816566 or chr1:816066-816566")
  }
  
  chromosomes <- sapply(parts, `[`, 1)
  starts <- as.numeric(sapply(parts, `[`, 2)) - 1  # Convert to 0-based
  ends <- as.numeric(sapply(parts, `[`, 3))
  
  # Create a GRanges object
  granges_obj <- GRanges(
    seqnames = Rle(chromosomes),
    ranges = IRanges(start = starts, end = ends)
  )
  return(granges_obj)
}


#' Convert data.frame to linkSet
#'
#' @param x A data.frame object
#' @param ... Additional arguments (not used)
#' @rdname Convert
#' @return A linkSet object
#' @export
setMethod("Convert", signature(x = "data.frame"), function(x, baitCol = "gene", oeCol = "peak", ...) {
  # Implementation for data.frame conversion
  # This is a placeholder and needs to be implemented based on your data.frame structure
  if (!baitCol %in% colnames(x) | !oeCol %in% colnames(x)) {
    stop("baitCol and oeCol must be columns in the data.frame")
  }
  bait <- x[[baitCol]]
  oe <- x[[oeCol]]
  oeGrange <- .convert_to_grange(oe)
  metadata <- x[ , !(colnames(x) %in% c(baitCol, oeCol)), drop = FALSE]
  ls <- linkSet(
    anchor1 = bait,
    anchor2 = oeGrange
  )
  mcols(ls) <- metadata
  return(ls)
})




#' Convert Pairs to linkSet
#' @param x A Pairs object
#' @param baitCol A character string specifying the column to use for bait naming
#' 
#' @rdname Convert
#'
#' @return A linkSet object
#' @export
setMethod("Convert", signature(x = "Pairs"), function(x,baitCol = NULL, ...) {
  # Extract first and second GRanges
  anchor1 <- x@first
  anchor2 <- x@second
  
  # Extract metadata
  metadata <- as.data.frame(x@elementMetadata)
  
  # Create specificCol (using 'symbol' if available, otherwise use ranges)
  if ("symbol" %in% colnames(metadata)) {
    specificCol <- metadata$symbol
  } else {
    specificCol <- paste0(anchor1)
  }

  if (is.null(baitCol)) {
    nameBait <- paste0(anchor1)
  } else {
    nameBait <- metadata[[baitCol]]
  }
  # Create linkSet object
  ls <- linkSet(
    anchor1 = anchor1,
    anchor2 = anchor2,
    specificCol = nameBait
  )
  
  # Add metadata to linkSet
  mcols(ls) <- metadata
  
  return(ls)
})

#' Default conversion method
#'
#' @param x An object of unsupported class
#' @param ... Additional arguments (not used)
#'
#' @rdname Convert
#' @return Nothing, throws an error
#' @export
setMethod("Convert", signature(x = "ANY"), function(x, ...) {
  stop(paste("Conversion from", class(x), "to linkSet is not supported"))
})




###############################################################
#' Convert GInteractions with bait range and oe ranges to linkSet
#' 
#' 
#'
#' @param gi A GInteractions object
#' @param geneGr A GRanges object representing genes
#' @param peakGr A GRanges object representing peaks
#' @param geneSymbol A character vector with same length as geneGr or column name in mcols(geneGr) for gene symbols
#'
#' @return A linkSet object
#' @export
setMethod("baitGInteractions", signature(x = "GInteractions", geneGr = "GRanges", peakGr = "GRanges"), function(x, geneGr, peakGr, geneSymbol=NULL) {
  gi <- x
  # validate geneSymbol in in mcols(geneGr)
  # validate length of geneSymbol is 1 or same length as geneGr
  if (length(geneSymbol) != 1 & length(geneSymbol) != length(geneGr)) {
    warning("geneSymbol must be a column name in mcols(geneGr) or a vector of the same length as geneGr, 
    using range to represent gene instead.")
    geneSymbol = NULL
  } else if (length(geneSymbol) == length(geneGr)) {
    geneSymbol <- geneSymbol
  } else if (!geneSymbol %in% colnames(mcols(geneGr))){
    warning("geneSymbol not found in geneGr, using range to represent gene instead.")
    geneSymbol = NULL
  } else {
    geneSymbol <- mcols(geneGr)[[geneSymbol]]
  }

  if (is.null(geneSymbol)) {
    geneSymbol <- paste0(geneGr)
  }
  out <- InteractionSet::linkOverlaps(gi, geneGr, peakGr)
  linkGene <- geneSymbol[out$subject1]
  linkGeneRange <- geneGr[out$subject1]
  linkPeak <- peakGr[out$subject2]
  linkSetObj <- linkSet(linkGeneRange, linkPeak, linkGene)
  
  # Add mcols from gi to linkSetObj
  mcols(linkSetObj) <- cbind(mcols(linkSetObj), mcols(gi)[out$query, ])
  
  return(linkSetObj)
})

