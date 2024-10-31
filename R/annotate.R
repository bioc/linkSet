# Annotation--------------

#' Annotate the link set with txDb. Give a gene list, and return a
#'
#' @aliases annotatePromoter
#' @param x linkSet
#' @param genome the genome you want to annotate
#' @param keyType the key type. You can check with AnnotationDbi::keytypes
#' @param upstream The upstream base from the gene
#' @param overwrite Whether to overwrite the regionsBait if it already exists
#'
#' @return linkSet object
#' @export
#'
#' @examples
#'   gr1 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                 ranges = IRanges(start = c(1000, 2000, 3000), width = 100),
#'                 strand = "+", symbol = c("BRCA1", "TP53", "NONEXISTENT"))
#'   gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                 ranges = IRanges(start = c(5000, 6000, 7000), width = 100),
#'                 strand = "+")
#'   ls <- linkSet(gr1, gr2, specificCol = "symbol")
#'
#'   # Test annotatePromoter
#'   annotated_ls <- suppressWarnings(annotatePromoter(ls, genome = "hg38", upstream = 500,overwrite = TRUE))
#'
#'
setMethod("annotatePromoter", "linkSet", function(x, genome = "hg38",
                                                  keyType = "symbol", upstream = 5000,
                                                  overwrite = FALSE) {
  if (!is.null(regionsBait(x)) && !overwrite) {
    warning("regionsBait already exists, set overwrite = TRUE to overwrite it")
    return(x)
  }
  if (!is.null(regionsBait(x)) && overwrite) {
    warning("regionsBait is being overwritten")
  }

  src <- NULL
  tryCatch({
    if (genome == "hg38") {
      src <- Organism.dplyr::src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
    } else if (genome == "mm10") {
      src <- Organism.dplyr::src_organism("TxDb.Mmusculus.UCSC.mm10.knownGene")
    } else {
      stop("Unsupported genome. Please use 'hg38' or 'mm10'.")
    }

    genes <- bait(x)
    geneGr <- Organism.dplyr::genes(src, filter = ~(symbol %in% genes))
    promoterGr <- IRanges::promoters(geneGr, upstream = upstream)

    index <- match(genes, geneGr$symbol)
    if (any(is.na(index))) {
      warning("Some genes are not found in the txDb, they will be set to chrNULL")
      newIndex <- which(!is.na(index))
      grMatch <- promoterGr[index[newIndex]]
      gr <- GRanges(seqnames = rep("chrNULL", length(genes)),
                    ranges = IRanges(rep(0, length(genes)), rep(0, length(genes))))
      mcols(grMatch) <- NULL
      gr[newIndex] <- grMatch
    } else {
      gr <- promoterGr[index]
    }

    regionsBait(x) <- gr
    return(x)
  }, error = function(e) {
    warning(paste("An error occurred:", e$message))
    return(x)
  }, finally = {
    if (!is.null(src) && inherits(src, "src_dbi")) {
      tryCatch({
        DBI::dbDisconnect(src$con)
      }, error = function(e) {
        warning(paste("Failed to disconnect from database:", e$message))
      })
    }
  })
})
