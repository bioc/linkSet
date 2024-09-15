# Annotation--------------

#' Annotate the link set with txDb. Give a gene list, and return a
#'
#'
#' @param x linkSet
#' @param genome the genome you want to annotate
#' @param keyType the key type. You can check with AnnotationDbi::keytypes
#' @param upstream The upstream base from the gene
#'
#'
#' @return linkSet object
#' @export
#'
#'
#'
#'
#'
setMethod("annotatePromoter", "linkSet", function(x,genome = "hg38",
                                                keyType = "symbol",upstream = 500,
                                                overwrite = FALSE) {
  if (!is.null(regionsBait(x)) && !overwrite){
    stop("regionsBait is already exist, set overwrite = TRUE to overwrite it")
  }
  if (!is.null(regionsBait(x)) && overwrite){
    warning("regionsBait is overwritten")
  }
  if (genome=="hg38"){
    src <- Organism.dplyr::src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
  } else if (genome == "mm10"){
    src <- Organism.dplyr::src_organism("TxDb.Mmusculus.UCSC.mm10.knownGene")
  }
  genes <- bait(x)
  geneGr <- Organism.dplyr::genes(src,filter = ~(symbol %in% genes))
  promoterGr <- IRanges::promoters(geneGr,upstream = 500)
  index <- match(genes,geneGr$symbol)
  if (any(is.na(index))){
    warning("Some genes are not found in the txDb, they will be set to chrNULL")
    newIndex <- which(!is.na(index))
    grMatch <- promoterGr[index[newIndex]]
    gr <- GRanges(seqnames = rep("chrNULL", length(genes)), ranges = IRanges(rep(0, length(genes)), rep(0, length(genes))))
    mcols(grMatch) <- NULL
    gr[newIndex] <- grMatch
  } else {
    gr <- promoterGr[index]
  }
  regionsBait(x) <- gr
  return(x)
})
