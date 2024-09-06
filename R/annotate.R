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
setMethod("annotatePromoter", "linkSet", function(x,genome = "hg38", keyType = "symbol",upstream = 500) {
  if (genome=="hg38"){
    src <- Organism.dplyr::src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
  }
  genes <- x@nameBait
  geneGr <- Organism.dplyr::genes(src,filter = ~(symbol %in% genes))
  promoterGr <- IRanges::promoters(geneGr,upstream = 500)
  index <- match(genes,geneGr$symbol)
  x@regionBait <- promoterGr[index]
  return(x)
})
