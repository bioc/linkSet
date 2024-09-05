# Annotation--------------

setMethod("annotatePromoter", "linkSet", function(x,genome = "hg38", keyType = "symbol",upstream = 500) {
  if (genome=="hg38"){
    src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
  }
  genes <- LK_test@nameBait
  geneGr <- genes(src,filter = ~(symbol %in% genes))
  promoterGr <- promoters(geneGr,upstream = 500)
  index <- match(genes,geneGr$symbol)
  x@regionBait <- promoterGr[index]
  return(x)
})
