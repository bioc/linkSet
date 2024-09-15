# getset.R
#' @export
setGeneric("anchors", function(x, ...) standardGeneric("anchors"))

#' @export
setGeneric("anchor1", function(x) standardGeneric("anchor1"))

#' @export
setGeneric("anchor2", function(x) standardGeneric("anchor2"))

#' @export
setGeneric("regions", function(x, ...) standardGeneric("regions"))

#' @export
setGeneric("bait", function(x) standardGeneric("bait"))

#' @export
setGeneric("oe", function(x) standardGeneric("oe"))

#' @export
setGeneric("anchorIds", function(x, ...) standardGeneric("anchorIds"))



## set functions for linkSet
#' @export
setGeneric("bait<-", function(x, value) standardGeneric("bait<-"))

#' @export
setGeneric("regions<-", function(x, value) standardGeneric("regions<-"))

#' @export
setGeneric("anchor1<-", function(x, value) standardGeneric("anchor1<-"))

#' @export
setGeneric("anchor2<-", function(x, value) standardGeneric("anchor2<-"))

setGeneric("unchecked_regions<-", function(x, value) standardGeneric("unchecked_regions<-"))
setGeneric("unchecked_anchor1<-", function(x, value) standardGeneric("unchecked_anchor1<-"))
setGeneric("unchecked_anchor2<-", function(x, value) standardGeneric("unchecked_anchor2<-"))

#' @export
setGeneric("subsetBait", function(x, ...) {
  standardGeneric("subsetBait")
})

#' @export
setGeneric("subsetBaitRegion", function(x, ...) {
  standardGeneric("subsetBaitRegion")
})

#' @export
setGeneric("subsetOE", function(x, ...) {
  standardGeneric("subsetOE")
})


# annotate.R
#' @export
setGeneric("annotatePromoter", function(x, ...) standardGeneric("annotatePromoter"))

#' @export
setGeneric("oe<-", function(x, value) standardGeneric("oe<-"))

#' @export
setGeneric("regionsBait<-", function(x, value) standardGeneric("regionsBait<-"))

# method.R
#' @export
setGeneric("linkSet", function(anchor1, anchor2, specificCol, ...) standardGeneric("linkSet"))

#' @export
setGeneric("showLinkSet", function(x, margin="", print.seqinfo=FALSE, print.classinfo=FALSE, baitRegion=FALSE,...) standardGeneric("showLinkSet"))

#' @export
setGeneric("regionsBait", function(x) standardGeneric("regionsBait"))

#' @export
setGeneric("clean_unused_regions", function(x) standardGeneric("clean_unused_regions"))

# distance.R
#' @export
setGeneric("annotateInter", function(x) standardGeneric("annotateInter"))

#' @export
setGeneric("pairdist", function(x, type="mid") standardGeneric("pairdist"))


# convert.R
#' Convert different object types to linkSet
#' @aliases Convert
#' @param x Object to be converted
#' @param ... Additional arguments passed to methods
#'
#' @return A linkSet object
#' @export
setGeneric("Convert", function(x, ...) {
  standardGeneric("Convert")
})

#' @export
setGeneric("baitGInteractions", function(x, geneGr, peakGr,...) {
  standardGeneric("baitGInteractions")
})


# grange_method.R

#' @export
setGeneric("resizeRegions", function(x,width = 1000, fix = "center", use.names = TRUE,region="both", ...) standardGeneric("resizeRegions"))

#' @export
setGeneric("narrowRegions", function(x,  start = 1L, end = 1000L, width = 1000L, use.names = TRUE,region="both",...) standardGeneric("narrowRegions"))

#' @export
setGeneric("shiftRegions", function(x, width = 1000, shift = 0L, use.names = TRUE,region="both",...) standardGeneric("shiftRegions"))

#' @export
setGeneric("flankRegions", function(x,width = 1000,start=TRUE, both=FALSE, use.names=TRUE,ignore.strand=FALSE,region="both", ...) standardGeneric("flankRegions"))

#' @export
setGeneric("promoterRegions", function(x,upstream=2000, downstream=200, use.names=TRUE,region="both", ...) standardGeneric("promoterRegions"))



# count.R
#' @export
setGeneric("countInteractions", function(x, baitRegions = TRUE) standardGeneric("countInteractions"))

#' @export
setGeneric("reduceRegions", function(x, region = "both", countInteractions = TRUE, ...) {
  standardGeneric("reduceRegions")
})

# plot.R
#' @export
setGeneric("geom_linkset", function(linkSet, score.col = "count", score.threshold = NULL, score.color = c("grey70", "#56B1F7", "#132B43"),
 scale.range = 10,  plot.space = 0.1, plot.height = 0.2, arrow.size = 0.2,  remove_x_axis = FALSE, 
 link_plot_on_top = FALSE,extend.base = 10000,show.rect = FALSE) {
  standardGeneric("geom_linkset")
})
