# getset.R
setGeneric("anchors", function(x, ...) standardGeneric("anchors"))
setGeneric("anchor1", function(x) standardGeneric("anchor1"))
setGeneric("anchor2", function(x) standardGeneric("anchor2"))
setGeneric("regions", function(x, ...) standardGeneric("regions"))
setGeneric("bait", function(x) standardGeneric("bait"))
setGeneric("oe", function(x) standardGeneric("oe"))

# annotate.R
setGeneric("annotatePromoter", function(x, ...) standardGeneric("annotatePromoter"))

setGeneric("linkSet", function(anchor1, anchor2, specificCol, ...) standardGeneric("linkSet"))
