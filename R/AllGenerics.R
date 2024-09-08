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

# annotate.R
#' @export
setGeneric("annotatePromoter", function(x, ...) standardGeneric("annotatePromoter"))

#' @export
setGeneric("linkSet", function(anchor1, anchor2, specificCol, ...) standardGeneric("linkSet"))

#' @export
setGeneric("regionsBait", function(x) standardGeneric("regionsBait"))

# distance.R
#' @export
setGeneric("annotateInter", function(x) standardGeneric("annotateInter"))

#' @export
setGeneric("pairdist", function(x, type="mid") standardGeneric("pairdist"))