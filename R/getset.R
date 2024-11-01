###############################################################
# Simple getters and setters (some of these are exported):

#' linkSet-accessors
#' @rdname linkSet-accessors
#' @aliases anchor1
#' @description
#' Methods to get and set fields in an linkSet object.
#' @return For the getters, values in various slots of x are returned, 
#' while for the setters, the slots of x are modified accordingly – see Details. 
#' @param x A linkSet object
#' @author Gilbert Han
#' @export
setMethod("anchor1", "linkSet", function(x) { x@anchor1 })

#' 
#' @rdname linkSet-accessors
#' @aliases anchor2
#' @param x A linkSet object
#' @export
setMethod("anchor2", "linkSet", function(x) { x@anchor2 })

#' 
#' @rdname linkSet-accessors
#' @aliases regions
#' @param x A linkSet object
#' @return A vector of the regions
#' @export
setMethod("regions", "linkSet", function(x) {x@regions})


###############################################################
# Seqinfo getting and setting.
#' @rdname linkSet-accessors
#' @aliases seqinfo
#' @param x A linkSet object
#' @export
setMethod("seqinfo", "linkSet", function(x) {
  seqinfo(regions(x))
})

###############################################################
# Exported getters for anchors.

#' @export
setMethod("anchorIds", "linkSet", function(x, type="both") {
  type <- match.arg(type, c("both", "bait", "oe"))
  if (type=="both") {
    out <- list(bait=anchor1(x), oe=anchor2(x))
    names(out$bait) <- names(out$oe) <- names(x)
  } else if (type=="bait") {
    out <- anchor1(x)
    names(out) <- names(x)
  } else {
    out <- anchor2(x)
    names(out) <- names(x)
  }
  return(out)
})


#' @export
setMethod("anchors", "linkSet", function(x, type="both", id=FALSE) {
  if (id) {
    return(anchorIds(x, type=type))
  }

  type <- match.arg(type, c("both", "bait", "oe"))
  if (type=="both") {
    out <- list(bait=x@nameBait, oe=regions(x)[anchor2(x)])
    names(out$bait) <- names(out$oe) <- names(x)
  } else if (type=="bait") {
    out <- x@nameBait
    names(out) <- names(x)
  } else {
    out <- regions(x)[anchor2(x)]
    names(out) <- names(x)
  }
  return(out)
})

# Defining some convenience methods.
#' @export
setMethod("first", "linkSet", function(x) { anchors(x, type="bait") })

#' @export
setMethod("second", "linkSet", function(x) { anchors(x, type="oe") })

#' @export
setMethod("bait", "linkSet", function(x) { anchors(x, type="bait") })

#' @export
setMethod("oe", "linkSet", function(x) { anchors(x, type="oe") })

#' @export
setMethod("regions","linkSet",function(x){
  return(x@regions)
})

#' @export
setMethod("regionsBait", "linkSet", function(x) {
  if (length(x@anchor1) == 0) {
    return(NULL)
  }
  regions(x)[anchor1(x)]
})

#' @export
# replace method for bait
setReplaceMethod("bait", "linkSet", function(x, value) {
  x@nameBait <- value
  return(x)
})

#' @export
setReplaceMethod("unchecked_regions", "linkSet", function(x, value) {
    x@regions <- value
    return(x)        
})

#' @export
setReplaceMethod("unchecked_anchor1", "linkSet", function(x, value) {
    x@anchor1 <- value 
    return(x)        
})

#' @export
setReplaceMethod("unchecked_anchor2", "linkSet", function(x, value) {
    x@anchor2 <- value 
    return(x)        
})

#' @export
setReplaceMethod("regions", "linkSet", function(x, value) {
  if (length(value)!=length(regions(x))) { 
      stop("assigned value must be of the same length as 'regions(x)'")
  }
  out <- .resort_regions(anchor1(x), anchor2(x), value)
  unchecked_anchor1(x) <- out$anchor1
  unchecked_anchor2(x) <- out$anchor2
  unchecked_regions(x) <- out$regions
  validObject(x)
  return(x)
})



#' @export
setReplaceMethod("regionsBait", "linkSet", function(x, value) {
  if (!is(value, "GRanges")) {
    stop("The 'value' must be a GRanges object")
  }
  if (length(value) != length(oe(x))) {
    stop("The length of 'value' must be equal to the number of bait regions")
  }
  metadata <- mcols(value)
  mcols(value) <- NULL
  if (length(anchor1(x)) == 0) {
    x@regions <- c(regions(x),value)
    x@anchor1 <- (length(regions(x))-length(value)+1) : length(regions(x))
  } else {
    regions(x)[anchor1(x)] <- value
  }
  mcols(x) <- cbind(mcols(x),metadata)
  return(x)
})

#' @export
setReplaceMethod("oe", "linkSet", function(x, value) {
  if (!is(value, "GRanges")) {
    stop("The 'value' must be a GRanges object")
  }
  if (length(value) != length(regions(x)[anchor2(x)])) {
    stop("The length of 'value' must be equal to the number of other end regions")
  }
  metadata <- mcols(value)
  mcols(value) <- NULL
  regions(x)[anchor2(x)] <- value
  mcols(x) <- cbind(mcols(x),metadata)
  return(x)
})

###############################################################
# Defining some other getters and setters.
#' @export
setMethod("$", "linkSet", function(x, name) {
    return(mcols(x)[[name]])
})

#' @export
setReplaceMethod("$", "linkSet", function(x, name, value) {
    mcols(x)[[name]] <- value
    return(x)
})


###############################################################
# Name getting and setting.

setMethod("names", "linkSet", function(x) { 
    x@NAMES 
})

setReplaceMethod("names", "linkSet", function(x, value) {
    if (!is.null(value) && !is.character(value)) { value <- as.character(value) }                
    x@NAMES <- value
    validObject(x)
    return(x)
})
