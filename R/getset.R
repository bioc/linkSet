###############################################################
# Simple getters and setters (some of these are exported):
#' @export
setMethod("anchor1", "linkSet", function(x) { x@anchor1 })

#' @export
setMethod("anchor2", "linkSet", function(x) { x@anchor2 })

#' @export
setMethod("regions", "linkSet", function(x) {x@regions})


###############################################################
# Seqinfo getting and setting.
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
