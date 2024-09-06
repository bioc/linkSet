###############################################################
# Simple getters and setters (some of these are exported):
#' @export
setMethod("anchor1", "linkSet", function(x) { x@anchor1 })

#' @export
setMethod("anchor2", "linkSet", function(x) { x@anchor2 })

#' @export
setMethod("regions", "linkSet", function(x,type = "both",combine = F) {
  type <- match.arg(type, c("both", "bait", "oe"))
  if (type=="both" & !combine) {
    out <- list(bait=x@regionBait, oe=x@regionOE)
  }else if(type=="both" & combine){
    out <- c(x@regionBait, x@regionOE)
  }else if (type=="bait") {
    out <- x@regionBait
  }else {
    out=x@regionOE
  }
  return(out)
  })


###############################################################
# Seqinfo getting and setting.
#' @export
setMethod("seqinfo", "linkSet", function(x) {
  seqinfo(regions(x,combine = T))
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
    out <- list(bait=x@nameBait[anchor1(x)], oe=regions(x,"oe")[anchor2(x)])
    names(out$bait) <- names(out$oe) <- names(x)
  } else if (type=="bait") {
    out <- x@nameBait[anchor1(x)]
    names(out) <- names(x)
  } else {
    out <- regions(x,"oe")[anchor2(x)]
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
