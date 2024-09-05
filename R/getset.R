###############################################################
# Simple getters and setters (some of these are exported):

setMethod("anchor1", "linkSet", function(x) { x@anchor1 })
setMethod("anchor2", "linkSet", function(x) { x@anchor2 })

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

setMethod("seqinfo", "linkSet", function(x) {
  seqinfo(regions(x,combine = T))
})

###############################################################
# Exported getters for anchors. 


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
setMethod("first", "linkSet", function(x) { anchors(x, type="bait") })
setMethod("second", "linkSet", function(x) { anchors(x, type="oe") })
setMethod("bait", "linkSet", function(x) { anchors(x, type="bait") })
setMethod("oe", "linkSet", function(x) { anchors(x, type="oe") })
