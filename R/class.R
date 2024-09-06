#== make linkObj-------------------------
#' @export
setClass("linkSet",
         contains="Vector",
         representation(
           anchor1="integer",
           anchor2="integer",
           nameBait = "character",
           region="GenomicRanges_OR_missing",
           NAMES="character_OR_NULL",
           elementMetadata="DataFrame"
         )
)

#' @export
setClassUnion("character_Or_missing", c("character", "missing"))
