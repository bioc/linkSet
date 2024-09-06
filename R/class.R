#== make linkObj-------------------------
#' @export
setClass("linkSet",
         contains="Vector",
         representation(
           anchor1="integer",
           anchor2="integer",
           nameBait = "character",
           regionBait="GenomicRanges_OR_missing",
           regionOE="GRanges",
           NAMES="character_OR_NULL",
           elementMetadata="DataFrame"
         )
)

#' @export
setClassUnion("character_Or_missing", c("character", "missing"))
