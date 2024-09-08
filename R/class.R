#== make linkObj-------------------------
#' @export
setClass("linkSet",
         contains="Vector",
         representation(
           nameBait = "character",
           anchor1="integer",
           anchor2="integer",
           regions="GenomicRanges_OR_missing",
           NAMES="character_OR_NULL",
           elementMetadata="DataFrame"
         )
)

#' @export
setClassUnion("character_Or_missing", c("character", "missing"))
# setClassUnion("integer_Or_missing", c("integer", "missing"))
