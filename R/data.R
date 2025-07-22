#' Example linkSet Object
#'
#' A dataset containing example genomic interactions in linkSet format. This example dataset
#' was created to demonstrate the functionality of the linkSet package for representing and
#' analyzing genomic interactions such as those from Hi-C or promoter-capture Hi-C experiments.
#'
#' The dataset represents simulated chromatin interactions between regulatory elements (enhancers)
#' and promoters across several chromosomes. It includes interaction counts, genomic coordinates
#' for both anchors of the interactions, and associated metadata.
#'
#' @format A linkSet object with example interactions. The object contains:
#'   \itemize{
#'     \item Bait regions (anchor1): GRanges object representing promoter regions
#'     \item Other end regions (anchor2): GRanges object representing potential enhancer regions
#'     \item Metadata columns including: count (interaction strength), baitID (unique identifiers
#'           for bait regions), and additional annotations
#'   }
#'
#'   The data was simulated to reflect typical patterns seen in chromatin interaction data,
#'   including distance-dependent interaction frequencies and varying interaction strengths.
#'
#' @source This is a synthetic dataset created specifically for the linkSet package
#'   to demonstrate various analysis workflows. The genomic coordinates are based on
#'   the human genome (hg38), but the interaction patterns were simulated.
#'
#' @usage data(linkExample)
#' @examples
#' data(linkExample)
#' show(linkExample)
#'
#' # Examine the structure
#' regions(linkExample)
#'
#' # View metadata
#' head(mcols(linkExample))
"linkExample"
