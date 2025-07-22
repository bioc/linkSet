# Helper function to create a sample linkSet object
#' Create sample linkSet object
#'
#' @description
#' Create a sample linkSet object for testing purposes
#'
#' @keywords internal
#' @return A linkSet object with sample data
createSampleLinkSet <- function() {
  gr1 <- GRanges(
    seqnames = c("chr1", "chr1", "chr2"),
    ranges = IRanges(start = c(1, 100, 200), width = 50),
    strand = "+",
    symbol = c("Gene1", "Gene2", "Gene3")
  )
  gr2 <- GRanges(
    seqnames = c("chr1", "chr2", "chr2"),
    ranges = IRanges(start = c(50, 150, 250), width = 50),
    strand = "+"
  )
  linkSet(gr1, gr2, specificCol = "symbol")
}
