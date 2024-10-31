# Helper function to create a sample linkSet object
create_sample_linkSet <- function() {
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