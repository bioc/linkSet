test_that("countInteractions works correctly", {
  # Create a sample linkSet object
  gr1 <- GRanges(seqnames = c("chr1", "chr1", "chr2", "chr2"),
                 ranges = IRanges(start = c(1, 100, 200, 300), width = 10),
                 strand = "+", symbol = c("Gene1", "Gene2", "Gene3", "Gene4"))
  gr2 <- GRanges(seqnames = c("chr1", "chr1", "chr2", "chr2"),
                 ranges = IRanges(start = c(50, 150, 250, 350), width = 10),
                 strand = "+")
  ls <- linkSet(gr1, gr2, specificCol = "symbol")
  
  # Duplicate some interactions
  ls_dup <- c(ls, ls[c(1, 2)])
  
  # Test countInteractions with baitRegions = TRUE
  result_with_bait <- countInteractions(ls_dup, baitRegions = TRUE)
  
  # Check if the result is a linkSet object
  expect_s4_class(result_with_bait, "linkSet")
  
  # Check if the number of unique interactions is correct
  expect_equal(length(result_with_bait), 4)
  
  # Check if the counts are correct
  expect_equal(mcols(result_with_bait)$count, c(2, 2, 1, 1))
  
  # Test countInteractions with baitRegions = FALSE
  result_without_bait <- countInteractions(ls_dup, baitRegions = FALSE)
  
  # Check if the result is a linkSet object
  expect_s4_class(result_without_bait, "linkSet")
  
  # Check if the number of unique interactions is correct
  expect_equal(length(result_without_bait), 4)
  
  # Check if the counts are correct
  expect_equal(mcols(result_without_bait)$count, c(2, 2, 1, 1))
  
  # Test with all unique interactions
  result_unique <- countInteractions(ls)
  expect_equal(length(result_unique), 4)
  expect_equal(mcols(result_unique)$count, c(1, 1, 1, 1))
})