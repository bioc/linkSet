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


##########################################

test_that("filterLinks works correctly", {
  # Create a sample linkSet object
  gr1 <- GRanges(
    seqnames = c("chr1", "chr1", "chr2", "chr2", "chrNULL"),
    ranges = IRanges(start = c(1, 100, 200, 300, 400), width = 10),
    strand = "+",
    symbol = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
  )
  gr2 <- GRanges(
    seqnames = c("chr1", "chr2", "chr2", "chr3", "chr3"),
    ranges = IRanges(start = c(50, 150, 300, 350, 450), width = 10),
    strand = "+"
  )
  ls <- suppressWarnings(linkSet(gr1, gr2, specificCol = "symbol"))

  # Annotate inter/intra and calculate distances
  ls <- annotateInter(ls)
  ls <- pairdist(ls)

  # Test filtering intrachromosomal interactions
  result_intra <- filterLinks(ls, filter_intra = TRUE, filter_unannotate = FALSE, distance = NULL)
  expect_equal(length(result_intra), 2)  # Should keep only inter-chromosomal interactions

  # Test filtering unannotated interactions
  result_unannotate <- filterLinks(ls, filter_intra = FALSE, filter_unannotate = TRUE, distance = NULL)
  expect_equal(length(result_unannotate), 4)  # Should remove the chrNULL interaction

  # Test filtering by distance
  result_distance <- filterLinks(ls, filter_intra = FALSE, filter_unannotate = FALSE, distance = 100)
  expect_true(all(result_distance$distance <= 100))

  # Test all filters combined
  result_all <- filterLinks(ls, filter_intra = TRUE, filter_unannotate = TRUE, distance = 200)
  expect_true(length(result_all) < length(ls))
  expect_true(all(result_all$inter_type == "inter"))
  expect_true(all(as.character(seqnames(regionsBait(result_all))) != "chrNULL"))
  expect_true(all(result_all$distance <= 200))

  # Test with no filtering
  result_none <- filterLinks(ls, filter_intra = FALSE, filter_unannotate = FALSE, distance = NULL)
  expect_equal(length(result_none), length(ls))
})

test_that("filterLinks handles edge cases", {

  # Create a linkSet with only intrachromosomal interactions
  gr1 <- GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges(start = c(1, 100), width = 10),
    strand = "+",
    symbol = c("Gene1", "Gene2")
  )
  gr2 <- GRanges(
    seqnames = c("chr2", "chr2"),
    ranges = IRanges(start = c(50, 150), width = 10),
    strand = "+"
  )
  intra_ls <- suppressWarnings(linkSet(gr1, gr2, specificCol = "symbol"))
  intra_ls <- annotateInter(intra_ls)
  intra_ls <- pairdist(intra_ls)

  # Test filtering intrachromosomal interactions on intra-only linkSet
  result_intra_only <- filterLinks(intra_ls, filter_intra = TRUE)
  expect_equal(length(result_intra_only), 0)
})
