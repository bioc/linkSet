context("Distance calculations")

test_that("pairdist works correctly for linkSet objects", {
  # Create a mock linkSet object
  set.seed(42)
  gr1 <- GRanges(seqnames = c("chr1", "chr1", "chr2", "chr3", "chr4"),
                 ranges = IRanges(start = c(1, 100, 200, 300, 400), width = 10),
                 strand = "+")
  gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr2", "chr4", "chr4"),
                 ranges = IRanges(start = c(50, 150, 250, 350, 450), width = 10),
                 strand = "+")
  x <- linkSet(gr1, gr2)

  # Annotate interactions
  x <- annotateInter(x)

  # Test mid distance
  x_mid <- pairdist(x, type="mid")
  distances_mid <- mcols(x_mid)$distance
  expect_equal(length(distances_mid), length(x))
  expect_equal(distances_mid[1], 49)    # inter-chromosomal: abs((1 + 11) - (50 + 60)) / 2 = 49
  expect_true(is.na(distances_mid[2]))  # intra-chromosomal
  expect_equal(distances_mid[3], 50)  # intra-chromosomal
  expect_true(is.na(distances_mid[4]))  # intra-chromosomal
  expect_equal(distances_mid[5], 50)    # inter-chromosomal: abs((400 + 410) - (450 + 460)) / 2 = 50

  # Test gap distance
  x_gap <- pairdist(x, type="gap")
  distances_gap <- mcols(x_gap)$distance
  expect_equal(length(distances_gap), length(x))
  expect_equal(distances_gap[1], 39)    # inter-chromosomal: 50 - 11 - 1 = 38
  expect_true(is.na(distances_gap[2]))  # intra-chromosomal
  expect_equal(distances_gap[3], 40)  # intra-chromosomal
  expect_true(is.na(distances_gap[4]))  # intra-chromosomal
  expect_equal(distances_gap[5], 40)    # inter-chromosomal: 450 - 410 - 1 = 39

  # Test span distance
  x_span <- pairdist(x, type="span")
  distances_span <- mcols(x_span)$distance
  expect_equal(length(distances_span), length(x))
  expect_equal(distances_span[1], 59)   # inter-chromosomal: 60 - 1 + 1 = 60
  expect_true(is.na(distances_span[2])) # intra-chromosomal
  expect_equal(distances_span[3], 60) # intra-chromosomal
  expect_true(is.na(distances_span[4])) # intra-chromosomal
  expect_equal(distances_span[5], 60)   # inter-chromosomal: 460 - 400 + 1 = 61

  # Check if inter_type column exists
  expect_true("inter_type" %in% colnames(mcols(x)))

  # Check interaction types
  expect_equal(mcols(x)$inter_type, c("inter", "intra", "inter", "intra", "inter"))

  # Test with non-annotated regionsBait
  x_no_bait <- linkSet(paste0(gr1), gr2)
  expect_error(pairdist(x_no_bait, type="mid"), "No regionsBait found. Please annotate regionsBait first.")

  # Test error for invalid distance type
  expect_error(pairdist(x, type="invalid"), "should be one of")
})

test_that("annotateInteractions works correctly", {
  # Create a mock linkSet object
  set.seed(42)
  gr1 <- GRanges(seqnames = c("chr1", "chr1", "chr2"),
                 ranges = IRanges(start = c(1, 100, 200), width = 10),
                 strand = "+")
  gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr2"),
                 ranges = IRanges(start = c(50, 150, 250), width = 10),
                 strand = "+")
  x <- linkSet(gr1, gr2)



  # Run annotateInteractions
  x_annotated <- annotateInter(x)

  # Check if the new column is added
  expect_true("inter_type" %in% colnames(mcols(x_annotated)))

  # Check if the annotations are correct
  expected_annotations <- c("inter", "intra", "inter")
  expect_equal(mcols(x_annotated)$inter_type, expected_annotations)

  # Test error when regionsBait is not annotated
  x_no_bait <- linkSet(paste(gr1), gr2)
  expect_error(annotateInter(x_no_bait),
               "regionsBait is not annotated. Please run annotatePromoter first.")
})
