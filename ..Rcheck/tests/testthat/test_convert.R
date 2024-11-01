test_that("Convert function works correctly for GInteractions", {
  # Create a sample GInteractions object
  gi <- InteractionSet::GInteractions(
    anchor1 = c(1, 3, 5),
    anchor2 = c(2, 4, 6),
    regions = GRanges(
      seqnames = c("chr1", "chr1", "chr2", "chr2", "chr3", "chr3"),
      ranges = IRanges(start = c(1, 50, 100, 150, 200, 250), width = 10)
    ))
  mcols(gi) <- DataFrame(symbol = c("Gene1", "Gene2", "Gene3"), score = c(0.1, 0.2, 0.3))

  # Convert GInteractions to linkSet
  ls <- Convert(gi,baitCol="symbol")

  # Test that the conversion was successful
  expect_s4_class(ls, "linkSet")
  expect_equal(length(ls), 3)
  expect_equal(bait(ls), c("Gene1", "Gene2", "Gene3"))
  expect_equal(mcols(ls)$symbol, c("Gene1", "Gene2", "Gene3"))
  expect_equal(mcols(ls)$score, c(0.1, 0.2, 0.3))
})

test_that("Convert method for data.frame works correctly", {
  # Setup test data
  test_df <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    peak = c("chr1:1000-2000", "chr2:3000-4000", "chr3:5000-6000"),
    score = c(0.5, 0.7, 0.9),
    stringsAsFactors = FALSE
  )

  # Test successful conversion
  result <- Convert(test_df)
  expect_s4_class(result, "linkSet")
  expect_equal(length(result), 3)
  expect_equal(mcols(result)$score, c(0.5, 0.7, 0.9))

  # Test with custom column names
  result_custom <- Convert(test_df, baitCol = "gene", oeCol = "peak")
  expect_s4_class(result_custom, "linkSet")
  expect_equal(length(result_custom), 3)

  # Test error when columns are missing
  expect_error(Convert(test_df, baitCol = "nonexistent", oeCol = "peak"),
               "baitCol and oeCol must be columns in the data.frame")

  # Test with single row data
  single_row_df <- test_df[1, , drop = FALSE]
  result_single <- Convert(single_row_df)
  expect_s4_class(result_single, "linkSet")
  expect_equal(length(result_single), 1)

  # Test with no metadata columns
  no_metadata_df <- test_df[, c("gene", "peak")]
  result_no_metadata <- Convert(no_metadata_df)
  expect_s4_class(result_no_metadata, "linkSet")
  expect_equal(ncol(mcols(result_no_metadata)), 0)

  # Test with invalid peak format
  invalid_df <- test_df
  invalid_df$peak[1] <- "invalid_format"
  expect_error(Convert(invalid_df), "Please input peak format like chr1.816066.816566 or chr1:816066-816566")
})

test_that("Convert function throws error for unsupported types", {
  unsupported_object <- matrix(1:9, 3, 3)
  expect_error(Convert(unsupported_object), "Conversion from matrix to linkSet is not supported")
})

test_that("baitGInteractions function works correctly", {
  # Create sample data
  gi <- InteractionSet::GInteractions(
    anchor1 = c(1, 3, 5),
    anchor2 = c(2, 4, 6),
    regions = GRanges(
      seqnames = c("chr1", "chr1", "chr2", "chr2", "chr3", "chr3"),
      ranges = IRanges(start = c(1, 50, 100, 150, 200, 250), width = 10)
    )
  )

  geneGr <- GRanges(
    seqnames = c("chr1", "chr2", "chr3"),
    ranges = IRanges(start = c(5, 105, 205), width = 20),
    symbol = c("Gene1", "Gene2", "Gene3")
  )

  peakGr <- GRanges(
    seqnames = c("chr1", "chr2", "chr3"),
    ranges = IRanges(start = c(45, 145, 245), width = 10)
  )

  # Run baitGInteractions
  result <- baitGInteractions(gi, geneGr, peakGr, geneSymbol = "symbol")

  # Test the result
  expect_s4_class(result, "linkSet")
  expect_equal(length(result), 3)
  expect_equal(bait(result), c("Gene1", "Gene2", "Gene3"))
  expect_equal(as.character(seqnames(oe(result))), c("chr1", "chr2", "chr3"))
  expect_equal(start(oe(result)), c(45, 145, 245))
})
