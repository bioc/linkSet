
test_that("resizeRegions works correctly", {
  ls <- create_sample_linkSet()

  # Test resizing both regions
  result_both <- resizeRegions(ls, width = 100, fix = "center", region = "both")
  expect_equal(width(regions(result_both)), rep(100, 6))

  # Test resizing only bait region
  result_bait <- resizeRegions(ls, width = 75, fix = "start", region = "bait")
  expect_equal(width(regionsBait(result_bait)), rep(75, 3))
  expect_equal(width(oe(result_bait)), c(50, 50,50))

  # Test resizing only oe region
  result_oe <- resizeRegions(ls, width = 60, fix = "end", region = "oe")
  expect_equal(width(regionsBait(result_oe)), c(50, 50, 50))
  expect_equal(width(oe(result_oe)), rep(60, 3))
})

test_that("narrowRegions works correctly", {
  ls <- create_sample_linkSet()

  # Test narrowing both regions
  result_both <- narrowRegions(ls, start = 10, end = 40, region = "both")
  expect_equal(width(regions(result_both)), rep(31, 6))

  # Test narrowing only bait region
  result_bait <- narrowRegions(ls, start = 5, width = 20, region = "bait")
  expect_equal(width(regionsBait(result_bait)), rep(20, 3))
  expect_equal(width(oe(result_bait)), rep(50, 3))

  # Test narrowing only oe region
  result_oe <- narrowRegions(ls, end = 30, width = 25, region = "oe")
  expect_equal(width(regionsBait(result_oe)), rep(50, 3))
  expect_equal(width(oe(result_oe)), rep(25, 3))
})

test_that("shiftRegions works correctly", {
  ls <- create_sample_linkSet()

  # Test shifting both regions
  result_both <- shiftRegions(ls, shift = 10, region = "both")
  expect_equal(start(regions(result_both)), c(11L, 60L, 110L, 160L, 210L, 260L))

  # Test shifting only bait region
  result_bait <- shiftRegions(ls, shift = -5, region = "bait")
  expect_equal(start(regions(result_bait)), c(-4L, 50L, 95L, 150L, 195L, 250L))

  # Test shifting only oe region
  result_oe <- shiftRegions(ls, shift = 15, region = "oe")
  expect_equal(start(regions(result_oe)), c(1L, 65L, 100L, 165L, 200L, 265L))
})

test_that("flankRegions works correctly", {
  ls <- create_sample_linkSet()

  # Test flanking both regions
  result_both <- flankRegions(ls, width = 10, start = TRUE, region = "both")
  expect_equal(width(regions(result_both)), rep(10, 6))
  expect_equal(end(regions(result_both)), c(0L, 49L, 99L, 149L, 199L, 249L))
  expect_equal(start(regions(result_both)), c(-9L, 40L, 90L, 140L, 190L, 240L))


  # Test flanking only bait region
  result_bait <- flankRegions(ls, width = 15, start = FALSE, region = "bait")
  expect_equal(width(regionsBait(result_bait)), rep(15, 3))
  expect_equal(start(regionsBait(result_bait)), c(51, 150, 250))
  expect_equal(oe(result_bait), oe(ls))

  # Test flanking only oe region
  result_oe <- flankRegions(ls, width = 20, both = TRUE, region = "oe")
  expect_equal(width(oe(result_oe)), rep(20, 3))
  expect_equal(start(oe(result_oe)), c(30, 130, 230))
  expect_equal(regionsBait(result_oe), regionsBait(ls))
})

test_that("promoterRegions works correctly", {
  ls <- create_sample_linkSet()

  # Test creating promoter regions for both
  result_both <- promoterRegions(ls, upstream = 1000, downstream = 100, region = "both")
  expect_equal(width(regions(result_both)), rep(1100, 6))
  expect_equal(start(regions(result_both)), c(-999L, -950L, -900L, -850L, -800L, -750L))

  # Test creating promoter regions for bait only
  result_bait <- promoterRegions(ls, upstream = 500, downstream = 50, region = "bait")
  expect_equal(width(regionsBait(result_bait)), rep(550, 3))
  expect_equal(start(regionsBait(result_bait)), c(-499, -400, -300))
  expect_equal(oe(result_bait), oe(ls))

  # Test creating promoter regions for oe only
  result_oe <- promoterRegions(ls, upstream = 200, downstream = 20, region = "oe")
  expect_equal(width(oe(result_oe)), rep(220, 3))
  expect_equal(start(oe(result_oe)), c(-150, -50, 50))
  expect_equal(regionsBait(result_oe), regionsBait(ls))
})

test_that("reduceRegions works correctly", {
  # Create a sample linkSet object
  gr1 <- GRanges(
    seqnames = c("chr1", "chr1", "chr1", "chr2"),
    ranges = IRanges(start = c(1, 10, 20, 100), width = 15),
    strand = "+",
    symbol = c("Gene1", "Gene1", "Gene2", "Gene3")
  )
  gr2 <- GRanges(
    seqnames = c("chr1", "chr1", "chr2", "chr2"),
    ranges = IRanges(start = c(30, 40, 110, 120), width = 15),
    strand = "+"
  )
  ls <- linkSet(gr1, gr2, specificCol = "symbol")

  # Test reducing both regions
  result_both <- reduceRegions(ls, region = "both", countInteractions = FALSE)

  expect_equal(length(result_both), 4)
  expect_equal(length(unique(regionsBait(result_both))), 2)
  expect_equal(length(unique(oe(result_both))), 2)

  # Test reducing only bait regions
  result_bait <- reduceRegions(ls, region = "bait", countInteractions = FALSE)

  expect_equal(length(result_bait), 4)
  expect_equal(length(unique(regionsBait(result_bait))), 2)
  expect_equal(length(unique(oe(result_bait))), 4)

  # Test reducing only oe regions
  result_oe <- reduceRegions(ls, region = "oe", countInteractions = FALSE)

  expect_equal(length(result_oe), 4)
  expect_equal(length(unique(regionsBait(result_oe))), 4)
  expect_equal(length(unique(oe(result_oe))), 2)

  # Test with countInteractions = TRUE
  result_count <- reduceRegions(ls, region = "both", countInteractions = TRUE)

  expect_true("count" %in% colnames(mcols(result_count)))
  expect_equal(sum(mcols(result_count)$count), 4)

  # Test with additional arguments
  result_args <- reduceRegions(ls, region = "both", countInteractions = FALSE, ignore.strand = TRUE)

  expect_equal(length(result_args), 4)
  expect_equal(length(unique(regionsBait(result_args))), 2)
  expect_equal(length(unique(oe(result_args))), 2)
})
