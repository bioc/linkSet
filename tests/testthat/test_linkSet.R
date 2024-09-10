# Tests the construction and manipulation of linkSet objects.
# library(InteractionSet); library(testthat); source("test-GI.R")

# construct with methods 1 and 2
set.seed(7000)
N <- 40
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))
genes = c(rep("SP7",4),rep("ASPN",10),rep("XBP1",6))


Np <- 20
all.anchor1 <- sample(N, Np)
gr1 <- all.regions[all.anchor1]
gr1$symbol <- genes
all.anchor2 <- setdiff(1:40,all.anchor1)
gr2 <- all.regions[all.anchor2]
x <- linkSet(gr1, gr2,specificCol = "symbol")
x2 <- linkSet(genes, gr2)


test_that("show methods work for linkSet objects", {
  expect_output(sub("[0-9]", " ", show(x)), "linkSet object with 20 interactions and 1 metadata column:
              bait     seqnames_oe ranges_oe | anchor1.symbol
       <character>           <Rle> <IRanges> |    <character>
   [1]         SP7 ---        chrA     47-58 |            SP7
   [2]         SP7 ---        chrA     41-54 |            SP7
   [3]         SP7 ---        chrA     20-33 |            SP7
   [4]         SP7 ---        chrA    87-104 |            SP7
   [5]        ASPN ---        chrA     59-76 |           ASPN
   ...         ... ...         ...       ... .            ...
  [16]        XBP1 ---        chrB     54-61 |           XBP1
  [17]        XBP1 ---        chrB      1-18 |           XBP1
  [18]        XBP1 ---        chrB     83-91 |           XBP1
  [19]        XBP1 ---        chrB     33-45 |           XBP1
  [20]        XBP1 ---        chrB    96-114 |           XBP1
  -------
  regions: 20 ranges and 0 metadata columns
  seqinfo: 2 sequences from an unspecified genome; no seqlengths", fixed=TRUE)

  expect_output(sub("[0-9]", " ", show(x2)), "linkSet object with 20 interactions and 0 metadata columns:
              bait     seqnames_oe ranges_oe
       <character>           <Rle> <IRanges>
   [1]         SP7 ---        chrA     47-58
   [2]         SP7 ---        chrA     41-54
   [3]         SP7 ---        chrA     20-33
   [4]         SP7 ---        chrA    87-104
   [5]        ASPN ---        chrA     59-76
   ...         ... ...         ...       ...
  [16]        XBP1 ---        chrB     54-61
  [17]        XBP1 ---        chrB      1-18
  [18]        XBP1 ---        chrB     83-91
  [19]        XBP1 ---        chrB     33-45
  [20]        XBP1 ---        chrB    96-114
  -------
  regions: 20 ranges and 0 metadata columns
  seqinfo: 2 sequences from an unspecified genome; no seqlengths", fixed=TRUE)
})

######################################

o <- order(all.regions)
ref.regions <- all.regions[o]
new.pos <- integer(length(o))
new.pos[o] <- seq_along(new.pos)
ref.anchor1 <- new.pos[all.anchor1]
ref.anchor2 <- new.pos[all.anchor2]

test_that("slot access works in linkSet objects", {
  expect_that(x, is_a("linkSet"))
  expect_true(!is.unsorted(regions(x)))
  expect_identical(regions(x), ref.regions)

  expect_identical(anchors(x, id=TRUE, type="bait"), ref.anchor1)
  expect_identical(anchors(x, id=TRUE, type="oe"), ref.anchor2)
  expect_identical(anchors(x, id=TRUE), list(bait=ref.anchor1, oe=ref.anchor2))

  expect_identical(anchors(x, type="bait"), genes)
  expect_identical(anchors(x, type="oe"), ref.regions[ref.anchor2])
  expect_identical(anchors(x), list(bait=genes, oe=ref.regions[ref.anchor2]))
  expect_identical(anchors(x, type="bait"), first(x))
  expect_identical(anchors(x, type="oe"), second(x))
  expect_identical(anchors(x, type="bait"), bait(x))
  expect_identical(anchors(x, type="oe"), oe(x))
})

test_that("regionsBait works correctly", {
  # Test for linkSet with non-empty anchor1
  x_with_anchor1 <- linkSet(gr1, gr2, specificCol = "symbol")

  # Test regionsBait
  expect_identical(regionsBait(x_with_anchor1), regions(x_with_anchor1)[anchor1(x_with_anchor1)])

  # Check if the result is a GRanges object
  expect_s4_class(regionsBait(x_with_anchor1), "GRanges")

  # Check if the length of regionsBait matches the length of anchor1
  expect_equal(length(regionsBait(x_with_anchor1)), length(anchor1(x_with_anchor1)))

  # Test for linkSet with empty anchor1
  x_without_anchor1 <- linkSet(genes, gr2)

  # regionsBait should return NULL for linkSet with empty anchor1
  expect_null(regionsBait(x_without_anchor1))
})

test_that("setter functions work correctly", {
  # Create a sample linkSet object for testing
  gr1 <- GRanges(seqnames = c("chr1", "chr1", "chr2"),
                ranges = IRanges(start = c(1, 100, 200), width = 10),
                strand = "+", symbol = c("Gene1", "Gene2", "Gene3"))

  gr2 <- GRanges(seqnames = c("chr1", "chr1", "chr2"),
                ranges = IRanges(start = c(11, 110, 210), width = 10),
                strand = "+", symbol = c("Enh1", "Enh2", "Enh3"))

  ls <- linkSet(gr1, gr2, specificCol = "symbol")
  # Test setting bait
  new_bait <- c("NewGene1", "NewGene2", "NewGene3")
  bait(ls) <- new_bait
  expect_equal(bait(ls), new_bait)
  # Test setting names
  new_names <- c("Interaction1", "Interaction2", "Interaction3")
  names(ls) <- new_names
  expect_equal(names(ls), new_names)
  
  # Test setting metadata columns
  ls$new_score <- runif(length(ls))
  expect_true("new_score" %in% colnames(mcols(ls)))
})
######################################

test_that("parallel_slot_names works correctly", {
  # Create a linkSet object with non-empty anchor1
  x_with_anchor1 <- linkSet(gr1, gr2, specificCol = "symbol")

  # Create a linkSet object with empty anchor1
  x_without_anchor1 <- linkSet(genes, gr2)

  # Test for linkSet with non-empty anchor1
  expect_equal(
    parallel_slot_names(x_with_anchor1),
    c("anchor1", "anchor2", "nameBait", "NAMES", "elementMetadata")
  )

  # Test for linkSet with empty anchor1
  expect_equal(
    parallel_slot_names(x_without_anchor1),
    c("anchor2", "nameBait", "NAMES", "elementMetadata")
  )
})


######################################

