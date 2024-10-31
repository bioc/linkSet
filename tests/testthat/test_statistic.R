test_that("run_chicane processes linkSet correctly", {
  # Load example data
  # Prepare the linkSet object
  linkExample <- countInteractibility(linkExample)
  linkExample <- pairdist(linkExample)
  
  # Run the function
  result <- run_chicane(linkExample)
  
  # Check that the result is a linkSet object
  expect_s4_class(result, "linkSet")
  
  # Check that the result contains expected columns
  expected_columns <- c("expected", "p.value", "q.value")
  expect_true(all(expected_columns %in% colnames(mcols(result))))
  
  # Check that the expected values are calculated (not NA)
  expect_false(any(is.na(mcols(result)$expected)))
  
  # Check that p-values and q-values are calculated (not NA)
  expect_false(any(is.na(mcols(result)$p.value)))
  expect_false(any(is.na(mcols(result)$q.value)))
  
  # Check that the bait.id and bait.to.bait columns are added
  expect_true("bait.id" %in% colnames(mcols(result)))
  expect_true("bait.to.bait" %in% colnames(mcols(result)))
  
  # Check that bait.to.bait is set to FALSE by default
  expect_true(all(mcols(result)$bait.to.bait == FALSE))
})
