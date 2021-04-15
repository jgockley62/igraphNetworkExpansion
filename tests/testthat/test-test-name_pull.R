testthat::test_that("name_pull", {
  example_path <- list(
    c('319', "49", "23")
  )
  names(example_path[[1]]) <- c("GeneA", "GeneZ", "GeneAlpha")
  nam <- name_pull(example_path[[1]])
  exp_nam <- c("GeneA","GeneZ","GeneAlpha")
  testthat::expect_equal(nam, exp_nam)
})
