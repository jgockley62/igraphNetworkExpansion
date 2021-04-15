testthat::test_that("path_calc", {
  example_path <- list(
    c('319', "49", "23", "86", "690", "238", "102")
  )
  names(example_path[[1]]) <- c(
    "GeneA","GeneZ", "GeneAlpha",
    "GeneB", "GeneX", "GeneOmega"
  )
  sampweights <- c(1.45, 2.45, 0.89, .003, 1.3, 2.1)
  names(sampweights) <- c(
    "GeneA","GeneZ", "GeneAlpha",
    "GeneB", "GeneX", "GeneOmega"
  )
  exp_pc <- 1.3486
  pc <- path_calc(example_path[[1]], sampweights)
  testthat::expect_equal(exp_pc, pc)
})
