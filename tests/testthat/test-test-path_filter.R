testthat::test_that("path_filter Works", {
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
  pf_1 <- path_filter(example_path[[1]], sampweights, 1)
  pf_2 <- path_filter(example_path[[1]], sampweights, 1.5)
  
  exp_pf1 <- list( 
    Keep = c("GeneZ", "GeneAlpha", "GeneB", "GeneX","GeneOmega"),
    Genes = c("GeneZ", "GeneAlpha", "GeneB", "GeneX", "GeneOmega"),
    Used = 1,
    Passed = 1
  )
  exp_pf2 <- list( 
    Keep = NA,
    Genes = c("GeneZ", "GeneAlpha", "GeneB", "GeneX", "GeneOmega"),
    Used = 1,
    Passed = 0
  )
  
  testthat::expect_equal(pf_1, exp_pf1)
  testthat::expect_equal(pf_2, exp_pf2)
  
})
