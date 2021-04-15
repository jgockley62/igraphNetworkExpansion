testthat::test_that("find_limit works", {
  example_path <- list()
  example_path$res <- list(
    c('319', "49", "23", "86", "690", "238"),
    c('422', "899", "37", "240", "970", "28")
  )
  names(example_path$res[[1]]) <- c(
    "GeneA","GeneZ", "GeneAlpha",
    "GeneB", "GeneX", "GeneOmega"
  )
  names(example_path$res[[2]]) <- c(
    "Gene1","Gene2", "GeneUno",
    "Gene12", "Gene13", "GeneOcho"
  )
  example_path_b <-  example_path 
  sampweights <- c(1.45, 2.45, 0.89, .003, 1.3, 2.1, 0.02, 0, 0, 0.2, 0.6, .70)
  names(sampweights) <- c(
    "GeneA","GeneZ", "GeneAlpha",
    "GeneB", "GeneX", "GeneOmega",
    "Gene1","Gene2", "GeneUno",
    "Gene12", "Gene13", "GeneOcho"
  )
  
  limit <- find_limit(
    s_path = example_path,
    t_path =  example_path_b,
    weights = sampweights,
    cores = 2
  )
  exp_limit <- list()
  exp_limit$cutoff['Median'] <- 0.780375
  exp_limit$t_verts <- c(0, 0, 0)
  exp_limit$s_verts <- c(0, 0, 0)
    
  testthat::expect_equal(limit, exp_limit)
})
