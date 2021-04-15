testthat::test_that("path_obj_filter works", {
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
  
  sampweights <- c(1.45, 2.45, 0.89, .003, 1.3, 2.1, 0.02, 0, 0, 0.2, 0.6, .70)
  names(sampweights) <- c(
    "GeneA","GeneZ", "GeneAlpha",
    "GeneB", "GeneX", "GeneOmega",
    "Gene1","Gene2", "GeneUno",
    "Gene12", "Gene13", "GeneOcho"
  )
  filt_obj <- path_obj_filter(
    path_obj = example_path,
    path_name =  "test",
    vertices = c(0,0,0),
    weights = sampweights,
    lim = 1,
    cores = 2
  )
  exp_filtobj <- c( "GeneZ", "GeneAlpha", "GeneB", "GeneX" )
  testthat::expect_equal(filt_obj, exp_filtobj)
})
