testthat::test_that("trace_filter works", {
  obj <- list(list( 
    Inter = c(
      "GeneA","GeneZ", "GeneAlpha",
      "GeneB", "GeneX", "GeneOmega",
      "Gene1","Gene2", "GeneUno",
      "Gene12", "Gene13", "GeneOcho"
    ),
    Sentinal = c(
      "GeneEh","GeneZee", "GeneAlpha",
      "GeneB", "GeneXray", "GeneOmega",
      "Gene11","Gene21", "GeneUno",
      "Gene1", "Gene3", "GeneOcho")
  ))
  tf_obj <- trace_filter( obj )
  exp_tf <- c("GeneAlpha", "GeneB", "GeneOmega", "GeneUno", "Gene1", "GeneOcho") 
  testthat::expect_equal(tf_obj, exp_tf)
})
