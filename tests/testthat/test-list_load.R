testthat::test_that("list_load ", {
  data(slim_net, package = "igraphNetworkExpansion")
  data("genelists", package = "igraphNetworkExpansion")
  all_goterms <- list(
    c(
      'syn25185319',
      system.file(
        "extdata/inputlists/",
        "APP_Metabolism.txt",
        package = "igraphNetworkExpansion"
      ),
      "APP Metabolism"
      ),
    c(
      'syn25185320',
      system.file(
        "extdata/inputlists/",
        "Endolysosomal.txt",
        package = "igraphNetworkExpansion"
      ),
      "Endolysosomal"
      )
  )
  output <- list_load(
    all_goterms[[1]][2],
    network = slim_net
  )
  expected <- c(
    "GRIN2B", "TREM2", "IGF1R", "APOE", "CLU", "GRIN2A", "APP", "NAMPT", "INSR",
    "GJA1", "LGMN", "FYN", "GRIN1", "EPHA4", "ITGB2", "PICALM", "LRPAP1", 
    "GRIA2", "LRP1", "ABCA7", "BCL2L11", "GRIA3", "DLGAP3", "CD74"
  )
    testthat::expect_equal(expected, output)
})
