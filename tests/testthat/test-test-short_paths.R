testthat::test_that("short_paths works", {
  data(slim_net, package = "igraphNetworkExpansion")
  data(genelists, package = "igraphNetworkExpansion")
  example_path <- short_paths(
    tnet = slim_net,
    target = genelists$targets$APP_Metabolism[1],
    targets = genelists$targets$APP_Metabolism,
    sentinals = genelists$sentinals$Immune_Response,
    cores = 1
  )
  exp_path <- list(
    Inter = c( 
      "PRKCH", "PRKACB", "RAC1", "PRKCG", "PRKCE", "BDNF", "PRKCA", "FYN", 
      "PRKCD","CAMK2A", "GRIN2B", "GRIN1", "SYT1", "SNAP25", "STX1A", "SLC24A2",
      "VAMP2", "SRC", "PFN1", "VLDLR", "ACTN1", "FN1", "PEBP1", "HSP90AA1", 
      "SLC8A1", "SCO1", "SH3GL1", "APP", "CALM1", "CACNB2", "SLC8A2", "RXRA", 
      "CAMK2D", "VIM", "NMNAT1", "PTGS2", "RUNX1", "TGFB1", "GRIK1", "VAV3",
      "TYROBP", "VAV2", "PLCG1", "BTK", "SYK", "ERBB4", "PLCG2", "PRKD1", 
      "MAPK3", "NCOA1", "NRG1", "NCOA2", "PRKCB", "DAPK1", "MAPK1", "CNGB1",
      "DLG2", "PRKCZ", "ACTN2", "ACTB", "RTN4", "PKN1", "SLC4A10", "SLC24A4", 
      "GRIA3", "DLG4", "GRIN2A", "LRP1", "SDC3", "SDC4", "SLC3A2", "CSNK2A1",
      "DLG3", "SORT1", "SQSTM1", "PRKCI", "MYD88", "TLR7", "CAMK2B", "ATG7",
      "GABARAPL1", "YWHAZ", "SNCA", "BECN1", "FOXF2", "PFN2", "AKT1", "SRPK2",
      "AR", "RELN", "DAB1", "ABL1", "GPC1", "SLC1A6", "NOG"
    ),
    Sentinal = c(
      "FYN", "SYT1", "VAMP2", "GRIA3", "PRKCD", "NCOA1", "NCOA2", "KAT2B", 
      "FABP7", "RAC1", "YWHAZ", "HPCAL1", "PRKCA", "CALM1", "ACTB", "ACTN2", 
      "SRC", "SDHA", "CFTR", "AHCYL1", "GRIN1", "ATP2A2", "RUNX1", "MAPK3", 
      "ERBB4", "NRG1", "PRKCB", "DAPK1", "RXRA", "PRKCE", "MAP3K1", "BDNF", 
      "RAP1A", "FRS2", "PRKCZ", "PRKCH", "FRS3", "PRKCG", "GAB1", "PIK3R2", 
      "RAPGEF1", "ZAP70", "BCR", "DOCK1", "DAB1", "PDGFRB", "CRK", "ABL1", 
      "CBL", "GAB2", "SOS1", "DLG4", "SHC1", "PIK3R1", "DLG1", "IGF1R", "BCAR1",
      "PXN", "VAV1", "JAK2", "GRB2", "LYN", "STAT5A", "STAT5B", "WAS", "PTPN11",
      "HIF1A", "SLC6A7", "TRAF2", "CAMK2A", "MAPT", "PRKCI", "PKN1", "UBA52",
      "RPS27A", "UBB", "UBC", "TGFB1", "HEXIM1", "HDAC5", "SLC20A1", "TRIM25", 
      "HLA-DRB1", "KTN1", "CDH2", "SLC24A2", "JUP", "CTNNA1", "CTNNB1", "MAX", 
      "CSNK2A1", "DLG3", "SLC3A2", "STX1A", "SLC13A5", "AHR", "HSP90AB1", 
      "PSMC1", "DESI1", "SNCA", "CASP4", "ARHGEF6", "APC", "GIT1", "ARHGEF7", 
      "PRKACA", "SNAP25", "CYFIP1", "PCNA", "SPG21", "MSH6", "MSH2", "SUMO1", 
      "BRCA1", "XRCC6", "TOP1", "CDKN1A", "FLOT1", "CDK5", "GRIN2A", "EGFR",
      "PKN3", "EP300", "PRKACB", "CCND1", "FUS", "CPSF7", "SLC6A6", "ATP1B1",
      "TRIM14", "SLC4A8", "SLC25A4"
    )
  )
  testthat::expect_equal(example_path, exp_path)
})
