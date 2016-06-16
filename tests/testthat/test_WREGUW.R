context("WREG.UW test")

test_that("Run WREG.UW",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/wreg.uw.staticOut.rda"))
  )
  
  # Organizing input data
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
  transY <- "none"
  
  # Make simple weighting using inverse record lengths
  inverseRecLen <- diag(1 / diag(staticData_peakFQ$recLen))
  
  # Run user-weights regression
  result_test <- WREG.UW(Y, X, customWeight = inverseRecLen, transY)
  
  expect_equal(result_test,result)
})

