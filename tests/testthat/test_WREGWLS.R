context("WREG.WLS test")

test_that("Run WREG.WLS",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/wreg.wls.staticOut.rda"))
  )
  
  # Organizing input data
  lp3Data <- staticData_peakFQ$LP3f
  lp3Data$K <- staticData_peakFQ$LP3k$AEP_0.5
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
  recordLengths <- staticData_peakFQ$recLen
  transY <- "none"
  
  # Run WLS regression
  result_test <- WREG.WLS(Y, X, recordLengths, LP3 = lp3Data, transY)
  
  expect_equal(result_test,result)
})