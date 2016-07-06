context("Omega.GLS tests")

test_that("Run OMEGA GLS",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/omega.gls.out.rda"))
  )
  
  # Organizing input data
  lp3Data <- staticData_peakFQ$LP3f
  lp3Data$K <- staticData_peakFQ$LP3k$AEP_0.5
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
  
  # Compute weighting matrix
  weightingResult <- Omega.GLS(alpha = 0.01, theta = 0.98,
                               independent = staticData_peakFQ$BasChars, X = X,
                               Y = Y, recordLengths = staticData_peakFQ$recLen,
                               LP3 = lp3Data, MSEGR = NA, TY = 20, peak = TRUE, distMeth = 2)
  expect_equal(weightingResult,omega.gls.out)
})

test_that("Run OMEGA GLS with MSEGR",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/omega.gls.out.msegr.rda"))
  )
  
  # Organizing input data
  lp3Data <- staticData_peakFQ$LP3f
  lp3Data$K <- staticData_peakFQ$LP3k$AEP_0.5
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
  
  # Compute weighting matrix
  weightingResult <- Omega.GLS(alpha = 0.01, theta = 0.98,
                               independent = staticData_peakFQ$BasChars, X = X,
                               Y = Y, recordLengths = staticData_peakFQ$recLen,
                               LP3 = lp3Data, MSEGR = 1, TY = 20, peak = TRUE, distMeth = 2)
  expect_equal(weightingResult,omega.gls.out.msegr)
})
