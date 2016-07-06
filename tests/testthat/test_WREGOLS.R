context("WREG.OLS test")

test_that("Run WREG.OLS",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/wreg.ols.staticOut.rda"))
  )
  
  # Organizing input data
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
  transY <- "none"
  
  # Run WLS regression
  expect_silent(
    resultTest <- WREG.OLS(Y, X, transY)
  )
  
  expect_silent(print(resultTest))
  
  expect_equal(resultTest,wreg.ols.staticOut)
})
  
  
test_that("Run WREG.OLS with log10",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/wreg.ols.staticOut.log10.rda"))
  )
  
  # Organizing input data
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
  transY <- "log10"
  
  # Run WLS regression with log10
  expect_silent(
    resultTest <- WREG.OLS(Y, X, transY)
  )
  expect_equal(resultTest,wreg.ols.staticOut.log10)
})

test_that("Run WREG.OLS with ln",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/wreg.ols.staticOut.ln.rda"))
  )
  
  # Organizing input data
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
  transY <- "log10"
  
  # Run WLS regression with log10
  expect_silent(
    resultTest <- WREG.OLS(Y, X, transY)
  )
  expect_equal(resultTest,wreg.ols.staticOut.ln)
})

test_that("check errors",{
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )

  # Organizing input data
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
  transY <- "log10"
  
  expect_error(WREG.OLS(Y="jazandapus", X, transY),
               "Invalid inputs were provided. See warnings()."
  )
  
  expect_error(WREG.OLS(Y, X, transY="jazandapus"),
               "Invalid inputs were provided. See warnings()."
  )
  
})