context("Leverage")


test_that("leverage ROI = FALSE",{
  # Import some example data
  expect_silent({
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
    load(paste0(system.file("testData", package = "WREG"),"/leverageResult.out.rda"))
    
  })
  
  # Run a simple regression
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
  transY <- "none"
  result <- WREG.OLS(Y, X, transY)
  
  # calculate leverage of each point
  leverageResult <- Leverage(X = X,
                             Omega = result$Weighting)
  
  expect_equal(leverageResult,leverageResult.out)
  
})

test_that("leverage error checking",{
  # Import some example data

  expect_silent({
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))

  })
  
  # Run a simple regression
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
  transY <- "none"
  result <- WREG.OLS(Y, X, transY)
  
  expect_error(
    leverageResult <- Leverage(X = "jazandapus",
                               Omega = result$Weighting))
  expect_error(
    leverageResult <- Leverage(X,
                               Omega = "jazandapus")
  )
  
})