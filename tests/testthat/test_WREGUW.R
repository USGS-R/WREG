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
  resultTest <- WREG.UW(Y, X, customWeight = inverseRecLen, transY)
  
  expect_equal(resultTest,wreg.uw.staticOut)
})

test_that("Run WREG.UW with model variance list",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/wreg.uw.modVar.out.rda"))
  )
  
  # Organizing input data
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
  transY <- "none"
  
  # Make simple weighting using inverse record lengths
  inverseRecLen <- list(Omega=diag(1 / diag(staticData_peakFQ$recLen)),
                        var.modelerror.k = 2,
                        var.modelerror.0 = 1
  )
  
  # Run user-weights regression
  resultTest <- WREG.UW(Y, X, customWeight = inverseRecLen, transY)
  
  expect_equal(resultTest,wreg.uw.modvar.out)
})

test_that("Error checking",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )

  # Organizing input data
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("Sand", "OutletElev", "Slope")]
  transY <- "none"
  
  # Make simple weighting using inverse record lengths
  inverseRecLen <- list(Omega=diag(1 / diag(staticData_peakFQ$recLen)),
                        var.modelerror.k = 2,
                        var.modelerror.0 = 1
  )
  
  # Run user-weights regression
  expect_error(WREG.UW(Y="jazandapus", X, customWeight = inverseRecLen, transY),
               "Invalid inputs were provided. See warnings()"
  )
  expect_error(WREG.UW(Y, X="jazandapus", customWeight = inverseRecLen, transY),
               "missing value where TRUE/FALSE needed"
  )

  expect_error(WREG.UW(Y, X, customWeight = inverseRecLen, transY ="jazandapus"),
               "Invalid inputs were provided. See warnings()"
  )
})