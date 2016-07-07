context("WREG.GLS tests")



test_that("Run WREG.GLS with transY=none",{
  
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda")))
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/wreg.gls.staticOut.rda")))
  
  
  #Parameterize model.
  #all rows must be ordered identically between dataframes that contain site-specific information.
  ##Select Y variable
  Y <- staticData_peakFQ$Y$AEP_0.5
  ##Select X variables
  X <- staticData_peakFQ$X[c("Sand",
                             "OutletElev",
                             "Slope")]
  ##Define record lengths
  recordLengths <- staticData_peakFQ$recLen
  ##Build LP3 matrix
  ###Select LP3K and merge with LP3f
  LP3 <- data.frame(S=staticData_peakFQ$LP3f$S,
                    K=staticData_peakFQ$LP3k$AEP_0.5,
                    G=staticData_peakFQ$LP3f$G)
  ##Select basin characteristics
  basinChars <- staticData_peakFQ$BasChars
  transY <- "none"
  
  
  
  
  #Run WREG.GLS
  
  expect_silent({
    wreg.gls.out <- WREG.GLS(Y, X, recordLengths, LP3, basinChars, transY)
  })
  
  
  #Coefs
  expect_equal(wreg.gls.out$Coefs,wreg.gls.staticOut$Coefs)
  expect_equal(wreg.gls.out$ResLevInf,wreg.gls.staticOut$ResLevInf)
  expect_equal(wreg.gls.out$InflLim,wreg.gls.staticOut$InflLim)
  expect_equal(wreg.gls.out$LevInf.Sig,wreg.gls.staticOut$LevInf.Sig)
  expect_equal(wreg.gls.out$PerformanceMetrics,wreg.gls.staticOut$PerformanceMetrics)
  expect_equal(wreg.gls.out$X,wreg.gls.staticOut$X)
  expect_equal(wreg.gls.out$Y,wreg.gls.staticOut$Y)
  expect_equal(wreg.gls.out$fitted.values,wreg.gls.staticOut$fitted.values)
  expect_equal(wreg.gls.out$residuals,wreg.gls.staticOut$residuals)
  expect_equal(wreg.gls.out$Weighting,wreg.gls.staticOut$Weighting)
  expect_equal(wreg.gls.out$Inputs,wreg.gls.staticOut$Inputs)
  
})

test_that("Run WREG.GLS with transY=log10",{
  
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda")))
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/wreg.gls.staticOut.log10.rda")))
  
  #Parameterize model.
  #all rows must be ordered identically between dataframes that contain site-specific information.
  ##Select Y variable
  Y <- staticData_peakFQ$Y$AEP_0.5
  ##Select X variables
  X <- staticData_peakFQ$X[c("Sand",
                             "OutletElev",
                             "Slope")]
  ##Define record lengths
  recordLengths <- staticData_peakFQ$recLen
  ##Build LP3 matrix
  ###Select LP3K and merge with LP3f
  LP3 <- data.frame(S=staticData_peakFQ$LP3f$S,
                    K=staticData_peakFQ$LP3k$AEP_0.5,
                    G=staticData_peakFQ$LP3f$G)
  ##Select basin characteristics
  basinChars <- staticData_peakFQ$BasChars
  transY <- "log10"
  
  #Run WREG.GLS with log10
  expect_silent({
    wreg.gls.out <- WREG.GLS(Y, X, recordLengths, LP3, basinChars, transY)
  })
  
  
  #Coefs
  expect_equal(wreg.gls.out$Coefs,wreg.gls.staticOut$Coefs)
  expect_equal(wreg.gls.out$ResLevInf,wreg.gls.staticOut$ResLevInf)
  expect_equal(wreg.gls.out$InflLim,wreg.gls.staticOut$InflLim)
  expect_equal(wreg.gls.out$LevInf.Sig,wreg.gls.staticOut$LevInf.Sig)
  expect_equal(wreg.gls.out$PerformanceMetrics,wreg.gls.staticOut$PerformanceMetrics)
  expect_equal(wreg.gls.out$X,wreg.gls.staticOut$X)
  expect_equal(wreg.gls.out$Y,wreg.gls.staticOut$Y)
  expect_equal(wreg.gls.out$fitted.values,wreg.gls.staticOut$fitted.values)
  expect_equal(wreg.gls.out$residuals,wreg.gls.staticOut$residuals)
  expect_equal(wreg.gls.out$Weighting,wreg.gls.staticOut$Weighting)
  expect_equal(wreg.gls.out$Inputs,wreg.gls.staticOut$Inputs)
  
})

test_that("Run WREG.GLS with transY=ln",{
  
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda")))
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/wreg.gls.staticOut.log10.rda")))
  
  #Parameterize model.
  #all rows must be ordered identically between dataframes that contain site-specific information.
  ##Select Y variable
  Y <- staticData_peakFQ$Y$AEP_0.5
  ##Select X variables
  X <- staticData_peakFQ$X[c("Sand",
                             "OutletElev",
                             "Slope")]
  ##Define record lengths
  recordLengths <- staticData_peakFQ$recLen
  ##Build LP3 matrix
  ###Select LP3K and merge with LP3f
  LP3 <- data.frame(S=staticData_peakFQ$LP3f$S,
                    K=staticData_peakFQ$LP3k$AEP_0.5,
                    G=staticData_peakFQ$LP3f$G)
  ##Select basin characteristics
  basinChars <- staticData_peakFQ$BasChars
  transY <- "ln"
  
  #Run WREG.GLS with ln
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/wreg.gls.staticOut.ln.rda")))
  
  expect_silent({
    wreg.gls.out <- WREG.GLS(Y, X, recordLengths, LP3, basinChars, transY)
  })
  
  expect_warning(print(wreg.gls.out),regexp=NA)
  expect_error(print(wreg.gls.out),regexp=NA)
  
  #Coefs
  expect_equal(wreg.gls.out$Coefs,wreg.gls.staticOut$Coefs)
  expect_equal(wreg.gls.out$ResLevInf,wreg.gls.staticOut$ResLevInf)
  expect_equal(wreg.gls.out$InflLim,wreg.gls.staticOut$InflLim)
  expect_equal(wreg.gls.out$LevInf.Sig,wreg.gls.staticOut$LevInf.Sig)
  expect_equal(wreg.gls.out$PerformanceMetrics,wreg.gls.staticOut$PerformanceMetrics)
  expect_equal(wreg.gls.out$X,wreg.gls.staticOut$X)
  expect_equal(wreg.gls.out$Y,wreg.gls.staticOut$Y)
  expect_equal(wreg.gls.out$fitted.values,wreg.gls.staticOut$fitted.values)
  expect_equal(wreg.gls.out$residuals,wreg.gls.staticOut$residuals)
  expect_equal(wreg.gls.out$Weighting,wreg.gls.staticOut$Weighting)
  expect_equal(wreg.gls.out$Inputs,wreg.gls.staticOut$Inputs)
  
})

test_that("Error check tests",{
  
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda")))
  #Parameterize model.
  #all rows must be ordered identically between dataframes that contain site-specific information.
  ##Select Y variable
  Y <- staticData_peakFQ$Y$AEP_0.5
  ##Select X variables
  X <- staticData_peakFQ$X[c("Sand",
                             "OutletElev",
                             "Slope")]
  ##Define record lengths
  recordLengths <- staticData_peakFQ$recLen
  ##Build LP3 matrix
  ###Select LP3K and merge with LP3f
  LP3 <- data.frame(S=staticData_peakFQ$LP3f$S,
                    K=staticData_peakFQ$LP3k$AEP_0.5,
                    G=staticData_peakFQ$LP3f$G)
  ##Select basin characteristics
  basinChars <- staticData_peakFQ$BasChars
  transY <- "none"
  
  expect_error(
    WREG.GLS(Y="jazandapus", X, recordLengths, LP3, basinChars, transY),
    "Invalid inputs were provided. See warnings()."
  )
  expect_error(
    WREG.GLS(Y, X="jazandapus", recordLengths, LP3, basinChars, transY),
    "missing value where TRUE/FALSE needed"
  )
  expect_error(
    WREG.GLS(Y, X, recordLengths="jazandapus", LP3, basinChars, transY),
    "argument is of length zero"
  )
  expect_error(
    WREG.GLS(Y, X, recordLengths, LP3="jazandapus", basinChars, transY),
    "Invalid inputs were provided.  See warnings()."
  )
  expect_error(
    WREG.GLS(Y, X, recordLengths, LP3, basinChars="jazandapus", transY)
  )
  expect_error(
    WREG.GLS(Y, X, recordLengths, LP3, basinChars, transY="jazandapus"),
    "Invalid inputs were provided. See warnings()."
  )
})


