context("WREG.GLS tests")



test_that("Run WREG.GLS",{
  
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda")))
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/wreg.gls.staticOut.rda")))
  
  expect_silent({
    #Parameterize model.
    #all rows must be ordered identically between dataframes that contain site-specific information.
    ##Select Y variable
    Y <- staticData_peakFQ$Y$AEP_0.5
  })
  
  expect_silent({
    ##Select X variables
    X <- staticData_peakFQ$X[c("Sand",
                               "OutletElev",
                               "Slope")]
  })
  
  expect_silent({
    ##Define record lengths
    recordLengths <- staticData_peakFQ$recLen
  })
  
  expect_silent({
    ##Build LP3 matrix
    ###Select LP3K and merge with LP3f
    LP3 <- data.frame(S=staticData_peakFQ$LP3f$S,
                      K=staticData_peakFQ$LP3k$AEP_0.5,
                      G=staticData_peakFQ$LP3f$G)
  })
  
  expect_silent({
    ##Select basin characteristics
    basinChars <- staticData_peakFQ$BasChars
    ##Specify Y transformation if any
  })
  
  expect_silent({
    transY <- "none"
  })



  
  expect_silent({
    #Run WREG.GLS
    wreg.gls.out <- WREG.GLS(Y, X, recordLengths, LP3, basinChars, transY)
  })


  #Coefs
  expect_identical(wreg.gls.out$Coefs,wreg.gls.staticOut$Coefs)
  expect_identical(wreg.gls.out$ResLevInf,wreg.gls.staticOut$ResLevInf)
  expect_identical(wreg.gls.out$InflLim,wreg.gls.staticOut$InflLim)
  expect_identical(wreg.gls.out$LevInf.Sig,wreg.gls.staticOut$LevInf.Sig)
  expect_identical(wreg.gls.out$PerformanceMetrics,wreg.gls.staticOut$PerformanceMetrics)
  expect_identical(wreg.gls.out$X,wreg.gls.staticOut$X)
  expect_identical(wreg.gls.out$Y,wreg.gls.staticOut$Y)
  expect_identical(wreg.gls.out$fitted.values,wreg.gls.staticOut$fitted.values)
  expect_identical(wreg.gls.out$residuals,wreg.gls.staticOut$residuals)
  expect_identical(wreg.gls.out$Weighting,wreg.gls.staticOut$Weighting)
  expect_identical(wreg.gls.out$Inputs,wreg.gls.staticOut$Inputs)
  
})

  
  