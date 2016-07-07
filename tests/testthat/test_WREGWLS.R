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
  expect_silent(
    resultTest <- WREG.WLS(Y, X, recordLengths, LP3 = lp3Data, transY)
  )
  
  #wreg.wls.staticOut <- resultTest
  #save(wreg.wls.staticOut,file="inst/testData/wreg.wls.staticOut.rda")
  
  expect_warning(print(resultTest),regexp=NA)
  expect_error(print(resultTest),regexp=NA)
  
  expect_equal(resultTest,wreg.wls.staticOut)
  
})

test_that("Run WREG.WLS with log 10",{
  # Run WLS regression with log10
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/wreg.wls.staticOut.log10.rda"))
  )
  
  expect_silent(
    resultTest <- WREG.WLS(Y, X, recordLengths, LP3 = lp3Data, transY = "log10")
  )
  
  #wreg.wls.staticOut <- resultTest
  #save(wreg.wls.staticOut,file="inst/testData/wreg.wls.staticOut.log10.rda")
  
  expect_equal(resultTest,wreg.wls.staticOut)
  
})

test_that("Run WREG.WLS with ln",{
  # Run WLS regression with ln
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/wreg.wls.ln.staticOut.rda"))
  )
  
  expect_silent(
    resultTest <- WREG.WLS(Y, X, recordLengths, LP3 = lp3Data, transY = "ln")
  )
  
  #wreg.wls.staticOut <- resultTest
  #save(wreg.wls.staticOut,file="inst/testData/wreg.wls.staticOut.ln.rda")
  
  expect_equal(resultTest,wreg.wls.staticOut)
  
})


# #Check warnings
# expect_error(WREG.WLS(Y="jazandapus",X,recordLengths,LP3=lp3Data,transY),
#              "Invalid inputs were provided.  See warnings()."
#              )
# 
# expect_error(WREG.WLS(Y,X="jazandapus",recordLengths,LP3=lp3Data,transY)
# )
# expect_error(WREG.WLS(Y,X,recordLengths="jazandapus",LP3=lp3Data,transY),
#              "Invalid inputs were provided.  See warnings()."
# )
# expect_error(WREG.WLS(Y,X,recordLengths,LP3=lp3Data,transY="jazandapus"),
#              "Invalid inputs were provided.  See warnings()."
# )
