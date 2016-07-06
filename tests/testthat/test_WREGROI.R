context("WREG.RoI")

test_that("Run WREG.RoI with GROI",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/wreg.roi.out.rda"))
  )
  
# Run a simple regression
Y <- staticData_peakFQ$Y$AEP_0.5
X <- staticData_peakFQ$X[c("A")]
transY <- "none"
basinChars <- staticData_peakFQ$BasChars

result <- WREG.RoI(Y, X, Reg = "OLS", transY, BasinChars = basinChars,
                   ROI='GRoI', n = 10L)
expect_equal(result,wreg.roi.out)
})

test_that("Run WREG.RoI with legacy",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/wreg.roi.out.legacy.rda"))
  )
  
  # Run a simple regression
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("A")]
  transY <- "none"
  basinChars <- staticData_peakFQ$BasChars
  
  result <- WREG.RoI(Y, X, Reg = "OLS", transY, BasinChars = basinChars,
                     ROI='GRoI', n = 10L,Legacy=TRUE)
  expect_equal(result,wreg.roi.out.legacy)
})

test_that("Run WREG.RoI PROI",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/wreg.roi.out.proi.rda"))
  )
  
  # Run a simple regression
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("A")]
  transY <- "none"
  basinChars <- staticData_peakFQ$BasChars
  
  result <- WREG.RoI(Y, X, Reg = "OLS", transY, BasinChars = basinChars,
                     ROI='PRoI', n = 10L,Legacy=TRUE)
  expect_equal(result,wreg.roi.out.proi)
})

test_that("Run WREG.RoI HRoI",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/wreg.roi.out.hroi.rda"))
  )
  
  # Run a simple regression
  Y <- staticData_peakFQ$Y$AEP_0.5
  X <- staticData_peakFQ$X[c("A")]
  transY <- "none"
  basinChars <- staticData_peakFQ$BasChars
  
  result <- WREG.RoI(Y, X, Reg = "OLS", transY, BasinChars = basinChars,
                     ROI='HRoI', n = 10L,Legacy=TRUE)
  expect_equal(result,wreg.roi.out.hroi)
})

test_that("Run WREG.RoI GLS",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/wreg.roi.out.gls.rda"))
  )
  
  # Run a GLS regression
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
  
  result <- WREG.RoI(Y,X, Reg = "GLS", 
                     transY, 
                     BasinChars = basinChars,
                     recordLengths = recordLengths,
                     LP3 = LP3,
                     ROI='HRoI', n = 10L,Legacy=TRUE)
  expect_equal(result,wreg.roi.out.gls)
})