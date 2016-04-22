###importPeakFQ tests

#Import the data
test_that("importPeakFQ", {
  #Load comparison dataset
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda")))
  
  #Import the data
  expect_silent({
    peakFQdir <- paste0(system.file("exampleDirectory", package = "WREG"),"/pfqImport")
    gisFilePath <- paste0(peakFQdir,"/pfqSiteInfo.txt")
    importedData <- importPeakFQ(pfqPath = peakFQdir,
                                 gisFile = gisFilePath)
  })
  
  #Check that imported data match static data
  expect_identical(importedData$Y,staticData_peakFQ$Y)
  expect_identical(importedData$sites,staticData_peakFQ$sites)
  expect_identical(importedData$AEP,staticData_peakFQ$AEP)
  expect_identical(importedData$X,staticData_peakFQ$X)
  expect_identical(importedData$LP3f,staticData_peakFQ$LP3f)
  expect_identical(importedData$LP3k,staticData_peakFQ$LP3k)
  expect_identical(importedData$BasChars,staticData_peakFQ$BasChars)
  expect_identical(importedData$MSEGR,staticData_peakFQ$MSEGR)
  
  #Check recLen, needs additional checks for rownames and colnames because contain siteIDs
  expect_identical(importedData$recLen,staticData_peakFQ$recLen)
  expect_identical(row.names(importedData$recLen),row.names(staticData_peakFQ$recLen))
  expect_identical(colnames(importedData$recLen),colnames(staticData_peakFQ$recLen))
  
  #Check recCor, needs additional checks for rownames and colnames because contain siteIDs
  expect_identical(importedData$recCor,staticData_peakFQ$recCor)
  expect_identical(row.names(importedData$recCor),row.names(staticData_peakFQ$recCor))
  expect_identical(colnames(importedData$recCor),colnames(staticData_peakFQ$recCor))
  
  
})