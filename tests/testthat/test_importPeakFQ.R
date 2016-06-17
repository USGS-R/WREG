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
  expect_equal(importedData$Y,staticData_peakFQ$Y)
  expect_equal(importedData$sites,staticData_peakFQ$sites)
  expect_equal(importedData$AEP,staticData_peakFQ$AEP)
  expect_equal(importedData$X,staticData_peakFQ$X)
  expect_equal(importedData$LP3f,staticData_peakFQ$LP3f)
  expect_equal(importedData$LP3k,staticData_peakFQ$LP3k)
  expect_equal(importedData$BasChars,staticData_peakFQ$BasChars)
  expect_equal(importedData$MSEGR,staticData_peakFQ$MSEGR)
  
  #Check recLen, needs additional checks for rownames and colnames because contain siteIDs
  expect_equal(importedData$recLen,staticData_peakFQ$recLen)
  expect_equal(row.names(importedData$recLen),row.names(staticData_peakFQ$recLen))
  expect_equal(colnames(importedData$recLen),colnames(staticData_peakFQ$recLen))
  
  #Check recCor, needs additional checks for rownames and colnames because contain siteIDs
  expect_equal(importedData$recCor,staticData_peakFQ$recCor)
  expect_equal(row.names(importedData$recCor),row.names(staticData_peakFQ$recCor))
  expect_equal(colnames(importedData$recCor),colnames(staticData_peakFQ$recCor))
  
  
})