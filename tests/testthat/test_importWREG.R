context("Import data function tests")



###importWREG tests

#Import the data
test_that("importWREG", {
  #Load comparison dataset
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/staticData_wreg.rda")))
  
  #Import the data
  expect_silent({
    wregDir <- paste0(system.file("exampleDirectory", package = "WREG"),"/matlabImport")
    importedData <- importWREG(wregPath = wregDir)
  })
  
  #Check that imported data match static data
  expect_identical(importedData$Y,staticData_WREG$Y)
  expect_identical(importedData$sites,staticData_WREG$sites)
  expect_identical(importedData$AEP,staticData_WREG$AEP)
  expect_identical(importedData$X,staticData_WREG$X)
  expect_identical(importedData$LP3f,staticData_WREG$LP3f)
  expect_identical(importedData$LP3k,staticData_WREG$LP3k)
  expect_identical(importedData$BasChars,staticData_WREG$BasChars)

  #Check recLen, needs additional checks for rownames and colnames because contain siteIDs
  expect_equal(importedData$recLen,staticData_WREG$recLen)
  expect_identical(row.names(importedData$recLen),row.names(staticData_WREG$recLen))
  expect_identical(colnames(importedData$recLen),colnames(staticData_WREG$recLen))
  
  #Check recCor, needs additional checks for rownames and colnames because contain siteIDs
  expect_equal(importedData$recCor,staticData_WREG$recCor)
  expect_identical(row.names(importedData$recCor),row.names(staticData_WREG$recCor))
  expect_identical(colnames(importedData$recCor),colnames(staticData_WREG$recCor))
  
  
})
