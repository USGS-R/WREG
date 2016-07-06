context("Supporting function tests")



test_that("Dist.WREG",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  
  # For two sites, compute the inter-site distance
  # Use the haversine approximation
  expect_silent(intersiteDistance_meth1_test <- Dist.WREG(Lat1 = staticData_peakFQ$BasChars$Lat[1],
                                                          Long1 = staticData_peakFQ$BasChars$Lat[1],
                                                          Lat2 = staticData_peakFQ$BasChars$Lat[2],
                                                          Long2 = staticData_peakFQ$BasChars$Lat[2],
                                                          method = 1)
  )
  
  
  expect_silent(intersiteDistance_meth2_test <- Dist.WREG(Lat1 = staticData_peakFQ$BasChars$Lat[1],
                                                          Long1 = staticData_peakFQ$BasChars$Lat[1],
                                                          Lat2 = staticData_peakFQ$BasChars$Lat[2],
                                                          Long2 = staticData_peakFQ$BasChars$Lat[2],
                                                          method = 2)
  )
  
  expect_error(Dist.WREG(Lat1 = staticData_peakFQ$BasChars$Lat[1],
                           Long1 = staticData_peakFQ$BasChars$Lat[1],
                           Lat2 = staticData_peakFQ$BasChars$Lat[2],
                           Long2 = staticData_peakFQ$BasChars$Lat[2],
                           method = 3)
                 )
  
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/Dist.WREG.out.rda")))
  expect_identical(intersiteDistance_meth1_test,intersiteDistance_meth1)
  expect_identical(intersiteDistance_meth2_test,intersiteDistance_meth2)
}
)

test_that("xcorrplot",{
  
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  
  # Build cross-correlation plot
  expect_silent(
    p1_test <- xcorPlot(object = staticData_peakFQ, alpha = 0.01, theta = 0.98,
                        concurrentMin = 10)
  )
  
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/xcorPlot.out.rda")))
  
  #expect_equal(p1_test,p1)
  
  
  
})



