context("Supporting function tests")



test_that("Dist.WREG",{
  # Import some example data
  expect_silent(
    
    load(paste0(system.file("testData", package = "WREG"),"/staticData_peakFQ.rda"))
  )
  
  # For two sites, compute the inter-site distance
  # Use the haversine approximation
  intersiteDistance_meth1_test <- Dist.WREG(Lat1 = staticData_peakFQ$BasChars$Lat[1],
                                            Long1 = staticData_peakFQ$BasChars$Lat[1],
                                            Lat2 = staticData_peakFQ$BasChars$Lat[2],
                                            Long2 = staticData_peakFQ$BasChars$Lat[2],
                                            method = 1)
  
  
  intersiteDistance_meth2_test <- Dist.WREG(Lat1 = staticData_peakFQ$BasChars$Lat[1],
                                            Long1 = staticData_peakFQ$BasChars$Lat[1],
                                            Lat2 = staticData_peakFQ$BasChars$Lat[2],
                                            Long2 = staticData_peakFQ$BasChars$Lat[2],
                                            method = 2)
  
  
  expect_silent(load(paste0(system.file("testData", package = "WREG"),"/Dist.WREG.out.rda")))
  expect_identical(intersiteDistance_meth1_test,intersiteDistance_meth1)
  expect_identical(intersiteDistance_meth2_test,intersiteDistance_meth2)
}
)


